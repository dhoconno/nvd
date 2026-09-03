"""Shared LabKey list I/O helpers for the BLAST hits and FASTA upload scripts.

Both upload scripts dedup against their destination LabKey list (the list is its
own completion ledger). This module centralizes the two pieces of client logic
they would otherwise duplicate:

- ``rows_present`` — a presence check with the ``QueryFilter`` -> filter-string
  fallback so it still works against older labkey-api-python versions.
- ``insert_records`` — an insert that deliberately lets failures propagate so
  callers can hard-fail rather than log-and-continue, optionally splitting the
  payload into fixed-size batches.

Batching exists because unassembled reads became BLAST query classes by default
in v3.4.0. A ``(sample_id, query_class)`` unit used to mean contigs only; it now
also covers read-scale hit tables, and sending one of those as a single request
can hang the LabKey server.

The cost is atomicity. A failed batch leaves earlier batches committed, and the
presence-keyed ledger cannot tell partial state from complete state, so a retry
would skip the unit and silently drop the remainder. ``PartialInsertError``
exists to make callers report that instead of letting it pass quietly.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    from collections.abc import Mapping


class LabKeyQueryApi(Protocol):
    """The subset of the LabKey APIWrapper ``query`` interface used here."""

    def select_rows(self, **kwargs: object) -> dict[str, object]: ...

    def insert_rows(self, **kwargs: object) -> object: ...


def rows_present(
    query_api: LabKeyQueryApi,
    schema: str,
    table: str,
    filters: Mapping[str, object],
) -> bool:
    """Return whether the list holds at least one row matching every filter.

    ``filters`` maps column name to the value it must equal; all are ANDed. The
    destination list is its own completion ledger, so any matching row means the
    unit was already uploaded and the caller should skip re-inserting it.

    Uses ``QueryFilter`` objects when available, falling back to filter-string
    syntax (``col~eq=value``) on labkey-api-python versions that lack them.
    """
    try:
        # Imported lazily so the filter-string fallback can handle older
        # labkey-api-python versions that do not ship QueryFilter.
        from labkey.query import QueryFilter  # noqa: PLC0415

        result = query_api.select_rows(
            schema_name=schema,
            query_name=table,
            filter_array=[
                QueryFilter(column, value, "eq") for column, value in filters.items()
            ],
            max_rows=1,
        )
    except (ImportError, AttributeError):
        result = query_api.select_rows(
            schema_name=schema,
            query_name=table,
            filter_array=[f"{column}~eq={value}" for column, value in filters.items()],
        )
    return bool(result and result.get("rows"))


class PartialInsertError(RuntimeError):
    """A batched insert failed partway, leaving earlier batches committed.

    Carries ``rows_committed`` so the caller can tell the operator exactly how
    much survived. That matters because the destination list is its own
    completion ledger keyed on row presence: rows left behind by a failed run
    make a later retry treat the unit as already uploaded and skip the rest.
    """

    def __init__(self, rows_committed: int, total_rows: int, cause: object) -> None:
        self.rows_committed = rows_committed
        self.total_rows = total_rows
        super().__init__(
            f"LabKey insert failed after committing {rows_committed} of "
            f"{total_rows} rows: {cause}",
        )


def insert_records(
    query_api: LabKeyQueryApi,
    schema: str,
    table: str,
    records: list[dict[str, object]],
    *,
    batch_size: int | None = None,
) -> None:
    """Insert records, optionally splitting them across several calls.

    Without ``batch_size`` this is a single atomic call, as it has always been.
    With one, records are sent in fixed-size batches: read-derived query classes
    produce payloads large enough that a single call can hang the LabKey server.

    Batching trades atomicity for a bounded request size. When a later batch
    fails the earlier ones stay committed, and that raises
    ``PartialInsertError`` to force callers to report the partial state.
    Nothing is rolled back; recovering means clearing the unit's rows before
    re-running. A failure in the very first batch commits nothing, so the
    original exception propagates unwrapped and callers can still honestly
    report that no rows landed.
    """
    if batch_size is None:
        query_api.insert_rows(
            schema_name=schema,
            query_name=table,
            rows=records,
        )
        return

    if batch_size < 1:
        msg = f"batch_size must be at least 1, got {batch_size}"
        raise ValueError(msg)

    rows_committed = 0
    for start in range(0, len(records), batch_size):
        batch = records[start : start + batch_size]
        try:
            query_api.insert_rows(
                schema_name=schema,
                query_name=table,
                rows=batch,
            )
        except Exception as error:
            if rows_committed == 0:
                # Nothing landed, so the pre-batching contract still holds
                # exactly: the caller can report "no rows committed" and the
                # original LabKey error reaches the operator unwrapped.
                raise
            raise PartialInsertError(rows_committed, len(records), error) from error
        rows_committed += len(batch)
