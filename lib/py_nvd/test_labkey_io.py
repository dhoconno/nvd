"""Tests for the shared LabKey list I/O helpers."""

from __future__ import annotations

import pytest

from py_nvd.labkey_io import PartialInsertError, insert_records, rows_present


class _FakeQuery:
    def __init__(
        self,
        rows: list[dict[str, object]] | None = None,
        *,
        raise_on_insert: bool = False,
        raise_on_insert_call: int | None = None,
    ) -> None:
        self._rows = rows or []
        self.raise_on_insert = raise_on_insert
        self.raise_on_insert_call = raise_on_insert_call
        self.select_calls: list[dict[str, object]] = []
        self.insert_calls: list[dict[str, object]] = []

    def select_rows(self, **kwargs: object) -> dict[str, object]:
        self.select_calls.append(kwargs)
        return {"rows": self._rows}

    def insert_rows(self, **kwargs: object) -> dict[str, object]:
        self.insert_calls.append(kwargs)
        if self.raise_on_insert or len(self.insert_calls) == self.raise_on_insert_call:
            msg = "insert failed"
            raise RuntimeError(msg)
        return {"rows": kwargs["rows"]}


def test_rows_present_true_when_any_row_returned() -> None:
    query = _FakeQuery(rows=[{"Key": 1}])
    assert (
        rows_present(query, "lists", "t", {"experiment": 7, "sample_id": "s1"}) is True
    )


def test_rows_present_false_when_empty() -> None:
    query = _FakeQuery(rows=[])
    assert (
        rows_present(query, "lists", "t", {"experiment": 7, "sample_id": "s1"}) is False
    )


def test_rows_present_builds_one_eq_filter_per_column() -> None:
    query = _FakeQuery(rows=[])
    rows_present(
        query,
        "lists",
        "t",
        {"experiment": 7, "sample_id": "s1", "query_class": "single_read"},
    )
    (call,) = query.select_calls
    filter_array = call["filter_array"]
    assert isinstance(filter_array, list)
    assert len(filter_array) == 3
    rendered = " ".join(str(f) for f in filter_array)
    assert "experiment" in rendered
    assert "sample_id" in rendered
    assert "query_class" in rendered


def test_insert_records_inserts_all_rows_in_one_call() -> None:
    query = _FakeQuery()
    records = [{"a": 1}, {"a": 2}]
    insert_records(query, "lists", "t", records)
    (call,) = query.insert_calls
    assert call["rows"] == records
    assert call["schema_name"] == "lists"
    assert call["query_name"] == "t"


def test_insert_records_propagates_api_failure() -> None:
    query = _FakeQuery(raise_on_insert=True)
    with pytest.raises(RuntimeError):
        insert_records(query, "lists", "t", [{"a": 1}])


def test_insert_records_splits_large_payloads_into_batches() -> None:
    """Read-derived query classes make single-call inserts large enough to hang."""
    query = _FakeQuery()
    records = [{"a": index} for index in range(250)]

    insert_records(query, "lists", "t", records, batch_size=100)

    assert [len(call["rows"]) for call in query.insert_calls] == [100, 100, 50]
    sent = [row for call in query.insert_calls for row in call["rows"]]
    assert sent == records


def test_insert_records_without_batch_size_stays_one_call() -> None:
    """Callers that have not opted in keep the previous single-call behavior."""
    query = _FakeQuery()
    records = [{"a": index} for index in range(250)]

    insert_records(query, "lists", "t", records)

    assert len(query.insert_calls) == 1


def test_failed_batch_reports_how_many_rows_were_committed() -> None:
    """Partial state is not rolled back, so the count must reach the operator.

    The destination list is its own completion ledger keyed on presence, so rows
    left behind by a failed run make a retry skip the unit entirely. The caller
    can only warn about that if it knows how much was already written.
    """
    query = _FakeQuery(raise_on_insert_call=3)
    records = [{"a": index} for index in range(250)]

    with pytest.raises(PartialInsertError) as excinfo:
        insert_records(query, "lists", "t", records, batch_size=100)

    assert excinfo.value.rows_committed == 200
    assert excinfo.value.total_rows == 250
    assert isinstance(excinfo.value.__cause__, RuntimeError)
    # The underlying LabKey error is what an operator needs to act on, so it
    # must survive being wrapped.
    assert "insert failed" in str(excinfo.value)


def test_first_batch_failure_propagates_the_original_error() -> None:
    """Nothing committed means the pre-batching contract still holds exactly.

    Callers report "no rows committed" on this path, which stays true, and the
    raw LabKey error reaches the operator rather than a partial-state warning
    about rows that do not exist.
    """
    query = _FakeQuery(raise_on_insert_call=1)

    with pytest.raises(RuntimeError) as excinfo:
        insert_records(query, "lists", "t", [{"a": 1}], batch_size=100)

    assert not isinstance(excinfo.value, PartialInsertError)
    assert "insert failed" in str(excinfo.value)
