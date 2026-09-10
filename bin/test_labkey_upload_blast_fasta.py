"""Tests for deduping the query FASTA LabKey upload by (experiment, sample_id, query_class)."""

import sys

import labkey.api_wrapper
import pytest
from labkey_upload_blast_fasta import combo_already_uploaded, insert_records, main


class _FakeQuery:
    def __init__(self, rows):
        self._rows = rows
        self.calls = []

    def select_rows(self, **kwargs):
        self.calls.append(kwargs)
        return {"rows": self._rows}


def test_combo_present_is_detected() -> None:
    q = _FakeQuery([{"Key": 1}])
    assert combo_already_uploaded(q, "lists", "fasta", 7, "s1", "single_read") is True


def test_combo_absent_is_detected() -> None:
    q = _FakeQuery([])
    assert combo_already_uploaded(q, "lists", "fasta", 7, "s1", "single_read") is False
    (call,) = q.calls
    assert any("experiment" in str(f) for f in call["filter_array"])
    assert any("sample_id" in str(f) for f in call["filter_array"])


def test_guard_keys_on_query_class() -> None:
    """A sample's contigs being present must not mask an unuploaded read class.

    The list is its own completion ledger. Keyed on (experiment, sample_id)
    alone, the contig batch would make every later query class for that sample
    look already uploaded, and its rows would be dropped without an error.
    """
    q = _FakeQuery([])
    assert (
        combo_already_uploaded(q, "lists", "fasta", 7, "s1", "overlap_merged_pair")
        is False
    )
    (call,) = q.calls
    assert any("query_class" in str(f) for f in call["filter_array"])


def test_insert_records_propagates_api_failure() -> None:
    """A failed insert_rows call must not be swallowed by insert_records."""

    class _FailingQuery:
        def insert_rows(self, **kwargs):
            raise RuntimeError("insert boom")

    with pytest.raises(RuntimeError, match="insert boom"):
        insert_records(_FailingQuery(), "lists", "fasta", [{"experiment": 1}])


class _FakeQueryAPI:
    """Stand-in for the LabKey query API: configurable presence check + insert."""

    def __init__(self, rows, insert_error=None):
        self._rows = rows
        self._insert_error = insert_error
        self.inserted_rows = None

    def select_rows(self, **kwargs):
        return {"rows": self._rows}

    def insert_rows(self, **kwargs):
        if self._insert_error is not None:
            raise self._insert_error
        self.inserted_rows = kwargs.get("rows")
        return {"rows": [{"Key": 1}]}


def _fake_api_wrapper_class(fake_query: "_FakeQueryAPI"):
    """Build a fake APIWrapper class whose .query is the given fake, never the network."""

    class _FakeAPIWrapper:
        def __init__(self, *args, **kwargs):
            self.query = fake_query

    return _FakeAPIWrapper


def _cli_argv(table_name: str = "fasta") -> list[str]:
    return [
        "labkey_upload_blast_fasta.py",
        "--experiment-id",
        "7",
        "--sample-id",
        "s1",
        "--query-class",
        "single_read",
        "--labkey-server",
        "https://example.org",
        "--labkey-project-name",
        "proj",
        "--labkey-api-key",
        "key",
        "--labkey-schema",
        "lists",
        "--table-name",
        table_name,
    ]


def test_main_exits_nonzero_when_insert_fails(tmp_path, monkeypatch, capsys) -> None:
    """A failed atomic insert must hard-fail the run, not report success."""
    csv_path = tmp_path / "s1_fasta.csv"
    csv_path.write_text(
        "experiment,sample_id,contig_id,contig_sequence,notes,nextflow_run_id\n"
        "7,s1,contig_1,ACGT,,run1\n",
    )

    fake_query = _FakeQueryAPI(rows=[], insert_error=RuntimeError("insert boom"))
    monkeypatch.setattr(
        labkey.api_wrapper,
        "APIWrapper",
        _fake_api_wrapper_class(fake_query),
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", _cli_argv())

    with pytest.raises(SystemExit) as exc_info:
        main()

    assert exc_info.value.code == 1
    log_text = (tmp_path / "fasta_labkey_upload.log").read_text()
    assert "insert boom" in log_text
    assert "FASTA UPLOAD COMPLETE" not in log_text
    assert "insert boom" in capsys.readouterr().err


def test_main_skip_present_sample_exits_zero(tmp_path, monkeypatch) -> None:
    """A unit already present in the destination list is a no-op success, not a failure."""
    csv_path = tmp_path / "s1_fasta.csv"
    csv_path.write_text(
        "experiment,sample_id,contig_id,contig_sequence,notes,nextflow_run_id\n"
        "7,s1,contig_1,ACGT,,run1\n",
    )

    fake_query = _FakeQueryAPI(rows=[{"Key": 1}])
    monkeypatch.setattr(
        labkey.api_wrapper,
        "APIWrapper",
        _fake_api_wrapper_class(fake_query),
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", _cli_argv())

    main()  # must return normally (implicit exit 0), never raise SystemExit

    log_text = (tmp_path / "fasta_labkey_upload.log").read_text()
    assert "SKIP: combo already uploaded" in log_text
    assert "query_class=single_read" in log_text
    assert fake_query.inserted_rows is None


class _LedgerQueryAPI:
    """Query API backed by rows that actually honour the filters it is given.

    The destination list is its own completion ledger, and LabKey returns every
    row matching the filters supplied. A fake that ignored filters could not
    distinguish "this unit is uploaded" from "some other unit for this sample
    is", which is the exact distinction these tests exist to pin.
    """

    def __init__(self, rows):
        self._rows = rows
        self.inserted_rows = None

    def select_rows(self, **kwargs):
        wanted = {f.column_name: f.value for f in kwargs["filter_array"]}
        matches = [
            row
            for row in self._rows
            if all(row.get(column) == value for column, value in wanted.items())
        ]
        return {"rows": matches[:1]}

    def insert_rows(self, **kwargs):
        self.inserted_rows = kwargs.get("rows")
        return {"rows": [{"Key": 1}]}


def test_read_class_uploads_when_contigs_are_already_present(
    tmp_path,
    monkeypatch,
) -> None:
    """A sample's uploaded contigs must not suppress its read query classes."""
    csv_path = tmp_path / "s1_fasta.csv"
    csv_path.write_text(
        "experiment,sample_id,query_class,contig_id,contig_sequence,notes,"
        "nextflow_run_id\n"
        "7,s1,single_read,nvdReadQuery_s1_000001,ACGT,,run1\n",
    )

    fake_query = _LedgerQueryAPI(
        rows=[
            {
                "experiment": 7,
                "sample_id": "s1",
                "query_class": "short_assembly_contig",
            },
        ],
    )
    monkeypatch.setattr(
        labkey.api_wrapper,
        "APIWrapper",
        _fake_api_wrapper_class(fake_query),
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", _cli_argv())

    main()

    assert fake_query.inserted_rows is not None, (
        "single_read rows were dropped because this sample's contigs were "
        "already in the list"
    )


class _PartialInsertQueryAPI:
    """Commits the first batch, then fails — the partial-commit path."""

    def __init__(self):
        self.batches_committed = 0

    def select_rows(self, **kwargs):
        return {"rows": []}

    def insert_rows(self, **kwargs):
        if self.batches_committed >= 1:
            raise RuntimeError("insert boom")
        self.batches_committed += 1
        return {"rows": [{"Key": 1}]}


def test_partial_insert_guidance_names_the_query_class(
    tmp_path,
    monkeypatch,
    capsys,
) -> None:
    """Committed-but-incomplete rows must identify the unit to clear.

    The ledger is keyed on (experiment, sample_id, query_class), so an operator
    told only the sample would not know which of its batches to delete before
    re-running.
    """
    csv_path = tmp_path / "s1_fasta.csv"
    csv_path.write_text(
        "experiment,sample_id,query_class,contig_id,contig_sequence,notes,"
        "nextflow_run_id\n"
        "7,s1,single_read,nvdReadQuery_s1_000001,ACGT,,run1\n"
        "7,s1,single_read,nvdReadQuery_s1_000002,TTGC,,run1\n",
    )

    fake_query = _PartialInsertQueryAPI()
    monkeypatch.setattr(
        labkey.api_wrapper,
        "APIWrapper",
        _fake_api_wrapper_class(fake_query),
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [*_cli_argv(), "--insert-batch-size", "1"])

    with pytest.raises(SystemExit) as excinfo:
        main()
    assert excinfo.value.code == 1

    log_text = (tmp_path / "fasta_labkey_upload.log").read_text()
    assert "query_class=single_read" in log_text
    assert "query_class=single_read" in capsys.readouterr().err
