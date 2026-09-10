"""Tests for the FASTA list schema validate_labkey.py checks before upload."""

import ast
import re
from pathlib import Path

from prepare_blast_labkey import COLUMN_RENAMES, LABKEY_INPUT_COLUMNS
from validate_labkey import blast_fasta_fields

ROOT = Path(__file__).resolve().parents[1]


def _prepared_fasta_fieldnames() -> list[str]:
    """Columns LABKEY_PREPARE_FASTA writes, read out of the process script.

    The prep step is inline Python inside modules/labkey.nf, so there is no
    import seam. Reading the literal keeps the pipeline's writer and the
    validator's expectations pinned to each other.
    """
    text = (ROOT / "modules" / "labkey.nf").read_text(encoding="utf-8")
    match = re.search(r"fieldnames = (\[[^\]]*\])", text, re.DOTALL)
    assert match is not None, "LABKEY_PREPARE_FASTA fieldnames literal not found"
    return ast.literal_eval(re.sub(r"\s+", " ", match.group(1)))


def _dummy_row_keys() -> set[str]:
    """Keys of the blast_fasta dummy row validate_labkey.py inserts."""
    text = (ROOT / "bin" / "validate_labkey.py").read_text(encoding="utf-8")
    block = re.search(
        r'elif args\.type == "blast_fasta":\s*dummy = \{(.*?)\}',
        text,
        re.DOTALL,
    )
    assert block is not None, "blast_fasta dummy row not found"
    return set(re.findall(r'"([^"]+)":', block.group(1)))


def test_prepared_fasta_uses_query_scoped_column_names() -> None:
    """The list holds every queried class, so contig-scoped names are wrong."""
    fieldnames = _prepared_fasta_fieldnames()

    assert "qseqid" in fieldnames
    assert "query_sequence" in fieldnames
    assert "contig_id" not in fieldnames
    assert "contig_sequence" not in fieldnames


def test_validator_expects_the_same_query_scoped_columns() -> None:
    """A writer/validator mismatch fails the run, so pin them together."""
    assert "Qseqid" in blast_fasta_fields
    assert "Query Sequence" in blast_fasta_fields
    assert "Contig Id" not in blast_fasta_fields
    assert "Contig Sequence" not in blast_fasta_fields


def test_qseqid_matches_the_hits_list_column_name() -> None:
    """The two lists join on this column, so the name must agree exactly."""
    assert "qseqid" in LABKEY_INPUT_COLUMNS
    assert "qseqid" not in COLUMN_RENAMES
    assert "qseqid" in _prepared_fasta_fieldnames()


def test_dummy_row_covers_the_validated_schema() -> None:
    """The dummy insert proves the real schema accepts a row, so it must be complete."""
    assert _dummy_row_keys() == set(blast_fasta_fields)
