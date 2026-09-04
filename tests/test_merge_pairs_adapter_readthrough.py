"""Adapter read-through behavior of MERGE_PAIRS, which v3.4.0 enables by default.

Merged pairs became a default BLAST query class in v3.4.0, so whatever MERGE_PAIRS
emits is submitted to BLAST rather than being an opt-in experiment. That makes one
property load-bearing: a pair whose insert is shorter than the read length reads
through into adapter, and the merged product must not carry that adapter.

bbmerge already guarantees this without any adapter flag -- it trims sequence to
the right of the detected overlap. These tests pin that behavior so a future flag
or version change cannot silently start feeding adapter sequence to BLAST.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BBMAP_MODULE = ROOT / "modules" / "bbmap"

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None
    or shutil.which("bbmerge.sh") is None
    or shutil.which("reformat.sh") is None,
    reason="MERGE_PAIRS behavior requires Nextflow and locked BBTools",
)

ILLUMINA_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
INSERT = (
    "ACGTTGCAGATCGTAGCTAGGCTAACGATCGTACGGCATTCGATGCTAGCATCGTACGTAGGCATTCAGGATCCA"
)
READ_LENGTH = 100


def reverse_complement(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))


def write_read_through_pair(path: Path) -> None:
    """Write one interleaved pair whose insert is shorter than the read length."""
    read1 = (INSERT + ILLUMINA_ADAPTER)[:READ_LENGTH].ljust(READ_LENGTH, "A")
    read2 = (reverse_complement(INSERT) + ILLUMINA_ADAPTER)[:READ_LENGTH].ljust(
        READ_LENGTH,
        "A",
    )
    quality = "I" * READ_LENGTH
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"@pair_1 1:N:0\n{read1}\n+\n{quality}\n")
        handle.write(f"@pair_1 2:N:0\n{read2}\n+\n{quality}\n")


def run_merge_pairs(tmp_path: Path) -> tuple[str, str]:
    """Run the real MERGE_PAIRS process; return merged and unmerged FASTQ text."""
    reads = tmp_path / "sample_A.interleaved.fastq.gz"
    write_read_through_pair(reads)

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ MERGE_PAIRS }} from '{BBMAP_MODULE}'

workflow {{
    MERGE_PAIRS(
        Channel.of(tuple(
            'sample_A',
            'illumina',
            'interleaved',
            file('{reads}'),
        )),
    )
}}
""",
        encoding="utf-8",
    )

    completed = subprocess.run(  # noqa: S603
        [shutil.which("nextflow"), "-C", "/dev/null", "run", str(workflow)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        check=False,
    )
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics

    merged = sorted(tmp_path.glob("work/**/sample_A.overlap_merged_pair.fastq.gz"))
    unmerged = sorted(tmp_path.glob("work/**/sample_A.single_read.fastq.gz"))
    assert merged, diagnostics
    assert unmerged, diagnostics

    with gzip.open(merged[-1], "rt", encoding="utf-8") as handle:
        merged_text = handle.read()
    with gzip.open(unmerged[-1], "rt", encoding="utf-8") as handle:
        unmerged_text = handle.read()
    return merged_text, unmerged_text


def test_merged_pair_recovers_the_insert_without_adapter(tmp_path: Path) -> None:
    """A read-through pair merges back to exactly its insert, adapter trimmed."""
    merged_text, _ = run_merge_pairs(tmp_path)

    sequences = merged_text.splitlines()[1::4]
    assert sequences == [INSERT]
    assert ILLUMINA_ADAPTER[:15] not in merged_text


def test_read_through_pairs_do_not_survive_into_the_unmerged_batch(
    tmp_path: Path,
) -> None:
    """Read-through implies overlap, so such pairs never reach the single_read class.

    This is why unmerged reads need no separate adapter handling: an insert short
    enough to read into adapter is short enough for its mates to overlap.
    """
    _, unmerged_text = run_merge_pairs(tmp_path)

    assert unmerged_text.strip() == ""
