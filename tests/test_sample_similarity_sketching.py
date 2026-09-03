"""SAMPLE_SIMILARITY_QC owns its sketching after rapid screening was removed.

The sketch step used to live in RAPID_SCREENING, which passed its sketches to
SAMPLE_SIMILARITY_QC. Removing rapid screening leaves similarity QC as the only
consumer, so it now takes post-QC read batches directly and sketches them
itself.

These tests pin the channel contract of that move -- that the subworkflow
accepts the profiled_batches_by_sample shape and reaches the sketch step, and
that it still skips samples the profiler counted as empty. sourmash and seqkit
are faked so this stays a wiring test rather than a tool test.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SAMPLE_SIMILARITY_QC = ROOT / "subworkflows" / "sample_similarity_qc"

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None,
    reason="Subworkflow wiring requires Nextflow",
)


def write_executable(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)


def write_fakes(bin_dir: Path) -> None:
    """Stand in for seqkit and sourmash so only the wiring is under test."""
    write_executable(
        bin_dir / "seqkit",
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
out = args[args.index("--out-file") + 1]
Path(out).write_bytes(b"")
print(f"FAKE_SEQKIT {out}", file=sys.stderr)
""",
    )
    # The similarity chain downstream of sketching is unchanged by the rehoming;
    # fake it so a missing numpy/plotting stack cannot mask a wiring failure.
    for script, outputs in (
        ("sample_sketch_distances.py", ("sample_sketch_distances.{metric}.tsv",)),
        (
            "sample_ordination.py",
            (
                "sample_distance_matrix.{metric}.tsv",
                "sample_ordination.{metric}.tsv",
                "sample_ordination_variance.{metric}.tsv",
            ),
        ),
        (
            "report_possible_sample_mixups.py",
            (
                "nearest_sample_neighbors.{metric}.tsv",
                "possible_sample_mixup_candidates.{metric}.tsv",
            ),
        ),
        (
            "sample_similarity_candidate_evidence.py",
            ("sample_similarity_candidate_evidence.tsv",),
        ),
        (
            "plot_sample_ordination.py",
            ("sample_ordination.{metric}.html", "sample_ordination.{metric}.png"),
        ),
    ):
        rendered = "\n".join(f'    "{name}",' for name in outputs)
        write_executable(
            bin_dir / script,
            f"""#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
metric = args[args.index("--metric") + 1] if "--metric" in args else "abund"
for template in (
{rendered}
):
    Path(template.format(metric=metric)).write_text("", encoding="utf-8")
""",
        )

    write_executable(
        bin_dir / "sourmash",
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "-o" in args:
    out = args[args.index("-o") + 1]
    Path(out).write_bytes(b"")
if "--output" in args:
    Path(args[args.index("--output") + 1]).write_bytes(b"")
if "--csv" in args:
    Path(args[args.index("--csv") + 1]).write_bytes(b"")
print(f"FAKE_SOURMASH {' '.join(args[:2])}", file=sys.stderr)
""",
    )


def run_workflow(tmp_path: Path, batches: str) -> subprocess.CompletedProcess[str]:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fakes(bin_dir)

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ SAMPLE_SIMILARITY_QC }} from '{SAMPLE_SIMILARITY_QC}'

params.sourmash_ksize = 31
params.sourmash_scaled = 50

workflow {{
    SAMPLE_SIMILARITY_QC({batches})

    SAMPLE_SIMILARITY_QC.out.query_sketches.view {{ sketch ->
        "SKETCH: ${{sketch[0]}}"
    }}
}}
""",
        encoding="utf-8",
    )

    environment = os.environ.copy()
    environment["PATH"] = f"{bin_dir}{os.pathsep}{environment['PATH']}"
    nextflow = shutil.which("nextflow")
    assert nextflow is not None
    return subprocess.run(  # noqa: S603
        [nextflow, "-C", "/dev/null", "run", str(workflow)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        env=environment,
        check=False,
    )


def batch_tuple(tmp_path: Path, sample_id: str, sequence_count: int) -> str:
    reads = tmp_path / f"{sample_id}.fastq.gz"
    reads.write_bytes(b"")
    return (
        f"tuple("
        f"[id: '{sample_id}', platform: 'illumina', read_structure: 'single', "
        f"sequence_count: {sequence_count}], "
        f"[[meta: [query_class: 'single_read'], reads: file('{reads}')]])"
    )


def test_sketches_post_qc_batches_without_rapid_screening(tmp_path: Path) -> None:
    """The subworkflow takes profiled batches directly and sketches them."""
    batches = batch_tuple(tmp_path, "sample_A", 10)
    completed = run_workflow(tmp_path, f"Channel.of({batches})")
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"

    assert completed.returncode == 0, diagnostics
    assert "SKETCH: sample_A" in completed.stdout, diagnostics


def test_empty_samples_are_not_sketched(tmp_path: Path) -> None:
    """A sample the profiler counted as empty cannot produce a sketch."""
    empty = batch_tuple(tmp_path, "empty_sample", 0)
    completed = run_workflow(tmp_path, f"Channel.of({empty})")
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"

    assert completed.returncode == 0, diagnostics
    # No sketch is emitted. Nextflow lists every process in its progress display
    # whether or not it ran any task, so absence of the emit is the signal here,
    # not absence of the process name.
    assert "SKETCH:" not in completed.stdout, diagnostics
    assert "FAKE_SEQKIT" not in completed.stderr, diagnostics
