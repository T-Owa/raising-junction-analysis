#!/usr/bin/env python3
"""Actual external-aligner smoke test using a wholly synthetic host and reads.

Run explicitly after installing all requirements. This is not a biological
accuracy benchmark or a replay of the original full study FASTQs.
"""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def run(script, *arguments):
    subprocess.run([sys.executable, "-B", str(ROOT / script), *map(str, arguments)], check=True)


def main():
    run("tools/check_environment.py")
    with tempfile.TemporaryDirectory(prefix="raising_synthetic_validation_") as td:
        work = Path(td) / "demo"
        run("examples/make_demo.py", "--outdir", work)
        run("pipelines/fourway/raising_batch.py", "--root", work,
            "--host-fasta", work / "ref/synthetic_host.fa", "--threads", 2, "--allow-index-build")
        run("pipelines/junction/raising_junction_batch.py", "--root", work, "--expected-count", 2)
        batches = list((work / "results/junction_batch_v1").glob("*/batch_summary.json"))
        if len(batches) != 1:
            raise AssertionError("Expected one completed synthetic batch")
        figures = work / "figure-data"
        run("tools/prepare_figure_data.py", "--input-dir", batches[0].parent,
            "--output-dir", figures, "--metadata", work / "sample_metadata.csv")
        designed = json.loads((work / "expected_demo.json").read_text())
        with (work / "results/all_libraries_fourway_batch_v1/all_figure_counts.csv").open(encoding="utf-8-sig") as handle:
            counts = {(r["sample"], r["category"]): int(r["count"]) for r in csv.DictReader(handle)}
        for library in designed["libraries"]:
            for category, count in library["expected_fourway_counts"].items():
                if counts[(library["sample"], category)] != count:
                    raise AssertionError("Synthetic four-way category differs from its designed expectation")
        with (figures / "metrics.csv").open(encoding="utf-8-sig") as handle:
            rows = [r for r in csv.DictReader(handle) if r["subset"] == "nuclear_all"]
        if len(rows) != 4:
            raise AssertionError("Expected two libraries times two support thresholds")
        expected = {"N_patterns": 4, "N_direct": 3, "N_insertion_candidate": 1,
                    "N_direct_mh0": 2, "N_direct_mh1": 0, "N_direct_mh2plus": 1,
                    "N_direct_deletion_min_gt0": 1, "N_direct_deletion_zero_exact": 1,
                    "N_direct_deletion_ambiguous_zero": 1}
        for row in rows:
            for key, value in expected.items():
                if int(row[key]) != value:
                    raise AssertionError("Synthetic junction metric differs: " + key)
        with (figures / "candidate_annotations.csv").open(encoding="utf-8-sig") as handle:
            candidates = list(csv.DictReader(handle))
        insertions = [r for r in candidates if int(r["insertion_length_bp"]) > 0]
        if len(insertions) != 2 or any(r["mh_for_analysis_bp"] != "" for r in insertions):
            raise AssertionError("Insertion-interface MH must remain unassessed")
        print("External round trip passed: both sides, four categories, donor loss, MH, insertion-interface MH=NA.")


if __name__ == "__main__":
    main()
