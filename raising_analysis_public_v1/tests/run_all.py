#!/usr/bin/env python3
"""Run isolated synthetic test suites; optionally exercise real external tools."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SUITES = (
    "pipelines/fourway/tests/test_batch.py",
    "pipelines/junction/raising_junction_pilot_v1/tests/test_junction.py",
    "pipelines/junction/raising_junction_refine_v1_1/tests/test_refine.py",
    "pipelines/junction/tests/test_batch.py",
    "tests/test_figure_data.py",
    "tests/test_demo.py",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-external", action="store_true",
                        help="Also run the small synthetic FASTQ-to-Figure round trip with actual BWA/BLAST/samtools")
    args = parser.parse_args()
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    for suite in SUITES:
        print("Running " + suite, flush=True)
        path = ROOT / suite
        result = subprocess.run([sys.executable, "-B", "-m", "unittest", "discover",
                                 "-s", str(path.parent), "-p", path.name],
                                cwd=ROOT, env=environment)
        if result.returncode:
            return result.returncode
    if args.include_external:
        return subprocess.run([sys.executable, "-B", str(ROOT / "tests/test_external_roundtrip.py")],
                              cwd=ROOT, env=environment).returncode
    print("All synthetic/unit suites passed. External-aligner round trip was not requested.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
