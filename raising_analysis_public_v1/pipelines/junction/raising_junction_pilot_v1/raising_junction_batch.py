#!/usr/bin/env python3
"""Run the junction pilot sequentially on completed full-depth libraries."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--list", action="store_true", help="Show complete input libraries only")
    p.add_argument("--skip-complete", action="store_true", help="Reuse this exact pilot's completed runs")
    args = p.parse_args(argv)
    root = args.root.expanduser().resolve()
    here = Path(__file__).resolve().parent
    source = {name: sha256(here / name) for name in ("raising_junction_pilot.py", "junction_core.py")}
    runs = {}
    for suffix in ("_fourway_full_v1", "_fourway_full_v1_1"):
        for manifest in sorted((root / "results").glob("*" + suffix + "/run_summary.json")):
            s = json.loads(manifest.read_text(encoding="utf-8"))
            if s.get("status") != "complete":
                print("Skipping incomplete upstream run: " + str(manifest), flush=True)
                continue
            if s["sample"] in runs:
                raise ValueError("More than one complete result for sample " + s["sample"])
            runs[s["sample"]] = (manifest, s)
    if not runs:
        raise ValueError("No complete full-depth libraries found")
    for sample, (manifest, s) in sorted(runs.items()):
        print(f"{sample}\t{s['side']}\t{s['total_pairs']} input pairs", flush=True)
    if args.list:
        return 0
    for number, (sample, (manifest, s)) in enumerate(sorted(runs.items()), 1):
        output = root / "results" / (sample + "_junction_pilot_v1")
        previous = output / "summary.json"
        if args.skip_complete and previous.is_file():
            old = json.loads(previous.read_text(encoding="utf-8"))
            matches = (old.get("analysis_type") == "raising_exact_anchor_junction_pilot" and
                       old.get("status") == "complete" and old.get("sample") == sample and
                       old.get("source_sha256") == source and
                       old.get("upstream_input_sha256", {}).get("run_summary.json") == sha256(manifest) and
                       old.get("config", {}).get("host_window_padding_bp") == 500)
            if not matches:
                raise ValueError("Existing output differs from current source/input; use a separate output for re-analysis: " + str(output))
            for name, digest in old.get("output_sha256", {}).items():
                if not (output / name).is_file() or sha256(output / name) != digest:
                    raise ValueError("Completed pilot output was changed: " + str(output / name))
            print(f"[{number}/{len(runs)}] Reusing completed pilot: {sample}", flush=True)
            continue
        print(f"[{number}/{len(runs)}] Running: {sample}", flush=True)
        subprocess.run([sys.executable, str(here / "raising_junction_pilot.py"), "--root", str(root),
                        "--sample", sample, "--results-dir", str(manifest.parent)], check=True)
    subprocess.run([sys.executable, str(here / "summarize_junctions.py"), "--root", str(root)], check=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        raise SystemExit(2)
