#!/usr/bin/env python3
"""Check required programs and report versions without local executable paths."""
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys


def check():
    result = {"python": platform.python_version(), "platform": platform.system(),
              "architecture": platform.machine(), "programs": {}, "errors": [], "notes": []}
    if sys.version_info < (3, 12):
        result["errors"].append("Use Python 3.12 or newer; validation used Python 3.12.")
    if platform.system() == "Windows":
        result["errors"].append("Use Linux/Ubuntu in WSL; native Windows lacks the required POSIX file locking.")
    try:
        result["biopython"] = importlib.metadata.version("biopython")
        from Bio.Align import PairwiseAligner
        PairwiseAligner()
    except (ImportError, importlib.metadata.PackageNotFoundError) as exc:
        result["biopython"] = None
        result["errors"].append("Biopython with Bio.Align.PairwiseAligner is required.")
    if result["biopython"] and result["biopython"] != "1.88":
        result["notes"].append("The study refinement records report Biopython 1.88; retain the version used for this run.")
    for name, args in (("bwa", []), ("blastn", ["-version"]),
                       ("makeblastdb", ["-version"]), ("samtools", ["--version"])):
        if not shutil.which(name):
            result["programs"][name] = {"available": False, "version": None}
            result["errors"].append("Required program not found: " + name)
            continue
        try:
            process = subprocess.run([name, *args], capture_output=True,
                                     text=True, errors="replace", timeout=15)
            text = process.stdout + "\n" + process.stderr
            pattern = (r"(?m)^Version:\s*(\S+)" if name == "bwa" else
                       rf"(?m)^{re.escape(name)}:?\s+(\S+)")
            match = re.search(pattern, text)
            version = match.group(1) if match else None
            result["programs"][name] = {"available": version is not None, "version": version}
            if version is None:
                result["errors"].append("Program version could not be read: " + name)
        except (OSError, subprocess.TimeoutExpired):
            result["programs"][name] = {"available": False, "version": None}
            result["errors"].append("Program could not be executed: " + name)
    result["ready"] = not result["errors"]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional JSON version record; an existing file is not overwritten")
    args = parser.parse_args()
    result = check()
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8") as handle:
            handle.write(text)
    print(text, end="")
    return 0 if result["ready"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
