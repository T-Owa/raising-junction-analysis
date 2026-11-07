#!/usr/bin/env python3
"""Sequential, resumable exact-pilot and v0.1.1 junction reassessment.

Uses completed full-depth four-way results, never sampled three-way results.
All existing outputs are reused only after provenance/content checks. Nothing
is overwritten by this wrapper. Counts describe read pairs/candidate patterns,
not independent integration events. The summary contains no new MH calls.
"""
import argparse
import ast
import collections
import csv
import datetime
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import uuid

VERSION = "1.0.0"
HERE = Path(__file__).resolve().parent
SUFFIXES = ("_fourway_full_v1", "_fourway_full_v1_1")
PILOT_TYPE = "raising_exact_anchor_junction_pilot"
REFINE_TYPE = "raising_existing_template_gapped_reassessment"
STATUS_KEYS = ("template_supported", "template_ambiguous", "alignment_ambiguous",
               "junction_variant_review", "no_qualifying_template", "no_nearby_template")


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def atomic_json(path, value):
    temporary = Path(str(path) + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def static_value(path, name):
    for node in ast.parse(Path(path).read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise ValueError(f"Missing static {name}: {path}")


class HashCache:
    """One digest per unchanged path during this invocation; never trust mtime alone."""
    def __init__(self):
        self.cache = {}

    def __call__(self, path):
        path = Path(path).resolve()
        stat = path.stat()
        key = (str(path), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
        if key not in self.cache:
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                    digest.update(block)
            stable = lambda x: (x.st_dev, x.st_ino, x.st_size, x.st_mtime_ns, x.st_ctime_ns)
            if stable(path.stat()) != stable(stat):
                raise ValueError("File changed while hashing: " + str(path))
            self.cache[key] = digest.hexdigest()
        return self.cache[key]


def safe_sample(sample):
    if not isinstance(sample, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", sample) or sample in (".", ".."):
        raise ValueError("Unsafe sample identifier: " + repr(sample))


def discover(root, samples=None, preference=None):
    choices = collections.defaultdict(list)
    for path in sorted((root / "results").glob("*/run_summary.json")):
        suffix = next((s for s in reversed(SUFFIXES) if path.parent.name.endswith(s)), None)
        if suffix is None:
            continue
        manifest = load_json(path)
        sample = manifest.get("sample")
        safe_sample(sample)
        if path.parent.name != sample + suffix:
            raise ValueError("Full-run folder/sample mismatch: " + str(path))
        if samples is not None and sample not in samples:
            continue
        if manifest.get("status") != "complete" or not manifest.get("inputs", {}).get("lanes"):
            raise ValueError("Incomplete/non-full-depth result: " + str(path))
        if manifest.get("side") not in ("3p", "5p"):
            raise ValueError("Missing library side: " + str(path))
        validation = manifest.get("validation", {})
        if not validation or not all(v is True for v in validation.values()):
            raise ValueError("Full-run validation is absent or failed: " + str(path))
        choices[sample].append({"sample": sample, "side": manifest["side"],
            "directory": path.parent.resolve(), "manifest_path": path.resolve(),
            "manifest": manifest, "suffix": suffix})
    if samples is not None and set(samples) - set(choices):
        raise ValueError("No completed full-run results for: " + ", ".join(sorted(set(samples) - set(choices))))
    if not choices:
        raise ValueError("No completed *_fourway_full_v1[_1]/run_summary.json files found")
    libraries = []
    for sample, items in sorted(choices.items()):
        if len(items) > 1:
            preferred = [x for x in items if x["suffix"] == "_fourway_full_" + str(preference)] if preference else []
            if len(preferred) != 1:
                raise ValueError("Multiple full-run versions for " + sample + ": " +
                    ", ".join(str(x["directory"]) for x in items) +
                    ". Select deliberately with --prefer-fullrun v1 or v1_1.")
            chosen = preferred[0]
            chosen["excluded_fullrun_directories"] = [str(x["directory"]) for x in items if x is not chosen]
        else:
            chosen = items[0]
        libraries.append(chosen)
    return libraries


def find_engine(explicit, root, folder, script):
    candidates = [explicit] if explicit else [HERE / folder / script, root / "scripts" / folder / script]
    for candidate in candidates:
        path = Path(candidate).expanduser().resolve()
        if path.is_file():
            return path
    raise ValueError("Engine script missing; tried " + ", ".join(str(x) for x in candidates))


def engine_metadata(pilot, refine, sha):
    if static_value(pilot, "VERSION") != "0.1.0" or static_value(refine, "VERSION") != "0.1.1":
        raise ValueError("This batch requires exact pilot 0.1.0 and refinement 0.1.1")
    if static_value(pilot, "ANALYSIS_TYPE") != PILOT_TYPE or static_value(refine, "ANALYSIS_TYPE") != REFINE_TYPE:
        raise ValueError("Unexpected analysis engine type")
    pilot_sources = {n: sha(pilot.parent / n) for n in ("raising_junction_pilot.py", "junction_core.py")}
    refine_sources = {p.relative_to(refine.parent).as_posix(): sha(p)
                      for p in sorted(refine.parent.glob("*.py")) + sorted((refine.parent / "support").glob("*.py"))}
    pilot_config = {**static_value(pilot.parent / "junction_core.py", "DEFAULT_CONFIG"),
                    "host_window_padding_bp": 500, "complex_host_locus_distance_bp": 5000, "pair_merging": False}
    refine_config = {"cluster_total_span_bp": 10, "common_host_flank_bp": 120,
        "accepted_host_template_search_padding_bp": 20, "complex_host_distance_bp": 5000,
        "comparison_context_merge_padding_bp": 0, "old_template_retention_min_accepted_host_overlap_bp": 30,
        "alignment": static_value(refine.parent / "refine_alignment.py", "DEFAULT_CONFIG"),
        "pair_merging": False, "abundance_based_correction": False, "novel_template_discovery": False}
    return {"pilot_sources": pilot_sources, "refine_sources": refine_sources,
            "pilot_config": pilot_config, "refine_config": refine_config}


def evidence_hashes(library, sha):
    expected = library["manifest"].get("final_output_sha256", {})
    hashes = {"run_summary.json": sha(library["manifest_path"])}
    for stem in ("pair_classification", "read_evidence", "host_alignment_evidence"):
        found = [library["directory"] / (stem + suffix) for suffix in (".tsv.gz", ".tsv")
                 if (library["directory"] / (stem + suffix)).is_file()]
        if len(found) != 1:
            raise ValueError("Need exactly one " + stem + ".tsv[.gz] in " + str(library["directory"]))
        hashes[found[0].name] = sha(found[0])
        if expected.get(found[0].name) != hashes[found[0].name]:
            raise ValueError("Full-run evidence missing/changed checksum: " + str(found[0]))
    return hashes


def validate_raw(manifest, sha):
    for lane in manifest["inputs"]["lanes"]:
        for mate in ("R1", "R2"):
            entry = lane[mate]
            if sha(Path(entry["path"]).expanduser()) != entry.get("sha256"):
                raise ValueError("Raw FASTQ checksum differs: " + str(entry["path"]))


def validate_stage(directory, filename, kind, version, library, sources, config, inputs, sha,
                   reference_hash, raw_inputs, pilot_hashes=None):
    path = directory / filename
    if not path.is_file():
        raise ValueError("Existing output is incomplete: " + str(directory) + ". Preserve it and choose a clean output location manually.")
    summary = load_json(path)
    for key, value in {"analysis_type": kind, "analysis_version": version, "status": "complete",
                       "sample": library["sample"], "library_side": library["side"]}.items():
        if summary.get(key) != value:
            raise ValueError("Existing result " + key + " differs: " + str(path))
    if summary.get("source_sha256") != sources:
        raise ValueError("Existing result used different source code: " + str(path))
    for key, value in config.items():
        if summary.get("config", {}).get(key) != value:
            raise ValueError("Existing result configuration differs (" + key + "): " + str(path))
    validation = summary.get("validation", {})
    if not validation or not all(v is True for v in validation.values()):
        raise ValueError("Existing result validation failed: " + str(path))
    input_key = "upstream_input_sha256" if kind == PILOT_TYPE else "fullrun_input_sha256"
    raw_key = "upstream_inputs" if kind == PILOT_TYPE else "raw_inputs"
    if summary.get(input_key) != inputs or summary.get(raw_key) != raw_inputs:
        raise ValueError("Existing result input provenance differs: " + str(path))
    if summary.get("reference", {}).get("sha256") != reference_hash:
        raise ValueError("Existing result reference differs: " + str(path))
    if pilot_hashes is not None and summary.get("pilot_input_sha256") != pilot_hashes:
        raise ValueError("Existing reassessment pilot provenance differs: " + str(path))
    outputs = summary.get("output_sha256", {})
    required = {"junctions.tsv", "pair_junction_status.tsv.gz"} if kind == PILOT_TYPE else {
        "junction_reassessment.tsv", "identical_sequence_groups.tsv", "pair_refinement.tsv.gz", "read_realignment.tsv.gz"}
    if not required <= set(outputs):
        raise ValueError("Existing result lacks required output hashes: " + str(path))
    for name, expected in outputs.items():
        if Path(name).name != name or sha(directory / name) != expected:
            raise ValueError("Existing result output changed: " + str(directory / name))
    stats = summary["stats"]
    if int(stats["total_input_pairs"]) != int(library["manifest"]["total_pairs"]) or int(stats["host_supported_pairs"]) != int(library["manifest"]["host_supported_pairs"]):
        raise ValueError("Existing result denominators differ: " + str(path))
    if kind == REFINE_TYPE and sum(int(stats[k]) for k in STATUS_KEYS) != int(stats["host_supported_pairs"]):
        raise ValueError("Refinement states do not reconcile: " + str(path))
    return summary


def run_child(command, logfile, cwd):
    print("RUN " + " ".join(str(x) for x in command), flush=True)
    with logfile.open("w", encoding="utf-8") as log:
        log.write(json.dumps({"started_utc": now(), "command": command}) + "\n")
        log.flush()
        process = subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
        try:
            for line in process.stdout:
                log.write(line); log.flush()
                print(line, end="", flush=True)
            code = process.wait()
        except KeyboardInterrupt:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill(); process.wait()
            raise
        if code:
            raise RuntimeError(f"Analysis exited {code}; see {logfile}")


def write_batch(path, data):
    atomic_json(path / "batch_summary.json", data)
    fields = ["sample", "library_side", "status", "pilot_action", "refine_action", "total_input_pairs",
        "host_supported_pairs", *STATUS_KEYS, "supported_identical_sequence_groups", "refine_output", "error"]
    temporary = path / "batch_status.tsv.tmp"
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(data["libraries"])
    temporary.replace(path / "batch_status.tsv")


def merge_results(run_dir, records, sha):
    """Preserve per-library rows; append only explicit identifiers/MH NA labels.

    Empty engine tables can have abbreviated headers. Take the ordered union
    of schemas so a library with no candidates can still join the batch.
    """
    products = {}
    metadata = ["sample", "library_side", "condition", "biological_replicate"]
    for name in ("junction_reassessment.tsv", "identical_sequence_groups.tsv", "template_comparisons.tsv"):
        destination = run_dir / ("all_" + name + ".gz")
        fields = list(metadata)
        inputs = []
        populated_schema = None
        for record in records:
            directory = Path(record["refine_output"])
            summary = load_json(directory / "refinement_summary.json")
            source = directory / name
            if sha(source) != summary["output_sha256"][name]:
                raise ValueError("Source changed before merge: " + str(source))
            with source.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                names = reader.fieldnames or []
                if next(reader, None) is not None:
                    if populated_schema is None:
                        populated_schema = names
                    elif names != populated_schema:
                        raise ValueError("Different nonempty table schemas: " + name)
                fields.extend(key for key in names if key not in fields)
            inputs.append((source, summary))
        if name == "junction_reassessment.tsv":
            fields += ["mh_for_analysis_bp", "mh_assessment"]
        row_count = 0
        with gzip.open(destination, "wt", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fields, delimiter="\t")
            writer.writeheader()
            for source, summary in inputs:
                with source.open(encoding="utf-8", newline="") as handle:
                    for row in csv.DictReader(handle, delimiter="\t"):
                        if None in row:
                            raise ValueError("Malformed row: " + str(source))
                        for key in metadata:
                            value = summary.get(key, "")
                            if key in row and row[key] != str(value):
                                raise ValueError("Table/sample metadata mismatch: " + str(source))
                            row.setdefault(key, value)
                        if name == "junction_reassessment.tsv":
                            insertion = int(row["insertion_length_bp"]) > 0
                            row["mh_for_analysis_bp"] = "" if insertion else row["template_reference_mh_length_bp"]
                            row["mh_assessment"] = ("not_assessed_insertion_interfaces" if insertion
                                else "direct_junction_reference_annotation")
                        writer.writerow(row)
                        row_count += 1
        products[destination.name] = {"rows": row_count, "sha256": sha(destination)}
    return products


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--samples", nargs="+", help="Optional exact sample IDs")
    parser.add_argument("--expected-count", type=int, help="Stop if selected library count differs (e.g. 48)")
    parser.add_argument("--prefer-fullrun", choices=("v1", "v1_1"), help="Explicit choice only where duplicate versions exist")
    parser.add_argument("--pilot-script", type=Path)
    parser.add_argument("--refine-script", type=Path)
    parser.add_argument("--reference", type=Path, help="Same existing combined uncompressed FASTA if its path moved")
    parser.add_argument("--dry-run", action="store_true", help="Show inventory and planned actions; do not run or write outputs")
    parser.add_argument("--continue-on-error", action="store_true", help="Record failure and try remaining samples; exit nonzero if any failed")
    return parser.parse_args(argv)


def run(args):
    root = args.root.expanduser().resolve()
    if args.samples:
        for sample in args.samples:
            safe_sample(sample)
        if len(set(args.samples)) != len(args.samples):
            raise ValueError("Repeated --samples ID")
    libraries = discover(root, set(args.samples) if args.samples else None, args.prefer_fullrun)
    if args.expected_count is not None and (args.expected_count < 1 or len(libraries) != args.expected_count):
        raise ValueError(f"Expected {args.expected_count} libraries but found {len(libraries)}")
    sha = HashCache()
    pilot = find_engine(args.pilot_script, root, "raising_junction_pilot_v1", "raising_junction_pilot.py")
    refine = find_engine(args.refine_script, root, "raising_junction_refine_v1_1", "raising_junction_refine.py")
    engines = engine_metadata(pilot, refine, sha)
    print(f"Selected {len(libraries)} libraries: " + ", ".join(f"{side}={sum(x['side']==side for x in libraries)}" for side in ("3p", "5p")), flush=True)
    for library in libraries:
        sample = library["sample"]
        pdir = root / "results" / (sample + "_junction_pilot_v1")
        rdir = root / "results" / (sample + "_junction_refine_v1_1")
        print(f"{sample}\t{library['side']}\tpilot={'validate/reuse' if pdir.exists() else 'run'}\trefine={'validate/reuse' if rdir.exists() else 'run'}", flush=True)
    if args.dry_run:
        print("DRY RUN: inventory only. Existing outputs/raw data will be verified by SHA-256 during the actual run.", flush=True)
        return 0
    batch_root = root / "results" / "junction_batch_v1"
    batch_root.mkdir(parents=True, exist_ok=True)
    with (batch_root / ".running.lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError("Another junction batch is already running for this --root")
        run_id = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
        run_dir = batch_root / run_id
        run_dir.mkdir()
        data = {"batch_version": VERSION, "status": "running", "started_utc": now(), "root": str(root),
            "selected_library_count": len(libraries), "expected_count": args.expected_count,
            "scope": "selected completed full-depth libraries; not independent biological replicates",
            "pilot_script": str(pilot), "refine_script": str(refine), "engine_provenance": engines,
            "batch_source_sha256": sha(Path(__file__)),
            "libraries": [{"sample": x["sample"], "library_side": x["side"], "status": "pending",
                "fullrun_directory": str(x["directory"]), "excluded_fullrun_directories": x.get("excluded_fullrun_directories", [])} for x in libraries]}
        write_batch(run_dir, data)
        try:
            for number, (library, record) in enumerate(zip(libraries, data["libraries"]), 1):
                sample = library["sample"]
                print(f"\n[{number}/{len(libraries)}] {sample}", flush=True)
                record.update(status="running", started_utc=now())
                write_batch(run_dir, data)
                try:
                    inputs = evidence_hashes(library, sha)
                    validate_raw(library["manifest"], sha)
                    reference = (args.reference or Path(library["manifest"]["reference"]["combined_fasta"])).expanduser().resolve()
                    if reference.suffix == ".gz":
                        raise ValueError("Use the existing uncompressed combined FASTA")
                    reference_hash = sha(reference)
                    pdir = root / "results" / (sample + "_junction_pilot_v1")
                    rdir = root / "results" / (sample + "_junction_refine_v1_1")
                    common = ["--root", str(root), "--sample", sample, "--results-dir", str(library["directory"])]
                    if pdir.exists():
                        record["pilot_action"] = "reused_after_verification"
                    else:
                        record["pilot_action"] = "ran"
                        run_child([sys.executable, str(pilot), *common, "--reference", str(reference)], run_dir / (sample + ".pilot.log"), root)
                    ps = validate_stage(pdir, "summary.json", PILOT_TYPE, "0.1.0", library,
                        engines["pilot_sources"], engines["pilot_config"], inputs, sha, reference_hash, library["manifest"]["inputs"])
                    pilot_hashes = {name: sha(pdir / name) for name in ("summary.json", "junctions.tsv", "pair_junction_status.tsv.gz")}
                    if rdir.exists():
                        record["refine_action"] = "reused_after_verification"
                    else:
                        record["refine_action"] = "ran"
                        run_child([sys.executable, str(refine), *common, "--reference", str(reference), "--pilot-dir", str(pdir)], run_dir / (sample + ".refine.log"), root)
                    rs = validate_stage(rdir, "refinement_summary.json", REFINE_TYPE, "0.1.1", library,
                        engines["refine_sources"], engines["refine_config"], inputs, sha, reference_hash,
                        library["manifest"]["inputs"], pilot_hashes)
                    record.update(rs["stats"])
                    record.update(status="complete", completed_utc=now(), pilot_output=str(pdir), refine_output=str(rdir),
                        pilot_summary_sha256=sha(pdir / "summary.json"), refine_summary_sha256=sha(rdir / "refinement_summary.json"))
                except (OSError, ValueError, KeyError, RuntimeError) as exc:
                    record.update(status="failed", error=str(exc), completed_utc=now())
                    print("FAILED: " + str(exc), file=sys.stderr, flush=True)
                    write_batch(run_dir, data)
                    if not args.continue_on_error:
                        break
                write_batch(run_dir, data)
        except KeyboardInterrupt:
            data["status"] = "interrupted"
            for record in data["libraries"]:
                if record["status"] == "running":
                    record["status"] = "interrupted"
            data["completed_utc"] = now()
            write_batch(run_dir, data)
            print("Interrupted; completed verified samples can be reused on restart. Records: " + str(run_dir), file=sys.stderr)
            return 130
        complete = sum(x["status"] == "complete" for x in data["libraries"])
        if complete == len(libraries):
            try:
                data["merged_outputs"] = merge_results(run_dir, data["libraries"], sha)
            except (OSError, ValueError, KeyError) as exc:
                data.update(status="merge_failed", completed_utc=now(), completed_library_count=complete, error=str(exc))
                for path in run_dir.glob("all_*.tsv.gz"):
                    path.unlink()
                write_batch(run_dir, data)
                print("MERGE FAILED: " + str(exc), file=sys.stderr)
                return 2
        data.update(status="complete" if complete == len(libraries) else "incomplete", completed_utc=now(),
                    completed_library_count=complete)
        write_batch(run_dir, data)
        print(f"BATCH {data['status'].upper()}: {complete}/{len(libraries)} libraries. Records: {run_dir}", flush=True)
        return 0 if data["status"] == "complete" else 2


def main(argv=None):
    try:
        return run(parse_args(argv))
    except (OSError, ValueError, KeyError) as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
