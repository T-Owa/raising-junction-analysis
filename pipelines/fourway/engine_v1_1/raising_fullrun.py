#!/usr/bin/env python3
"""Full-depth paired RAISING classification with validation and resumable chunks.

Version 1.1.0 supports the supplied 800 bp donor and 3-prime or 5-prime libraries
named *_3prm_* or *_5prm_*. It preserves the bundled screening thresholds.
Run in the existing Ubuntu/WSL ``raising`` environment. Raw inputs are read-only.
"""
import argparse
import contextlib
import datetime
import fcntl
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import sys
import time
import traceback
import uuid

import raising_threeway_classify as core

ANALYSIS_VERSION = "1.1.0"
SOURCE_FILES = ("raising_fullrun.py", "raising_fullrun_report.py",
                "raising_threeway_classify.py", "raising_concatemer_screen.py")
REQUIRED_CHUNK_FILES = ("pair_classification.tsv", "read_evidence.tsv",
                        "donor_junction_evidence.tsv", "host_alignment_evidence.tsv",
                        "chunk_metrics.json")


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def atomic_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    tmp.replace(path)


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=True).encode("ascii")).hexdigest()


def file_metadata(path):
    path = Path(path).resolve()
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def update_stream_hash(digest, rid, raw1, raw2):
    # Length framing avoids ambiguous concatenations; headers are represented
    # by the paired normalized ID. Compressed-file hashes additionally preserve
    # every original byte, including comments and gzip metadata.
    for value in (rid, raw1[1], raw1[3], raw2[1], raw2[3]):
        encoded = value.encode("ascii")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)


def discover_fastq_pairs(root, sample, r1=None, r2=None):
    root = Path(root).expanduser().resolve()
    if bool(r1) != bool(r2):
        raise ValueError("Supply both --r1 and --r2, with matching files in the same order")
    if r1:
        if len(r1) != len(r2):
            raise ValueError("--r1 and --r2 must contain the same number of files")
        pairs = [(Path(a).expanduser().resolve(), Path(b).expanduser().resolve())
                 for a, b in zip(r1, r2)]
    else:
        rawdir = root / "vendor_results" / "Raw_fastq"
        firsts = sorted(rawdir.glob(sample + "_*_R1_*.fastq.gz"))
        pairs = [(p.resolve(), p.with_name(p.name.replace("_R1_", "_R2_", 1)).resolve())
                 for p in firsts]
        found_seconds = {p.resolve() for p in rawdir.glob(sample + "_*_R2_*.fastq.gz")}
        expected_seconds = {b for _, b in pairs}
        if found_seconds != expected_seconds:
            raise ValueError("R1/R2 lane files are incomplete or mismatched; use --r1/--r2 explicitly")
    if not pairs:
        raise ValueError("No full raw FASTQ pairs found in vendor_results/Raw_fastq")
    all_paths = [p for pair in pairs for p in pair]
    if len(set(all_paths)) != len(all_paths):
        raise ValueError("A FASTQ path was supplied more than once or used for both mates")
    for p in all_paths:
        if not p.is_file():
            raise ValueError("FASTQ file not found: " + str(p))
        if p.stat().st_size == 0:
            raise ValueError("Empty FASTQ file: " + str(p))
    return pairs


def validate_record(raw, path, record_number):
    _, sequence, _, quality = raw
    if set(sequence) - set("ACGTN"):
        raise ValueError(f"Non-ACGTN sequence in {path}, record {record_number}")
    if any(ord(c) < 33 or ord(c) > 126 for c in quality):
        raise ValueError(f"Invalid Phred+33 quality character in {path}, record {record_number}")


def iter_paired_fastq(pairs):
    """Yield (normalized_id, R1 four-tuple, R2 four-tuple); never drop mates."""
    for first, second in pairs:
        for number, (raw1, raw2) in enumerate(itertools.zip_longest(
                core.donor_screen.read_fastq(Path(first)),
                core.donor_screen.read_fastq(Path(second))), 1):
            if raw1 is None or raw2 is None:
                raise ValueError(f"Different R1/R2 record counts at record {number}: {first}; {second}")
            validate_record(raw1, first, number)
            validate_record(raw2, second, number)
            try:
                rid1 = core.donor_screen.normalized_id(raw1[0])
                rid2 = core.donor_screen.normalized_id(raw2[0])
            except IndexError as exc:
                raise ValueError(f"Empty FASTQ ID at record {number}: {first}; {second}") from exc
            if rid1 != rid2:
                raise ValueError(f"R1/R2 ID mismatch at record {number}: {rid1!r} != {rid2!r}; "
                                 "use raw paired files in their original matching order")
            if not rid1:
                raise ValueError("Empty normalized FASTQ ID")
            yield rid1, raw1, raw2


def scan_inputs(pairs, sqlite_path):
    """Validate every pair, disk-check global IDs, and hash all raw inputs."""
    sqlite_path = Path(sqlite_path)
    if sqlite_path.exists():
        raise ValueError("Validation database already exists: " + str(sqlite_path))
    sqlite_path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(sqlite_path)
    total, lanes = 0, []
    digest = hashlib.sha256()
    try:
        connection.execute("CREATE TABLE read_ids (read_id TEXT PRIMARY KEY) WITHOUT ROWID")
        for lane, pair in enumerate(pairs, 1):
            before = [file_metadata(p) for p in pair]
            count = 0
            with connection:
                for rid, raw1, raw2 in iter_paired_fastq([pair]):
                    try:
                        connection.execute("INSERT INTO read_ids VALUES (?)", (rid,))
                    except sqlite3.IntegrityError as exc:
                        raise ValueError("Duplicate normalized read ID across the inputs: " + rid) from exc
                    update_stream_hash(digest, rid, raw1, raw2)
                    count += 1
                    total += 1
                    if total % 100000 == 0:
                        print(f"Preflight: validated {total:,} paired records ...", flush=True)
            if not count:
                raise ValueError("FASTQ lane contains no records: " + str(pair[0]))
            data = {"lane_order": lane, "paired_records": count}
            for mate, path, original in zip(("R1", "R2"), pair, before):
                original["sha256"] = sha256_file(path)
                if file_metadata(path) != {k: original[k] for k in ("path", "size", "mtime_ns")}:
                    raise ValueError("Input changed during preflight: " + str(path))
                data[mate] = original
            lanes.append(data)
    finally:
        connection.close()
    if not total:
        raise ValueError("No input pairs")
    return {"total_pairs": total, "lanes": lanes,
            "normalized_paired_stream_sha256": digest.hexdigest()}


def inspect_reference(host, donor, root):
    """Check the original core's cache key without modifying its index."""
    host, root = Path(host).resolve(), Path(root).resolve()
    signature = {"host": core.file_metadata(host),
                 "donor_sha256": hashlib.sha256(donor.encode()).hexdigest(),
                 "donor_contig": core.DONOR_CONTIG, "format_version": 1}
    digest = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:20]
    directory = root / "ref" / ("threeway_" + digest)
    fasta = directory / "human_plus_donor.fa"
    marker = directory / "index_complete.json"
    result = {"complete": False, "fasta": str(fasta), "signature": signature,
              "directory": str(directory), "metadata": None,
              "reason": "Completed BWA index was not found for this exact host file"}
    indexes = [Path(str(fasta) + suffix) for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")]
    if marker.is_file() and fasta.is_file() and all(p.is_file() and p.stat().st_size for p in indexes):
        try:
            data = json.loads(marker.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            result["reason"] = "Index completion marker is unreadable"
        else:
            if data.get("signature") == signature:
                result.update(complete=True, metadata=data, reason="Completed existing BWA index will be reused")
            else:
                result["reason"] = "Index marker signature does not match the requested reference"
    return result


def reference_index_metadata(fasta):
    """Record large reference/index file identity without rehashing the index."""
    fasta = Path(fasta)
    paths = [fasta] + [Path(str(fasta) + suffix)
                      for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")]
    return [file_metadata(path) for path in paths if path.is_file()]


@contextlib.contextmanager
def output_lock(outdir):
    path = Path(outdir) / ".raising_fullrun.lock"
    with path.open("a+", encoding="utf-8") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("Another fullrun process is using this output directory") from exc
        handle.seek(0)
        handle.truncate()
        handle.write(json.dumps({"pid": os.getpid(), "locked_at": utc_now()}) + "\n")
        handle.flush()
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def source_hashes():
    directory = Path(__file__).resolve().parent
    return {name: sha256_file(directory / name) for name in SOURCE_FILES}


def check_no_input_output_overlap(outdir, input_paths):
    outdir = Path(outdir).resolve()
    for value in input_paths:
        path = Path(value).resolve()
        if path == outdir or outdir in path.parents:
            raise ValueError("An input is inside the output directory: " + str(path))


def validate_completed_chunk(chunk_dir, signature):
    chunk_dir = Path(chunk_dir)
    marker = chunk_dir / "complete.json"
    if not marker.is_file():
        raise ValueError("Chunk has no completion marker: " + str(chunk_dir))
    data = json.loads(marker.read_text(encoding="utf-8"))
    if data.get("status") != "complete" or data.get("signature") != signature:
        raise ValueError("Chunk signature/status mismatch: " + str(chunk_dir))
    hashes = data.get("output_sha256", {})
    if not set(REQUIRED_CHUNK_FILES).issubset(hashes):
        raise ValueError("Chunk manifest lacks required output hashes: " + str(chunk_dir))
    for relative, expected in hashes.items():
        target = chunk_dir / relative
        if target.is_symlink() or chunk_dir.resolve() not in target.resolve().parents:
            raise ValueError("Unsafe path in chunk manifest")
        if not target.is_file() or sha256_file(target) != expected:
            raise ValueError("Chunk output is missing or changed: " + str(target))
    return data


def make_chunk_records(raw_pairs, side):
    records, raw_r1 = {"R1": {}, "R2": {}}, {}
    for rid, raw1, raw2 in raw_pairs:
        records["R1"][rid] = core.clean_record(raw1, "R1", side)
        records["R2"][rid] = core.clean_record(raw2, "R2", side)
        raw_r1[rid] = raw1
    return records, raw_r1


def process_chunk(raw_pairs, index, chunk_root, run_signature, donor, sample,
                  side, threads, combined, resume):
    import raising_fullrun_report as report
    digest = hashlib.sha256()
    for rid, raw1, raw2 in raw_pairs:
        update_stream_hash(digest, rid, raw1, raw2)
    signature = {"run_signature": run_signature, "chunk_index": index,
                 "input_pairs": len(raw_pairs), "stream_sha256": digest.hexdigest()}
    final = chunk_root / f"chunk_{index:06d}"
    partial = chunk_root / f".chunk_{index:06d}.partial"
    if final.exists():
        if not resume:
            raise ValueError("Chunk already exists; use --resume: " + str(final))
        manifest = validate_completed_chunk(final, signature)
        return final, manifest, True
    if partial.exists():
        if not resume or partial.is_symlink() or partial.parent.resolve() != chunk_root.resolve():
            raise ValueError("Incomplete chunk exists; use --resume: " + str(partial))
        shutil.rmtree(partial)  # Only our own unfinished transaction is removed.
    partial.mkdir()
    work = partial / "work"
    work.mkdir()
    started = time.perf_counter()
    atomic_json(partial / "in_progress.json", {"signature": signature, "started_utc": utc_now()})
    records, raw_r1 = make_chunk_records(raw_pairs, side)
    fastqs = {mate: partial / f"cleaned_{mate}.fastq.gz" for mate in records}
    for mate, group in records.items():
        core.write_cleaned_fastq(fastqs[mate], group)
    commands = []
    core.align_donor(records, donor, sample, threads, work, commands)
    core.map_host(records, fastqs, combined, threads, commands, work)
    metrics = report.write_chunk_outputs(records, raw_r1, donor, sample, side, partial)
    if metrics.get("input_pairs") != len(raw_pairs) or sum(metrics.get("class_counts", {}).values()) != len(raw_pairs):
        raise ValueError("Chunk report counts do not reconcile with its paired inputs")
    if any(not (partial / name).is_file() for name in REQUIRED_CHUNK_FILES):
        raise ValueError("Chunk report did not create all required evidence files")
    elapsed = time.perf_counter() - started
    hashes = {}
    for p in sorted(partial.rglob("*")):
        if p.is_symlink():
            raise ValueError("Unexpected symlink in chunk output")
        if p.is_file():
            hashes[str(p.relative_to(partial))] = sha256_file(p)
    manifest = {"status": "complete", "signature": signature, "completed_utc": utc_now(),
                "input_pairs": len(raw_pairs), "metrics": metrics,
                "compute_seconds": elapsed, "output_sha256": hashes}
    atomic_json(partial / "complete.json", manifest)
    partial.rename(final)
    return final, manifest, False


def validate_input_metadata(inputs):
    for lane in inputs["lanes"]:
        for mate in ("R1", "R2"):
            old = lane[mate]
            if file_metadata(old["path"]) != {k: old[k] for k in ("path", "size", "mtime_ns")}:
                raise ValueError("Input changed during analysis: " + old["path"])


def human_time(seconds):
    seconds = max(0, int(seconds))
    hours, rest = divmod(seconds, 3600)
    minutes, seconds = divmod(rest, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path.home() / "raising")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--r1", nargs="+", type=Path)
    parser.add_argument("--r2", nargs="+", type=Path)
    parser.add_argument("--host-fasta", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--chunk-pairs", type=int, default=10000)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--preflight", action="store_true", help="Validate all inputs and index; do not classify")
    parser.add_argument("--allow-index-build", action="store_true",
                        help="Allow a new index only if the exact existing reference cache is unavailable")
    parser.add_argument("--condition", default="")
    parser.add_argument("--biological-replicate", default="")
    return parser


def sample_side(sample):
    """Infer one unambiguous direction, keeping the legacy core inference safe."""
    if (not re.fullmatch(r"[A-Za-z0-9._-]+", sample)
            or sample in (".", "..")):
        raise ValueError("Sample names may contain only letters, numbers, dots, underscores and hyphens")
    markers = re.findall(r"_(3|5)prm(?=_|$)", sample)
    if (len(markers) != 1 or ("_3prm" in sample and "_5prm" in sample)
            or re.search(r"_[35]prm(?!_|$)", sample)):
        raise ValueError("Version 1.1.0 requires exactly one _3prm or _5prm direction token in the sample name")
    return markers[0] + "p"


def run_analysis(args):
    overall_start = time.perf_counter()
    side = sample_side(args.sample)
    if args.threads < 1 or args.chunk_pairs < 1:
        raise ValueError("--threads and --chunk-pairs must be positive")
    if args.resume and args.preflight:
        raise ValueError("Use --preflight separately from --resume")
    root = args.root.expanduser().resolve()
    host = args.host_fasta.expanduser().resolve()
    if not host.is_file():
        raise ValueError("Host FASTA not found: " + str(host))
    host_initial_metadata = file_metadata(host)
    pairs = discover_fastq_pairs(root, args.sample, args.r1, args.r2)
    outdir = (args.outdir or root / "results" / (args.sample + "_fourway_full_v1_1")).expanduser().resolve()
    check_no_input_output_overlap(outdir, [p for pair in pairs for p in pair] + [host])
    bundled_donor = Path(__file__).resolve().parent / "RI_donor_800bp.fasta"
    donor = core.donor_screen.load_donor(bundled_donor)
    local_donor = root / "ref" / "RI_donor_800bp.fasta"
    if local_donor.exists() and core.donor_screen.load_donor(local_donor) != donor:
        raise ValueError("The root/ref donor differs from the bundled validated 800 bp donor")
    for binary in ("bwa", "blastn", "makeblastdb"):
        if shutil.which(binary) is None:
            raise ValueError(binary + " not found; run micromamba activate raising")
    ref = inspect_reference(host, donor, root)
    if not ref["complete"] and not args.allow_index_build:
        raise ValueError(ref["reason"] + ". Expected: " + ref["fasta"] +
                         ". Use the SAME --root and --host-fasta as your completed threeway run. "
                         "Only use --allow-index-build if you deliberately want a new index.")
    source = source_hashes()
    outdir.mkdir(parents=True, exist_ok=True)
    with output_lock(outdir):
        manifest_path = outdir / "run_manifest.json"
        previous = None
        if manifest_path.exists():
            if not args.resume and not args.preflight:
                raise ValueError("An analysis already exists here. Use --resume to validate/reuse it, or a new --outdir")
            previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        elif args.resume:
            raise ValueError("No run manifest exists; start this output directory without --resume")
        elif not args.preflight:
            unexpected = {p.name for p in outdir.iterdir()} - {
                ".raising_fullrun.lock", "preflight_summary.json", "preflight_runs"}
            if unexpected:
                raise ValueError("Output directory contains unrelated/incomplete files: " + ", ".join(sorted(unexpected)))
        run_id = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
        preflight_dir = outdir / "preflight_runs" / run_id
        preflight_dir.mkdir(parents=True)
        commands = []
        versions = core.tool_versions(commands, preflight_dir)
        print("Full input validation and SHA-256 calculation ...", flush=True)
        inputs = scan_inputs(pairs, preflight_dir / "input_ids.sqlite")
        print("Recording host FASTA SHA-256 (existing BWA index is unchanged) ...", flush=True)
        host_sha256 = sha256_file(host)
        if file_metadata(host) != host_initial_metadata:
            raise ValueError("Host FASTA changed during preflight: " + str(host))
        ref["host_fasta_sha256"] = host_sha256
        ref["index_file_metadata"] = reference_index_metadata(ref["fasta"])
        preflight_elapsed = time.perf_counter() - overall_start
        preflight = {"status": "preflight_complete", "sample": args.sample, "side": side,
                     "analysis_version": ANALYSIS_VERSION, "total_pairs": inputs["total_pairs"],
                     "inputs": inputs, "reference": ref, "source_sha256": source,
                     "tool_versions": versions, "elapsed_seconds": preflight_elapsed,
                     "completed_utc": utc_now()}
        atomic_json(preflight_dir / "preflight.json", preflight)
        atomic_json(outdir / "preflight_summary.json", preflight)
        print(f"Validated {inputs['total_pairs']:,} pairs across {len(pairs)} lane file pair(s). "
              f"Preflight elapsed {human_time(preflight_elapsed)}.", flush=True)
        print(ref["reason"] + ": " + ref["fasta"], flush=True)
        if args.preflight:
            print("Preflight complete; classification was not run.", flush=True)
            return preflight
        config = {"sample": args.sample, "side": side, "threads": args.threads,
                  "chunk_pairs": args.chunk_pairs, "condition": args.condition,
                  "biological_replicate": args.biological_replicate,
                  "analysis_version": ANALYSIS_VERSION,
                  "donor_anchor_independent_matches": 30, "donor_min_identity_percent": 95,
                  "donor_max_gap_bp": 50, "donor_max_overlap_bp": 20, "donor_min_rewind_bp": 20,
                  "host_independent_matches": 30, "host_min_identity_percent": 95,
                  "host_min_mapq": 20, "host_min_base_quality": 20,
                  "donor_base_quality_filter": "none_in_preserved_core",
                  "short_pattern": f"exact_{side}_donor_prefix_then_A15; adapter_confirmed_flag_separate"}
        signature_data = {"config": config, "inputs": inputs, "source_sha256": source,
                          "tool_versions": versions, "reference_signature": ref["signature"],
                          "host_fasta_sha256": host_sha256}
        signature = canonical_hash(signature_data)
        if previous is not None and previous.get("run_signature") != signature:
            raise ValueError("Resume conditions changed (inputs, code, tools, reference, or arguments). "
                             "Keep the original conditions or choose a new --outdir")
        if previous is not None and (outdir / "run_summary.json").is_file():
            old_summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            if old_summary.get("status") == "complete":
                print("--resume requested for a completed run: validating every chunk and "
                      "regenerating the final count tables from those verified chunks.", flush=True)
        manifest = {"run_signature": signature, "signature_data": signature_data,
                    "created_utc": previous.get("created_utc") if previous else utc_now(),
                    "initial_run_id": previous.get("initial_run_id") if previous else run_id,
                    "input_pairs": inputs["total_pairs"], "analysis_version": ANALYSIS_VERSION}
        if previous is None:
            atomic_json(manifest_path, manifest)
        state = {"status": "running", "classification_available": False, "sample": args.sample,
                 "run_signature": signature, "run_id": run_id, "analysis_version": ANALYSIS_VERSION,
                 "total_pairs": inputs["total_pairs"], "completed_pairs": 0,
                 "note": "Use final count tables only when this summary reports status complete."}
        atomic_json(outdir / "run_summary.json", state)
        try:
            index_start = time.perf_counter()
            if ref["complete"]:
                combined, reference_metadata = Path(ref["fasta"]), ref["metadata"]
            else:
                combined, reference_metadata = core.prepare_reference(host, donor, root, commands, preflight_dir)
            reference_metadata = {**(reference_metadata or {}), "host_fasta_sha256": host_sha256,
                                  "index_file_metadata": reference_index_metadata(combined)}
            index_elapsed = time.perf_counter() - index_start
            chunk_root = outdir / "chunks"
            chunk_root.mkdir(exist_ok=True)
            expected_chunks = math.ceil(inputs["total_pairs"] / args.chunk_pairs)
            existing_names = {p.name for p in chunk_root.iterdir() if p.name.startswith("chunk_")}
            expected_names = {f"chunk_{i:06d}" for i in range(1, expected_chunks + 1)}
            if existing_names - expected_names:
                raise ValueError("Unexpected extra completed chunks in output directory")
            processed, new_pairs, new_seconds, reused = 0, 0, 0.0, 0
            chunk_dirs, chunk_details, buffer = [], [], []
            all_digest = hashlib.sha256()
            process_start = time.perf_counter()

            def finish_buffer(rows):
                nonlocal processed, new_pairs, new_seconds, reused
                number = len(chunk_dirs) + 1
                start = time.perf_counter()
                directory, completed, did_reuse = process_chunk(
                    rows, number, chunk_root, signature, donor, args.sample, side,
                    args.threads, combined, args.resume)
                wall = time.perf_counter() - start
                chunk_dirs.append(directory)
                processed += len(rows)
                reused += int(did_reuse)
                if not did_reuse:
                    new_pairs += len(rows)
                    new_seconds += wall
                detail = {"chunk": number, "pairs": len(rows), "reused": did_reuse,
                          "compute_seconds": completed["compute_seconds"], "current_wall_seconds": wall}
                chunk_details.append(detail)
                eta = ((inputs["total_pairs"] - processed) * new_seconds / new_pairs) if new_pairs else None
                state.update(completed_pairs=processed, completed_chunks=number,
                             elapsed_seconds=time.perf_counter() - overall_start,
                             eta_seconds=eta, last_chunk=detail)
                atomic_json(outdir / "run_summary.json", state)
                atomic_json(outdir / "chunk_timings.json", chunk_details)
                print(f"Chunk {number}/{expected_chunks}: {'reused' if did_reuse else 'completed'}; "
                      f"{processed:,}/{inputs['total_pairs']:,} pairs; "
                      f"elapsed {human_time(state['elapsed_seconds'])}; "
                      f"ETA {human_time(eta) if eta is not None else 'pending first newly computed chunk'}", flush=True)

            for rid, raw1, raw2 in iter_paired_fastq(pairs):
                update_stream_hash(all_digest, rid, raw1, raw2)
                buffer.append((rid, raw1, raw2))
                if len(buffer) == args.chunk_pairs:
                    finish_buffer(buffer)
                    buffer = []
            if buffer:
                finish_buffer(buffer)
            if processed != inputs["total_pairs"] or all_digest.hexdigest() != inputs["normalized_paired_stream_sha256"]:
                raise ValueError("Paired input stream changed after preflight")
            validate_input_metadata(inputs)
            if file_metadata(host) != host_initial_metadata:
                raise ValueError("Host FASTA changed during analysis: " + str(host))
            if source_hashes() != source:
                raise ValueError("Analysis source files changed during execution; start a new output directory")
            processing_elapsed = time.perf_counter() - process_start
            import raising_fullrun_report as report
            aggregate_start = time.perf_counter()
            metadata = {"sample": args.sample, "side": side, "total_pairs": processed,
                        "condition": args.condition, "biological_replicate": args.biological_replicate,
                        "analysis_version": ANALYSIS_VERSION, "config": config, "inputs": inputs,
                        "tool_versions": versions, "source_sha256": source,
                        "run_signature": signature, "reference": reference_metadata,
                        "timing_seconds": {"preflight": preflight_elapsed, "index": index_elapsed,
                            "processing_current_invocation": processing_elapsed,
                            "chunk_compute_total": sum(c["compute_seconds"] for c in chunk_details)}}
            aggregate = report.aggregate_outputs(chunk_dirs, outdir, metadata)
            if source_hashes() != source:
                raise ValueError("Analysis source files changed during final aggregation")
            total = aggregate.get("total_pairs", aggregate.get("input_pairs"))
            if total != processed or sum(aggregate.get("class_counts", {}).values()) != processed:
                raise ValueError("Final count tables do not reconcile with all validated input pairs")
            metadata["timing_seconds"].update(aggregation=time.perf_counter() - aggregate_start,
                                             total_wall_current_invocation=time.perf_counter() - overall_start)
            final_hashes = {p.name: sha256_file(p) for p in sorted(outdir.iterdir())
                            if p.is_file() and (p.suffix == ".csv" or p.name.endswith(".tsv.gz") or p.name == "methods.txt")}
            final = {**metadata, **aggregate, "status": "complete", "classification_available": True,
                     "completed_utc": utc_now(), "completed_chunks": len(chunk_dirs),
                     "reused_chunks": reused, "final_output_sha256": final_hashes}
            atomic_json(outdir / "run_summary.json", final)
            print(f"COMPLETE: {processed:,} pairs classified. Results: {outdir}", flush=True)
            print("Main Figure table: " + str(outdir / "figure_counts.csv"), flush=True)
            return final
        except BaseException as exc:
            state.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
                         classification_available=False, error_type=type(exc).__name__,
                         error=str(exc), elapsed_seconds=time.perf_counter() - overall_start)
            atomic_json(outdir / "run_summary.json", state)
            with (preflight_dir / "exception.txt").open("w", encoding="utf-8") as handle:
                traceback.print_exc(file=handle)
            raise


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        run_analysis(args)
        return 0
    except KeyboardInterrupt:
        print("Interrupted. Re-run the same command with --resume to reuse completed chunks.", file=sys.stderr)
        return 130
    except (OSError, ValueError, EOFError, UnicodeError, sqlite3.Error) as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
