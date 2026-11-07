#!/usr/bin/env python3
"""Conservative exact-anchor pilot for existing full-depth RAISING evidence.

This is a post-processing caller, not a reimplementation of the four-way screen.
Only same-read donor--human junctions satisfying junction_core are called. All
host-supported pairs remain in the resolution denominator. Raw data and the
four-way result directory are never modified.
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

import junction_core as core

VERSION = "0.1.0"
ANALYSIS_TYPE = "raising_exact_anchor_junction_pilot"
DONOR_CONTIG = "__RAISING_DONOR_800__"
STATUSES = ("called_simple_candidate", "called_prior_donor_rearrangement",
            "called_complex_host_evidence", "ambiguous", "unresolved")
PAIR_FIELDS = ["sample", "library_side", "read_id", "status", "event_id",
    "prior_classification", "prior_donor_rearrangement_supported", "supporting_mates",
    "distinct_candidate_events", "complex_host_evidence", "diagnostics_json"]
EVENT_FIELDS = ["sample", "library_side", "event_id", "event_key_json", "chrom",
    "host_strand", "host_boundary0", "host_boundary_min0", "host_boundary_max0",
    "donor_cut_oriented_min0", "donor_cut_oriented_max0",
    "donor_cut_original_min0", "donor_cut_original_max0",
    "donor_deletion_min_bp", "donor_deletion_max_bp", "mh_length_bp", "mh_sequence",
    "insertion_length_bp", "insertion_sequence", "pair_support", "simple_pair_support",
    "prior_dd_pair_support", "complex_host_pair_support", "example_read_id",
    "example_mate", "example_candidate_json"]


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def truth(value):
    return str(value).strip().lower() in ("1", "true", "yes")


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_key(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def text_open(path, mode="rt"):
    return gzip.open(path, mode, encoding="ascii", newline="") if str(path).endswith(".gz") else open(path, mode, encoding="ascii", newline="")


def tsv_rows(path, required):
    with text_open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        for row in reader:
            if None in row:
                raise ValueError(f"Malformed TSV row in {path}")
            yield row


def evidence_path(directory, stem):
    found = [directory / (stem + suffix) for suffix in (".tsv.gz", ".tsv")
             if (directory / (stem + suffix)).is_file()]
    if len(found) != 1:
        raise ValueError(f"Expected exactly one {stem}.tsv[.gz] in {directory}")
    return found[0]


def discover_results(root, sample, explicit=None):
    if explicit:
        directory = Path(explicit).expanduser().resolve()
    else:
        found = [root / "results" / (sample + suffix)
                 for suffix in ("_fourway_full_v1", "_fourway_full_v1_1")
                 if (root / "results" / (sample + suffix) / "run_summary.json").is_file()]
        if len(found) != 1:
            raise ValueError("Expected one full-run result directory; use --results-dir explicitly. Found: " + repr([str(p) for p in found]))
        directory = found[0]
    manifest_path = directory / "run_summary.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete" or manifest.get("sample") != sample:
        raise ValueError("Four-way run_summary must be complete and match --sample")
    if manifest.get("validation") and not all(manifest["validation"].values()):
        raise ValueError("Upstream full-run validation is not all true")
    if not manifest.get("inputs", {}).get("lanes"):
        raise ValueError("Full-run input manifest is absent; sampled three-way output is not supported")
    return directory, manifest_path, manifest


class FastaReader:
    """Read an uncompressed FASTA through its samtools-compatible .fai."""
    def __init__(self, path):
        self.path = Path(path)
        if self.path.suffix == ".gz":
            raise ValueError("Use the existing uncompressed human_plus_donor.fa, not the downloaded .fa.gz")
        fai = Path(str(self.path) + ".fai")
        if not fai.is_file() or fai.stat().st_mtime_ns < self.path.stat().st_mtime_ns:
            if not shutil.which("samtools"):
                raise ValueError("samtools is required only to make the missing/stale .fai index")
            print("Creating small FASTA .fai index with samtools faidx (no BWA indexing) ...", flush=True)
            subprocess.run(["samtools", "faidx", str(self.path)], check=True)
        self.entries = {}
        for line in fai.read_text(encoding="ascii").splitlines():
            fields = line.split("\t")
            if len(fields) < 5 or fields[0] in self.entries:
                raise ValueError("Invalid FASTA .fai entry")
            length, offset, bases, width = map(int, fields[1:5])
            if length < 0 or bases <= 0 or width < bases:
                raise ValueError("Invalid FASTA .fai dimensions")
            self.entries[fields[0]] = (length, offset, bases, width)
        self.handle = self.path.open("rb")

    def close(self):
        self.handle.close()

    def length(self, chrom):
        if chrom not in self.entries:
            raise ValueError("Contig absent from FASTA: " + chrom)
        return self.entries[chrom][0]

    def fetch(self, chrom, start0, end0):
        length, offset, bases, width = self.entries.get(chrom, (None, None, None, None))
        if length is None or not 0 <= start0 <= end0 <= length:
            raise ValueError(f"Invalid reference interval: {chrom}:{start0}-{end0}")
        if start0 == end0:
            return ""
        begin = offset + (start0 // bases) * width + start0 % bases
        last = end0 - 1
        finish = offset + (last // bases) * width + last % bases + 1
        self.handle.seek(begin)
        seq = self.handle.read(finish - begin).replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
        if len(seq) != end0 - start0 or ">" in seq:
            raise ValueError("FASTA/.fai mismatch; recreate the .fai")
        return seq


def cigar_intervals(position_1based, cigar):
    tokens = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not tokens or "".join(a + b for a, b in tokens) != cigar:
        raise ValueError("Invalid CIGAR: " + cigar)
    position = int(position_1based) - 1
    if position < 0:
        raise ValueError("Invalid 1-based host position")
    intervals = []
    for count, operation in tokens:
        count = int(count)
        if count <= 0:
            raise ValueError("Zero-length CIGAR operation")
        if operation in "M=X":
            intervals.append((position, position + count))
            position += count
        elif operation in "DN":
            position += count
    if not intervals:
        raise ValueError("Accepted host CIGAR has no aligned bases")
    return intervals


def normalized_id(header):
    fields = header.lstrip("@").split()
    if not fields:
        raise ValueError("Empty FASTQ identifier")
    rid = fields[0]
    return rid[:-2] if rid.endswith(("/1", "/2")) else rid


def fastq_records(path):
    with text_open(path) as handle:
        number = 0
        while True:
            header = handle.readline()
            if not header:
                return
            record = [header.rstrip("\r\n")] + [handle.readline().rstrip("\r\n") for _ in range(3)]
            number += 1
            if not record[0].startswith("@") or not record[2].startswith("+") or not record[1] or len(record[1]) != len(record[3]):
                raise ValueError(f"Invalid FASTQ record {number} in {path}")
            if set(record[1].upper()) - set("ACGTN") or any(not 33 <= ord(c) <= 126 for c in record[3]):
                raise ValueError(f"Invalid sequence/quality at record {number} in {path}")
            yield tuple(record)


def update_stream_hash(digest, rid, raw1, raw2):
    for value in (rid, raw1[1], raw1[3], raw2[1], raw2[3]):
        encoded = value.encode("ascii")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)


def cleaned_record(raw, evidence):
    sequence, quality = raw[1], raw[3]
    if len(sequence) != int(evidence["raw_length"]):
        raise ValueError("Raw read length differs from full-run evidence")
    length = int(evidence["cleaned_length"])
    if not length:
        return "", ""
    start = int(evidence["cleaned_raw_start_1based"]) - 1
    end = int(evidence["cleaned_raw_end_1based"])
    if not 0 <= start < end <= len(sequence) or end - start != length:
        raise ValueError("Invalid cleaned-read coordinates in full-run evidence")
    return sequence[start:end].upper(), quality[start:end]


def candidate_key(candidate):
    return json_key(candidate["event_key"])


def event_id(key):
    return "J" + hashlib.sha256(key.encode("ascii")).hexdigest()[:16]


def select_pair_result(calls, prior_dd, host_hits):
    """Collapse identical event keys; no best-score resolution of distinct joins."""
    by_event = collections.defaultdict(list)
    for call in calls:
        by_event[candidate_key(call["candidate"])].append(call)
    if not by_event:
        return "unresolved", None, [], False
    if len(by_event) != 1:
        return "ambiguous", None, [], False
    key, support = next(iter(by_event.items()))
    # More independent anchor bases is only an exemplar choice within one event.
    support.sort(key=lambda c: (-int(c["candidate"].get("host_independent_bases", 0)),
                                -int(c["candidate"].get("donor_independent_bases", 0)),
                                c["mate"], c["query_orientation"]))
    example = support[0]
    c = example["candidate"]
    chrom = c.get("chrom", c.get("host_chrom", example["chrom"]))
    boundary = int(c["host_boundary0"])
    complex_host = any(h["reference"] != chrom or
                       min(abs(boundary - a) if boundary < a else abs(boundary - b) if boundary > b else 0
                           for a, b in h["intervals0"]) > 5000
                       for h in host_hits)
    status = "called_complex_host_evidence" if complex_host else (
        "called_prior_donor_rearrangement" if prior_dd else "called_simple_candidate")
    return status, example, support, complex_host


def representative_rows(ordered):
    """High-support examples plus examples spanning important call strata."""
    chosen = ordered[:20]
    predicates = (
        lambda r: int(r["mh_length_bp"]) == 0 and int(r["insertion_length_bp"]) == 0,
        lambda r: int(r["mh_length_bp"]) > 0,
        lambda r: int(r["insertion_length_bp"]) > 0,
        lambda r: r["prior_dd_pair_support"] > 0,
        lambda r: r["complex_host_pair_support"] > 0,
        lambda r: r["pair_support"] == 1,
        lambda r: int(r["donor_deletion_min_bp"]) > 0,
    )
    identifiers = {r["event_id"] for r in chosen}
    for predicate in predicates:
        example = next((r for r in ordered if predicate(r)), None)
        if example is not None and example["event_id"] not in identifiers:
            chosen.append(example)
            identifiers.add(example["event_id"])
    for row in ordered:
        if len(chosen) >= 30:
            break
        if row["event_id"] not in identifiers:
            chosen.append(row)
            identifiers.add(row["event_id"])
    return chosen


def write_representatives(path, ordered, donor, side):
    """Show exact reference blocks on the normalized query's coordinate axis."""
    oriented_donor = donor if side == "3p" else core.reverse_complement(donor)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("Representative exact-anchor candidates: 20 highest-support patterns plus\n"
                     "examples of blunt/MH/insertion/prior-DD/complex/singleton calls, up to 30.\n"
                     "The read is oriented donor -> human. D and H show exact reference matches.\n"
                     "'^' = shared MH bases; '+' = intervening bases; '|' = blunt boundary.\n"
                     "All coordinates are 0-based half-open intervals or between-base boundaries.\n"
                     "These are review records, not independent molecule or insertion counts.\n\n")
        for row in representative_rows(ordered):
            example = row["_example"]
            c = example["candidate"]
            query = example["normalized_query"]
            donor_block, host_block = c["donor_block"], c["host_block"]
            d_start, d_end = donor_block["query_start0"], donor_block["query_end0"]
            h_start, h_end = host_block["query_start0"], host_block["query_end0"]
            d_track, h_track, join_track = [" "] * len(query), [" "] * len(query), [" "] * (len(query) + 1)
            d_track[d_start:d_end] = oriented_donor[donor_block["target_start0"]:donor_block["target_end0"]]
            h_track[h_start:h_end] = example["oriented_host_window"][host_block["target_start0"]:host_block["target_end0"]]
            if d_end > h_start:
                join_track[h_start:d_end] = "^" * (d_end - h_start)
            elif d_end < h_start:
                join_track[d_end:h_start] = "+" * (h_start - d_end)
            else:
                join_track[d_end] = "|"
            handle.write(f"{row['event_id']}  pairs={row['pair_support']}  simple={row['simple_pair_support']}  prior_DD={row['prior_dd_pair_support']}  complex={row['complex_host_pair_support']}\n")
            handle.write(f"Read {example['read_id']} {example['mate']} ({example['query_orientation']})\n")
            handle.write(f"MH={c['mh_length_bp']} bp ({c['mh_sequence'] or '-'})  intervening={c['insertion_length_bp']} bp ({c['insertion_sequence'] or '-'})  donor deletion={c['donor_deletion_min_bp']}..{c['donor_deletion_max_bp']} bp\n")
            handle.write("QUERY  " + query + "\n")
            handle.write("DONOR  " + "".join(d_track) + "\n")
            handle.write("HUMAN  " + "".join(h_track) + "\n")
            handle.write("JOIN   " + "".join(join_track).rstrip() + "\n")
            handle.write(f"Donor reference ({side} outward orientation): [{donor_block['target_start0']},{donor_block['target_end0']}); query [{d_start},{d_end})\n")
            handle.write(f"Human: {example['chrom']} strand={example['host_strand']} boundary0={c['host_boundary0']}; query [{h_start},{h_end})\n")
            handle.write(f"Independent exact Q20 bases: donor={c['donor_independent_bases']}, human={c['host_independent_bases']}; global host overlap={c['host_global_overlap_bases']}\n")
            handle.write("Flags: " + json_key(c.get("flags", [])) + "\n")
            handle.write("Equivalent boundaries: " + json_key(c.get("equivalent_boundaries", [])) + "\n\n")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--results-dir", type=Path)
    p.add_argument("--reference", type=Path, help="Existing uncompressed combined FASTA, if its path has moved")
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--side", choices=("auto", "3p", "5p"), default="auto")
    p.add_argument("--window", type=int, default=500)
    p.add_argument("--overwrite", action="store_true", help="Replace this pilot's own previously completed outputs")
    return p.parse_args(argv)


def run(args):
    started = time.monotonic()
    root = args.root.expanduser().resolve()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.sample):
        raise ValueError("Unsafe sample name")
    if args.window < 100:
        raise ValueError("--window must be >=100")
    directory, manifest_path, manifest = discover_results(root, args.sample, args.results_dir)
    side = manifest.get("side") if args.side == "auto" else args.side
    if side not in ("3p", "5p") or side != manifest.get("side"):
        raise ValueError("Library side must match the upstream full-run manifest")
    total_expected = int(manifest["total_pairs"])
    outdir = (args.output_dir or root / "results" / (args.sample + "_junction_pilot_v1")).expanduser().resolve()
    reference = Path(args.reference or manifest["reference"]["combined_fasta"]).expanduser().resolve()
    if not reference.is_file():
        raise ValueError("Combined FASTA absent; supply --reference with the same uncompressed human_plus_donor.fa")
    if outdir == directory or outdir in directory.parents or outdir == reference.parent or outdir in reference.parents:
        raise ValueError("Output cannot contain existing classification/reference inputs")
    if outdir.exists():
        previous = outdir / "summary.json"
        if not args.overwrite or not previous.is_file() or json.loads(previous.read_text()).get("analysis_type") != ANALYSIS_TYPE:
            raise ValueError("Output exists; choose a new --output-dir or use --overwrite for this pilot's previous output")
    paths = {stem: evidence_path(directory, stem) for stem in
             ("pair_classification", "read_evidence", "host_alignment_evidence")}
    hashes = {"run_summary.json": sha256_file(manifest_path)}
    print("Validating existing full-run evidence ...", flush=True)
    for stem, path in paths.items():
        hashes[path.name] = sha256_file(path)
        expected = manifest.get("final_output_sha256", {}).get(path.name)
        if expected and expected != hashes[path.name]:
            raise ValueError("Evidence SHA-256 differs from run_summary: " + str(path))
    selected, all_ids = {}, set()
    for row in tsv_rows(paths["pair_classification"], ("read_id", "classification", "host_mapping_supported")):
        rid = row["read_id"]
        if not rid or rid in all_ids:
            raise ValueError("Duplicate/empty pair ID in full-run classification")
        all_ids.add(rid)
        if truth(row["host_mapping_supported"]):
            selected[rid] = row
    if len(all_ids) != total_expected or len(selected) != int(manifest["host_supported_pairs"]):
        raise ValueError("Evidence pair counts differ from complete full-run summary")
    reads, hits = {}, collections.defaultdict(list)
    for row in tsv_rows(paths["read_evidence"], ("read_id", "mate", "raw_length", "cleaned_length", "cleaned_raw_start_1based", "cleaned_raw_end_1based")):
        if row["read_id"] in selected:
            key = (row["read_id"], row["mate"])
            if row["mate"] not in ("R1", "R2") or key in reads:
                raise ValueError("Duplicate or invalid read-evidence mate")
            reads[key] = row
    if set(reads) != {(rid, mate) for rid in selected for mate in ("R1", "R2")}:
        raise ValueError("Every selected pair must have both read-evidence rows")
    for row in tsv_rows(paths["host_alignment_evidence"], ("read_id", "mate", "reference", "position_1based", "cigar", "strand", "accepted")):
        if row["read_id"] in selected and truth(row["accepted"]):
            if row["mate"] not in ("R1", "R2") or row["strand"] not in ("+", "-") or row["reference"] == DONOR_CONTIG:
                raise ValueError("Invalid accepted host-evidence row")
            row["intervals0"] = cigar_intervals(row["position_1based"], row["cigar"])
            hits[(row["read_id"], row["mate"])].append(row)
    if any(not hits.get((rid, "R1")) and not hits.get((rid, "R2")) for rid in selected):
        raise ValueError("Host-supported pair lacks an accepted host-alignment row")
    print(f"Input pairs: {total_expected:,}; host-supported pairs (including both-supported): {len(selected):,}", flush=True)
    fasta = FastaReader(reference)
    donor_contig = manifest["reference"].get("signature", {}).get("donor_contig", DONOR_CONTIG)
    donor = fasta.fetch(donor_contig, 0, fasta.length(donor_contig))
    donor_sha = hashlib.sha256(donor.encode("ascii")).hexdigest()
    expected_donor = manifest["reference"].get("signature", {}).get("donor_sha256")
    if len(donor) != 800 or (expected_donor and donor_sha != expected_donor):
        raise ValueError("Reference donor must match the original 800-bp donor")
    # Source reference metadata is checked when available; a content digest is
    # recorded here independently for reproducibility of the local re-analysis.
    for item in manifest["reference"].get("index_file_metadata", []):
        if Path(item["path"]).resolve() == reference:
            if reference.stat().st_size != item["size"] or reference.stat().st_mtime_ns != item["mtime_ns"]:
                raise ValueError("Combined reference metadata changed since the full run")
    print("Recording reference SHA-256 ...", flush=True)
    reference_sha = sha256_file(reference)
    lanes = []
    for lane in manifest["inputs"]["lanes"]:
        pair = []
        for mate in ("R1", "R2"):
            path = Path(lane[mate]["path"]).expanduser().resolve()
            if outdir == path.parent or outdir in path.parents:
                raise ValueError("Output contains a raw input")
            if not path.is_file():
                raise ValueError("Original FASTQ is missing: " + str(path))
            current_sha = sha256_file(path)
            if lane[mate].get("sha256") and current_sha != lane[mate]["sha256"]:
                raise ValueError("Raw FASTQ SHA-256 differs from full-run manifest: " + str(path))
            pair.append(path)
        lanes.append(tuple(pair))
    outdir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="." + outdir.name + ".partial_", dir=outdir.parent))
    counts, events, seen = collections.Counter(), {}, set()
    digest = hashlib.sha256()
    attempted_reads = 0
    calls_path = staging / "read_candidate_evidence.tsv.gz"
    pair_path = staging / "pair_junction_status.tsv.gz"
    window_cache = {}
    with gzip.open(pair_path, "wt", encoding="utf-8", newline="") as ph, gzip.open(calls_path, "wt", encoding="utf-8", newline="") as ch:
        pw = csv.DictWriter(ph, PAIR_FIELDS, delimiter="\t"); pw.writeheader()
        cw = csv.DictWriter(ch, ["read_id", "mate", "query_orientation", "chrom", "host_strand", "core_status", "candidate_json", "diagnostics_json"], delimiter="\t"); cw.writeheader()
        for first, second in lanes:
            print("Reading paired FASTQ: " + first.name, flush=True)
            for raw1, raw2 in itertools.zip_longest(fastq_records(first), fastq_records(second)):
                if raw1 is None or raw2 is None:
                    raise ValueError("R1/R2 FASTQ counts differ")
                rid = normalized_id(raw1[0])
                if rid != normalized_id(raw2[0]) or rid in seen or rid not in all_ids:
                    raise ValueError("FASTQ pair ID mismatch, duplicate, or absence from full-run evidence: " + rid)
                seen.add(rid)
                update_stream_hash(digest, rid, raw1, raw2)
                if len(seen) % 100000 == 0:
                    print(f"Scanned {len(seen):,}/{total_expected:,} pairs; evaluated {sum(counts.values()):,} host-supported pairs ...", flush=True)
                if rid not in selected:
                    continue
                pair_calls, diagnostics = [], []
                for mate, raw in (("R1", raw1), ("R2", raw2)):
                    query, quality = cleaned_record(raw, reads[(rid, mate)])
                    if not query or not hits.get((rid, mate)):
                        continue
                    attempted_reads += 1
                    for hit in hits[(rid, mate)]:
                        chrom = hit["reference"]
                        intervals = hit["intervals0"]
                        start = max(0, min(a for a, b in intervals) - args.window)
                        end = min(fasta.length(chrom), max(b for a, b in intervals) + args.window)
                        cache_key = (chrom, start, end)
                        if cache_key not in window_cache:
                            if len(window_cache) >= 10000:
                                window_cache.clear()
                            window_cache[cache_key] = fasta.fetch(chrom, start, end)
                        human = window_cache[cache_key]
                        for flip in (False, True):
                            q = core.reverse_complement(query) if flip else query
                            qual = quality[::-1] if flip else quality
                            strand = ("-" if hit["strand"] == "+" else "+") if flip else hit["strand"]
                            host = core.reverse_complement(human) if strand == "-" else human
                            result = core.call_junctions(q, qual, donor, side, host, chrom, start, strand,
                                                         accepted_host_interval0=intervals)
                            diagnostics.append({"mate": mate, "query_orientation": "reverse_complement" if flip else "forward", "chrom": chrom, "status": result["status"], "diagnostics": result.get("diagnostics", {})})
                            base = {"read_id": rid, "mate": mate, "query_orientation": "reverse_complement" if flip else "forward", "chrom": chrom, "host_strand": strand,
                                    "core_status": result["status"], "diagnostics_json": json_key(result.get("diagnostics", {}))}
                            if not result.get("candidates"):
                                cw.writerow({**base, "candidate_json": ""})
                            for candidate in result.get("candidates", []):
                                cw.writerow({**base, "candidate_json": json_key(candidate)})
                                pair_calls.append({"candidate": candidate, "read_id": rid, "mate": mate,
                                    "query_orientation": base["query_orientation"], "chrom": chrom, "host_strand": strand,
                                    "normalized_query": q, "oriented_host_window": host, "host_window_start0": start})
                prior_dd = truth(selected[rid].get("donor_rearrangement_supported", selected[rid].get("both_supported", "0")))
                status, example, support, complex_host = select_pair_result(pair_calls, prior_dd, hits.get((rid, "R1"), []) + hits.get((rid, "R2"), []))
                counts[status] += 1
                key = candidate_key(example["candidate"]) if example else None
                identifier = event_id(key) if key else ""
                pw.writerow({"sample": args.sample, "library_side": side, "read_id": rid, "status": status,
                    "event_id": identifier, "prior_classification": selected[rid]["classification"],
                    "prior_donor_rearrangement_supported": int(prior_dd),
                    "supporting_mates": ";".join(sorted({c["mate"] for c in support})),
                    "distinct_candidate_events": len({candidate_key(c["candidate"]) for c in pair_calls}),
                    "complex_host_evidence": int(complex_host), "diagnostics_json": json_key(diagnostics)})
                if example:
                    if key not in events:
                        c = example["candidate"]
                        events[key] = {**c, "sample": args.sample, "library_side": side, "event_id": identifier,
                            "event_key_json": key, "chrom": c.get("chrom", c.get("host_chrom", example["chrom"])),
                            "host_strand": c.get("host_strand", example["host_strand"]),
                            "pair_support": 0, "simple_pair_support": 0, "prior_dd_pair_support": 0,
                            "complex_host_pair_support": 0, "example_read_id": rid, "example_mate": example["mate"],
                            "example_candidate_json": json_key(c), "_example": example}
                    events[key]["pair_support"] += 1
                    field = {"called_simple_candidate": "simple_pair_support", "called_prior_donor_rearrangement": "prior_dd_pair_support", "called_complex_host_evidence": "complex_host_pair_support"}[status]
                    events[key][field] += 1
    fasta.close()
    if seen != all_ids or len(seen) != total_expected or sum(counts.values()) != len(selected):
        raise ValueError("Final FASTQ/evidence/pair counts do not reconcile")
    expected_stream = manifest["inputs"].get("normalized_paired_stream_sha256")
    if expected_stream and digest.hexdigest() != expected_stream:
        raise ValueError("Normalized FASTQ stream differs from upstream full run")
    ordered = sorted(events.values(), key=lambda r: (-r["pair_support"], r["event_id"]))
    with (staging / "junctions.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, EVENT_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(ordered)
    write_representatives(staging / "representative_alignments.txt", ordered, donor, side)
    called = sum(counts[s] for s in STATUSES if s.startswith("called_"))
    stats = {"total_input_pairs": total_expected, "host_supported_pairs": len(selected),
             **{s: counts[s] for s in STATUSES}, "called_pairs_total": called,
             "observed_junction_patterns_total": len(events),
             "simple_candidate_patterns": sum(r["simple_pair_support"] > 0 for r in ordered),
             "attempted_host_supported_reads": attempted_reads,
             "called_fraction_of_host_supported_pairs": called / len(selected) if selected else None}
    summary = {"analysis_type": ANALYSIS_TYPE, "analysis_version": VERSION, "status": "complete",
        "sample": args.sample, "library_side": side, "condition": manifest.get("condition", ""),
        "biological_replicate": manifest.get("biological_replicate", ""), "stats": stats,
        "python_version": sys.version,
        "completed_utc": utc_now(), "elapsed_seconds": time.monotonic() - started,
        "config": {**core.DEFAULT_CONFIG, "host_window_padding_bp": args.window,
                   "complex_host_locus_distance_bp": 5000, "pair_merging": False},
        "upstream_results_directory": str(directory), "upstream_analysis_version": manifest.get("analysis_version"),
        "upstream_input_sha256": hashes, "upstream_inputs": manifest["inputs"],
        "reference": {"path": str(reference), "sha256": reference_sha, "donor_sha256": donor_sha,
                      "upstream_reference": manifest["reference"], "explicit_override": args.reference is not None},
        "source_sha256": {name: sha256_file(Path(__file__).resolve().parent / name)
                          for name in ("raising_junction_pilot.py", "junction_core.py")},
        "validation": {"all_input_ids_reconciled": seen == all_ids,
                       "host_pair_statuses_reconcile": sum(counts.values()) == len(selected),
                       "event_support_reconciles": sum(r["pair_support"] for r in ordered) == called,
                       "event_strata_reconcile": all(r["pair_support"] == r["simple_pair_support"] + r["prior_dd_pair_support"] + r["complex_host_pair_support"] for r in ordered)},
        "interpretation": "Pilot exact-anchor same-read donor-human junction candidates; not all host-mapped pairs, independent molecules, integration events, or human-side deletion calls.",
        "limitations": ["No gapped realignment, overlapping-mate consensus, or mate-only junction inference.",
                        "Primer-distal exact donor anchors impose deletion-dependent ascertainment.",
                        "Unresolved and ambiguous cases have no MH/deletion value; not zero.",
                        "Insertion-positive joins must be stratified separately from blunt/MH joins.",
                        "Prior donor-rearrangement and discordant host-locus evidence are separate from simple candidates."]}
    if not all(summary["validation"].values()):
        raise ValueError("Internal output count reconciliation failed")
    summary["output_sha256"] = {p.name: sha256_file(p) for p in staging.iterdir() if p.is_file()}
    (staging / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if outdir.exists():
        for path in staging.iterdir():
            os.replace(path, outdir / path.name)
        staging.rmdir()
    else:
        staging.rename(outdir)
    print(f"COMPLETE: {called:,}/{len(selected):,} host-supported pairs called; {len(events):,} observed patterns.", flush=True)
    print("Status counts: " + json_key({s: counts[s] for s in STATUSES}), flush=True)
    print("Results: " + str(outdir), flush=True)
    return summary


def main(argv=None):
    try:
        run(parse_args(argv))
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
