#!/usr/bin/env python3
"""Conservative RAISING read-pair classification after technical trimming.

Requires raising_concatemer_screen.py beside this file. Prepare-only needs
Python's standard library; final analysis also needs blastn, makeblastdb, bwa
and a local human genomic FASTA. Existing inputs and earlier results are kept.
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import raising_concatemer_screen as donor_screen

R1_ADAPTER = "TACTGCTTGACTGACCTGCTGTC"
R2_ADAPTER = donor_screen.revcomp(R1_ADAPTER)
DONOR_CONTIG = "__RAISING_DONOR_800__"
FINAL_CLASSES = ("donor_donor_supported", "host_mappable_other", "unresolved")


def write_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def file_metadata(path):
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else 10**9


def homopolymer_positions(sequence, minimum=10):
    positions = set()
    for match in re.finditer(r"([ACGT])\1{%d,}" % (minimum - 1), sequence):
        positions.update(range(*match.span()))
    return positions


def confirmed_tail(sequence, adapter=R1_ADAPTER, max_mismatches=2):
    """Return polyA start only when its adjacent complete adapter is recognized.

    An A-only tract is not sufficient. Allow a one-base boundary ambiguity when
    a leading adapter T was sequenced as A and absorbed into the polyA run.
    """
    candidates = []
    for run in re.finditer("A{15,}", sequence):
        for end in (run.end(), run.end() - 1):
            if end - run.start() >= 15 and hamming(sequence[end:end + len(adapter)], adapter) <= max_mismatches:
                candidates.append((run.start(), end))
    return min(candidates)[0] if candidates else None


def clean_record(record, mate, side):
    header, raw, plus, quality = record
    flags = set()
    left, right = 0, len(raw)
    if mate == "R1":
        _, left, _ = donor_screen.trim_leading_primer(raw, side)
        if left:
            flags.add("leading_P3_trimmed")
        else:
            flags.add("leading_P3_not_recognized")
        tail = confirmed_tail(raw)
        if tail is not None:
            right = max(left, tail)
            flags.add("confirmed_polyA_adapter_trimmed")
        elif re.search("A{15,}", raw):
            flags.add("polyA_without_confirmed_adapter")
    else:
        # This adapter is a library-end feature; require it at the very start.
        if hamming(raw[:len(R2_ADAPTER)], R2_ADAPTER) <= 2:
            left = len(R2_ADAPTER)
            run = re.match("T{10,}", raw[left:])
            if run:
                left += run.end()
                flags.add("R2_adapter_polyT_trimmed")
            else:
                flags.add("R2_adapter_trimmed_polyT_not_recognized")
        else:
            flags.add("R2_adapter_not_recognized")
        # Read-through can reach the other physical PCR primer. Keep internal
        # primer motifs; donor alignments mask them for host evidence below.
    technical_right = right
    # Trim only a low-quality 3' edge, not isolated internal poor bases.
    while right > left:
        window = quality[max(left, right - 5):right]
        if sum(ord(c) - 33 for c in window) / len(window) >= 20:
            break
        right -= 1
    if right < technical_right:
        flags.add("low_quality_3prime_trimmed")
    sequence = raw[left:right]
    cleaned_quality = quality[left:right]
    if homopolymer_positions(sequence):
        flags.add("low_complexity_bases_masked_for_host")
    if len(sequence) < 30:
        flags.add("cleaned_length_below_30")
    return {"header": header, "read_id": donor_screen.normalized_id(header),
            "mate": mate, "raw_length": len(raw), "raw_start0": left,
            "raw_end0": right, "sequence": sequence, "quality": cleaned_quality,
            "flags": flags, "technical_suffix_bp": len(raw) - technical_right,
            "quality_suffix_bp": technical_right - right}


def load_records(r1path, r2path, side):
    records = {"R1": {}, "R2": {}}
    extra_r2 = 0
    for mate, path in (("R1", r1path), ("R2", r2path)):
        if path is None:
            continue
        seen = set()
        for raw in donor_screen.read_fastq(path):
            rid = donor_screen.normalized_id(raw[0])
            if rid in seen:
                raise ValueError(f"Duplicate {mate} read ID: {rid}")
            seen.add(rid)
            if mate == "R2" and rid not in records["R1"]:
                extra_r2 += 1
                continue
            records[mate][rid] = clean_record(raw, mate, side)
    if not records["R1"]:
        raise ValueError("R1 input is empty")
    return records, extra_r2


def write_cleaned_fastq(path, records):
    count = 0
    with gzip.open(path, "wt", encoding="ascii") as handle:
        for rec in records.values():
            if rec["sequence"]:
                handle.write(rec["header"] + "\n" + rec["sequence"] + "\n+\n" + rec["quality"] + "\n")
                count += 1
    return count


def qc_row(rec):
    return {"read_id": rec["read_id"], "mate": rec["mate"], "raw_length": rec["raw_length"],
            "cleaned_length": len(rec["sequence"]), "cleaned_raw_start_1based": rec["raw_start0"] + 1 if rec["sequence"] else "",
            "cleaned_raw_end_1based": rec["raw_end0"] if rec["sequence"] else "",
            "trimmed_prefix_bp": rec["raw_start0"], "technical_suffix_bp": rec["technical_suffix_bp"],
            "quality_suffix_bp": rec["quality_suffix_bp"], "qc_flags": ";".join(sorted(rec["flags"]))}


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_command(command, commands, work, stdout_path=None):
    commands.append([str(x) for x in command])
    write_json(work / "commands.json", commands)
    with (work / "tools.stderr.log").open("a", encoding="utf-8") as log:
        log.write("\nCOMMAND " + json.dumps(commands[-1]) + "\n")
        log.flush()
        if stdout_path is not None:
            with stdout_path.open("w", encoding="utf-8") as out:
                result = subprocess.run(command, stdout=out, stderr=log, text=True)
        else:
            result = subprocess.run(command, stdout=log, stderr=log, text=True)
    if result.returncode:
        raise ValueError(f"Command failed ({result.returncode}); see {work / 'tools.stderr.log'}")


def tool_versions(commands, work):
    versions = {}
    for command in (["bwa"], ["blastn", "-version"], ["makeblastdb", "-version"]):
        commands.append(command)
        result = subprocess.run(command, capture_output=True, text=True)
        lines = (result.stdout + result.stderr).strip().splitlines()
        selected = [line for line in lines if "Version:" in line] if command[0] == "bwa" else lines[:3]
        versions[command[0]] = {"executable": shutil.which(command[0]), "version": "\n".join(selected),
                               "version_command_exit_code": result.returncode}
    write_json(work / "commands.json", commands)
    return versions


def prepare_reference(host, donor, root, commands, work):
    signature = {"host": file_metadata(host), "donor_sha256": hashlib.sha256(donor.encode()).hexdigest(),
                 "donor_contig": DONOR_CONTIG, "format_version": 1}
    digest = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:20]
    directory = root / "ref" / ("threeway_" + digest)
    directory.mkdir(parents=True, exist_ok=True)
    fasta, marker = directory / "human_plus_donor.fa", directory / "index_complete.json"
    indexes = [Path(str(fasta) + suffix) for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")]
    if marker.is_file() and fasta.is_file() and all(p.is_file() and p.stat().st_size for p in indexes):
        data = json.loads(marker.read_text())
        if data.get("signature") == signature:
            print("Reusing completed BWA reference index.", flush=True)
            return fasta, data
    # A marker is written only after bwa index succeeds. Incomplete indexes
    # from an interrupted run are never reused.
    marker.unlink(missing_ok=True)
    print("Preparing human + donor FASTA; initial BWA indexing may take hours.", flush=True)
    temporary = directory / "human_plus_donor.fa.tmp"
    opener = gzip.open if host.suffix == ".gz" else open
    sequences, names, total_bases = 0, set(), 0
    with opener(host, "rt", encoding="ascii") as source, temporary.open("w", encoding="ascii") as target:
        for line in source:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                name = stripped[1:].split()[0]
                if name == DONOR_CONTIG or name in names:
                    raise ValueError("Duplicate/reserved reference name: " + name)
                names.add(name)
                sequences += 1
            elif not sequences or set(stripped.upper()) - set("ACGTNRYKMSWBDHV"):
                raise ValueError("Human reference must be a nucleotide genomic FASTA")
            else:
                total_bases += len(stripped)
            target.write(stripped + "\n")
        target.write(">" + DONOR_CONTIG + "\n" + donor + "\n")
    if not sequences or not total_bases:
        raise ValueError("Empty human reference")
    temporary.replace(fasta)
    run_command(["bwa", "index", "-a", "bwtsw", str(fasta)], commands, work)
    if not all(p.is_file() and p.stat().st_size for p in indexes):
        raise ValueError("BWA index did not produce every required index file")
    data = {"signature": signature, "combined_fasta": str(fasta), "host_sequences": sequences,
            "host_bases": total_bases, "completed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    write_json(marker, data)
    return fasta, data


def align_donor(records, donor, sample, threads, work, commands):
    seq_to_uid, queries = {}, {}
    for group in records.values():
        for rec in group.values():
            sequence = rec["sequence"]
            rec["uid"] = ""
            if len(sequence) >= 20:
                uid = seq_to_uid.setdefault(sequence, "U%08d" % (len(seq_to_uid) + 1))
                queries[uid] = sequence
                rec["uid"] = uid
    query_path, donor_path = work / "unique_cleaned.fasta", work / "donor.fasta"
    query_path.write_text("".join(">" + uid + "\n" + seq + "\n" for uid, seq in queries.items()), encoding="ascii")
    donor_path.write_text(">donor\n" + donor + "\n", encoding="ascii")
    hits = collections.defaultdict(list)
    if queries:
        database = work / "donor_db"
        run_command(["makeblastdb", "-in", str(donor_path), "-dbtype", "nucl", "-out", str(database)], commands, work)
        blast_path = work / "donor_hits.tsv"
        print(f"Checking donor junctions in {len(queries):,} distinct cleaned sequences ...", flush=True)
        run_command(["blastn", "-query", str(query_path), "-db", str(database), "-task", "blastn-short",
                     "-word_size", "7", "-dust", "no", "-soft_masking", "false", "-evalue", "100",
                     "-max_hsps", "50", "-max_target_seqs", "1", "-num_threads", str(threads),
                     "-outfmt", "6 " + donor_screen.BLAST_FIELDS], commands, work, blast_path)
        primer_mask = donor_screen.primer_reference_mask(donor, sample)
        with blast_path.open() as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                if float(fields[5]) >= 95 and int(fields[7]) >= 20:
                    hits[fields[0]].append(donor_screen.parse_hit(fields, queries, donor, primer_mask))
    params = SimpleNamespace(min_anchor=30, min_identity=95.0, max_gap=50, max_overlap=20, min_rewind=20)
    for group in records.values():
        for rec in group.values():
            hsps = hits.get(rec["uid"], [])
            rec["donor_mask"] = set().union(*(set(range(h.q0, h.q1)) for h in hsps)) if hsps else set()
            # All donor-covered bases, INCLUDING P3 coordinates, are excluded
            # from host support. P3 exclusion only affects independent DD anchors.
            result = donor_screen.classify_hits(hsps, len(rec["sequence"]), params)
            rec["dd"] = result["classification"] == donor_screen.CLASSES[0]
            rec["donor_class"] = result["classification"]
            rec["donor_details"] = donor_screen.best_fields(result, rec["raw_start0"])


def cigar_events(cigar, query_length, reverse=False):
    """Reference-column events for MD, with positions in original cleaned read.

    H contributes to full-query coordinates even though absent from SAM SEQ.
    S/I consume query; D consumes a reference column represented by None.
    """
    operations = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]
    if not operations or "".join(str(n) + op for n, op in operations) != cigar:
        raise ValueError("Invalid CIGAR: " + cigar)
    qpos, events, aligned, inserted = 0, [], 0, 0
    for length, op in operations:
        if op in "M=X":
            events.extend((query_length - 1 - i if reverse else i) for i in range(qpos, qpos + length))
            qpos += length
            aligned += length
        elif op in "ISH":
            qpos += length
            inserted += length if op == "I" else 0
        elif op == "D":
            events.extend([None] * length)
        elif op == "N":
            raise ValueError("Spliced CIGAR not supported for this genomic workflow")
    if qpos != query_length:
        raise ValueError(f"CIGAR/query length mismatch: {cigar}, {query_length}")
    return events, aligned, inserted


def md_matches(md, events):
    tokens = re.findall(r"\d+|\^[A-Za-z]+|[A-Za-z]", md)
    if "".join(tokens) != md:
        raise ValueError("Invalid MD tag")
    cursor, matches = 0, set()
    for token in tokens:
        if token.isdigit():
            length = int(token)
            part = events[cursor:cursor + length]
            if len(part) != length or any(x is None for x in part):
                raise ValueError("MD match run inconsistent with CIGAR")
            matches.update(part)
            cursor += length
        elif token.startswith("^"):
            length = len(token) - 1
            if events[cursor:cursor + length] != [None] * length:
                raise ValueError("MD deletion inconsistent with CIGAR")
            cursor += length
        else:
            if cursor >= len(events) or events[cursor] is None:
                raise ValueError("MD mismatch inconsistent with CIGAR")
            cursor += 1
    if cursor != len(events):
        raise ValueError("MD does not cover CIGAR alignment")
    return matches


def host_evidence(fields, rec):
    flag, reference, mapq, cigar = int(fields[1]), fields[2], int(fields[4]), fields[5]
    if flag & (4 | 256) or reference == DONOR_CONTIG:
        return None
    tags = {}
    for text in fields[11:]:
        parts = text.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    evidence = {"reference": reference, "position_1based": int(fields[3]), "mapq": mapq,
                "cigar": cigar, "strand": "-" if flag & 16 else "+", "supplementary": bool(flag & 2048),
                "accepted": False, "reasons": []}
    if "NM" not in tags:
        evidence["reasons"].append("host_NM_unavailable")
        return evidence
    events, aligned, inserted = cigar_events(cigar, len(rec["sequence"]), bool(flag & 16))
    nm = int(tags["NM"])
    # Conservative identity: NM includes deletions, while the denominator is
    # aligned query bases (M/= /X + I), not all alignment columns including D.
    identity = max(0.0, 100 * (1 - nm / max(1, aligned + inserted)))
    low_complexity = homopolymer_positions(rec["sequence"])
    allowed = {i for i in events if i is not None and i not in rec["donor_mask"]
               and i not in low_complexity
               and ord(rec["quality"][i]) - 33 >= 20 and rec["sequence"][i] in "ACGT"}
    if "MD" in tags:
        matched = allowed & md_matches(tags["MD"], events)
        independent = len(matched)
        match_method = "MD_exact_matches"
    else:
        matched = allowed
        independent = max(0, len(allowed) - nm)
        match_method = "conservative_all_NM_subtracted"
    evidence.update({"identity_percent": identity, "independent_host_matches": independent,
                     "match_method": match_method, "query_start_1based": min(matched) + 1 if matched else "",
                     "query_end_1based": max(matched) + 1 if matched else "",
                     "raw_query_start_1based": min(matched) + rec["raw_start0"] + 1 if matched else "",
                     "raw_query_end_1based": max(matched) + rec["raw_start0"] + 1 if matched else ""})
    if mapq == 255:
        evidence["reasons"].append("host_MAPQ_unavailable")
    elif mapq < 20:
        evidence["reasons"].append("host_low_MAPQ")
    if identity < 95:
        evidence["reasons"].append("host_identity_below_95")
    if independent < 30:
        evidence["reasons"].append("host_18_29_independent_matches" if independent >= 18 else "host_insufficient_independent_matches")
    evidence["accepted"] = not evidence["reasons"]
    return evidence


def map_host(records, fastqs, combined, threads, commands, work):
    for mate, group in records.items():
        for rec in group.values():
            rec["host"] = False
            rec["host_hits"] = []
        if not any(rec["sequence"] for rec in group.values()):
            continue
        print(f"Mapping {mate} independently to human + donor ...", flush=True)
        sam = work / (mate + ".sam")
        run_command(["bwa", "mem", "-a", "-Y", "-T", "20", "-t", str(threads),
                     str(combined), str(fastqs[mate])], commands, work, sam)
        with sam.open() as handle:
            for line in handle:
                if line.startswith("@"):
                    continue
                fields = line.rstrip().split("\t")
                rid = donor_screen.normalized_id(fields[0])
                if rid not in group:
                    raise ValueError("SAM read ID absent from input: " + rid)
                evidence = host_evidence(fields, group[rid])
                if evidence is not None:
                    group[rid]["host_hits"].append(evidence)
                    group[rid]["host"] |= evidence["accepted"]


def assign_pair(ends):
    dd, host = any(r["dd"] for r in ends), any(r["host"] for r in ends)
    label = FINAL_CLASSES[0] if dd else FINAL_CLASSES[1] if host else FINAL_CLASSES[2]
    reasons = set()
    if label == "unresolved":
        for rec in ends:
            if len(rec["sequence"]) < 30:
                reasons.add("short_informative_sequence")
            if rec["donor_mask"]:
                reasons.add("donor_signal_without_supported_junction")
            if "confirmed_polyA_adapter_trimmed" in rec["flags"]:
                reasons.add("confirmed_technical_tail")
            if "low_quality_3prime_trimmed" in rec["flags"]:
                reasons.add("low_quality_tail")
            for evidence in rec["host_hits"]:
                reasons.update(evidence["reasons"])
        if not reasons:
            reasons.add("no_supported_donor_junction_or_host_mapping")
    return label, dd, host, reasons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.home() / "raising")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--r1", type=Path)
    parser.add_argument("--r2", type=Path)
    parser.add_argument("--host-fasta", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if Path(args.sample).name != args.sample or args.sample in (".", "..") or args.threads < 1:
        parser.error("Invalid sample name or thread count")
    if not args.prepare_only and not args.host_fasta:
        parser.error("Final classification requires --host-fasta; use --prepare-only for technical QC only")
    root = args.root.expanduser().resolve()
    outdir = root / "results" / (args.sample + ("_threeway_prepare" if args.prepare_only else "_threeway"))
    work = outdir / "work"
    r1path = (args.r1 or root / "results" / (args.sample + "_first_check") / "sample_10000_R1.fastq.gz").expanduser().resolve()
    r2path = (args.r2 or root / "results" / (args.sample + "_concatemer_screen") / "matched_sample_R2.fastq.gz").expanduser().resolve()
    if args.host_fasta:
        args.host_fasta = args.host_fasta.expanduser().resolve()
    for input_path in (r1path, r2path, args.host_fasta):
        if input_path is not None and (input_path == outdir or outdir in input_path.parents):
            raise ValueError("Input files must be outside the selected output directory: " + str(input_path))
    work.mkdir(parents=True, exist_ok=True)
    run_id = datetime.datetime.now(datetime.timezone.utc).isoformat() + "_pid" + str(os.getpid())
    if not args.prepare_only:
        write_json(outdir / "threeway_summary.json", {"status": "running", "run_id": run_id,
                   "sample": args.sample, "classification_available": False,
                   "note": "This run is incomplete. Existing TSVs may belong to an earlier run; use only after status complete."})
    if not r1path.is_file() or not r2path.is_file():
        raise ValueError("Both sampled FASTQ files are required; supply --r1 and --r2 if stored elsewhere")
    if not args.prepare_only:
        for binary in ("bwa", "blastn", "makeblastdb"):
            if shutil.which(binary) is None:
                raise ValueError(binary + " not found; activate the raising environment")
        if not args.host_fasta.is_file():
            raise ValueError("Human genomic FASTA not found: " + str(args.host_fasta))
    donor = donor_screen.load_donor(root / "ref" / "RI_donor_800bp.fasta")
    side = "3p" if "_3prm" in args.sample else "5p" if "_5prm" in args.sample else None
    records, extra_r2 = load_records(r1path, r2path, side)
    fastqs = {mate: outdir / ("cleaned_" + mate + ".fastq.gz") for mate in records}
    exported = {mate: write_cleaned_fastq(fastqs[mate], group) for mate, group in records.items()}
    qc = [qc_row(rec) for group in records.values() for rec in group.values()]
    write_tsv(outdir / "technical_qc.tsv", qc, list(qc[0]))
    commands = []
    summary = {"sample": args.sample, "run_id": run_id, "status": "complete" if args.prepare_only else "running",
               "classification_available": False,
               "mode": "prepare_only_no_final_classification" if args.prepare_only else "final_threeway",
               "inputs": {"R1": file_metadata(r1path), "R2": file_metadata(r2path)},
               "input_R1_reads": len(records["R1"]), "matching_R2_reads": len(records["R2"]),
               "missing_R2_partners": len(records["R1"]) - len(records["R2"]), "ignored_extra_R2_reads": extra_r2,
               "nonempty_cleaned_fastq_reads": exported,
               "qc_flags_by_mate": {m: dict(collections.Counter(flag for r in group.values() for flag in r["flags"])) for m, group in records.items()},
               "thresholds": {"donor_junction_independent_matches_each_side": 30, "donor_identity_percent": 95,
                              "host_independent_matches": 30, "host_identity_percent": 95, "host_min_MAPQ": 20,
                              "host_min_base_quality": 20, "host_low_complexity_homopolymer_length": 10,
                              "donor_mask_min_HSP_matches": 20, "technical_adapter_max_substitutions": 2},
               "notes": ["Read/pair proportions describe the amplified sequencing library, not independent cells or original molecule frequencies.",
                         "Donor-donor and host evidence are independent; exclusive categories prioritize donor-donor, then other host-mappable, then unresolved.",
                         "A donor-donor join may be a concatemer, partial rearrangement, circular donor, or PCR artifact; genomic integration is not established.",
                         "Host mapping alone, particularly R2-only mapping, does not prove a donor-human junction.",
                         "Unresolved is not proof of non-concatemer or absence of a genomic match.",
                         "The chosen human FASTA may differ from the vendor reference. Primary assembly excludes ALT/patch sequences and affects uniqueness.",
                         "BWA-MEM -T20 is a conservative screen; 18-29 bp hits are flagged only when emitted, not exhaustively searched.",
                         "Host identity is conservatively defined as 100*(1-NM/(aligned_query_bases+insertions)); deletions contribute to NM but not the denominator.",
                         "No original-input sequence or earlier analysis file is modified."]}
    print(f"R1: {len(records['R1']):,}; matching R2: {len(records['R2']):,}", flush=True)
    if args.prepare_only:
        write_json(outdir / "technical_qc.json", summary)
        print("Technical preparation complete. No human mapping or final three-way classification was performed.")
        print("Results: " + str(outdir))
        return 0
    summary["tool_versions"] = tool_versions(commands, work)
    align_donor(records, donor, args.sample, args.threads, work, commands)
    combined, reference_metadata = prepare_reference(args.host_fasta, donor, root, commands, work)
    map_host(records, fastqs, combined, args.threads, commands, work)
    read_rows, evidence_rows, donor_rows = [], [], []
    for mate, group in records.items():
        for rec in group.values():
            row = qc_row(rec)
            row.update({"donor_donor_supported": int(rec["dd"]), "host_mapping_supported": int(rec["host"]),
                        "both_supported": int(rec["dd"] and rec["host"]), "donor_class": rec["donor_class"],
                        "junction_type": rec["donor_details"].get("junction_type", ""),
                        "donor_masked_bases": len(rec["donor_mask"]),
                        "host_best_independent_matches": max((h.get("independent_host_matches", 0) for h in rec["host_hits"]), default=0)})
            read_rows.append(row)
            if rec["dd"]:
                donor_rows.append({"read_id": rec["read_id"], "mate": mate,
                                   "cleaned_raw_start_1based": rec["raw_start0"] + 1,
                                   "cleaned_raw_end_1based": rec["raw_end0"], **rec["donor_details"]})
            for hit in rec["host_hits"]:
                item = {"read_id": rec["read_id"], "mate": mate, **hit}
                item["reasons"] = ";".join(item["reasons"])
                evidence_rows.append(item)
    pair_rows = []
    for rid, first in records["R1"].items():
        ends = [first] + ([records["R2"][rid]] if rid in records["R2"] else [])
        label, dd, host, reasons = assign_pair(ends)
        pair_rows.append({"read_id": rid, "classification": label, "donor_donor_supported": int(dd),
                          "host_mapping_supported": int(host), "both_supported": int(dd and host),
                          "R2_available": int(len(ends) == 2), "unresolved_reasons": ";".join(sorted(reasons))})
    write_tsv(outdir / "read_evidence.tsv", read_rows, list(read_rows[0]))
    write_tsv(outdir / "pair_classification.tsv", pair_rows, list(pair_rows[0]))
    donor_fields = ["read_id", "mate", "cleaned_raw_start_1based", "cleaned_raw_end_1based", "junction_type",
                    "candidate_alignments", "left_independent_matches", "right_independent_matches", "unaligned_gap_bp",
                    "alignment_overlap_bp", "reference_rewind_bp"]
    for side_name in ("left", "right"):
        donor_fields.extend(side_name + suffix for suffix in ("_read_start_1based", "_read_end_1based", "_donor_start_1based",
                            "_donor_end_1based", "_strand", "_identity_percent", "_aligned_read", "_aligned_donor"))
    write_tsv(outdir / "donor_junction_evidence.tsv", donor_rows, donor_fields)
    host_fields = ["read_id", "mate", "reference", "position_1based", "mapq", "cigar", "strand", "supplementary", "accepted",
                   "identity_percent", "independent_host_matches", "match_method", "query_start_1based", "query_end_1based",
                   "raw_query_start_1based", "raw_query_end_1based", "reasons"]
    write_tsv(outdir / "host_alignment_evidence.tsv", evidence_rows, host_fields)
    counts = collections.Counter(row["classification"] for row in pair_rows)
    summary.update({"status": "complete", "classification_available": True, "reference": reference_metadata, "commands": commands,
                    "pair_counts": {label: counts[label] for label in FINAL_CLASSES},
                    "pair_percentages": {label: 100 * counts[label] / len(pair_rows) for label in FINAL_CLASSES},
                    "pairs_with_both_supports": sum(row["both_supported"] for row in pair_rows),
                    "unresolved_reason_counts_nonexclusive": dict(collections.Counter(reason for row in pair_rows for reason in row["unresolved_reasons"].split(";") if reason))})
    write_json(outdir / "threeway_summary.json", summary)
    for label in FINAL_CLASSES:
        print(f"{label}: {counts[label]:,} / {len(pair_rows):,} ({100 * counts[label] / len(pair_rows):.2f}%)")
    print("Results: " + str(outdir))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, csv.Error) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        raise SystemExit(1)
