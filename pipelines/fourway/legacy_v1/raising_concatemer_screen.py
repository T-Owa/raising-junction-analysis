#!/usr/bin/env python3
"""Conservative donor-to-donor junction screen for RAISING R1/R2 FASTQ.

Optional standalone diagnostic screen (the full-depth runner imports this module):
  python raising_concatemer_screen.py --root ./work --sample sample001_3prm_S1 --all

Default: screen the previously sampled 10,000 R1 reads and retrieve their R2
partners if available. Add --all for all raw R1 reads, or --r1-only to skip R2.
Requires installed blastn and makeblastdb; Python dependencies are standard-only.

A positive means a sequenced end supports a donor-to-donor junction or donor
rearrangement. It does not establish genomic integration, molecule counts, or
exclude a circular donor/PCR chimera. A negative is not a non-concatemer call.
MH and exact deletion sizes are not reported by this initial screen.
"""

import argparse
import collections
import csv
import gzip
import itertools
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

P3 = {
    "3p": "GCTACTCTGGGCTCGATCATCCA",
    "5p": "GCACGCCGATCCGACACACTGAC",
}
BLAST_FIELDS = "qseqid qstart qend sstart send pident length nident mismatch gapopen evalue bitscore qseq sseq"
CLASSES = (
    "donor_donor_junction_supported",
    "donor_signal_without_supported_junction",
    "no_qualifying_donor_alignment",
    "insufficient_sequence_for_two_anchors",
)


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn-", "TGCANtgcan-"))[::-1]


def normalized_id(header):
    text = header.lstrip("@").split()[0]
    return text[:-2] if text.endswith(("/1", "/2")) else text


def read_fastq(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as handle:
        record_no = 0
        while True:
            header = handle.readline()
            if not header:
                break
            rest = [handle.readline() for _ in range(3)]
            record_no += 1
            if not all(rest):
                raise ValueError(f"Truncated FASTQ at record {record_no}: {path}")
            header, seq, plus, qual = [s.rstrip("\r\n") for s in [header] + rest]
            if not header.startswith("@") or not plus.startswith("+") or not seq or len(seq) != len(qual):
                raise ValueError(f"Invalid four-line FASTQ at record {record_no}: {path}")
            yield header, seq.upper(), plus, qual


def load_donor(path):
    records = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if line.startswith(">"):
            records.append("")
        elif line:
            if not records:
                raise ValueError("FASTA sequence appears before a header")
            records[-1] += line.upper()
    if len(records) != 1 or len(records[0]) != 800 or set(records[0]) - set("ACGT"):
        raise ValueError("Expected the supplied single 800-bp A/C/G/T donor FASTA")
    return records[0]


def primer_reference_mask(donor, sample):
    """Primer-aligned bases never contribute to the independent-anchor count."""
    side = "3p" if "_3prm" in sample else "5p" if "_5prm" in sample else None
    selected = {side: P3[side]} if side else P3
    mask = set()
    for primer in selected.values():
        patterns = [primer, revcomp(primer)]
        hits = [(donor.index(p), len(p)) for p in patterns if donor.count(p) == 1]
        if len(hits) != 1:
            raise ValueError("Final primer is not uniquely identified in donor")
        start, length = hits[0]
        mask.update(range(start, start + length))
    return mask


def trim_leading_primer(sequence, side=None):
    """Trim a leading P3 only; retain biological internal/trailing motifs.

    Reference-coordinate masking separately prevents any P3-aligned bases
    from contributing to independent junction anchors, including R2 tails.
    """
    primers = [P3[side]] if side else list(P3.values())
    patterns = list(dict.fromkeys(primers + [revcomp(p) for p in primers]))
    front = []
    for primer in patterns:
        for start in range(min(12, len(sequence) - len(primer)) + 1):
            differences = sum(a != b for a, b in zip(sequence[start:start + len(primer)], primer))
            if differences <= 2:
                front.append((differences, start, len(primer)))
    left = 0
    if front:
        _, start, length = min(front)
        left = start + length
    return sequence[left:], left, 0


@dataclass
class Hit:
    q0: int
    q1: int
    s0: int
    s1: int
    strand: str
    oriented_start: int
    oriented_end: int
    identity: float
    bitscore: float
    matched_query: frozenset
    independent_query: frozenset
    aligned_query: str
    aligned_donor: str


def parse_hit(fields, query_sequences, donor, mask):
    uid = fields[0]
    qs, qe, ss, se = map(int, fields[1:5])
    identity, bitscore = float(fields[5]), float(fields[11])
    qseq, sseq = fields[12].upper(), fields[13].upper()
    if qs > qe:
        qs, qe, ss, se = qe, qs, se, ss
        qseq, sseq = revcomp(qseq), revcomp(sseq)
    q0, q1 = qs - 1, qe
    if qseq.replace("-", "") != query_sequences[uid][q0:q1]:
        raise ValueError("BLAST query alignment disagrees with input sequence: " + uid)
    strand = "+" if ss <= se else "-"
    step = 1 if strand == "+" else -1
    qpos, spos = q0, ss - 1
    matched, independent = set(), set()
    for a, b in zip(qseq, sseq):
        if b != "-":
            expected = donor[spos] if step == 1 else revcomp(donor[spos])
            if b != expected:
                raise ValueError("BLAST subject alignment disagrees with supplied donor")
        if a in "ACGT" and a == b:
            matched.add(qpos)
            if spos not in mask:
                independent.add(qpos)
        qpos += a != "-"
        spos += step * (b != "-")
    if qpos != q1 or len(qseq) != len(sseq):
        raise ValueError("Invalid BLAST alignment coordinates")
    n = len(donor)
    oriented_start, oriented_end = (ss - 1, se) if step == 1 else (n - ss, n - se + 1)
    return Hit(q0, q1, min(ss, se) - 1, max(ss, se), strand,
               oriented_start, oriented_end, identity, bitscore,
               frozenset(matched), frozenset(independent), qseq, sseq)


def classify_hits(hits, query_length, params):
    if query_length < 2 * params.min_anchor:
        return {"classification": CLASSES[3], "best": None, "candidates": 0}
    hits = [h for h in hits if h.identity >= params.min_identity and len(h.independent_query) >= params.min_anchor]
    hits.sort(key=lambda h: (h.q0, h.q1, h.strand, h.s0))
    candidates = []
    for left, right in itertools.combinations(hits, 2):
        if right.q0 <= left.q0 or right.q1 <= left.q1:
            continue
        overlap = max(0, left.q1 - right.q0)
        gap = max(0, right.q0 - left.q1)
        if overlap > params.max_overlap or gap > params.max_gap:
            continue
        left_unique = sum(q < right.q0 for q in left.independent_query)
        right_unique = sum(q >= left.q1 for q in right.independent_query)
        if min(left_unique, right_unique) < params.min_anchor:
            continue
        if left.strand == right.strand:
            reference_gap = right.oriented_start - left.oriented_end
            # An unaligned insertion in a monomer is not a donor-copy rewind.
            # Shared query overlap is also insufficient evidence of duplication.
            rewind = -reference_gap - overlap
            if rewind < params.min_rewind:
                continue
            kind = "H-T_like"
        else:
            rewind = None
            kind = "T-T_like" if left.strand == "+" else "H-H_like"
        union = left.matched_query | right.matched_query
        # Prefer a continuous monomer explanation when it fits the same region.
        if any(h.q0 <= left.q0 + 2 and h.q1 >= right.q1 - 2 and
               len(h.matched_query & union) >= len(union) - 4 for h in hits):
            continue
        candidates.append({"kind": kind, "left": left, "right": right,
                           "left_independent_matches": left_unique,
                           "right_independent_matches": right_unique,
                           "query_gap": gap, "query_overlap": overlap,
                           "reference_rewind": rewind,
                           "union_matches": len(union),
                           "score": left.bitscore + right.bitscore})
    if candidates:
        best = max(candidates, key=lambda c: (c["union_matches"], c["score"]))
        return {"classification": CLASSES[0], "best": best, "candidates": len(candidates)}
    return {"classification": CLASSES[1] if hits else CLASSES[2], "best": None, "candidates": 0}


def run_command(command, logpath):
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with logpath.open("a", encoding="utf-8") as log:
        log.write(json.dumps(command) + "\n" + result.stdout + result.stderr + "\n")
    if result.returncode:
        raise ValueError("Command failed: " + command[0] + "\n" + result.stderr[-4000:])
    return result.stdout


def collect_reads(root, sample, use_all, r1_only):
    rawdir = root / "vendor_results" / "Raw_fastq"
    if use_all:
        r1paths = sorted(rawdir.glob(sample + "_L*_R1_*.fastq.gz"))
    else:
        r1paths = [root / "results" / (sample + "_first_check") / "sample_10000_R1.fastq.gz"]
    if not r1paths or any(not p.is_file() for p in r1paths):
        raise ValueError("R1 input missing. Supply the documented sampled FASTQ artifact, or use --all.")
    r1 = {}
    for path in r1paths:
        for rec in read_fastq(path):
            rid = normalized_id(rec[0])
            if rid in r1:
                raise ValueError("Duplicate R1 read ID: " + rid)
            r1[rid] = rec
    if not r1:
        raise ValueError("R1 input is empty")
    r2paths = [] if r1_only else sorted(rawdir.glob(sample + "_L*_R2_*.fastq.gz"))
    r2 = {}
    if r2paths:
        print("Retrieving R2 partners by read ID ...", flush=True)
        for path in r2paths:
            for rec in read_fastq(path):
                rid = normalized_id(rec[0])
                if rid in r1:
                    if rid in r2:
                        raise ValueError("Duplicate R2 partner ID: " + rid)
                    r2[rid] = rec
    return r1, r2, r1paths, r2paths


def best_fields(result, trim_start):
    best = result["best"]
    if best is None:
        return {"junction_type": "", "candidate_alignments": result["candidates"]}
    left, right = best["left"], best["right"]
    fields = {"junction_type": best["kind"], "candidate_alignments": result["candidates"],
              "left_independent_matches": best["left_independent_matches"],
              "right_independent_matches": best["right_independent_matches"],
              "unaligned_gap_bp": best["query_gap"], "alignment_overlap_bp": best["query_overlap"],
              "reference_rewind_bp": best["reference_rewind"]}
    for name, hit in [("left", left), ("right", right)]:
        fields.update({name + "_read_start_1based": hit.q0 + trim_start + 1,
                       name + "_read_end_1based": hit.q1 + trim_start,
                       name + "_donor_start_1based": hit.s0 + 1,
                       name + "_donor_end_1based": hit.s1,
                       name + "_strand": hit.strand,
                       name + "_identity_percent": hit.identity,
                       name + "_aligned_read": hit.aligned_query,
                       name + "_aligned_donor": hit.aligned_donor})
    return fields


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path.home() / "raising")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--all", action="store_true", help="Screen every raw R1 read instead of the saved sample")
    parser.add_argument("--r1-only", action="store_true", help="Do not retrieve or screen R2 partners")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-anchor", type=int, default=30, help="Independent matching bases required on EACH side")
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--max-gap", type=int, default=50)
    parser.add_argument("--max-overlap", type=int, default=20)
    parser.add_argument("--min-rewind", type=int, default=20)
    args = parser.parse_args()
    if Path(args.sample).name != args.sample or args.sample in (".", ".."):
        parser.error("--sample must be a name, not a path")
    if args.min_anchor < 20 or not 0 < args.min_identity <= 100 or min(args.threads, args.min_rewind) < 1 or min(args.max_gap, args.max_overlap) < 0:
        parser.error("Invalid alignment settings")
    for app in ["blastn", "makeblastdb"]:
        if shutil.which(app) is None:
            raise ValueError(app + " was not found. Run 'micromamba activate raising' first.")
    root = args.root.expanduser().resolve()
    donor = load_donor(root / "ref" / "RI_donor_800bp.fasta")
    mask = primer_reference_mask(donor, args.sample)
    side = "3p" if "_3prm" in args.sample else "5p" if "_5prm" in args.sample else None
    r1, r2, r1paths, r2paths = collect_reads(root, args.sample, args.all, args.r1_only)
    suffix = "_concatemer_screen_all" if args.all else "_concatemer_screen"
    outdir = root / "results" / (args.sample + suffix)
    work = outdir / "work"
    work.mkdir(parents=True, exist_ok=True)
    print(f"R1 reads: {len(r1):,}; matching R2 reads: {len(r2):,}", flush=True)
    seq_to_uid = {}
    queries = {}
    read_queries = {}
    for mate, records in [("R1", r1), ("R2", r2)]:
        for rid, record in records.items():
            trimmed, front, back = trim_leading_primer(record[1], side)
            uid = ""
            if len(trimmed) >= 2 * args.min_anchor:
                if trimmed not in seq_to_uid:
                    uid = f"U{len(seq_to_uid) + 1:08d}"
                    seq_to_uid[trimmed] = uid
                    queries[uid] = trimmed
                else:
                    uid = seq_to_uid[trimmed]
            read_queries[(mate, rid)] = (uid, front, back, len(trimmed))
    query_path = work / "unique_trimmed_reads.fasta"
    with query_path.open("w", encoding="ascii") as handle:
        for uid, sequence in queries.items():
            handle.write(f">{uid}\n{sequence}\n")
    donor_path = work / "donor.fasta"
    donor_path.write_text(">RI_donor_800bp\n" + donor + "\n", encoding="ascii")
    logpath = work / "commands.log"
    logpath.write_text("", encoding="utf-8")
    version = run_command(["blastn", "-version"], logpath).strip()
    db = work / "donor_db"
    hits_path = work / "blast_hits.tsv"
    blast_command = ["blastn", "-query", str(query_path), "-db", str(db),
                     "-task", "blastn-short", "-word_size", "7", "-dust", "no",
                     "-soft_masking", "false", "-evalue", "100", "-max_hsps", "50",
                     "-max_target_seqs", "1", "-num_threads", str(args.threads),
                     "-outfmt", "6 " + BLAST_FIELDS, "-out", str(hits_path)]
    hits = collections.defaultdict(list)
    if queries:
        run_command(["makeblastdb", "-in", str(donor_path), "-dbtype", "nucl", "-out", str(db)], logpath)
        print(f"Aligning {len(queries):,} unique trimmed sequences to the 800-bp donor ...", flush=True)
        run_command(blast_command, logpath)
        with hits_path.open(encoding="ascii") as handle:
            for line in handle:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 14:
                    raise ValueError("Unexpected BLAST output column count")
                if float(fields[5]) >= args.min_identity and int(fields[7]) >= args.min_anchor:
                    hit = parse_hit(fields, queries, donor, mask)
                    hits[fields[0]].append(hit)
    classified = {uid: classify_hits(hits.get(uid, []), len(seq), args) for uid, seq in queries.items()}
    fieldnames = ["read_id", "mate", "classification", "junction_type", "raw_read_length",
                  "trimmed_prefix_bp", "trimmed_suffix_bp", "left_independent_matches",
                  "right_independent_matches", "unaligned_gap_bp", "alignment_overlap_bp",
                  "reference_rewind_bp", "left_read_start_1based", "left_read_end_1based",
                  "left_donor_start_1based", "left_donor_end_1based", "left_strand", "left_identity_percent",
                  "right_read_start_1based", "right_read_end_1based", "right_donor_start_1based",
                  "right_donor_end_1based", "right_strand", "right_identity_percent",
                  "candidate_alignments", "left_aligned_read", "left_aligned_donor",
                  "right_aligned_read", "right_aligned_donor", "read_sequence"]
    summary_by_mate = {}
    positive_ids = set()
    paired_positive_ids = set()
    candidate_path = outdir / "concatemer_candidates.tsv"
    classification_path = outdir / "read_classification.tsv"
    with candidate_path.open("w", encoding="utf-8", newline="") as candidate_handle, classification_path.open("w", encoding="utf-8", newline="") as class_handle:
        candidates_writer = csv.DictWriter(candidate_handle, fieldnames=fieldnames, delimiter="\t")
        candidates_writer.writeheader()
        class_writer = csv.writer(class_handle, delimiter="\t")
        class_writer.writerow(["read_id", "mate", "classification", "junction_type"])
        for mate, records in [("R1", r1), ("R2", r2)]:
            counts = collections.Counter({name: 0 for name in CLASSES})
            types = collections.Counter()
            for rid, record in records.items():
                uid, front, back, length = read_queries[(mate, rid)]
                result = classified[uid] if uid else {"classification": CLASSES[3], "best": None, "candidates": 0}
                label = result["classification"]
                counts[label] += 1
                details = best_fields(result, front)
                class_writer.writerow([rid, mate, label, details["junction_type"]])
                if label == CLASSES[0]:
                    types[details["junction_type"]] += 1
                    positive_ids.add(rid)
                    if rid in r2:
                        paired_positive_ids.add(rid)
                    row = {"read_id": rid, "mate": mate, "classification": label,
                           "raw_read_length": len(record[1]), "trimmed_prefix_bp": front,
                           "trimmed_suffix_bp": back, "read_sequence": record[1]}
                    row.update(details)
                    candidates_writer.writerow(row)
            n = len(records)
            summary_by_mate[mate] = {"reads": n, "classes": dict(counts), "best_junction_type_counts": dict(types),
                                     "supported_percent_of_all_screened_reads": 100 * counts[CLASSES[0]] / n if n else None}
    if r2 and not args.all:
        with gzip.open(outdir / "matched_sample_R2.fastq.gz", "wt", encoding="ascii") as handle:
            for rid in r1:
                if rid in r2:
                    handle.write("\n".join(r2[rid]) + "\n")
    report = {
        "sample": args.sample, "scope": "all raw R1 reads" if args.all else "saved uniform random R1 sample",
        "input_R1_files": [str(p) for p in r1paths], "input_R2_files": [str(p) for p in r2paths],
        "blastn_version": version, "blastn_command": blast_command,
        "parameters": {"min_independent_matches_each_side": args.min_anchor,
                       "min_HSP_identity_percent": args.min_identity,
                       "max_unaligned_gap_bp": args.max_gap, "max_alignment_overlap_bp": args.max_overlap,
                       "min_same_strand_reference_rewind_bp": args.min_rewind,
                       "excluded_donor_reference_positions_1based": [p + 1 for p in sorted(mask)]},
        "by_mate": summary_by_mate,
        "read_pairs": {"input_R1_reads": len(r1), "matching_R2_reads": len(r2),
                       "missing_R2_partners": len(r1) - len(r2),
                       "R1_ids_with_a_supported_junction_in_an_available_end": len(positive_ids),
                       "supported_percent_of_input_R1_ids": 100 * len(positive_ids) / len(r1),
                       "complete_pairs_with_a_supported_junction": len(paired_positive_ids),
                       "supported_percent_of_complete_pairs": 100 * len(paired_positive_ids) / len(r2) if r2 else None},
        "notes": [
            "Positive: two donor alignments with independent anchors, inconsistent with a continuous monomer alignment.",
            "Final-primer reference bases and shared query overlap do not count toward the 30-base independent anchors.",
            "H-T_like includes reference rewinds/partial donor duplications; a full tandem-copy architecture is not established.",
            "H-H_like and T-T_like denote inverted donor junction candidates; internal donor rearrangements are also possible.",
            "Terminal deletions are allowed if both independent anchors remain. Short anchors or large inserted gaps may be missed.",
            "Alignment overlap is not a finalized microhomology length, and reference rewind is not a deletion length.",
            "No supported junction does not rule out concatemers when the junction is outside the read or its anchor is too short.",
            "This does not determine human mapping, genomic integration, original molecule frequencies, or PCR artifacts.",
            "Head-to-tail junctions can also be produced by circular donor molecules.",
            "Pair support requires at least one individual end to span a supported donor-to-donor junction; two donor-only ends are not enough.",
        ],
    }
    report_path = outdir / "concatemer_screen.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    for mate in ["R1", "R2"]:
        data = summary_by_mate[mate]
        if data["reads"]:
            print(f"{mate} junction-supported reads: {data['classes'][CLASSES[0]]:,} / {data['reads']:,} ({data['supported_percent_of_all_screened_reads']:.2f}%)")
            print(mate + " junction types: " + json.dumps(data["best_junction_type_counts"]))
    if r2:
        print(f"Complete pairs with >=1 supported end: {len(paired_positive_ids):,} / {len(r2):,} ({100 * len(paired_positive_ids) / len(r2):.2f}%)")
    print("Report: " + str(report_path))
    print("Candidates: " + str(candidate_path))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, csv.Error) as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        raise SystemExit(1)
