#!/usr/bin/env python3
"""Reassess existing RAISING junction templates with small-variant alignments.

Nearby boundaries form review groups, not independent integration events.
Only identical complete template sequences can share an equivalence group.
Minor alleles are never replaced by an abundant template merely because they
have low read support. No new junction templates are discovered by this tool.
"""
import argparse
import bisect
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
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "support"))
import pilot_io as io
import junction_core as exact
import refine_templates as templates_module
import refine_alignment as alignment_module

VERSION = "0.1.1"
ANALYSIS_TYPE = "raising_existing_template_gapped_reassessment"
STATUSES = ("template_supported", "template_ambiguous", "alignment_ambiguous",
            "junction_variant_review", "no_qualifying_template", "no_nearby_template")
PAIR_FIELDS = ["sample", "library_side", "read_id", "old_status", "old_event_id",
    "status", "change_type", "selected_equivalent_group", "selected_template_ids",
    "candidate_group_ids_json", "supporting_mates", "reference_compartment",
    "prior_donor_rearrangement_supported", "complex_host_evidence",
    "attempted_templates", "example_template_id", "example_mate",
    "example_substitution_bases", "example_gap_bases", "boundary_variant",
    "flags_json", "rejection_reasons_json"]
READ_FIELDS = ["read_id", "mate", "query_orientation", "template_id", "cluster_id",
    "comparison_context_id", "comparison_common_host_endpoint0",
    "effective_host_interval_start0", "effective_host_interval_end0",
    "equivalent_template_group", "chrom", "host_strand", "status", "score",
    "comparison_context_json", "best_alignment_json", "diagnostics_json"]


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def json_text(value):
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def compartment(chrom):
    return "mitochondrial_reference" if str(chrom).upper() in ("MT", "M", "CHRM", "CHRMT") else "nuclear_or_unplaced_reference"


def write_table(path, rows, fields=None, delimiter="\t"):
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json_text(v) if isinstance(v, (dict, list, tuple, set)) else v
                             for k, v in row.items()})


def summarize_upstream_reasons(old_rows):
    """Nonexclusive per-pair flags, not counts of orientation/alignment attempts."""
    selected = [r for r in old_rows if r["status"] in ("unresolved", "ambiguous")]
    denominators = collections.Counter(r["status"] for r in selected)
    counts = collections.Counter()
    for row in selected:
        reasons = set()
        details = row.get("diagnostics_json", "")
        details = json.loads(details) if isinstance(details, str) and details else details or []
        for attempt in details:
            diagnostic = attempt.get("diagnostics", {})
            for reason, count in diagnostic.get("rejected_counts", {}).items():
                if int(count) > 0:
                    reasons.add(reason)
        if row["status"] == "ambiguous":
            reasons.add("competing_exact_candidate_events")
        if not reasons:
            reasons.add("no_detailed_rejection_flag_recorded")
        for reason in reasons:
            counts[(row["status"], reason)] += 1
    return [{"old_status": status, "reason": reason, "affected_pairs": count,
             "denominator_pairs": denominators[status],
             "fraction": count / denominators[status], "nonexclusive": 1}
            for (status, reason), count in sorted(counts.items())]


def resolve_pair(read_results, group_members):
    """Rank comparable templates within each mate/accepted-host context.

    Support abundance is deliberately absent from this decision. Two mates
    support one pair. The intersection of their qualifying template groups can
    disambiguate ties; disjoint groups remain ambiguous. Alignment uncertainty
    and junction-proximal variants are retained even for a unique template.
    Display/review cluster IDs do NOT determine score competition. Every
    compared template must use the common human endpoint of its context.
    """
    qualifying = [r for r in read_results if r["result"]["status"] != "rejected"]
    if not qualifying:
        return {"status": "no_qualifying_template", "group": "", "template_ids": [],
                "candidate_groups": [], "supporting_mates": [], "selected_results": []}
    by_read_context = collections.defaultdict(list)
    for r in qualifying:
        by_read_context[(r["mate"], r["comparison_context_id"])].append(r)
    top = []
    for records in by_read_context.values():
        endpoints = {r["comparison_common_host_endpoint0"] for r in records
                     if "comparison_common_host_endpoint0" in r}
        if len(endpoints) > 1:
            raise ValueError("Competing templates do not share a human endpoint")
        best_score = max(float(r["result"].get("score") or 0) for r in records)
        top.extend(r for r in records if abs(float(r["result"].get("score") or 0) - best_score) < 1e-8)
    by_mate = collections.defaultdict(set)
    for r in top:
        by_mate[r["mate"]].add(r["equivalent_template_group"])
    sets = list(by_mate.values())
    union = set().union(*sets)
    intersection = set.intersection(*sets)
    chosen = intersection if intersection else union
    if len(chosen) != 1 or (len(sets) > 1 and not intersection):
        return {"status": "template_ambiguous", "group": "", "template_ids": [],
                "candidate_groups": sorted(union), "supporting_mates": sorted(by_mate),
                "selected_results": top}
    group = next(iter(chosen))
    selected = [r for r in top if r["equivalent_template_group"] == group]
    if any(r["result"]["status"] == "ambiguous_alignment" for r in selected):
        status = "alignment_ambiguous"
    elif any(r["result"]["status"] == "review" for r in selected):
        status = "junction_variant_review"
    else:
        status = "template_supported"
    return {"status": status, "group": group, "template_ids": sorted(group_members[group]),
            "candidate_groups": sorted(union),
            "supporting_mates": sorted({r["mate"] for r in selected}), "selected_results": selected}


def read_rows(path, required):
    # Input tables use ASCII DNA and normalized IDs, as in the frozen pilot.
    return list(io.tsv_rows(path, required))


def load_pilot(root, sample, explicit=None):
    directory = Path(explicit).expanduser().resolve() if explicit else root / "results" / (sample + "_junction_pilot_v1")
    path = directory / "summary.json"
    summary = json.loads(path.read_text(encoding="utf-8"))
    if summary.get("status") != "complete" or summary.get("sample") != sample or summary.get("analysis_type") != "raising_exact_anchor_junction_pilot":
        raise ValueError("Expected a complete matching exact-anchor pilot")
    if not summary.get("validation") or not all(summary["validation"].values()):
        raise ValueError("Pilot validation is not all true")
    paths = {"summary.json": path, "junctions.tsv": directory / "junctions.tsv",
             "pair_junction_status.tsv.gz": directory / "pair_junction_status.tsv.gz"}
    hashes = {}
    for name, file in paths.items():
        hashes[name] = io.sha256_file(file)
        if name != "summary.json" and hashes[name] != summary.get("output_sha256", {}).get(name):
            raise ValueError("Pilot output hash mismatch: " + str(file))
    rows = read_rows(paths["junctions.tsv"], ("event_id", "chrom", "host_strand", "host_boundary0", "example_candidate_json"))
    old_rows = read_rows(paths["pair_junction_status.tsv.gz"], ("read_id", "status", "event_id"))
    if len(rows) != int(summary["stats"]["observed_junction_patterns_total"]):
        raise ValueError("Pilot junction row count differs from its summary")
    if len(old_rows) != int(summary["stats"]["host_supported_pairs"]):
        raise ValueError("Pilot host-pair count differs from its summary")
    if len({r["event_id"] for r in rows}) != len(rows) or len({r["read_id"] for r in old_rows}) != len(old_rows):
        raise ValueError("Duplicate pilot pattern/pair identifier")
    return directory, summary, rows, old_rows, hashes


def template_index(templates):
    groups = collections.defaultdict(list)
    for t in templates:
        groups[(t["chrom"], t["host_strand"])].append(t)
    result = {}
    for key, members in groups.items():
        members.sort(key=lambda t: (int(t["pilot_row"]["host_boundary0"]), t["template_id"]))
        result[key] = ([int(t["pilot_row"]["host_boundary0"]) for t in members], members)
    return result


def nearby_templates(index, chrom, strand, intervals, padding=20):
    coordinates, members = index.get((chrom, strand), ([], []))
    selected = {}
    for start, end in intervals:
        lo = bisect.bisect_left(coordinates, start - padding)
        hi = bisect.bisect_right(coordinates, end + padding)
        for t in members[lo:hi]:
            selected[t["template_id"]] = t
    return selected


def assign_sequence_equivalence(templates, donor, side, fetch, contig_length):
    """Assign stable reconstructed-allele identity independently of clusters.

    A candidate is a fixed primer-distal donor prefix plus its intervening
    sequence, followed by the reference genome from its host boundary. Remove
    the longest suffix of that prefix identical to oriented reference bases
    immediately preceding the boundary, shifting the boundary accordingly.
    The remaining prefix and normalized boundary define the same sequence at
    every shared distal host endpoint. No substitutions or indels are erased.
    Original candidate IDs, boundary/MH descriptions and review groups remain.
    """
    donor_digest = hashlib.sha256(donor.upper().encode("ascii")).hexdigest()
    members = collections.defaultdict(list)
    for template in templates:
        chrom, strand = template["chrom"], template["host_strand"]
        boundary = int(template["pilot_row"]["host_boundary0"])
        prefix = template["sequence"][:int(template["host_start_t0"])].upper()
        length = int(contig_length(chrom))
        if not 0 <= boundary <= length:
            raise ValueError("Template boundary lies outside reference contig")
        if strand == "+":
            preceding = fetch(chrom, max(0, boundary - len(prefix)), boundary).upper()
        elif strand == "-":
            preceding = exact.reverse_complement(fetch(chrom, boundary, min(length, boundary + len(prefix))))
        else:
            raise ValueError("Invalid template host orientation")
        shared = 0
        for left, right in zip(reversed(prefix), reversed(preceding)):
            if left not in "ACGT" or left != right:
                break
            shared += 1
        normalized_boundary = boundary - shared if strand == "+" else boundary + shared
        remaining = prefix[:-shared] if shared else prefix
        key = ["reconstructed_allele_v1", donor_digest, side, chrom, strand,
               int(template["donor_start0"]), normalized_boundary, remaining]
        group = "E" + hashlib.sha256(json_text(key).encode("ascii")).hexdigest()[:16]
        template.setdefault("original_review_equivalent_group", template["equivalent_template_group"])
        template["equivalent_template_group"] = group
        template["equivalence_normalized_host_boundary0"] = normalized_boundary
        template["equivalence_remaining_prefix"] = remaining
        template["equivalence_reference_suffix_bp"] = shared
        template["equivalence_key_json"] = json_text(key)
        members[group].append(template["template_id"])
    for template in templates:
        template["equivalent_template_members"] = sorted(members[template["equivalent_template_group"]])
    return templates


def merge_intervals(intervals):
    """Union overlapping/touching genomic half-open alignment intervals."""
    merged = []
    for a, b in sorted(set((int(a), int(b)) for a, b in intervals)):
        if a < 0 or b <= a:
            raise ValueError("Invalid accepted genomic interval")
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
        else:
            merged.append((a, b))
    return merged


def accepted_host_contexts(host_hits, flip=False, merge_padding=0):
    """Build scoring contexts from overlapping accepted host-alignment spans.

    Contexts are distinct across contigs, oriented strands and query
    orientations. The default padding is zero: only overlapping/touching
    accepted alignment spans merge. Candidate-search padding is a separate
    parameter and never joins distinct accepted genomic loci for scoring.
    Actual M/=X intervals, not padded envelopes, validate independent support.
    """
    if not isinstance(merge_padding, int) or merge_padding < 0:
        raise ValueError("Context merge padding must be a nonnegative integer")
    grouped = collections.defaultdict(list)
    for number, hit in enumerate(host_hits):
        strand = ("-" if hit["strand"] == "+" else "+") if flip else hit["strand"]
        intervals = merge_intervals(hit["intervals0"])
        if not intervals:
            raise ValueError("Accepted host alignment has no reference intervals")
        grouped[(hit["reference"], strand)].append({
            "span_start0": intervals[0][0], "span_end0": intervals[-1][1],
            "accepted_intervals0": intervals, "accepted_alignment_indices": [number]})
    contexts = []
    for (chrom, strand), hits in sorted(grouped.items()):
        chunks = []
        for hit in sorted(hits, key=lambda x: (x["span_start0"], x["span_end0"])):
            if chunks and hit["span_start0"] <= chunks[-1]["span_end0"] + merge_padding:
                previous = chunks[-1]
                previous["span_end0"] = max(previous["span_end0"], hit["span_end0"])
                previous["accepted_intervals0"].extend(hit["accepted_intervals0"])
                previous["accepted_alignment_indices"].extend(hit["accepted_alignment_indices"])
            else:
                chunks.append({k: list(v) if isinstance(v, list) else v for k, v in hit.items()})
        for chunk in chunks:
            intervals = merge_intervals(chunk["accepted_intervals0"])
            identity = [chrom, strand, bool(flip), chunk["span_start0"], chunk["span_end0"], intervals]
            context_id = "H" + hashlib.sha256(json_text(identity).encode("ascii")).hexdigest()[:16]
            contexts.append({**chunk, "comparison_context_id": context_id,
                "chrom": chrom, "host_strand": strand, "flip": bool(flip),
                "query_orientation": "reverse_complement" if flip else "forward",
                "accepted_intervals0": intervals,
                "accepted_alignment_indices": sorted(chunk["accepted_alignment_indices"]),
                "context_merge_padding_bp": merge_padding})
    return contexts


def build_comparison_templates(context, index, by_id, old_event_id, donor, side,
                               fetch, contig_length, host_flank=120, padding=20):
    """Rebuild all competitors at this accepted locus over shared endpoints.

    The bounded display clusters and their original template identities remain
    unchanged. A context can include several display clusters. Its donor start
    is fixed by P3 and its distal host endpoint is shared by every candidate.
    The original assigned template can be retained beyond search padding only
    when its existing host-template interval overlaps >=30 accepted reference
    bases here. It is never injected into an unrelated same-contig locus.
    """
    selected = nearby_templates(index, context["chrom"], context["host_strand"],
                                context["accepted_intervals0"], padding=padding)
    old = by_id.get(old_event_id)
    if old and old["chrom"] == context["chrom"] and old["host_strand"] == context["host_strand"]:
        a, b = int(old["host_interval_start0"]), int(old["host_interval_end0"])
        overlap = sum(max(0, min(b, end) - max(a, start))
                      for start, end in context["accepted_intervals0"])
        if overlap >= 30:
            selected[old_event_id] = old
    if not selected:
        return []
    rows = [selected[key]["pilot_row"] for key in sorted(selected)]
    low = min(int(r["host_boundary0"]) for r in rows)
    high = max(int(r["host_boundary0"]) for r in rows)
    synthetic_group = {"cluster_id": context["comparison_context_id"],
        "chrom": context["chrom"], "host_strand": context["host_strand"],
        "boundary_min0": low, "boundary_max0": high, "rows": rows}
    rebuilt = templates_module.build_templates([synthetic_group], donor, side,
                                                fetch, contig_length, host_flank=host_flank)
    comparison_context = {**context, "candidate_template_ids": sorted(selected),
        "candidate_boundary_min0": low, "candidate_boundary_max0": high,
        "candidate_search_padding_bp": padding,
        "common_host_endpoint0": rebuilt[0]["common_host_endpoint0"],
        "common_donor_start0": rebuilt[0]["donor_start0"]}
    for template in rebuilt:
        original = selected[template["template_id"]]
        template["cluster_id"] = original["cluster_id"]
        template["equivalent_template_group"] = original["equivalent_template_group"]
        template["equivalent_template_members"] = original["equivalent_template_members"]
        template["comparison_context_id"] = context["comparison_context_id"]
        template["comparison_context"] = comparison_context
    return rebuilt


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--pilot-dir", type=Path)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--cluster-span", type=int, default=10)
    parser.add_argument("--host-flank", type=int, default=120)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def run(args):
    started = time.monotonic()
    root = args.root.expanduser().resolve()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.sample):
        raise ValueError("Unsafe sample name")
    if not 0 <= args.cluster_span <= 100 or args.host_flank < 60:
        raise ValueError("Require 0 <= cluster span <=100 and host flank >=60")
    print("Validating complete pilot and original full-run inputs ...", flush=True)
    pilot_dir, pilot, junction_rows, old_rows, pilot_hashes = load_pilot(root, args.sample, args.pilot_dir)
    old = {r["read_id"]: r for r in old_rows}
    directory, manifest_path, manifest = io.discover_results(root, args.sample, args.results_dir)
    full_summary_hash = io.sha256_file(manifest_path)
    if full_summary_hash != pilot.get("upstream_input_sha256", {}).get("run_summary.json"):
        raise ValueError("Original full-run summary differs from the one used by the pilot")
    side = manifest["side"]
    if side != pilot["library_side"]:
        raise ValueError("Pilot/full-run side mismatch")
    reference = Path(args.reference or pilot["reference"]["path"]).expanduser().resolve()
    outdir = (args.output_dir or root / "results" / (args.sample + "_junction_refine_v1_1")).expanduser().resolve()
    for inp in (pilot_dir, directory, reference):
        if outdir == inp or outdir in inp.parents:
            raise ValueError("Output cannot contain or replace existing inputs")
    if outdir.exists():
        previous = outdir / "refinement_summary.json"
        if not args.overwrite or not previous.is_file() or json.loads(previous.read_text()).get("analysis_type") != ANALYSIS_TYPE:
            raise ValueError("Output already exists; use a new --output-dir or --overwrite for this refinement's own outputs")
    paths = {stem: io.evidence_path(directory, stem) for stem in ("pair_classification", "read_evidence", "host_alignment_evidence")}
    input_hashes = {"run_summary.json": full_summary_hash}
    for stem, path in paths.items():
        input_hashes[path.name] = io.sha256_file(path)
        if input_hashes[path.name] != pilot.get("upstream_input_sha256", {}).get(path.name):
            raise ValueError("Full-run evidence changed since the pilot: " + str(path))
    selected, all_ids = {}, set()
    for row in io.tsv_rows(paths["pair_classification"], ("read_id", "classification", "host_mapping_supported")):
        rid = row["read_id"]
        if not rid or rid in all_ids:
            raise ValueError("Duplicate/empty original pair ID")
        all_ids.add(rid)
        if io.truth(row["host_mapping_supported"]):
            selected[rid] = row
    if set(selected) != set(old) or len(all_ids) != int(manifest["total_pairs"]):
        raise ValueError("Pilot/full-run host IDs or total denominator differ")
    read_evidence, hits = {}, collections.defaultdict(list)
    for row in io.tsv_rows(paths["read_evidence"], ("read_id", "mate", "raw_length", "cleaned_length", "cleaned_raw_start_1based", "cleaned_raw_end_1based")):
        if row["read_id"] in selected:
            key = (row["read_id"], row["mate"])
            if row["mate"] not in ("R1", "R2") or key in read_evidence:
                raise ValueError("Invalid/duplicate read-evidence mate")
            read_evidence[key] = row
    if set(read_evidence) != {(rid, mate) for rid in selected for mate in ("R1", "R2")}:
        raise ValueError("Selected host pairs lack matching read-evidence rows")
    for row in io.tsv_rows(paths["host_alignment_evidence"], ("read_id", "mate", "reference", "position_1based", "cigar", "strand", "accepted")):
        if row["read_id"] in selected and io.truth(row["accepted"]):
            if row["mate"] not in ("R1", "R2") or row["strand"] not in ("+", "-"):
                raise ValueError("Invalid accepted host alignment")
            row["intervals0"] = io.cigar_intervals(row["position_1based"], row["cigar"])
            hits[(row["read_id"], row["mate"])].append(row)
    if any(not hits.get((rid, "R1")) and not hits.get((rid, "R2")) for rid in selected):
        raise ValueError("Host-supported pair lacks accepted host evidence")
    print(f"Input {len(all_ids):,} pairs; reassessing all {len(selected):,} host-supported pairs.", flush=True)
    print("Verifying the existing combined reference (no BWA indexing) ...", flush=True)
    reference_hash = io.sha256_file(reference)
    if reference_hash != pilot["reference"]["sha256"]:
        raise ValueError("Combined reference sequence differs from the pilot reference")
    fasta = io.FastaReader(reference)
    donor_name = manifest["reference"].get("signature", {}).get("donor_contig", io.DONOR_CONTIG)
    donor = fasta.fetch(donor_name, 0, fasta.length(donor_name))
    if len(donor) != 800 or hashlib.sha256(donor.encode("ascii")).hexdigest() != pilot["reference"]["donor_sha256"]:
        raise ValueError("Donor reference differs from the pilot")
    clusters = templates_module.cluster_candidates(junction_rows, max_span=args.cluster_span)
    templates = templates_module.build_templates(clusters, donor, side, fasta.fetch, fasta.length, host_flank=args.host_flank)
    templates = assign_sequence_equivalence(templates, donor, side, fasta.fetch, fasta.length)
    by_id = {t["template_id"]: t for t in templates}
    if set(by_id) != {r["event_id"] for r in junction_rows}:
        raise ValueError("Template construction did not preserve every pilot junction")
    group_members = collections.defaultdict(list)
    for t in templates:
        group_members[t["equivalent_template_group"]].append(t["template_id"])
    index = template_index(templates)
    comparisons = templates_module.compare_templates(templates)
    print(f"Prepared {len(templates):,} candidate templates in {len(clusters):,} bounded review groups.", flush=True)
    lanes = []
    for lane in manifest["inputs"]["lanes"]:
        pair = []
        for mate in ("R1", "R2"):
            path = Path(lane[mate]["path"]).expanduser().resolve()
            if outdir == path or outdir in path.parents:
                raise ValueError("Output cannot contain raw input")
            if io.sha256_file(path) != lane[mate]["sha256"]:
                raise ValueError("Original FASTQ changed: " + str(path))
            pair.append(path)
        lanes.append(tuple(pair))
    outdir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="." + outdir.name + ".partial_", dir=outdir.parent))
    counts, transitions, changes = collections.Counter(), collections.Counter(), collections.Counter()
    group_support, group_review, group_mito = collections.Counter(), collections.Counter(), collections.Counter()
    group_strata = collections.defaultdict(collections.Counter)
    candidate_attempts, comparison_context_attempts, seen = 0, 0, set()
    digest = hashlib.sha256()
    best_examples = {}
    new_reasons = collections.Counter()
    with gzip.open(staging / "pair_refinement.tsv.gz", "wt", encoding="utf-8", newline="") as ph, gzip.open(staging / "read_realignment.tsv.gz", "wt", encoding="utf-8", newline="") as rh:
        pair_writer = csv.DictWriter(ph, PAIR_FIELDS, delimiter="\t"); pair_writer.writeheader()
        read_writer = csv.DictWriter(rh, READ_FIELDS, delimiter="\t"); read_writer.writeheader()
        for first, second in lanes:
            print("Reading paired FASTQ: " + first.name, flush=True)
            for raw1, raw2 in itertools.zip_longest(io.fastq_records(first), io.fastq_records(second)):
                if raw1 is None or raw2 is None:
                    raise ValueError("Different R1/R2 FASTQ counts")
                rid = io.normalized_id(raw1[0])
                if rid != io.normalized_id(raw2[0]) or rid in seen or rid not in all_ids:
                    raise ValueError("FASTQ pair IDs differ, are duplicate, or are absent from original evidence")
                seen.add(rid)
                io.update_stream_hash(digest, rid, raw1, raw2)
                if len(seen) % 100000 == 0:
                    print(f"Scanned {len(seen):,}/{len(all_ids):,} input pairs; reassessed {sum(counts.values()):,} host pairs ...", flush=True)
                if rid not in selected:
                    continue
                prior = old[rid]
                records, reasons, attempts = [], set(), 0
                for mate, raw in (("R1", raw1), ("R2", raw2)):
                    query, quality = io.cleaned_record(raw, read_evidence[(rid, mate)])
                    if not query or not hits.get((rid, mate)):
                        continue
                    for flip in (False, True):
                        q = exact.reverse_complement(query) if flip else query
                        qual = quality[::-1] if flip else quality
                        for context in accepted_host_contexts(hits[(rid, mate)], flip=flip):
                            context_templates = build_comparison_templates(
                                context, index, by_id, prior.get("event_id", ""), donor, side,
                                fasta.fetch, fasta.length, host_flank=args.host_flank)
                            if context_templates:
                                comparison_context_attempts += 1
                            for template in context_templates:
                                identifier = template["template_id"]
                                result = alignment_module.align_to_template(
                                    q, qual, template, donor, context["accepted_intervals0"])
                                attempts += 1; candidate_attempts += 1
                                comparison = template["comparison_context"]
                                provenance = {
                                    "comparison_context_id": context["comparison_context_id"],
                                    "comparison_common_host_endpoint0": template["common_host_endpoint0"],
                                    "comparison_common_donor_start0": template["donor_start0"],
                                    "effective_host_interval_start0": template["host_interval_start0"],
                                    "effective_host_interval_end0": template["host_interval_end0"]}
                                for aligned in result.get("alignments", []):
                                    aligned.update(provenance)
                                if result.get("best_alignment") is not None:
                                    result["best_alignment"].update(provenance)
                                for reason, value in result.get("diagnostics", {}).get("reason_counts", {}).items():
                                    if value:
                                        reasons.add(reason)
                                record = {"read_id": rid, "mate": mate,
                                    "query_orientation": context["query_orientation"],
                                    "template_id": identifier, "cluster_id": template["cluster_id"],
                                    "comparison_context_id": context["comparison_context_id"],
                                    "comparison_common_host_endpoint0": template["common_host_endpoint0"],
                                    "effective_host_interval_start0": template["host_interval_start0"],
                                    "effective_host_interval_end0": template["host_interval_end0"],
                                    "equivalent_template_group": template["equivalent_template_group"],
                                    "chrom": template["chrom"], "host_strand": template["host_strand"],
                                    "comparison_context_json": json_text(comparison), "result": result}
                                read_writer.writerow({**{k: v for k, v in record.items() if k != "result"},
                                    "status": result["status"], "score": result.get("score", ""),
                                    "best_alignment_json": json_text(result.get("best_alignment")) if result.get("best_alignment") is not None else "",
                                    "diagnostics_json": json_text(result.get("diagnostics", {}))})
                                if result["status"] != "rejected":
                                    records.append(record)
                decision = resolve_pair(records, group_members)
                if not attempts:
                    decision["status"] = "no_nearby_template"
                    reasons.add("no_existing_junction_template_near_accepted_host_alignment")
                status, group = decision["status"], decision["group"]
                if status not in STATUSES:
                    raise ValueError("Unexpected refinement status: " + status)
                if status == "template_supported":
                    if prior["status"] in ("unresolved", "ambiguous"):
                        change = "rescued_to_existing_template"
                    elif prior.get("event_id") in decision["template_ids"]:
                        change = "same_template_or_identical_sequence_representation"
                    else:
                        change = "different_existing_template"
                else:
                    change = "not_resolved_by_refinement"
                examples = decision["selected_results"]
                examples.sort(key=lambda r: (-(float(r["result"].get("score") or 0)), r["mate"], r["template_id"]))
                example = next((r for r in examples if r["result"].get("best_alignment") is not None), None)
                details = example["result"]["best_alignment"] if example else {}
                flags = sorted({flag for r in examples for flag in (r["result"].get("best_alignment") or {}).get("flags", [])})
                boundary_variant = any((r["result"].get("best_alignment") or {}).get("boundary_variant", False) for r in examples)
                prior_dd = io.truth(selected[rid].get("donor_rearrangement_supported", prior.get("prior_donor_rearrangement_supported", "0")))
                complex_host = io.truth(prior.get("complex_host_evidence", "0"))
                reference_compartment = ""
                if group:
                    selected_template = by_id[decision["template_ids"][0]]
                    reference_compartment = compartment(selected_template["chrom"])
                    boundary = int(selected_template["pilot_row"]["host_boundary0"])
                    for hit in hits.get((rid, "R1"), []) + hits.get((rid, "R2"), []):
                        distance = min(max(a - boundary, boundary - b, 0) for a, b in hit["intervals0"])
                        if hit["reference"] != selected_template["chrom"] or distance > 5000:
                            complex_host = True
                counts[status] += 1; transitions[(prior["status"], status)] += 1; changes[change] += 1
                for reason in reasons:
                    new_reasons[(status, reason)] += 1
                if status == "template_supported":
                    group_support[group] += 1
                    stratum = "complex_host" if complex_host else "prior_dd_without_complex_host" if prior_dd else "without_prior_dd_or_complex_host"
                    group_strata[group][stratum] += 1
                    group_mito[reference_compartment] += 1
                elif group:
                    group_review[group] += 1
                if example and group and group not in best_examples:
                    best_examples[group] = {**example, "pair_status": status}
                operations = details.get("operations", [])
                substitution_bases = sum(int(op.get("length", 1)) for op in operations if op.get("operation", op.get("op", op.get("type", ""))) in ("X", "SNV", "substitution"))
                pair_writer.writerow({"sample": args.sample, "library_side": side, "read_id": rid,
                    "old_status": prior["status"], "old_event_id": prior.get("event_id", ""),
                    "status": status, "change_type": change, "selected_equivalent_group": group,
                    "selected_template_ids": ";".join(decision["template_ids"]),
                    "candidate_group_ids_json": json_text(decision["candidate_groups"]),
                    "supporting_mates": ";".join(decision["supporting_mates"]),
                    "reference_compartment": reference_compartment,
                    "prior_donor_rearrangement_supported": int(prior_dd), "complex_host_evidence": int(complex_host),
                    "attempted_templates": attempts, "example_template_id": example["template_id"] if example else "",
                    "example_mate": example["mate"] if example else "", "example_substitution_bases": substitution_bases if example else "",
                    "example_gap_bases": details.get("gap_bases", ""), "boundary_variant": int(boundary_variant),
                    "flags_json": json_text(flags), "rejection_reasons_json": json_text(sorted(reasons))})
                if sum(counts.values()) % 500 == 0:
                    print(f"Reassessed {sum(counts.values()):,}/{len(selected):,} host pairs; supported {counts['template_supported']:,} ...", flush=True)
    fasta.close()
    if seen != all_ids or sum(counts.values()) != len(selected):
        raise ValueError("Final input/host pair reconciliation failed")
    if digest.hexdigest() != manifest["inputs"]["normalized_paired_stream_sha256"]:
        raise ValueError("Normalized raw read stream differs from the original full run")

    group_rows = []
    for cluster in clusters:
        members = cluster["rows"]
        group_rows.append({"sample": args.sample, "library_side": side, "cluster_id": cluster["cluster_id"],
            "chrom": cluster["chrom"], "host_strand": cluster["host_strand"],
            "boundary_min0": cluster["boundary_min0"], "boundary_max0": cluster["boundary_max0"],
            "span_bp": cluster["boundary_max0"] - cluster["boundary_min0"], "candidate_patterns": len(members),
            "original_support_pairs": sum(int(r["pair_support"]) for r in members),
            "template_ids": ";".join(sorted(r["event_id"] for r in members)),
            "reference_compartment": compartment(cluster["chrom"]),
            "interpretation": "bounded review group; not a validated insertion event"})
    write_table(staging / "boundary_groups.tsv", group_rows,
        ["sample", "library_side", "cluster_id", "chrom", "host_strand", "boundary_min0", "boundary_max0", "span_bp", "candidate_patterns", "original_support_pairs", "template_ids", "reference_compartment", "interpretation"])
    write_table(staging / "template_comparisons.tsv", comparisons, list(comparisons[0]) if comparisons else ["template_id", "relationship"])
    upstream_reasons = summarize_upstream_reasons(old_rows)
    write_table(staging / "upstream_unresolved_reasons.csv", upstream_reasons,
        ["old_status", "reason", "affected_pairs", "denominator_pairs", "fraction", "nonexclusive"], delimiter=",")
    write_table(staging / "refinement_reasons.csv",
        [{"status": status, "reason": reason, "affected_pairs": n, "denominator_pairs": counts[status],
          "fraction": n / counts[status], "nonexclusive": 1} for (status, reason), n in sorted(new_reasons.items())],
        ["status", "reason", "affected_pairs", "denominator_pairs", "fraction", "nonexclusive"], delimiter=",")
    reassessment = []
    eq_rows = []
    for group, identifiers in sorted(group_members.items()):
        representative = min(identifiers)
        t = by_id[representative]
        eq_rows.append({"equivalent_template_group": group, "representative_template_id": representative,
            "template_ids": ";".join(sorted(identifiers)),
            "review_cluster_ids": ";".join(sorted({by_id[identifier]["cluster_id"] for identifier in identifiers})),
            "equivalence_normalized_host_boundary0": t["equivalence_normalized_host_boundary0"],
            "equivalence_remaining_prefix": t["equivalence_remaining_prefix"],
            "reference_compartment": compartment(t["chrom"]),
            "supported_pairs": group_support[group], "review_pairs": group_review[group],
            **{k + "_pairs": group_strata[group][k] for k in ("without_prior_dd_or_complex_host", "prior_dd_without_complex_host", "complex_host")}})
        for identifier in sorted(identifiers):
            own_template = by_id[identifier]
            row = own_template["pilot_row"]
            reassessment.append({**row, "review_cluster_id": own_template["cluster_id"], "equivalent_template_group": group,
                "is_group_representative": int(identifier == representative),
                "equivalence_normalized_host_boundary0": own_template["equivalence_normalized_host_boundary0"],
                "equivalence_remaining_prefix": own_template["equivalence_remaining_prefix"],
                "equivalence_reference_suffix_bp": own_template["equivalence_reference_suffix_bp"],
                "supported_pair_count_for_identical_sequence_group": group_support[group],
                "supported_pair_count_once_per_group": group_support[group] if identifier == representative else 0,
                "review_pair_count_for_identical_sequence_group": group_review[group],
                "reference_compartment": compartment(t["chrom"]),
                "template_reference_mh_length_bp": row["mh_length_bp"],
                "template_reference_donor_deletion_min_bp": row["donor_deletion_min_bp"],
                "template_reference_donor_deletion_max_bp": row["donor_deletion_max_bp"],
                "interpretation": "Original candidate preserved; reference-template boundary values are not a new independent-event or consensus-MH call."})
    write_table(staging / "junction_reassessment.tsv", reassessment,
        list(reassessment[0]) if reassessment else ["event_id", "review_cluster_id", "supported_pair_count_once_per_group"])
    write_table(staging / "identical_sequence_groups.tsv", eq_rows,
        list(eq_rows[0]) if eq_rows else ["equivalent_template_group", "supported_pairs"])
    with (staging / "representative_realignments.txt").open("w", encoding="utf-8") as handle:
        handle.write("Small-variant realignment to pre-existing junction templates.\n"
                     "Templates are candidate explanations, not certified insertion events.\n"
                     "Minor variant templates were not abundance-corrected.\n\n")
        ordered_examples = sorted(best_examples.items(), key=lambda kv: (-group_support[kv[0]], kv[0]))[:30]
        for group, record in ordered_examples:
            result = record["result"]; d = result["best_alignment"]
            handle.write(f"{group} supported_pairs={group_support[group]} review_pairs={group_review[group]}\n")
            handle.write(f"Read {record['read_id']} {record['mate']} {record['query_orientation']} template={record['template_id']} status={result['status']}\n")
            handle.write("Comparison context: " + record["comparison_context_json"] + "\n")
            handle.write("TEMPLATE " + d.get("aligned_template", "") + "\n")
            handle.write("QUERY    " + d.get("aligned_query", "") + "\n")
            handle.write("Operations: " + json_text(d.get("operations", [])) + "\n")
            handle.write("Flags: " + json_text(d.get("flags", [])) + "\n\n")
    old_counts = collections.Counter(r["status"] for r in old_rows)
    stats = {"total_input_pairs": len(all_ids), "host_supported_pairs": len(selected),
             "upstream_called_pairs": int(pilot["stats"]["called_pairs_total"]),
             "upstream_unresolved_pairs": old_counts["unresolved"], "upstream_ambiguous_pairs": old_counts["ambiguous"],
             **{status: counts[status] for status in STATUSES},
             "rescued_to_existing_template_pairs": changes["rescued_to_existing_template"],
             "different_existing_template_pairs": changes["different_existing_template"],
             "same_template_or_identical_sequence_pairs": changes["same_template_or_identical_sequence_representation"],
             "original_candidate_patterns": len(junction_rows), "bounded_review_groups": len(clusters),
             "identical_sequence_groups": len(group_members),
             "supported_identical_sequence_groups": sum(n > 0 for n in group_support.values()),
             "candidate_read_template_alignments_attempted": candidate_attempts,
             "accepted_host_comparison_contexts_with_templates": comparison_context_attempts,
             "supported_fraction_of_host_pairs": counts["template_supported"] / len(selected) if selected else None,
             "supported_pairs_by_reference_compartment": dict(group_mito)}
    import Bio
    summary = {"analysis_type": ANALYSIS_TYPE, "analysis_version": VERSION, "status": "complete",
        "sample": args.sample, "library_side": side, "condition": pilot.get("condition", ""),
        "biological_replicate": pilot.get("biological_replicate", ""), "completed_utc": utc_now(),
        "elapsed_seconds": time.monotonic() - started, "stats": stats,
        "old_to_new_status_counts": [{"old_status": a, "new_status": b, "pairs": n} for (a, b), n in sorted(transitions.items())],
        "config": {"cluster_total_span_bp": args.cluster_span, "common_host_flank_bp": args.host_flank,
                   "accepted_host_template_search_padding_bp": 20, "complex_host_distance_bp": 5000,
                   "comparison_context_source": "overlapping accepted same-mate host-alignment spans, same contig/oriented strand/query orientation",
                   "comparison_context_merge_padding_bp": 0,
                   "score_competition_unit": "mate and accepted-host comparison context; independent of display cluster",
                   "shared_comparison_endpoints": "fixed primer-distal donor start and common distal host endpoint for every template in a comparison context",
                   "sequence_equivalence": "reference-suffix-normalized full reconstructed allele; independent of review cluster and comparison endpoint; no variant correction",
                   "old_template_retention_min_accepted_host_overlap_bp": 30,
                   "alignment": alignment_module.DEFAULT_CONFIG,
                   "pair_merging": False, "abundance_based_correction": False, "novel_template_discovery": False},
        "pilot_directory": str(pilot_dir), "pilot_input_sha256": pilot_hashes,
        "fullrun_directory": str(directory), "fullrun_input_sha256": input_hashes,
        "reference": {"path": str(reference), "sha256": reference_hash, "donor_sha256": pilot["reference"]["donor_sha256"]},
        "raw_inputs": manifest["inputs"], "python_version": sys.version, "biopython_version": Bio.__version__,
        "source_sha256": {p.relative_to(HERE).as_posix(): io.sha256_file(p) for p in sorted(HERE.glob("*.py")) + sorted((HERE / "support").glob("*.py"))},
        "validation": {"all_input_ids_reconciled": seen == all_ids,
                       "all_host_pairs_have_one_status": sum(counts.values()) == len(selected),
                       "supported_group_pair_counts_reconcile": sum(group_support.values()) == counts["template_supported"],
                       "status_transitions_reconcile": sum(transitions.values()) == len(selected),
                       "original_candidates_retained": len(reassessment) == len(junction_rows)},
        "limitations": ["Review groups are spatial neighborhoods, not validated independent integrations.",
            "Only previously observed candidate templates are tested; absence of a nearby template cannot exclude or discover a new junction.",
            "A higher-scoring minor-variant template is retained; PCR errors and true rare alleles cannot be distinguished from sequence alone.",
            "Template reference MH/deletion values are inherited descriptions, not newly inferred consensus or human-side deletion calls.",
            "Junction-proximal variants and alignment ambiguity remain review categories.",
            "Reason counts are nonexclusive per-pair flags across tested orientations; they are not exclusive biological causes.",
            "Mitochondrial-reference mappings are a separate compartment, not proof of mitochondrial integration."]}
    if not all(summary["validation"].values()):
        raise ValueError("Internal refinement count validation failed")
    summary["output_sha256"] = {p.name: io.sha256_file(p) for p in staging.iterdir() if p.is_file()}
    (staging / "refinement_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if outdir.exists():
        for path in staging.iterdir():
            os.replace(path, outdir / path.name)
        staging.rmdir()
    else:
        staging.rename(outdir)
    print(f"COMPLETE: {counts['template_supported']:,}/{len(selected):,} pairs supported an existing template; {changes['rescued_to_existing_template']:,} rescued from old unresolved/ambiguous.", flush=True)
    print("Results: " + str(outdir), flush=True)
    return summary


def main(argv=None):
    try:
        run(parse_args(argv))
    except (ValueError, KeyError, OSError) as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
