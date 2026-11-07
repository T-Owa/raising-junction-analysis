"""Small-variant-tolerant validation against *existing* junction templates.

This module does not discover new integration structures or decide whether a
minor allele is a PCR error.  Template-derived MH/deletion descriptions remain
explicitly separate from observed query-to-template substitutions and indels.
Reads with a variant near the junction are retained for review rather than
assigned the template's MH as a newly measured value.

Queries and templates are oriented donor -> host. Coordinates are zero-based,
half-open, and genomic boundary coordinates are between bases. Biopython's
local PairwiseAligner uses explicit scores, independent of version defaults.
"""
from collections import Counter
from itertools import islice

from Bio import Align
from support import junction_core as exact


DEFAULT_CONFIG = {
    "match_score": 2.0,
    "mismatch_score": -4.0,
    "open_gap_score": -6.0,
    "extend_gap_score": -1.0,
    "min_identity_percent": 95.0,
    "max_edit_bases": 3,
    "max_gap_bases": 3,
    "min_donor_anchor": 20,
    "min_host_anchor": 30,
    "min_global_overlap": 30,
    "min_quality": 20,
    "junction_flank": 10,
    "boundary_guard": 3,
    "donor_mask_seed": 20,
    "homopolymer_min": 10,
    "max_optimal_alignments": 16,
}
_ACGT = set("ACGT")


def _configuration(config):
    values = dict(DEFAULT_CONFIG)
    if config:
        unknown = set(config) - set(values)
        if unknown:
            raise ValueError("Unknown alignment configuration: " + ", ".join(sorted(unknown)))
        values.update(config)
    for name in ("max_edit_bases", "max_gap_bases", "min_donor_anchor", "min_host_anchor",
                 "min_global_overlap", "min_quality", "junction_flank", "boundary_guard",
                 "donor_mask_seed", "homopolymer_min", "max_optimal_alignments"):
        if not isinstance(values[name], int) or isinstance(values[name], bool) or values[name] < 0:
            raise ValueError("Invalid nonnegative integer configuration: " + name)
    for name in ("min_donor_anchor", "min_host_anchor", "min_global_overlap",
                 "donor_mask_seed", "homopolymer_min", "max_optimal_alignments"):
        if values[name] < 1:
            raise ValueError("Configuration must be positive: " + name)
    if not 0 < float(values["min_identity_percent"]) <= 100:
        raise ValueError("min_identity_percent must lie in (0,100]")
    if float(values["match_score"]) <= 0 or any(float(values[k]) >= 0 for k in
                ("mismatch_score", "open_gap_score", "extend_gap_score")):
        raise ValueError("Match score must be positive and mismatch/gap scores negative")
    return values


def _quality_scores(quality, length):
    scores = [ord(x) - 33 for x in quality] if isinstance(quality, str) else [int(x) for x in quality]
    if len(scores) != length or any(x < 0 or x > 93 for x in scores):
        raise ValueError("Query and Phred+33 qualities must have matching lengths")
    return scores


def _accepted_intervals(intervals):
    if intervals is None:
        return []
    if len(intervals) == 2 and isinstance(intervals[0], (int, float)):
        intervals = [intervals]
    result = []
    for a, b in intervals:
        a, b = int(a), int(b)
        if a < 0 or b <= a:
            raise ValueError("Accepted genomic intervals must be nonempty [start,end)")
        result.append((a, b))
    return result


def _region(start, end, donor_end, host_start):
    """Name a target span; a point exactly at a join is junction-associated."""
    if start == end:
        if start in (donor_end, host_start):
            return "junction"
        return "donor" if start < donor_end else "intervening" if start < host_start else "host"
    if end <= donor_end:
        return "donor"
    if start >= host_start:
        return "host"
    if start >= donor_end and end <= host_start:
        return "intervening"
    return "junction"


def _touches(start, end, low, high):
    return low <= start <= high if start == end else start <= high and end > low


def _donor_mask(query, donor, config):
    mask = set()
    for reference in (donor, exact.reverse_complement(donor)):
        for match in exact.maximal_exact_matches(query, reference, config["donor_mask_seed"]):
            mask.update(range(match["query_start0"], match["query_end0"]))
    return mask


def _parse_alignment(alignment, query, quality, template, accepted, donor_mask, poly_mask, cfg):
    target = template["sequence"]
    donor_end, host_start = int(template["donor_end_t0"]), int(template["host_start_t0"])
    row = template["pilot_row"]
    mh = int(row["mh_length_bp"])
    coords = alignment.coordinates.tolist()
    t_start, t_end = int(coords[0][0]), int(coords[0][-1])
    q_start, q_end = int(coords[1][0]), int(coords[1][-1])
    operations, matches, target_to_query = [], [], {}
    query_insertions, t_text, q_text = [], [], []
    mismatch_bases = gap_bases = max_gap = columns = 0

    for i in range(len(coords[0]) - 1):
        a, b = int(coords[0][i]), int(coords[0][i + 1])
        u, v = int(coords[1][i]), int(coords[1][i + 1])
        dt, dq = b - a, v - u
        if dt > 0 and dq > 0:
            if dt != dq:
                raise ValueError("Unexpected PairwiseAligner diagonal segment")
            t_text.append(target[a:b]); q_text.append(query[u:v]); columns += dt
            for k in range(dt):
                t, q = a + k, u + k
                target_to_query[t] = q
                if target[t] == query[q] and query[q] in _ACGT:
                    matches.append((t, q))
                else:
                    mismatch_bases += 1
                    operations.append({"operation": "SNV", "target_start0": t,
                        "target_end0": t + 1, "query_start0": q, "query_end0": q + 1,
                        "template_bases": target[t], "query_bases": query[q],
                        "region": _region(t, t + 1, donor_end, host_start),
                        "min_query_quality": quality[q]})
        elif dq > 0 and dt == 0:
            t_text.append("-" * dq); q_text.append(query[u:v]); columns += dq
            gap_bases += dq; max_gap = max(max_gap, dq)
            query_insertions.append((a, u, v))
            operations.append({"operation": "I", "target_start0": a,
                "target_end0": a, "query_start0": u, "query_end0": v,
                "template_bases": "", "query_bases": query[u:v],
                "region": _region(a, a, donor_end, host_start),
                "min_query_quality": min(quality[u:v])})
        elif dt > 0 and dq == 0:
            t_text.append(target[a:b]); q_text.append("-" * dt); columns += dt
            gap_bases += dt; max_gap = max(max_gap, dt)
            operations.append({"operation": "D", "target_start0": a,
                "target_end0": b, "query_start0": u, "query_end0": u,
                "template_bases": target[a:b], "query_bases": "",
                "region": _region(a, b, donor_end, host_start),
                "min_query_quality": None})
        else:
            raise ValueError("Invalid empty/reversed alignment segment")

    identity = 100.0 * len(matches) / columns if columns else 0.0
    reasons = []
    if mismatch_bases + gap_bases > cfg["max_edit_bases"]:
        reasons.append("too_many_edit_bases")
    if max_gap > cfg["max_gap_bases"]:
        reasons.append("gap_too_long")
    if identity < cfg["min_identity_percent"]:
        reasons.append("identity_below_threshold")

    # Canonical templates assign shared MH to donor. Those bases do not count
    # as independent donor anchors; the canonical human arm begins after them.
    d_positions = [(t, q) for t, q in matches if t < donor_end - mh and quality[q] >= cfg["min_quality"]]
    h_positions = [(t, q) for t, q in matches if t >= host_start
                   and quality[q] >= cfg["min_quality"] and q not in donor_mask and q not in poly_mask]
    if len(d_positions) < cfg["min_donor_anchor"]:
        reasons.append("insufficient_independent_donor_anchor")
    if len(h_positions) < cfg["min_host_anchor"]:
        reasons.append("insufficient_independent_host_anchor")

    def genomic_base(t):
        offset = t - host_start
        if template["host_strand"] == "+":
            return int(template["host_interval_start0"]) + offset
        return int(template["host_interval_end0"]) - 1 - offset

    global_overlap = sum(any(a <= genomic_base(t) < b for a, b in accepted) for t, q in h_positions)
    if global_overlap < cfg["min_global_overlap"]:
        reasons.append("insufficient_accepted_global_overlap")

    # Read sequence must bridge the entire shared-MH/intervening region and
    # both flanks. A local host-only alignment cannot pass this criterion.
    required_start = donor_end - mh - cfg["junction_flank"]
    required_end = host_start + cfg["junction_flank"]
    if required_start < 0 or required_end > len(target) or t_start > required_start or t_end < required_end:
        reasons.append("junction_and_flanks_not_spanned")
    observed = [q for t, q in target_to_query.items() if required_start <= t < required_end]
    for t, u, v in query_insertions:
        if required_start <= t <= required_end:
            observed.extend(range(u, v))
    if not observed:
        reasons.append("junction_not_observed")
    else:
        lo, hi = min(observed), max(observed) + 1
        if any(quality[q] < cfg["min_quality"] or query[q] not in _ACGT for q in range(lo, hi)):
            reasons.append("junction_quality_below_threshold")

    guard = cfg["boundary_guard"]
    review_low, review_high = donor_end - mh - guard, host_start + guard
    boundary_ops = [op for op in operations if _touches(op["target_start0"], op["target_end0"], review_low, review_high)]
    flags = []
    if boundary_ops:
        flags.append("variant_at_or_near_junction_do_not_transfer_template_MH")
    if operations:
        flags.append("query_differs_from_existing_template")
    if q_start:
        flags.append("unmodeled_query_prefix")
    if q_end < len(query):
        flags.append("unmodeled_query_suffix")
    if any(op["region"] == "donor" and op["operation"] in ("I", "D") for op in operations):
        flags.append("donor_arm_indel_not_terminal_deletion_measurement")
    if int(row["insertion_length_bp"]) > 0:
        flags.append("intervening_sequence_interfaces_MH_not_assessed")

    # Distinct optimal paths are equivalent for the local junction only when
    # their query projection around the boundary and local edits agree.
    projection = [[t, target_to_query.get(t)] for t in range(max(0, review_low), min(len(target), review_high + 1))]
    local_ops = [(op["operation"], op["target_start0"], op["target_end0"],
                  op["query_start0"], op["query_end0"]) for op in boundary_ops]
    signature = (tuple((a, b) for a, b in projection), tuple(local_ops))
    data = {
        "template_id": template["template_id"], "cluster_id": template["cluster_id"],
        "chrom": template["chrom"], "host_strand": template["host_strand"],
        "score": float(alignment.score), "identity_percent": round(identity, 6),
        "edit_bases": mismatch_bases + gap_bases, "substitution_bases": mismatch_bases,
        "gap_bases": gap_bases, "max_gap": max_gap,
        "donor_independent_bases": len(d_positions), "host_independent_bases": len(h_positions),
        "host_global_overlap_bases": global_overlap,
        "alignment_target_start0": t_start, "alignment_target_end0": t_end,
        "alignment_query_start0": q_start, "alignment_query_end0": q_end,
        "unmodeled_query_prefix_bp": q_start, "unmodeled_query_suffix_bp": len(query) - q_end,
        "operations": operations, "flags": flags, "boundary_variant": bool(boundary_ops),
        "template_reference_mh_length_bp": mh,
        "template_reference_mh_sequence": row.get("mh_sequence", ""),
        "template_reference_donor_deletion_min_bp": int(row["donor_deletion_min_bp"]),
        "template_reference_donor_deletion_max_bp": int(row["donor_deletion_max_bp"]),
        "template_reference_intervening_length_bp": int(row["insertion_length_bp"]),
        "template_reference_intervening_sequence": row.get("insertion_sequence", ""),
        "junction_projection": projection,
        "aligned_template": "".join(t_text), "aligned_query": "".join(q_text),
    }
    return data, reasons, signature


def align_to_template(query, quality, template, donor, accepted_intervals0, config=None):
    """Re-align one cleaned, D->H-oriented read to one existing template.

    A ``supported`` result supports compatibility with this existing template;
    it does not validate the integration event, correct a minor allele, or
    infer an original molecule. ``review`` retains adequate anchors but has an
    observed small variant at/near the junction. Such a result must not inherit
    the template's MH as a newly measured value. ``ambiguous_alignment`` means
    optimal alignment paths leave different junction projections or the
    conservative path-enumeration cap was exceeded.
    """
    cfg = _configuration(config)
    query, donor = str(query).upper(), str(donor).upper()
    quality = _quality_scores(quality, len(query))
    accepted = _accepted_intervals(accepted_intervals0)
    sequence = str(template["sequence"]).upper()
    if sequence != template["sequence"]:
        template = {**template, "sequence": sequence}
    d_end, h_start = int(template["donor_end_t0"]), int(template["host_start_t0"])
    mh = int(template["pilot_row"]["mh_length_bp"])
    if not 0 < d_end <= h_start < len(sequence) or not 0 <= mh <= d_end:
        raise ValueError("Invalid template arm/MH coordinates")
    if template["host_strand"] not in ("+", "-"):
        raise ValueError("Template host_strand must be + or -")
    if int(template["host_interval_end0"]) - int(template["host_interval_start0"]) != len(sequence) - h_start:
        raise ValueError("Host template interval length mismatch")
    diagnostics = {"reason_counts": {}, "optimal_paths_examined": 0,
                   "truncated": False, "distinct_junction_projections": 0}

    def early(reason):
        diagnostics["reason_counts"][reason] = 1
        return {"status": "rejected", "score": None, "best_alignment": None,
                "alignments": [], "diagnostics": diagnostics}

    if not accepted:
        return early("no_accepted_global_host_interval")
    if not query:
        return early("empty_cleaned_query")
    donor_mask = _donor_mask(query, donor, cfg)
    poly_mask = exact._homopolymer_mask(query, cfg["homopolymer_min"])
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = cfg["match_score"]
    aligner.mismatch_score = cfg["mismatch_score"]
    aligner.open_gap_score = cfg["open_gap_score"]
    aligner.extend_gap_score = cfg["extend_gap_score"]
    optimal = aligner.align(sequence, query)
    score = float(optimal.score)
    seen, rejections = {}, Counter()
    for index, alignment in enumerate(islice(optimal, cfg["max_optimal_alignments"] + 1)):
        if index == cfg["max_optimal_alignments"]:
            diagnostics["truncated"] = True
            rejections["optimal_alignment_path_limit"] += 1
            break
        diagnostics["optimal_paths_examined"] += 1
        candidate, reasons, signature = _parse_alignment(
            alignment, query, quality, template, accepted, donor_mask, poly_mask, cfg)
        if reasons:
            rejections.update(reasons)
            continue
        if signature not in seen:
            candidate["equivalent_optimal_paths"] = 1
            seen[signature] = candidate
        else:
            seen[signature]["equivalent_optimal_paths"] += 1
            if "distal_optimal_path_ambiguity" not in seen[signature]["flags"]:
                seen[signature]["flags"].append("distal_optimal_path_ambiguity")
    if not diagnostics["optimal_paths_examined"]:
        rejections["no_positive_score_alignment"] += 1
    diagnostics["reason_counts"] = dict(sorted(rejections.items()))
    diagnostics["distinct_junction_projections"] = len(seen)
    passing = sorted(seen.values(), key=lambda c: (c["boundary_variant"], c["edit_bases"],
                         -min(c["donor_independent_bases"], c["host_independent_bases"]),
                         c["alignment_query_start0"], c["alignment_target_start0"]))
    best = passing[0] if passing else None
    if diagnostics["truncated"] or len(passing) > 1:
        status = "ambiguous_alignment"
    elif not passing:
        status = "rejected"
    else:
        status = "review" if best["boundary_variant"] else "supported"
    return {"status": status, "score": score, "best_alignment": best,
            "alignments": passing, "diagnostics": diagnostics}
