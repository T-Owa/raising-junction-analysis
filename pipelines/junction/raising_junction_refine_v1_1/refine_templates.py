"""Review groups and candidate templates for RAISING junction refinement.

Nearby boundaries are grouped for comparison, never asserted to be one
integration event. Only identical sequences over shared comparison endpoints
receive the same template-equivalence identifier. Small substitutions or indels
remain distinct alternatives, irrespective of their abundance.
"""
import collections
import hashlib
import json


NUMERIC_FIELDS = (
    "host_boundary0", "host_boundary_min0", "host_boundary_max0",
    "donor_cut_oriented_min0", "donor_cut_oriented_max0",
    "donor_cut_original_min0", "donor_cut_original_max0",
    "donor_deletion_min_bp", "donor_deletion_max_bp", "mh_length_bp",
    "insertion_length_bp", "pair_support", "simple_pair_support",
    "prior_dd_pair_support", "complex_host_pair_support",
)


def _identifier(prefix, value):
    text = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return prefix + hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def normalize_row(row):
    """Copy a pilot row and normalize known numeric fields without mutation."""
    result = dict(row)
    for name in NUMERIC_FIELDS:
        if name in result and result[name] not in (None, ""):
            result[name] = int(result[name])
            if result[name] < 0:
                raise ValueError("Negative pilot coordinate/count: " + name)
    for name in ("event_id", "chrom", "host_strand", "host_boundary0",
                 "donor_cut_oriented_max0"):
        if name not in result or result[name] in (None, ""):
            raise ValueError("Missing pilot field: " + name)
    if result["host_strand"] not in ("+", "-"):
        raise ValueError("Invalid host strand")
    result["chrom"] = str(result["chrom"])
    result["insertion_sequence"] = str(result.get("insertion_sequence", "")).upper()
    if set(result["insertion_sequence"]) - set("ACGTN"):
        raise ValueError("Invalid intervening sequence")
    if ("insertion_length_bp" in result and
            result["insertion_length_bp"] != len(result["insertion_sequence"])):
        raise ValueError("Intervening sequence length mismatch")
    result.setdefault("pair_support", 0)
    return result


def cluster_candidates(rows, max_span=10):
    """Bounded-span review groups on the same contig and host orientation.

    Each group's maximum boundary minus minimum boundary is <=max_span.
    Groups are built in ascending boundary order; there is no single-linkage
    chaining beyond that total span. Boundaries are 0-based between bases.
    """
    if not isinstance(max_span, int) or max_span < 0:
        raise ValueError("max_span must be a nonnegative integer")
    groups = collections.defaultdict(list)
    seen = set()
    for raw in rows:
        row = normalize_row(raw)
        if row["event_id"] in seen:
            raise ValueError("Duplicate event_id: " + row["event_id"])
        seen.add(row["event_id"])
        groups[(row["chrom"], row["host_strand"])].append(row)
    result = []
    for (chrom, strand), members in sorted(groups.items()):
        members.sort(key=lambda r: (r["host_boundary0"], r["event_id"]))
        chunks = []
        for row in members:
            if not chunks or row["host_boundary0"] - chunks[-1][0]["host_boundary0"] > max_span:
                chunks.append([])
            chunks[-1].append(row)
        for chunk in chunks:
            low, high = chunk[0]["host_boundary0"], chunk[-1]["host_boundary0"]
            result.append({
                "cluster_id": _identifier("C", [chrom, strand, low, high,
                                                  [r["event_id"] for r in chunk]]),
                "chrom": chrom, "host_strand": strand,
                "boundary_min0": low, "boundary_max0": high,
                "rows": chunk,
            })
    return result


def _primer_info(donor, side):
    from support.junction_core import primer_info
    return primer_info(donor, side)


def _revcomp(sequence):
    from support.junction_core import reverse_complement
    return reverse_complement(sequence)


def build_templates(clusters, donor, side, fetch, contig_length, host_flank=120):
    """Build D->H templates over shared endpoints within each review group.

    fetch(chrom,start0,end0) returns forward-reference sequence. The same
    primer-distal donor start and distal human endpoint are used for every
    template in a group. Host intervals retain forward-genome coordinates even
    when their sequence is reverse-complemented for a minus-strand junction.
    """
    if not isinstance(host_flank, int) or host_flank < 30:
        raise ValueError("host_flank must be at least 30 bp")
    info = _primer_info(donor, side)
    oriented = info["oriented_sequence"]
    donor_start = info["primer_end0"]
    result = []
    for group in clusters:
        chrom, strand = group["chrom"], group["host_strand"]
        length = int(contig_length(chrom))
        if not 0 <= group["boundary_min0"] <= group["boundary_max0"] <= length:
            raise ValueError("Candidate boundary outside reference contig")
        endpoint = (min(length, group["boundary_max0"] + host_flank)
                    if strand == "+" else max(0, group["boundary_min0"] - host_flank))
        for raw in group["rows"]:
            row = normalize_row(raw)
            boundary = row["host_boundary0"]
            cut = row["donor_cut_oriented_max0"]
            if not donor_start <= cut <= len(oriented):
                raise ValueError("Candidate donor cut is outside primer-distal reference")
            human_start, human_end = ((boundary, endpoint) if strand == "+" else (endpoint, boundary))
            human = str(fetch(chrom, human_start, human_end)).upper()
            if len(human) != human_end - human_start or set(human) - set("ACGTN"):
                raise ValueError("Reference fetch has invalid length/bases")
            if strand == "-":
                human = _revcomp(human)
            d = oriented[donor_start:cut]
            ins = row["insertion_sequence"]
            sequence = d + ins + human
            template = {
                "template_id": row["event_id"], "cluster_id": group["cluster_id"],
                "chrom": chrom, "host_strand": strand, "sequence": sequence,
                "donor_start0": donor_start, "donor_cut_oriented0": cut,
                "donor_end_t0": len(d), "lenD": len(d),
                "host_start_t0": len(d) + len(ins),
                "host_interval_start0": human_start, "host_interval_end0": human_end,
                "common_host_endpoint0": endpoint,
                "pilot_row": row,
                "equivalent_template_group": _identifier("E", [group["cluster_id"], sequence]),
            }
            result.append(template)
    equivalences = collections.defaultdict(list)
    for template in result:
        equivalences[template["equivalent_template_group"]].append(template["template_id"])
    for template in result:
        template["equivalent_template_members"] = sorted(equivalences[template["equivalent_template_group"]])
    return result


def _global_affine(reference, query):
    """Global affine comparison; report one optimal path and possible ties.

    The comparison is descriptive. It does not assign error probabilities or
    infer whether a minor allele is a PCR artifact. Terminal gaps are penalized
    like internal gaps so that true terminal-donor deletions remain visible.
    """
    if reference == query:
        return {"score": 2 * len(reference), "alignment_ambiguous": False,
                "alignment_edit_bases": 0, "substitutions": 0, "indel_bases": 0,
                "gap_events": 0, "operations": []}
    if not reference or not query:
        operation = {"type": "insertion" if query else "deletion",
                     "reference_start0": 0, "reference_end0": len(reference),
                     "query_start0": 0, "query_end0": len(query),
                     "reference_bases": reference, "query_bases": query}
        n = max(len(reference), len(query))
        return {"score": -6 - (n - 1), "alignment_ambiguous": False,
                "alignment_edit_bases": n, "substitutions": 0, "indel_bases": n,
                "gap_events": 1, "operations": [operation]}
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -4
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -1
    alignments = iter(aligner.align(reference, query))
    alignment = next(alignments)
    ambiguous = next(alignments, None) is not None
    coordinates = alignment.coordinates
    operations = []
    for i in range(coordinates.shape[1] - 1):
        a, b = int(coordinates[0, i]), int(coordinates[0, i + 1])
        c, d = int(coordinates[1, i]), int(coordinates[1, i + 1])
        if b > a and d > c:
            if b - a != d - c:
                raise ValueError("Unexpected non-unit diagonal in pairwise alignment")
            for offset in range(b - a):
                if reference[a + offset] != query[c + offset]:
                    operations.append({"type": "substitution",
                        "reference_start0": a + offset, "reference_end0": a + offset + 1,
                        "query_start0": c + offset, "query_end0": c + offset + 1,
                        "reference_bases": reference[a + offset], "query_bases": query[c + offset]})
        else:
            operations.append({"type": "deletion" if b > a else "insertion",
                "reference_start0": a, "reference_end0": b,
                "query_start0": c, "query_end0": d,
                "reference_bases": reference[a:b], "query_bases": query[c:d]})
    substitutions = sum(op["type"] == "substitution" for op in operations)
    indel_bases = sum(max(len(op["reference_bases"]), len(op["query_bases"]))
                      for op in operations if op["type"] != "substitution")
    return {"score": float(alignment.score), "alignment_ambiguous": ambiguous,
            "alignment_edit_bases": substitutions + indel_bases,
            "substitutions": substitutions, "indel_bases": indel_bases,
            "gap_events": sum(op["type"] != "substitution" for op in operations),
            "operations": operations}


def _region(operation, reference_template, guard=3):
    a, b = operation["reference_start0"], operation["reference_end0"]
    donor_end, host_start = reference_template["donor_end_t0"], reference_template["host_start_t0"]
    mh = int(reference_template["pilot_row"].get("mh_length_bp", 0))
    low, high = donor_end - mh - guard, host_start + guard
    if (a == b and low <= a <= high) or (a < high and b > low):
        return "junction_proximal"
    if b <= donor_end:
        return "donor_internal"
    if a >= host_start:
        return "human_internal"
    return "intervening_sequence"


def _review_row(template, reference, max_edit_distance, limited_host_context=False):
    comparison = _global_affine(reference["sequence"], template["sequence"])
    operations = comparison.pop("operations")
    for op in operations:
        op["reference_region"] = _region(op, reference)
    if template["template_id"] == reference["template_id"]:
        relationship = "strongest_template"
    elif not comparison["alignment_edit_bases"]:
        relationship = "sequence_equivalent"
    elif comparison["alignment_edit_bases"] <= max_edit_distance:
        relationship = ("substitution_only_alternative" if not comparison["indel_bases"]
                        else "small_indel_alternative")
    else:
        relationship = "sequence_different"
    ref_support = int(reference["pilot_row"].get("pair_support", 0))
    support = int(template["pilot_row"].get("pair_support", 0))
    return {
        "cluster_id": template["cluster_id"], "chrom": template["chrom"],
        "host_strand": template["host_strand"],
        "template_id": template["template_id"], "reference_template_id": reference["template_id"],
        "host_boundary0": template["pilot_row"]["host_boundary0"],
        "reference_host_boundary0": reference["pilot_row"]["host_boundary0"],
        "pair_support": support, "reference_pair_support": ref_support,
        "reference_to_candidate_support_ratio": ref_support / support if support else None,
        "relationship": relationship,
        "limited_host_context": limited_host_context,
        "automatic_merge": False,
        "has_junction_proximal_difference": any(op["reference_region"] == "junction_proximal" for op in operations),
        "differences_only_donor_internal_substitutions": bool(operations) and all(
            op["type"] == "substitution" and op["reference_region"] == "donor_internal" for op in operations),
        **comparison,
        "operations_json": json.dumps(operations, sort_keys=True, separators=(",", ":")),
    }


def compare_templates(templates, max_edit_distance=3):
    """Compare every candidate with the strongest template in its review group.

    Abundance chooses a displayed comparison reference only. It does not select
    the true biological allele, correct reads, or merge distinct sequences.
    """
    if not isinstance(max_edit_distance, int) or max_edit_distance < 0:
        raise ValueError("max_edit_distance must be a nonnegative integer")
    groups = collections.defaultdict(list)
    for template in templates:
        groups[template["cluster_id"]].append(template)
    result = []
    for _, members in sorted(groups.items()):
        reference = min(members, key=lambda t: (-int(t["pilot_row"].get("pair_support", 0)), t["template_id"]))
        for template in sorted(members, key=lambda t: t["template_id"]):
            result.append(_review_row(template, reference, max_edit_distance))
    return result


def build_donor_only_review(rows, donor, side, max_edit_distance=3):
    """Audit donor-side alleles at the EXACT same host boundary without FASTA.

    This is intentionally a partial, attachment-only audit: no nearby-boundary
    comparison, host sequence verification, quality analysis or raw-read
    reassignment is performed. Every distinct original pattern is retained.
    Only exact-boundary groups containing multiple patterns are returned.
    """
    info = _primer_info(donor, side)
    oriented, start = info["oriented_sequence"], info["primer_end0"]
    result = []
    for group in cluster_candidates(rows, max_span=0):
        if len(group["rows"]) < 2:
            continue
        templates = []
        for row in group["rows"]:
            cut = row["donor_cut_oriented_max0"]
            if not start <= cut <= len(oriented):
                raise ValueError("Candidate donor cut outside primer-distal reference")
            d = oriented[start:cut]
            sequence = d + row["insertion_sequence"]
            templates.append({"template_id": row["event_id"], "cluster_id": group["cluster_id"],
                "chrom": group["chrom"], "host_strand": group["host_strand"],
                "sequence": sequence, "donor_end_t0": len(d), "host_start_t0": len(sequence),
                "pilot_row": row})
        reference = min(templates, key=lambda t: (-t["pilot_row"]["pair_support"], t["template_id"]))
        for template in sorted(templates, key=lambda t: t["template_id"]):
            record = _review_row(template, reference, max_edit_distance, limited_host_context=True)
            record.update(candidate_donor_side_sequence=template["sequence"],
                          reference_donor_side_sequence=reference["sequence"])
            result.append(record)
    return result
