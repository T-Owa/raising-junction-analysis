"""Conservative exact-anchor pilot for directly sequenced donor--human joins.

This is NOT a general structural-variant caller or a gapped realigner.  It
recognizes a subset of joins supported by two maximal exact matches in one
cleaned read.  Failure to call is unresolved, never evidence for zero deletion
or zero microhomology (MH).  Exact anchors, primer recovery, local reference
windows and the pre-existing global host-alignment filter impose selection
biases.  MH denotes sequence-compatible breakpoint ambiguity, not a repair
pathway.  An intervening sequence is not necessarily a nontemplated insertion.
Human deletion length cannot be obtained from this one-sided analysis.

All intervals are zero-based, half-open.  All boundaries are zero-based BETWEEN
bases.  ``donor`` is the original forward donor; ``query`` and ``host_sequence``
must already be oriented from donor toward host.  ``host_strand`` describes
that oriented host sequence relative to the forward genomic reference.

No third-party dependencies.  Values returned by public functions are JSON
serializable.  This module never changes the existing four-way classifier.
"""

from collections import Counter, defaultdict


P3 = {
    "3p": "GCTACTCTGGGCTCGATCATCCA",
    "5p": "GCACGCCGATCCGACACACTGAC",
}

DEFAULT_CONFIG = {
    "donor_seed": 20,
    "host_seed": 30,
    "min_donor_anchor": 20,
    "min_host_anchor": 30,
    "min_global_overlap": 30,
    "min_quality": 20,
    "max_insertion": 20,
    "junction_flank": 10,
    "homopolymer_min": 10,
}

_COMPLEMENT = str.maketrans("ACGTRYMKSWBDHVNacgtrymkswbdhvn",
                           "TGCAYRKMSWVHDBNtgcayrkmswvhdbn")
_ACGT = frozenset("ACGT")


def reverse_complement(sequence):
    """Return an upper-case reverse complement, including IUPAC ambiguity."""
    return str(sequence).translate(_COMPLEMENT)[::-1].upper()


def _side(side):
    value = str(side).lower().replace("'", "").replace("\u2032", "")
    if value in ("3", "3p", "3prm"):
        return "3p"
    if value in ("5", "5p", "5prm"):
        return "5p"
    raise ValueError("side must be '3p' or '5p'")


def primer_info(donor, side):
    """Locate the unique supplied P3 in the outward-oriented donor.

    The independent donor anchor must lie AFTER ``primer_end0``.  Primer bases
    cannot support an independent donor anchor even when present in the read.
    """
    side = _side(side)
    donor = str(donor).upper()
    oriented = donor if side == "3p" else reverse_complement(donor)
    primer = P3[side]
    start = oriented.find(primer)
    if start < 0 or oriented.find(primer, start + 1) >= 0:
        raise ValueError("The selected P3 must occur exactly once in the oriented donor")
    end = start + len(primer)
    original_start, original_end = ((start, end) if side == "3p"
                                    else (len(donor) - end, len(donor) - start))
    return {
        "side": side, "oriented_sequence": oriented,
        "primer_sequence": primer, "primer_start0": start, "primer_end0": end,
        "original_primer_start0": original_start,
        "original_primer_end0": original_end,
        "outward_bases": len(donor) - end,
        "donor_length": len(donor),
    }


def maximal_exact_matches(query, target, seed_length=20):
    """Return all maximal exact A/C/G/T matches of at least seed_length.

    Both sequences are used in the supplied orientation.  A seed index avoids
    an all-by-all base comparison.  Each matching seed extends to a mismatch,
    non-ACGT base, query edge or target edge.  Distinct diagonals are retained;
    no best-hit-only resolution or arbitrary upper MH limit is applied.
    """
    query, target = str(query).upper(), str(target).upper()
    seed_length = int(seed_length)
    if seed_length < 1:
        raise ValueError("seed_length must be positive")
    if len(query) < seed_length or len(target) < seed_length:
        return []
    index = defaultdict(list)
    for t in range(len(target) - seed_length + 1):
        word = target[t:t + seed_length]
        if set(word) <= _ACGT:
            index[word].append(t)
    # On a diagonal a previous maximal extension contains later seeds.
    diagonal_ranges = defaultdict(list)
    found = set()
    for q in range(len(query) - seed_length + 1):
        word = query[q:q + seed_length]
        if not set(word) <= _ACGT:
            continue
        for t in index.get(word, ()):
            diagonal = t - q
            if any(a <= q and q + seed_length <= b
                   for a, b in diagonal_ranges[diagonal]):
                continue
            q0, t0, q1, t1 = q, t, q + seed_length, t + seed_length
            while (q0 > 0 and t0 > 0 and query[q0 - 1] in _ACGT
                   and query[q0 - 1] == target[t0 - 1]):
                q0 -= 1
                t0 -= 1
            while (q1 < len(query) and t1 < len(target)
                   and query[q1] in _ACGT and query[q1] == target[t1]):
                q1 += 1
                t1 += 1
            found.add((q0, q1, t0, t1))
            diagonal_ranges[diagonal].append((q0, q1))
    return [
        {"query_start0": q0, "query_end0": q1,
         "target_start0": t0, "target_end0": t1, "length": q1 - q0,
         "target_left_edge": t0 == 0,
         "target_right_edge": t1 == len(target)}
        for q0, q1, t0, t1 in sorted(found)
    ]


def _phred(qualities, length):
    if isinstance(qualities, str):
        result = [ord(c) - 33 for c in qualities]
    else:
        result = [int(q) for q in qualities]
    if len(result) != length:
        raise ValueError("query and qualities must have the same length")
    if any(q < 0 or q > 93 for q in result):
        raise ValueError("qualities must be Phred+33 text or scores in 0..93")
    return result


def _homopolymer_mask(query, minimum):
    mask = set()
    start = 0
    while start < len(query):
        end = start + 1
        while end < len(query) and query[end] == query[start]:
            end += 1
        if query[start] in _ACGT and end - start >= minimum:
            mask.update(range(start, end))
        start = end
    return mask


def _intervals(value):
    if value is None:
        return []
    if (len(value) == 2 and isinstance(value[0], (int, float))
            and isinstance(value[1], (int, float))):
        value = [value]
    intervals = []
    for start, end in value:
        start, end = int(start), int(end)
        if start < 0 or end <= start:
            raise ValueError("accepted host intervals must be nonempty genomic [start,end)")
        intervals.append((start, end))
    return intervals


def _config(config):
    result = dict(DEFAULT_CONFIG)
    if config:
        unknown = set(config) - set(result)
        if unknown:
            raise ValueError("Unknown config keys: " + ", ".join(sorted(unknown)))
        result.update(config)
    for key, value in result.items():
        if not isinstance(value, int) or isinstance(value, bool):
            raise ValueError("Configuration values must be integers: " + key)
        minimum = 0 if key in ("min_quality", "max_insertion", "junction_flank") else 1
        if value < minimum:
            raise ValueError("Configuration value out of range: " + key)
    return result


def call_junctions(query, qualities, donor, side, host_sequence, chrom,
                   host_window_start0, host_strand,
                   accepted_host_interval0=None, config=None):
    """Call exact-anchor D->H explanations within one accepted host context.

    ``donor`` is the original forward sequence, not the side-oriented sequence.
    ``query``/qualities are cleaned and already oriented D->H.  ``host_sequence``
    is the forward genomic window for '+' or its reverse complement for '-'.
    ``host_window_start0`` is ALWAYS the low genomic coordinate of that window.
    ``accepted_host_interval0`` is either one genomic [start,end) or a list of
    M/= /X CIGAR reference blocks from an accepted, same-read global alignment.
    At least min_global_overlap independent Q20 exact host bases must overlap
    these blocks.  Omitting this provenance produces an unresolved result.

    Donor and host anchors are maximal exact matches, not gapped alignments.
    The required independent anchor counts are counted A/C/G/T bases with
    Phred>=min_quality, not necessarily one contiguous high-quality run.
    Donor support excludes overlap and everything through the selected P3.
    Host support excludes overlap, ALL donor MEM intervals on either donor
    orientation, and homopolymers.  Reverse-donor MEMs are used for masking only;
    the supplied library side still determines which D->H joins are tested.
    Every junction also requires high-quality bridge/MH and flanking bases.
    Any host MEM touching a supplied window edge is rejected conservatively;
    enlarge the window and rerun rather than treating truncation as MH=0.

    Returns status='called' (one distinct event), 'ambiguous' (>1), or
    'unresolved' (zero), candidates, and diagnostics.  Candidates retain all
    equivalent MH boundary assignments.  The canonical assignment gives all
    MH bases to donor, hence the minimum donor-terminal deletion.  An insertion
    has mh_length_bp=0 under this operational simple two-arm model; its two
    separate interfaces have NOT been analyzed for MH.
    """
    query, host_sequence = str(query).upper(), str(host_sequence).upper()
    donor, side = str(donor).upper(), _side(side)
    cfg = _config(config)
    phred = _phred(qualities, len(query))
    p3 = primer_info(donor, side)
    oriented_donor = p3["oriented_sequence"]
    if host_strand not in ("+", "-"):
        raise ValueError("host_strand must be '+' or '-'")
    host_window_start0 = int(host_window_start0)
    if host_window_start0 < 0:
        raise ValueError("host_window_start0 must be nonnegative")
    accepted = _intervals(accepted_host_interval0)
    diagnostics = {
        "query_length": len(query), "donor_length": len(donor),
        "host_window_length": len(host_sequence), "config": cfg,
        "primer_end_oriented0": p3["primer_end0"],
        "donor_outward_length": p3["outward_bases"],
        "donor_mem_count": 0, "reverse_donor_mem_count": 0,
        "host_mem_count": 0,
        "ordered_block_pairs": 0, "rejected_counts": {},
    }
    rejected = Counter()

    def result(candidates):
        diagnostics["rejected_counts"] = dict(sorted(rejected.items()))
        diagnostics["distinct_candidate_count"] = len(candidates)
        return {"status": ("called" if len(candidates) == 1 else
                           "ambiguous" if candidates else "unresolved"),
                "candidates": candidates, "diagnostics": diagnostics}

    if not accepted:
        rejected["no_accepted_global_host_interval"] += 1
        return result([])

    donor_mems = maximal_exact_matches(query, oriented_donor, cfg["donor_seed"])
    diagnostics["donor_mem_count"] = len(donor_mems)
    if not donor_mems:
        rejected["no_exact_donor_anchor"] += 1
        return result([])
    host_mems = maximal_exact_matches(query, host_sequence, cfg["host_seed"])
    diagnostics["host_mem_count"] = len(host_mems)
    if not host_mems:
        rejected["no_exact_host_anchor"] += 1
        return result([])

    reverse_donor_mems = maximal_exact_matches(
        query, reverse_complement(oriented_donor), cfg["donor_seed"])
    diagnostics["reverse_donor_mem_count"] = len(reverse_donor_mems)
    donor_mask = set()
    for mem in donor_mems + reverse_donor_mems:
        donor_mask.update(range(mem["query_start0"], mem["query_end0"]))
    homopolymer_mask = _homopolymer_mask(query, cfg["homopolymer_min"])

    def genomic_boundary(oriented_boundary):
        if host_strand == "+":
            return host_window_start0 + oriented_boundary
        return host_window_start0 + len(host_sequence) - oriented_boundary

    def genomic_base(oriented_base):
        if host_strand == "+":
            return host_window_start0 + oriented_base
        return host_window_start0 + len(host_sequence) - 1 - oriented_base

    def original_donor_boundary(oriented_boundary):
        return oriented_boundary if side == "3p" else len(donor) - oriented_boundary

    candidates_by_key = {}
    for d in donor_mems:
        for h in host_mems:
            dq0, dq1 = d["query_start0"], d["query_end0"]
            hq0, hq1 = h["query_start0"], h["query_end0"]
            if not (dq0 < hq0 and dq1 < hq1):
                continue
            diagnostics["ordered_block_pairs"] += 1
            gap = hq0 - dq1
            if gap > cfg["max_insertion"]:
                rejected["intervening_sequence_too_long"] += 1
                continue
            if h["target_left_edge"] or h["target_right_edge"]:
                rejected["host_window_edge_truncated"] += 1
                continue
            mh = max(0, -gap)
            d_positions = [
                q for q in range(dq0, min(dq1, hq0))
                if d["target_start0"] + q - dq0 >= p3["primer_end0"]
                and phred[q] >= cfg["min_quality"]
            ]
            if len(d_positions) < cfg["min_donor_anchor"]:
                rejected["insufficient_outward_donor_anchor"] += 1
                continue
            h_positions = [
                q for q in range(max(hq0, dq1), hq1)
                if q not in donor_mask and q not in homopolymer_mask
                and phred[q] >= cfg["min_quality"]
            ]
            if len(h_positions) < cfg["min_host_anchor"]:
                rejected["insufficient_independent_host_anchor"] += 1
                continue
            global_overlap = sum(
                any(a <= genomic_base(h["target_start0"] + q - hq0) < b
                    for a, b in accepted)
                for q in h_positions
            )
            if global_overlap < cfg["min_global_overlap"]:
                rejected["insufficient_accepted_alignment_overlap"] += 1
                continue
            bridge0, bridge1 = min(dq1, hq0), max(dq1, hq0)
            quality_start = bridge0 - cfg["junction_flank"]
            quality_end = bridge1 + cfg["junction_flank"]
            if quality_start < dq0 or quality_end > hq1:
                rejected["junction_flank_not_sequenced"] += 1
                continue
            if (any(phred[q] < cfg["min_quality"]
                    for q in range(quality_start, quality_end))
                    or any(query[q] not in _ACGT
                           for q in range(quality_start, quality_end))):
                rejected["low_quality_or_ambiguous_junction"] += 1
                continue

            donor_max = d["target_end0"]
            donor_min = donor_max - mh
            host_local_min = h["target_start0"]
            host_local_max = host_local_min + mh
            host_canonical = genomic_boundary(host_local_max)
            host_boundaries = sorted((genomic_boundary(host_local_min), host_canonical))
            original_boundaries = sorted((original_donor_boundary(donor_min),
                                          original_donor_boundary(donor_max)))
            insertion = query[dq1:hq0] if gap > 0 else ""
            mh_sequence = query[hq0:dq1] if mh else ""
            equivalent = [
                {"donor_cut_oriented0": donor_min + shift,
                 "donor_cut_original0": original_donor_boundary(donor_min + shift),
                 "host_boundary0": genomic_boundary(host_local_min + shift),
                 "donor_deletion_bp": len(donor) - donor_min - shift}
                for shift in range(mh + 1)
            ]
            flags = []
            if insertion:
                flags.append("intervening_sequence_origin_not_assessed")
            if dq0:
                flags.append("query_prefix_outside_two_arm_model")
            if hq1 < len(query):
                flags.append("query_suffix_outside_two_arm_model")
            candidate = {
                "side": side, "chrom": str(chrom), "host_strand": host_strand,
                "host_boundary0": host_canonical,
                "host_boundary_min0": host_boundaries[0],
                "host_boundary_max0": host_boundaries[1],
                "host_window_start0": host_window_start0,
                "host_window_end0": host_window_start0 + len(host_sequence),
                "donor_cut_oriented_min0": donor_min,
                "donor_cut_oriented_max0": donor_max,
                "donor_cut_original_min0": original_boundaries[0],
                "donor_cut_original_max0": original_boundaries[1],
                "donor_deletion_min_bp": len(donor) - donor_max,
                "donor_deletion_max_bp": len(donor) - donor_min,
                "mh_length_bp": mh, "mh_sequence": mh_sequence,
                "insertion_length_bp": len(insertion), "insertion_sequence": insertion,
                "donor_independent_bases": len(d_positions),
                "host_independent_bases": len(h_positions),
                "host_global_overlap_bases": global_overlap,
                "donor_block": dict(d), "host_block": dict(h),
                "unmodeled_query_prefix_bp": dq0,
                "unmodeled_query_suffix_bp": len(query) - hq1,
                "junction_query_start0": bridge0, "junction_query_end0": bridge1,
                "junction_sequence": query[quality_start:quality_end],
                "equivalent_boundaries": equivalent,
                "equivalent_alignment_count": 1,
                "flags": flags,
            }
            key = (side, donor_max, str(chrom), host_canonical, host_strand, insertion)
            candidate["event_key"] = list(key)
            previous = candidates_by_key.get(key)
            if previous is None:
                candidates_by_key[key] = candidate
            else:
                # This only chooses between explanations of the SAME event.
                # Distinct boundaries/loci/insertions are never score-resolved.
                def support(c):
                    return (min(c["donor_independent_bases"], c["host_independent_bases"]),
                            c["donor_independent_bases"] + c["host_independent_bases"],
                            c["donor_block"]["length"] + c["host_block"]["length"],
                            c["mh_length_bp"])
                count = previous["equivalent_alignment_count"] + 1
                chosen = candidate if support(candidate) > support(previous) else previous
                chosen["equivalent_alignment_count"] = count
                candidates_by_key[key] = chosen

    if not diagnostics["ordered_block_pairs"]:
        rejected["no_ordered_donor_to_host_blocks"] += 1
    return result([candidates_by_key[key] for key in sorted(candidates_by_key)])
