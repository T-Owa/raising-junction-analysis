#!/usr/bin/env python3
"""Combine completed pilot outputs without treating PCR reads as events."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import sys

VERSION = "0.1.0"
STRATA = {
    "without_prior_dd_or_complex_host_support": "simple_pair_support",
    "with_prior_dd_without_complex_host_support": "prior_dd_pair_support",
    "with_complex_host_support": "complex_host_pair_support",
    "all_called": "pair_support",
}
BASE = ["sample", "library_side", "condition", "biological_replicate"]
STATUSES = ["called_simple_candidate", "called_prior_donor_rearrangement",
            "called_complex_host_evidence", "ambiguous", "unresolved"]


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def write_csv(path, fields, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize(root, samples=None, output=None):
    root = Path(root).expanduser().resolve()
    paths = sorted((root / "results").glob("*_junction_pilot_v1/summary.json"))
    wanted = set(samples or [])
    runs, seen = [], set()
    signature = None
    for path in paths:
        s = json.loads(path.read_text(encoding="utf-8"))
        if wanted and s.get("sample") not in wanted:
            continue
        if s.get("analysis_type") != "raising_exact_anchor_junction_pilot" or s.get("status") != "complete":
            raise ValueError("Not a complete pilot output: " + str(path))
        if s["sample"] in seen:
            raise ValueError("Duplicate sample: " + s["sample"])
        seen.add(s["sample"])
        current = {"analysis_version": s["analysis_version"], "config": s["config"],
                   "source_sha256": s["source_sha256"], "reference_sha256": s["reference"]["sha256"]}
        if signature is None:
            signature = current
        elif signature != current:
            raise ValueError("Do not pool different caller/config/reference versions in one table")
        if not s.get("validation") or not all(s["validation"].values()):
            raise ValueError("Pilot validation failed: " + str(path))
        junctions = path.parent / "junctions.tsv"
        if sha256(junctions) != s["output_sha256"].get("junctions.tsv"):
            raise ValueError("junctions.tsv differs from recorded SHA-256: " + str(path))
        with junctions.open(encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            required = {"event_id", "mh_length_bp", "insertion_length_bp",
                        "donor_deletion_min_bp", "donor_deletion_max_bp", *STRATA.values()}
            if not required.issubset(set(reader.fieldnames or [])):
                raise ValueError("Missing junction fields in " + str(junctions))
            rows = list(reader)
        event_ids = set()
        for r in rows:
            if r["event_id"] in event_ids or r["sample"] != s["sample"] or r["library_side"] != s["library_side"]:
                raise ValueError("Duplicate/mismatched junction identity")
            event_ids.add(r["event_id"])
            for key in required - {"event_id"}:
                r[key] = int(r[key])
                if r[key] < 0:
                    raise ValueError("Negative count/length in junction table")
            if r["donor_deletion_min_bp"] > r["donor_deletion_max_bp"]:
                raise ValueError("Inverted deletion range")
            if r["pair_support"] != sum(r[k] for k in ("simple_pair_support", "prior_dd_pair_support", "complex_host_pair_support")):
                raise ValueError("Junction support strata do not reconcile")
        stats = s["stats"]
        if (len(rows) != stats["observed_junction_patterns_total"] or
            sum(r["pair_support"] for r in rows) != stats["called_pairs_total"] or
            sum(stats[k] for k in STATUSES) != stats["host_supported_pairs"]):
            raise ValueError("Junction/status counts disagree with summary")
        for field, status in (("simple_pair_support", "called_simple_candidate"),
                              ("prior_dd_pair_support", "called_prior_donor_rearrangement"),
                              ("complex_host_pair_support", "called_complex_host_evidence")):
            if sum(r[field] for r in rows) != stats[status]:
                raise ValueError("Junction support strata disagree with summary: " + field)
        runs.append((path, s, rows))
    if wanted - seen:
        raise ValueError("Missing completed pilot samples: " + ", ".join(sorted(wanted - seen)))
    if not runs:
        raise ValueError("No completed *_junction_pilot_v1 outputs found")

    qc_rows, dist_rows, all_events, metadata, inputs = [], [], [], [], []
    for path, s, rows in runs:
        base = {k: s.get(k, "") for k in BASE}
        metadata.append(base)
        inputs.append({"summary": str(path), "sha256": sha256(path)})
        stats = s["stats"]
        for stratum, support in STRATA.items():
            subset = [r for r in rows if r[support] > 0]
            qc = {**base, "stratum": stratum, **stats,
                  "stratum_called_pairs": sum(r[support] for r in subset),
                  "stratum_observed_patterns": len(subset),
                  "patterns_with_insertion": sum(r["insertion_length_bp"] > 0 for r in subset),
                  "patterns_with_ambiguous_deletion_boundary": sum(r["donor_deletion_min_bp"] != r["donor_deletion_max_bp"] for r in subset),
                  "patterns_supported_by_at_least_two_pairs_in_stratum": sum(r[support] >= 2 for r in subset)}
            qc_rows.append(qc)
            for metric, field in (("MH_length_no_intervening_sequence", "mh_length_bp"),
                                  ("intervening_sequence_length", "insertion_length_bp"),
                                  ("donor_terminal_deletion_min", "donor_deletion_min_bp"),
                                  ("donor_terminal_deletion_max", "donor_deletion_max_bp")):
                eligible = [r for r in subset if metric != "MH_length_no_intervening_sequence" or r["insertion_length_bp"] == 0]
                for unit in ("observed_junction_patterns", "supporting_read_pairs"):
                    counts = collections.Counter()
                    for r in eligible:
                        counts[r[field]] += 1 if unit == "observed_junction_patterns" else r[support]
                    denominator = sum(counts.values())
                    # Include every integer from zero through the observed maximum;
                    # empty eligible populations keep an explicit undefined row.
                    values = range(max(counts, default=0) + 1) if denominator else [None]
                    for value in values:
                        count = counts[value] if value is not None else 0
                        dist_rows.append({**base, "stratum": stratum, "metric": metric,
                                          "unit": unit, "value_bp": value, "count": count,
                                          "denominator": denominator,
                                          "fraction": count / denominator if denominator else None})
        all_events.extend({**r, **base} for r in rows)

    out = Path(output).expanduser().resolve() if output else root / "results" / "junction_pilot_summary_v1"
    out.mkdir(parents=True, exist_ok=True)
    # Only this program's output filenames are replaced. It never edits inputs.
    write_csv(out / "library_qc.csv", list(qc_rows[0]), qc_rows)
    write_csv(out / "figure_distributions.csv", BASE + ["stratum", "metric", "unit", "value_bp", "count", "denominator", "fraction"], dist_rows)
    write_csv(out / "sample_metadata_template.csv", BASE, metadata)
    event_fields = BASE + list(dict.fromkeys(k for r in all_events for k in r if k not in BASE))
    if not all_events:
        event_fields += ["event_id", "mh_length_bp", "donor_deletion_min_bp", "donor_deletion_max_bp", "insertion_length_bp", *dict.fromkeys(STRATA.values())]
    write_csv(out / "all_junction_patterns.csv", event_fields, all_events)
    report = {"analysis_version": VERSION, "libraries": len(runs), "shared_configuration": signature,
              "input_summaries": inputs,
              "notes": ["All-called unique patterns are not the sum of stratum-specific unique patterns: the same boundary can have different pair support strata.",
                        "All-called contains complex and prior-DD evidence; use without_prior_dd_or_complex_host_support for the initial primary descriptive panel.",
                        "MH distribution excludes insertion-positive junctions. Unresolved/ambiguous pairs are never MH zero.",
                        "Deletion min/max distributions are paired alternative boundary representations, not independent observations or confidence bounds.",
                        "At least two supporting pairs is a count threshold, not proof of independent molecules.",
                        "No biological group inference or statistical significance testing is performed."]}
    (out / "aggregation_summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Combined {len(runs)} libraries. Tables: {out}")
    return report


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True)
    p.add_argument("--sample", action="append", dest="samples", help="Optional sample restriction; repeatable")
    p.add_argument("--output-dir")
    a = p.parse_args(argv)
    try:
        summarize(a.root, a.samples, a.output_dir)
    except (ValueError, KeyError, OSError) as e:
        print("ERROR: " + str(e), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
