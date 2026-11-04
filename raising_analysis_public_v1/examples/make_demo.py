#!/usr/bin/env python3
"""Generate a deterministic, entirely synthetic RAISING demonstration dataset.

The toy host is randomly generated, not human sequence. The bundled experimental
donor is used only as an analysis reference. No participant or experimental reads
are included. R2 is the full reverse complement of the same synthetic amplicon as
R1; each designed template is represented by two separately named read pairs.
These repeats are software fixtures, not independent molecules or replicates.
"""

import argparse
import csv
import gzip
import io
import json
from pathlib import Path
import random


P3 = {"3p": "GCTACTCTGGGCTCGATCATCCA", "5p": "GCACGCCGATCCGACACACTGAC"}
ADAPTER = "TACTGCTTGACTGACCTGCTGTC"
CATEGORIES = (
    "donor_rearrangement_supported",
    "host_mappable_other",
    "short_donor_primer_candidate",
    "other_unresolved",
)


def revcomp(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def load_donor():
    path = (Path(__file__).resolve().parents[1] / "pipelines" / "fourway" /
            "legacy_v1" / "RI_donor_800bp.fasta")
    donor = "".join(line.strip() for line in path.read_text().splitlines()
                    if line and not line.startswith(">"))
    if len(donor) != 800 or set(donor) - set("ACGT"):
        raise ValueError("Expected the bundled single 800 bp donor")
    return donor


def different_base(base):
    return next(nucleotide for nucleotide in "ACGT" if nucleotide != base)


def write_fastq(path, reads):
    # Empty gzip filename and mtime=0 prevent machine/time-specific gzip headers.
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="ascii", newline="\n") as handle:
                for identifier, sequence in reads:
                    handle.write(f"@{identifier}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def generate(output, seed=20261024):
    output = Path(output).expanduser().absolute()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Output must be absent or an empty directory; existing files are not overwritten")
    donor = load_donor()
    rng = random.Random(seed)
    host = list("".join(rng.choice("ACGT") for _ in range(100000)))
    designs = []

    for side_index, side in enumerate(("3p", "5p")):
        oriented = donor if side == "3p" else revcomp(donor)
        start = oriented.index(P3[side])
        sample = "sample001_3prm_S1" if side == "3p" else "sample001_5prm_S2"
        base_position = 10000 + side_index * 40000
        for index, kind in enumerate(("direct", "donor_deletion_5", "mh_3", "insertion_3")):
            position = base_position + index * 5000
            cut = 795 if kind == "donor_deletion_5" else 800
            mh = 3 if kind == "mh_3" else 0
            inserted = "GCA" if kind == "insertion_3" else ""
            # Force the designed boundary: exclude accidental extension at the
            # donor/host edges and cap the deliberately introduced MH at 3 bp.
            host[position - 1] = different_base(
                inserted[-1] if inserted else oriented[cut - mh - 1])
            if mh:
                host[position:position + mh] = list(oriented[cut - mh:cut])
            if cut < len(oriented):
                host[position + mh] = different_base(oriented[cut])
            fragment = "".join(host[position + mh:position + 200])
            sequence = oriented[start:cut] + inserted + fragment
            designs.append({
                "sample": sample, "library_side": side, "fixture": kind,
                "category": CATEGORIES[1], "sequence": sequence,
                "host_contig": "synthetic_contig_1",
                "designed_host_match_start0": position,
                "designed_donor_loss_min_bp": 800 - cut,
                "designed_donor_loss_max_bp": 800 - cut + mh,
                "designed_mh_bp": None if inserted else mh,
                "designed_intervening_sequence": inserted,
            })
        # These three templates complete the exclusive four-way categories.
        extras = (
            ("donor_donor", CATEGORIES[0], oriented[start:] + oriented[:180]),
            ("short_donor_primer", CATEGORIES[2], oriented[start:start + 45] + "A" * 20 + ADAPTER),
            ("unresolved", CATEGORIES[3], "".join(rng.choice("ACGT") for _ in range(130))),
        )
        for kind, category, sequence in extras:
            designs.append({"sample": sample, "library_side": side, "fixture": kind,
                            "category": category, "sequence": sequence})

    raw = output / "vendor_results" / "Raw_fastq"
    reference = output / "ref"
    raw.mkdir(parents=True, exist_ok=True)
    reference.mkdir(parents=True, exist_ok=True)
    host_text = "".join(host)
    (reference / "synthetic_host.fa").write_text(
        ">synthetic_contig_1 entirely_synthetic_toy_reference_not_human\n" +
        "\n".join(host_text[i:i + 80] for i in range(0, len(host_text), 80)) + "\n",
        encoding="ascii")
    (reference / "RI_donor_800bp.fasta").write_text(
        ">RI_donor_800bp study_donor_reference\n" + donor + "\n", encoding="ascii")

    libraries = []
    for sample in sorted({design["sample"] for design in designs}):
        selected = [design for design in designs if design["sample"] == sample]
        first, second = [], []
        for design in selected:
            for repeat in (1, 2):
                identifier = f"{sample}:{design['fixture']}:{repeat}"
                first.append((identifier + "/1", design["sequence"]))
                second.append((identifier + "/2", revcomp(design["sequence"])))
        write_fastq(raw / f"{sample}_L001_R1_001.fastq.gz", first)
        write_fastq(raw / f"{sample}_L001_R2_001.fastq.gz", second)
        libraries.append({"sample": sample, "library_side": selected[0]["library_side"],
                          "total_pairs": len(first),
                          "expected_fourway_counts": dict(zip(CATEGORIES, (2, 8, 2, 2))),
                          "expected_supported_host_junction_groups": 4,
                          "expected_supported_pairs_per_host_group": 2})

    metadata = output / "sample_metadata.csv"
    with metadata.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "library_side", "condition", "biological_replicate"])
        writer.writeheader()
        for library in libraries:
            writer.writerow({"sample": library["sample"], "library_side": library["library_side"],
                             "condition": "synthetic_demo", "biological_replicate": ""})
    manifest = {
        "dataset_type": "entirely_synthetic_software_fixture_not_biological_data",
        "seed": seed, "toy_host_length": len(host),
        "counting_warning": "Repeated read pairs are software fixtures, not independent molecules or biological replicates.",
        "expected_results_status": "Designed expectations; verify with the installed aligner versions.",
        "libraries": libraries,
        "designs": [{key: value for key, value in design.items() if key != "sequence"} for design in designs],
    }
    (output / "expected_demo.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True, help="Absent or empty destination for the toy analysis root")
    parser.add_argument("--seed", type=int, default=20261024)
    args = parser.parse_args()
    try:
        result = generate(args.outdir, args.seed)
    except (ValueError, OSError) as error:
        parser.exit(1, str(error) + "\n")
    print(f"Created {len(result['libraries'])} synthetic libraries, 14 read pairs each, in {args.outdir}")
    print("Toy host: ref/synthetic_host.fa. Designed expectations: expected_demo.json.")


if __name__ == "__main__":
    main()
