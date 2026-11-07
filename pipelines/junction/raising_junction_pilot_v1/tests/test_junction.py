"""Synthetic junction truth tests; these do not validate real sample accuracy."""
import csv
import gzip
import hashlib
import json
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE))
import junction_core as core
import raising_junction_pilot as runner

DONOR = (
    "TAGAGATGTTAAGCTGTTGCTTACCTTCTGAATCACCGCAACATCGGTTGAGCAGAATAGCTCCGGTCAGTGTGTCGGAT"
    "CGGCGTGCCCTGAATTAGCATTAGGGCTAGTAACCCTGTACAATATGGGCCGCTCGTCCGCGTACCGGAAACCGTCTTGA"
    "ACGACCAGAGACATAGTTCTAAGGTGGCTCAAGAAGGCGAGTTAAATGCACCCGCTTAGCCCTCAGAATCTGCCAATATC"
    "AATGCGACGCTGATACGCGGCTAAAGGATAGCTGTTATACAGTGCGTGCACACTCCCACGTAGGACTCCGCCCTCTTTGA"
    "TGTGACCTAAAACCAAGTTGTTTCTGCTACTTAGCTGAGGGAGTTCAACCACGTGTATCATGAGTAGTAGAGCCCTCTGG"
    "ATAACGGTCATCTATTGCGGGTACTAACACCCGAGTACCTTCTTGAAGTATGCTAGGCGCCCCGAAGTGTCCGAATACCC"
    "CGATTTCAACATAGCGGTGAGCGGAAAGCATCGTGATTCAGACCCAACATGCTGAGCTCAGTGATGGGTCTCCCCGTAGG"
    "TATCTATGCTGACGTCTGATCGCTCTTAGGTTTAAGGAATTTTGCAGGCTGGCTCCTGTTTATAACGGCCTGCTCTCGGA"
    "ATGACAACCCAAGATACTATCGCACACGTCGTTAAGGAAGAGTAACGATGATCCTTAGGAGCTAAGTGTCAATATAGAGC"
    "TACTCTGGGCTCGATCATCCAGTTTACTGACGCGCGTTGCCCGCCGAGGACAGATTCACAGTAATGTGACTCTGTCGAAT"
)
PRIMERS = {"3p": "GCTACTCTGGGCTCGATCATCCA", "5p": "GCACGCCGATCCGACACACTGAC"}


def rc(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def fixture(side="3p", mh=0, insertion="", deletion=10, strand="+"):
    """Build an exact junction with independently specified breakpoint truth.

    Donor dcut is the retained endpoint AFTER shared bases; host hstart is
    BEFORE shared bases. Flanking mismatches prevent accidental MH extension.
    """
    if mh and insertion:
        raise ValueError("fixture represents either MH or insertion")
    donor = DONOR if side == "3p" else rc(DONOR)
    p_end = donor.index(PRIMERS[side]) + len(PRIMERS[side])
    dcut = len(donor) - deletion
    dpart = donor[p_end:dcut]
    rng = random.Random(3109 + (side == "5p"))
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    shared = donor[dcut-mh:dcut] if mh else ""
    # The first base after donor retention cannot extend the donor MEM.
    after = rand(70)
    disallowed = donor[dcut] if dcut < len(donor) else "N"
    after = next(c for c in "ACGT" if c != disallowed) + after[1:]
    if insertion:
        assert insertion[0] != disallowed
    hpart = shared + after
    prefix = rand(80)
    preceding = insertion[-1] if insertion else donor[dcut-mh-1]
    prefix = prefix[:-1] + next(c for c in "ACGT" if c != preceding)
    suffix = rand(80)
    sequence = prefix + hpart + suffix
    query = dpart + insertion + after
    window_start = 1000
    hstart, hend = len(prefix), len(prefix) + len(hpart)
    if strand == "+":
        interval = (window_start + hstart, window_start + hend)
        canonical_boundary = window_start + hstart + mh
    else:
        interval = (window_start + len(sequence) - hend,
                    window_start + len(sequence) - hstart)
        canonical_boundary = window_start + len(sequence) - hstart - mh
    return dict(query=query, qualities="I" * len(query), donor=DONOR,
                side=side, host_sequence=sequence, chrom="chrFixture",
                host_window_start0=window_start, host_strand=strand,
                accepted_host_interval0=interval), dict(
                    p_end=p_end, donor_cut=dcut, hstart=hstart, hend=hend,
                    dpart=dpart, hpart=hpart, prefix=prefix, suffix=suffix,
                    canonical_boundary=canonical_boundary, shared=shared)


class JunctionTruthTests(unittest.TestCase):
    def one_call(self, args):
        result = core.call_junctions(**args)
        self.assertEqual(result["status"], "called", result)
        self.assertEqual(len(result["candidates"]), 1, result)
        return result["candidates"][0]

    def test_fixture_donor_and_primer_positions(self):
        self.assertEqual(len(DONOR), 800)
        for side, outward in (("3p", 59), ("5p", 65)):
            seq = DONOR if side == "3p" else rc(DONOR)
            self.assertEqual(len(seq) - seq.index(PRIMERS[side]) - len(PRIMERS[side]), outward)
            self.assertEqual(core.primer_info(DONOR, side)["outward_bases"], outward)

    def test_three_prime_blunt_truth(self):
        args, truth = fixture()
        call = self.one_call(args)
        self.assertEqual(call["mh_length_bp"], 0)
        self.assertEqual(call["mh_sequence"], "")
        self.assertEqual(call["insertion_length_bp"], 0)
        self.assertEqual(call["donor_deletion_min_bp"], 10)
        self.assertEqual(call["donor_deletion_max_bp"], 10)
        self.assertEqual(call["donor_cut_original_min0"], 790)
        self.assertEqual(call["donor_cut_original_max0"], 790)
        self.assertEqual(call["host_boundary0"], truth["canonical_boundary"])

    def test_five_prime_boundary_is_original_donor_left_end(self):
        args, truth = fixture(side="5p", deletion=12)
        call = self.one_call(args)
        self.assertEqual(call["donor_deletion_min_bp"], 12)
        self.assertEqual(call["donor_deletion_max_bp"], 12)
        self.assertEqual(call["donor_cut_oriented_min0"], 788)
        self.assertEqual(call["donor_cut_original_min0"], 12)
        self.assertEqual(call["donor_cut_original_max0"], 12)
        self.assertEqual(call["host_boundary0"], truth["canonical_boundary"])

    def test_shared_sequence_yields_range_not_arbitrary_deletion(self):
        args, truth = fixture(mh=5)
        call = self.one_call(args)
        self.assertEqual(call["mh_length_bp"], 5)
        self.assertEqual(call["mh_sequence"], truth["shared"])
        self.assertEqual((call["donor_deletion_min_bp"], call["donor_deletion_max_bp"]), (10, 15))
        self.assertEqual((call["donor_cut_oriented_min0"], call["donor_cut_oriented_max0"]), (785, 790))
        self.assertEqual(call["host_boundary0"], truth["canonical_boundary"])
        self.assertEqual(len(call["equivalent_boundaries"]), 6)

    def test_five_prime_mh_and_host_reverse_coordinates(self):
        args, truth = fixture(side="5p", mh=3, deletion=9, strand="-")
        call = self.one_call(args)
        self.assertEqual(call["mh_length_bp"], 3)
        self.assertEqual((call["donor_deletion_min_bp"], call["donor_deletion_max_bp"]), (9, 12))
        self.assertEqual((call["donor_cut_original_min0"], call["donor_cut_original_max0"]), (9, 12))
        self.assertEqual(call["host_boundary0"], truth["canonical_boundary"])
        self.assertEqual(call["host_boundary_max0"] - call["host_boundary_min0"], 3)

    def test_intervening_sequence_is_reported_separately(self):
        # donor[790] is T; leading A prevents extension of donor anchor.
        args, _ = fixture(insertion="AGC")
        call = self.one_call(args)
        self.assertEqual(call["insertion_sequence"], "AGC")
        self.assertEqual(call["insertion_length_bp"], 3)
        self.assertEqual(call["mh_length_bp"], 0)
        self.assertEqual(call["donor_deletion_min_bp"], 10)

    def test_weak_donor_anchor_not_mh_zero_call(self):
        args, truth = fixture(deletion=40)
        self.assertEqual(len(truth["dpart"]), 19)
        result = core.call_junctions(**args)
        self.assertEqual(result["status"], "unresolved", result)
        self.assertFalse(result["candidates"])

    def test_primer_only_donor_does_not_count_as_donor_anchor(self):
        args, truth = fixture()
        args["query"] = PRIMERS["3p"] + truth["hpart"]
        args["qualities"] = "I" * len(args["query"])
        result = core.call_junctions(**args)
        self.assertEqual(result["status"], "unresolved", result)

    def test_host_anchor_which_is_also_donor_is_not_independent(self):
        for hpart in (DONOR[100:170], rc(DONOR[100:170])):
            with self.subTest(reverse=hpart != DONOR[100:170]):
                args, truth = fixture()
                args["query"] = truth["dpart"] + hpart
                args["qualities"] = "I" * len(args["query"])
                args["host_sequence"] = truth["prefix"] + hpart + truth["suffix"]
                args["accepted_host_interval0"] = (1080, 1150)
                result = core.call_junctions(**args)
                self.assertEqual(result["status"], "unresolved", result)

    def test_low_quality_at_junction_prevents_call(self):
        for mh in (0, 5):
            with self.subTest(mh=mh):
                args, truth = fixture(mh=mh)
                q = list(args["qualities"])
                q[len(truth["dpart"]) - 1] = "!"
                args["qualities"] = "".join(q)
                result = core.call_junctions(**args)
                self.assertEqual(result["status"], "unresolved", result)

    def test_two_valid_host_locations_remain_ambiguous(self):
        args, truth = fixture()
        # Both copies have controlled preceding mismatches and ample flank.
        args["host_sequence"] = (truth["prefix"] + truth["hpart"] +
                                 truth["prefix"] + truth["hpart"] + truth["suffix"])
        start2 = 80 + len(truth["hpart"]) + 80
        args["accepted_host_interval0"] = [(1080, 1080 + len(truth["hpart"])),
                                            (1000 + start2, 1000 + start2 + len(truth["hpart"]))]
        result = core.call_junctions(**args)
        self.assertEqual(result["status"], "ambiguous", result)
        self.assertEqual(len(result["candidates"]), 2, result)

    def test_global_alignment_locus_required(self):
        args, _ = fixture()
        args["accepted_host_interval0"] = (1200, 1230)
        self.assertEqual(core.call_junctions(**args)["status"], "unresolved")

    def test_window_edge_is_not_trusted_as_true_match_end(self):
        args, truth = fixture()
        args["host_sequence"] = truth["prefix"] + truth["hpart"]
        result = core.call_junctions(**args)
        self.assertEqual(result["status"], "unresolved", result)

    def test_mate_only_evidence_cannot_supply_missing_junction(self):
        args, truth = fixture()
        for query in (truth["dpart"], truth["hpart"]):
            args["query"] = query
            args["qualities"] = "I" * len(query)
            result = core.call_junctions(**args)
            self.assertEqual(result["status"], "unresolved", result)


def write_table(path, rows):
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_table(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_fasta_and_fai(path, records):
    entries = []
    with path.open("wb") as handle:
        for name, sequence in records:
            handle.write((">" + name + "\n").encode())
            offset = handle.tell()
            handle.write((sequence + "\n").encode())
            entries.append(f"{name}\t{len(sequence)}\t{offset}\t{len(sequence)}\t{len(sequence)+1}\n")
    Path(str(path) + ".fai").write_text("".join(entries))


def make_runner_fixture(root, insertion=""):
    sample = "Synthetic_3prm_S1"
    evidence = root / "results" / (sample + "_fourway_full_v1")
    raw = root / "vendor_results" / "Raw_fastq"
    reference = root / "ref" / "synthetic" / "human_plus_donor.fa"
    for path in (evidence, raw, reference.parent):
        path.mkdir(parents=True, exist_ok=True)
    args, truth = fixture(insertion=insertion)
    chrseq = "N" * 1000 + args["host_sequence"] + "N" * 1000
    write_fasta_and_fai(reference, [("chrFixture", chrseq), ("__RAISING_DONOR_800__", DONOR)])
    queries = {
        "pairA": {"R1": args["query"], "R2": rc(args["query"])},
        "pairB": {"R1": args["query"], "R2": truth["hpart"]},
        "pairC": {"R1": truth["dpart"], "R2": truth["hpart"]},
        "pairD": {"R1": truth["dpart"], "R2": rc(truth["dpart"])},
    }
    paths = {}
    for mate in ("R1", "R2"):
        path = raw / f"{sample}_L001_{mate}_001.fastq.gz"
        with gzip.open(path, "wt") as handle:
            for name, ends in queries.items():
                sequence = "TGT" + ends[mate] + "AC"
                handle.write(f"@{name}/{mate[-1]}\n{sequence}\n+\n{'I'*len(sequence)}\n")
        paths[mate] = path
    pairs, reads, hits = [], [], []
    for name, ends in queries.items():
        host = name != "pairD"
        dd = name == "pairB"
        pairs.append(dict(read_id=name,
                          classification="donor_rearrangement_supported" if dd else
                          "host_mappable_other" if host else "other_unresolved",
                          host_mapping_supported=int(host), donor_rearrangement_supported=int(dd)))
        for mate, query in ends.items():
            reads.append(dict(read_id=name, mate=mate, raw_length=len(query)+5,
                              cleaned_length=len(query), cleaned_raw_start_1based=4,
                              cleaned_raw_end_1based=len(query)+3))
    for name, mate, cigar, strand in (
        ("pairA", "R1", f"{len(truth['dpart'])+len(insertion)}S{len(truth['hpart'])}M", "+"),
        ("pairA", "R2", f"{len(truth['hpart'])}M{len(truth['dpart'])+len(insertion)}S", "-"),
        ("pairB", "R1", f"{len(truth['dpart'])+len(insertion)}S{len(truth['hpart'])}M", "+"),
        ("pairB", "R2", f"{len(truth['hpart'])}M", "+"),
        ("pairC", "R2", f"{len(truth['hpart'])}M", "+"),
    ):
        hits.append(dict(read_id=name, mate=mate, reference="chrFixture", position_1based=1081,
                         cigar=cigar, strand=strand, accepted=1, mapq=60,
                         supplementary=0))
    write_table(evidence / "pair_classification.tsv.gz", pairs)
    write_table(evidence / "read_evidence.tsv.gz", reads)
    write_table(evidence / "host_alignment_evidence.tsv.gz", hits)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = dict(sample=sample, side="3p", status="complete", total_pairs=4,
                    host_supported_pairs=3, analysis_version="1.0.0",
                    inputs=dict(lanes=[{mate: dict(path=str(path), sha256=sha(path))
                                        for mate, path in paths.items()}]),
                    reference=dict(combined_fasta=str(reference),
                                   signature=dict(donor_sha256=hashlib.sha256(DONOR.encode()).hexdigest(),
                                                  donor_contig="__RAISING_DONOR_800__")),
                    final_output_sha256={p.name: sha(p) for p in evidence.glob("*.tsv.gz")})
    (evidence / "run_summary.json").write_text(json.dumps(manifest))
    return sample, evidence, paths


class RunnerIntegrationTests(unittest.TestCase):
    def invoke(self, root, sample, *extra):
        return subprocess.run([sys.executable, str(PACKAGE / "raising_junction_pilot.py"),
                               "--root", str(root), "--sample", sample, *extra],
                              capture_output=True, text=True, timeout=30)

    def test_full_input_selection_pair_support_and_missing_junction(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample, evidence, _ = make_runner_fixture(root)
            output = root / "pilot"
            proc = self.invoke(root, sample, "--output-dir", str(output))
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            rows = read_table(output / "pair_junction_status.tsv.gz")
            self.assertEqual(len(rows), 3, rows)
            self.assertEqual({r["read_id"] for r in rows}, {"pairA", "pairB", "pairC"})
            by_id = {r["read_id"]: r for r in rows}
            self.assertEqual(by_id["pairA"]["status"], "called_simple_candidate")
            self.assertEqual(by_id["pairB"]["status"], "called_prior_donor_rearrangement")
            self.assertEqual(by_id["pairC"]["status"], "unresolved")
            # Shared sequence from A's two mates must count as ONE pair.
            patterns = read_table(output / "junctions.tsv")
            self.assertEqual(len(patterns), 1)
            self.assertEqual(int(patterns[0]["pair_support"]), 2)
            self.assertEqual(int(patterns[0]["simple_pair_support"]), 1)
            self.assertEqual(int(patterns[0]["prior_dd_pair_support"]), 1)
            summary = json.loads((output / "summary.json").read_text())
            self.assertEqual(summary["stats"]["total_input_pairs"], 4)
            self.assertEqual(summary["stats"]["host_supported_pairs"], 3)
            self.assertEqual(summary["stats"]["called_pairs_total"], 2)
            self.assertTrue(all(summary["validation"].values()))

    def test_changed_raw_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample, evidence, paths = make_runner_fixture(root)
            with paths["R1"].open("ab") as handle:
                handle.write(b"changed")
            proc = self.invoke(root, sample, "--output-dir", str(root / "pilot"))
            self.assertNotEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn("Raw FASTQ SHA-256 differs", proc.stderr)

    def test_changed_evidence_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample, evidence, _ = make_runner_fixture(root)
            table = evidence / "read_evidence.tsv.gz"
            rows = read_table(table)
            rows[0]["cleaned_length"] = "999"
            write_table(table, rows)
            proc = self.invoke(root, sample, "--output-dir", str(root / "pilot"))
            self.assertNotEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn("Evidence SHA-256 differs", proc.stderr)

    def test_inconsistent_cleaned_coordinates_rejected_after_hash_verification(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample, evidence, _ = make_runner_fixture(root)
            table = evidence / "read_evidence.tsv.gz"
            rows = read_table(table)
            rows[0]["cleaned_raw_start_1based"] = "5"
            write_table(table, rows)
            manifest_path = evidence / "run_summary.json"
            manifest = json.loads(manifest_path.read_text())
            manifest["final_output_sha256"][table.name] = hashlib.sha256(table.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest))
            proc = self.invoke(root, sample, "--output-dir", str(root / "pilot"))
            self.assertNotEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn("Invalid cleaned-read coordinates", proc.stderr)

    def test_distinct_mate_junctions_are_not_best_score_resolved(self):
        args, _ = fixture()
        first = core.call_junctions(**args)["candidates"][0]
        second = json.loads(json.dumps(first))
        second["host_boundary0"] += 50
        second["event_key"][3] += 50
        calls = [dict(candidate=first, mate="R1", query_orientation="forward", chrom="chrFixture"),
                 dict(candidate=second, mate="R2", query_orientation="reverse_complement", chrom="chrFixture")]
        status, example, support, complex_host = runner.select_pair_result(calls, False, [])
        self.assertEqual(status, "ambiguous")
        self.assertIsNone(example)
        self.assertFalse(support)

    def test_other_mate_discordant_chromosome_is_separate_complex_group(self):
        args, _ = fixture()
        candidate = core.call_junctions(**args)["candidates"][0]
        calls = [dict(candidate=candidate, mate="R1", query_orientation="forward", chrom="chrFixture")]
        hits = [dict(reference="otherChr", intervals0=[(1080, 1150)])]
        status, _, _, complex_host = runner.select_pair_result(calls, False, hits)
        self.assertEqual(status, "called_complex_host_evidence")
        self.assertTrue(complex_host)

    def test_fasta_reader_fetch_across_wrapped_lines(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wrapped.fa"
            sequence = DONOR[:173]
            path.write_bytes(b">test\n" + b"\n".join(sequence[i:i+60].encode() for i in range(0, len(sequence), 60)) + b"\n")
            Path(str(path)+".fai").write_text("test\t173\t6\t60\t61\n")
            reference = runner.FastaReader(path)
            try:
                self.assertEqual(reference.fetch("test", 58, 133), sequence[58:133])
                self.assertEqual(reference.fetch("test", 0, 173), sequence)
                self.assertEqual(reference.fetch("test", 173, 173), "")
            finally:
                reference.close()


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


class AggregationIntegrationTests(unittest.TestCase):
    def make_run(self, root, insertion=""):
        sample, _, _ = make_runner_fixture(root, insertion=insertion)
        proc = subprocess.run([sys.executable, str(PACKAGE / "raising_junction_pilot.py"),
                               "--root", str(root), "--sample", sample],
                              capture_output=True, text=True, timeout=30)
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        return root / "results" / (sample + "_junction_pilot_v1")

    def aggregate(self, root):
        return subprocess.run([sys.executable, str(PACKAGE / "summarize_junctions.py"),
                               "--root", str(root)], capture_output=True, text=True, timeout=30)

    def test_runner_to_distribution_preserves_denominators_and_strata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.make_run(root)
            proc = self.aggregate(root)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            output = root / "results" / "junction_pilot_summary_v1"
            qc = read_csv(output / "library_qc.csv")
            self.assertEqual(len(qc), 4)
            by_stratum = {r["stratum"]: r for r in qc}
            self.assertEqual(by_stratum["all_called"]["stratum_called_pairs"], "2")
            self.assertEqual(by_stratum["all_called"]["stratum_observed_patterns"], "1")
            self.assertEqual(by_stratum["without_prior_dd_or_complex_host_support"]["stratum_called_pairs"], "1")
            self.assertEqual(by_stratum["with_prior_dd_without_complex_host_support"]["stratum_called_pairs"], "1")
            self.assertTrue(all(r["host_supported_pairs"] == "3" for r in qc))
            distributions = read_csv(output / "figure_distributions.csv")
            mh = [r for r in distributions if r["stratum"] == "all_called" and
                  r["metric"] == "MH_length_no_intervening_sequence"]
            counts = {r["unit"]: (r["value_bp"], r["count"], r["denominator"], r["fraction"]) for r in mh}
            self.assertEqual(counts["observed_junction_patterns"], ("0", "1", "1", "1.0"))
            self.assertEqual(counts["supporting_read_pairs"], ("0", "2", "2", "1.0"))

    def test_insertion_positive_calls_are_not_mh_zero_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.make_run(root, insertion="AGC")
            proc = self.aggregate(root)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            rows = read_csv(root / "results" / "junction_pilot_summary_v1" / "figure_distributions.csv")
            mh = [r for r in rows if r["metric"] == "MH_length_no_intervening_sequence"]
            self.assertTrue(mh)
            self.assertTrue(all((r["denominator"], r["fraction"], r["value_bp"]) == ("0", "", "") for r in mh))
            insertion = [r for r in rows if r["metric"] == "intervening_sequence_length" and
                         r["stratum"] == "all_called" and r["unit"] == "supporting_read_pairs" and r["value_bp"] == "3"]
            self.assertEqual(len(insertion), 1)
            self.assertEqual((insertion[0]["count"], insertion[0]["denominator"]), ("2", "2"))

    def test_different_configuration_or_reference_cannot_be_pooled(self):
        for different in ("config", "reference"):
            with self.subTest(different=different), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                original = self.make_run(root)
                second = root / "results" / "SyntheticOther_3prm_S2_junction_pilot_v1"
                shutil.copytree(original, second)
                summary_path = second / "summary.json"
                summary = json.loads(summary_path.read_text())
                summary["sample"] = "SyntheticOther_3prm_S2"
                table = second / "junctions.tsv"
                table.write_text(table.read_text().replace("Synthetic_3prm_S1", "SyntheticOther_3prm_S2"))
                summary["output_sha256"]["junctions.tsv"] = hashlib.sha256(table.read_bytes()).hexdigest()
                if different == "config":
                    summary["config"]["min_donor_anchor"] += 1
                else:
                    summary["reference"]["sha256"] = "0" * 64
                summary_path.write_text(json.dumps(summary))
                proc = self.aggregate(root)
                self.assertNotEqual(proc.returncode, 0, proc.stdout + proc.stderr)
                self.assertIn("different caller/config/reference", proc.stderr)

    def test_duplicate_event_rejected_even_if_file_hash_is_updated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = self.make_run(root)
            table = output / "junctions.tsv"
            lines = table.read_text().splitlines()
            table.write_text("\n".join(lines + [lines[1]]) + "\n")
            summary_path = output / "summary.json"
            summary = json.loads(summary_path.read_text())
            summary["output_sha256"]["junctions.tsv"] = hashlib.sha256(table.read_bytes()).hexdigest()
            summary_path.write_text(json.dumps(summary))
            proc = self.aggregate(root)
            self.assertNotEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn("Duplicate/mismatched junction", proc.stderr)


class BatchIntegrationTests(unittest.TestCase):
    def invoke(self, root, *extra):
        return subprocess.run([sys.executable, str(PACKAGE / "raising_junction_batch.py"),
                               "--root", str(root), *extra], capture_output=True,
                              text=True, timeout=30)

    def test_list_then_process_then_reuse_exact_complete_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample, _, _ = make_runner_fixture(root)
            output = root / "results" / (sample + "_junction_pilot_v1")
            proc = self.invoke(root, "--list")
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn(sample, proc.stdout)
            self.assertFalse(output.exists())
            proc = self.invoke(root)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            original = (output / "summary.json").read_bytes()
            proc = self.invoke(root, "--skip-complete")
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn("Reusing completed pilot", proc.stdout)
            self.assertEqual((output / "summary.json").read_bytes(), original)
            aggregate = root / "results" / "junction_pilot_summary_v1" / "aggregation_summary.json"
            self.assertEqual(json.loads(aggregate.read_text())["libraries"], 1)
            summary = json.loads(original)
            summary["source_sha256"]["junction_core.py"] = "0" * 64
            (output / "summary.json").write_text(json.dumps(summary))
            proc = self.invoke(root, "--skip-complete")
            self.assertNotEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn("differs from current source/input", proc.stderr)


if __name__ == "__main__":
    unittest.main()
