"""Synthetic fixture generation tests; these tests do not invoke aligners."""
import gzip
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("raising_demo_generator", ROOT / "examples" / "make_demo.py")
DEMO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DEMO)


def fastq(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    return [lines[index:index + 4] for index in range(0, len(lines), 4)]


class DemoTests(unittest.TestCase):
    def test_generation_is_byte_reproducible_and_anonymous(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "first", Path(directory) / "second"
            a, b = DEMO.generate(first), DEMO.generate(second)
            self.assertEqual(a, b)
            files = [path.relative_to(first) for path in first.rglob("*") if path.is_file()]
            self.assertEqual(len(files), 8)
            for relative in files:
                self.assertEqual(hashlib.sha256((first / relative).read_bytes()).digest(),
                                 hashlib.sha256((second / relative).read_bytes()).digest())
            self.assertIn("not_human", (first / "ref/synthetic_host.fa").read_text().splitlines()[0])

    def test_pairs_are_consistent_and_four_classes_are_designed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "demo"
            manifest = DEMO.generate(root)
            for library in manifest["libraries"]:
                self.assertEqual(list(library["expected_fourway_counts"].values()), [2, 8, 2, 2])
                base = root / "vendor_results/Raw_fastq"
                first = fastq(base / f"{library['sample']}_L001_R1_001.fastq.gz")
                second = fastq(base / f"{library['sample']}_L001_R2_001.fastq.gz")
                self.assertEqual(len(first), 14)
                self.assertEqual(len(second), 14)
                self.assertEqual(len({record[0][:-2] for record in first}), 14)
                for r1, r2 in zip(first, second):
                    self.assertEqual(r1[0][:-2], r2[0][:-2])
                    self.assertEqual(r2[1], DEMO.revcomp(r1[1]))
                    self.assertEqual(len(r1[1]), len(r1[3]))
                    self.assertEqual(set(r1[3]), {"I"})

    def test_refuses_nonempty_output_without_changing_contents(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sentinel = root / "keep.txt"
            sentinel.write_text("original")
            with self.assertRaisesRegex(ValueError, "not overwritten"):
                DEMO.generate(root)
            self.assertEqual(list(root.iterdir()), [sentinel])
            self.assertEqual(sentinel.read_text(), "original")

    def test_each_side_has_deletion_mh_and_insertion_designs(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = DEMO.generate(Path(directory) / "demo")
            for side in ("3p", "5p"):
                designs = {row["fixture"]: row for row in manifest["designs"] if row["library_side"] == side}
                self.assertEqual(designs["donor_deletion_5"]["designed_donor_loss_min_bp"], 5)
                self.assertEqual(designs["mh_3"]["designed_mh_bp"], 3)
                self.assertEqual(designs["mh_3"]["designed_donor_loss_max_bp"], 3)
                self.assertIsNone(designs["insertion_3"]["designed_mh_bp"])
                self.assertEqual(len(designs["insertion_3"]["designed_intervening_sequence"]), 3)


if __name__ == "__main__":
    unittest.main()
