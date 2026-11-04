import contextlib
import csv
import gzip
import importlib.util
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import raising_junction_batch as batch


def engine_dir(name):
    return HERE / name if (HERE / name).is_dir() else HERE.parent / name


class BatchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def manifest(self, sample="Synthetic_3prm_S1", suffix="_fourway_full_v1"):
        directory = self.root / "results" / (sample + suffix)
        directory.mkdir(parents=True)
        path = directory / "run_summary.json"
        path.write_text(json.dumps({"sample": sample, "side": "3p", "status": "complete",
            "inputs": {"lanes": [{"R1": {}, "R2": {}}]}, "validation": {"passed": True}}))
        return path

    def test_duplicate_fullrun_requires_explicit_preference(self):
        self.manifest()
        self.manifest(suffix="_fourway_full_v1_1")
        with self.assertRaisesRegex(ValueError, "Multiple full-run"):
            batch.discover(self.root)
        chosen = batch.discover(self.root, preference="v1_1")
        self.assertEqual(chosen[0]["suffix"], "_fourway_full_v1_1")
        self.assertEqual(len(chosen[0]["excluded_fullrun_directories"]), 1)

    def test_missing_requested_library_and_incomplete_manifest_fail(self):
        path = self.manifest()
        with self.assertRaisesRegex(ValueError, "No completed"):
            batch.discover(self.root, {"Missing_3prm_S1"})
        data = json.loads(path.read_text()); data["status"] = "running"
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            batch.discover(self.root)

    def test_hash_cache_detects_replaced_bytes(self):
        path = self.root / "data"
        path.write_bytes(b"old")
        cache = batch.HashCache()
        before = cache(path)
        path.write_bytes(b"new")
        self.assertNotEqual(before, cache(path))

    def test_merge_keeps_original_mh_and_marks_insertion_na(self):
        records = []
        for i, insertion in enumerate((0, 2)):
            directory = self.root / str(i); directory.mkdir()
            sample = f"Sample{i}_3prm_S1"
            tables = {
                "junction_reassessment.tsv": [{"event_id": "same_event", "insertion_length_bp": insertion,
                    "mh_length_bp": 0, "template_reference_mh_length_bp": 0,
                    "supported_pair_count_once_per_group": 1}],
                "identical_sequence_groups.tsv": [{"equivalent_template_group": "same_group", "supported_pairs": 1}],
                "template_comparisons.tsv": [{"template_id": "same_event", "relationship": "strongest_template"}]}
            for name, rows in tables.items():
                with (directory / name).open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
                    writer.writeheader(); writer.writerows(rows)
            sha = batch.HashCache()
            (directory / "refinement_summary.json").write_text(json.dumps({"sample": sample,
                "library_side": "3p", "condition": "", "biological_replicate": "",
                "output_sha256": {name: sha(directory / name) for name in tables}}))
            records.append({"refine_output": str(directory)})
        output = self.root / "merged"; output.mkdir()
        products = batch.merge_results(output, records, batch.HashCache())
        with gzip.open(output / "all_junction_reassessment.tsv.gz", "rt") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([r["mh_length_bp"] for r in rows], ["0", "0"])
        self.assertEqual([r["mh_for_analysis_bp"] for r in rows], ["0", ""])
        self.assertEqual(rows[1]["mh_assessment"], "not_assessed_insertion_interfaces")
        self.assertNotEqual(rows[0]["sample"], rows[1]["sample"])
        self.assertEqual(len(products), 3)

    def test_merge_accepts_empty_fallback_headers(self):
        record_dir = self.root / "zero"; record_dir.mkdir()
        names = {"junction_reassessment.tsv": "event_id\treview_cluster_id\tsupported_pair_count_once_per_group\n",
                 "identical_sequence_groups.tsv": "equivalent_template_group\tsupported_pairs\n",
                 "template_comparisons.tsv": "template_id\trelationship\n"}
        for name, text in names.items():
            (record_dir / name).write_text(text)
        sha = batch.HashCache()
        (record_dir / "refinement_summary.json").write_text(json.dumps({"sample": "Zero_3prm_S1", "library_side": "3p",
            "output_sha256": {name: sha(record_dir / name) for name in names}}))
        out = self.root / "merge"; out.mkdir()
        products = batch.merge_results(out, [{"refine_output": str(record_dir)}], sha)
        self.assertTrue(all(p["rows"] == 0 for p in products.values()))
        with gzip.open(out / "all_junction_reassessment.tsv.gz", "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertIn("mh_for_analysis_bp", reader.fieldnames)
            self.assertEqual(list(reader), [])

    def test_merge_rejects_nonempty_schema_mismatch(self):
        records = []
        for i in range(2):
            directory = self.root / str(i); directory.mkdir()
            path = directory / "junction_reassessment.tsv"
            path.write_text("event_id\tinsertion_length_bp\ttemplate_reference_mh_length_bp" + ("\textra" if i else "") +
                            "\nj1\t0\t0" + ("\tvalue" if i else "") + "\n")
            (directory / "refinement_summary.json").write_text(json.dumps({"sample": "sample" + str(i), "library_side": "3p",
                "output_sha256": {path.name: batch.HashCache()(path)}}))
            records.append({"refine_output": str(directory)})
        out = self.root / "merge"; out.mkdir()
        with self.assertRaisesRegex(ValueError, "Different nonempty table schemas"):
            batch.merge_results(out, records, batch.HashCache())

    def test_real_small_pipeline_resume_and_changed_output_rejection(self):
        # Reuse the engine's independently tested FASTQ/FASTA/evidence fixture,
        # then run both real command-line engines rather than mocked processes.
        refine = engine_dir("raising_junction_refine_v1_1")
        pilot = engine_dir("raising_junction_pilot_v1")
        sys.path.insert(0, str(refine))
        spec = importlib.util.spec_from_file_location("refine_fixture", refine / "tests" / "test_refine.py")
        fixture = importlib.util.module_from_spec(spec); spec.loader.exec_module(fixture)
        made = fixture.make_runner_fixture(self.root)
        sample = "Synthetic_3prm_S1"
        pdir = self.root / "results" / (sample + "_junction_pilot_v1")
        shutil.rmtree(pdir)  # Remove only the fixture's mock pilot, not any user data.
        manifest_path = self.root / "results" / (sample + "_fourway_full_v1") / "run_summary.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["validation"] = {"synthetic_fixture": True}
        manifest_path.write_text(json.dumps(manifest))
        argv = ["--root", str(self.root), "--expected-count", "1", "--pilot-script", str(pilot / "raising_junction_pilot.py"),
                "--refine-script", str(refine / "raising_junction_refine.py")]
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(batch.main(argv + ["--dry-run"]), 0)
        self.assertFalse((self.root / "results" / "junction_batch_v1").exists())
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            code = batch.main(argv)
        summaries = sorted((self.root / "results" / "junction_batch_v1").glob("*/batch_summary.json"))
        self.assertEqual(code, 0, [p.read_text() for p in summaries])
        original_summary = load = json.loads(summaries[-1].read_text())
        self.assertEqual(load["libraries"][0]["pilot_action"], "ran")
        self.assertEqual(load["libraries"][0]["refine_action"], "ran")
        self.assertEqual(len(load["merged_outputs"]), 3)
        rdir = self.root / "results" / (sample + "_junction_refine_v1_1")
        digest = batch.HashCache()(rdir / "refinement_summary.json")
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(batch.main(argv), 0)
        self.assertEqual(digest, batch.HashCache()(rdir / "refinement_summary.json"))
        summaries2 = set((self.root / "results" / "junction_batch_v1").glob("*/batch_summary.json")) - set(summaries)
        self.assertEqual(len(summaries2), 1)
        reused = json.loads(next(iter(summaries2)).read_text())
        self.assertEqual(reused["libraries"][0]["pilot_action"], "reused_after_verification")
        self.assertEqual(reused["libraries"][0]["refine_action"], "reused_after_verification")
        with (rdir / "junction_reassessment.tsv").open("a") as handle:
            handle.write("changed\n")
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(batch.main(argv), 2)
        self.assertTrue((rdir / "junction_reassessment.tsv").read_text().endswith("changed\n"))


if __name__ == "__main__":
    unittest.main()
