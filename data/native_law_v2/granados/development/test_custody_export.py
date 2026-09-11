"""Synthetic-only tests of the custodian's data adapter and release boundary.

No native measurement file, fitter, prediction, or score is loaded or executed.
"""
import copy
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[4]
SPEC = importlib.util.spec_from_file_location("granados_custody", ROOT / "scripts/export_granados_training.py")
CUSTODY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CUSTODY)


def synthetic_source():
    groups = {}
    for rep, rows in CUSTODY.DEV_ROWS.items():
        groups[rep] = {
            "origin": 49,
            "max5": [[2.0 + cell for _ in range(97)] for cell in range(rows)],
            "median": [[1.0 + cell for _ in range(97)] for cell in range(rows)],
            "times": [[2.5 * (frame - 49) for frame in range(1, 98)] for _ in range(rows)],
        }
    groups["rep1"]["max5"][0][27] = None
    return groups


class CustodyAdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.protocol = CUSTODY.dev_strict_json((CUSTODY.DEV_ROOT / "contracts/development_protocol.json").read_bytes())
        cls.source = synthetic_source()
        cls.train = CUSTODY.dev_projection_payload(cls.source, "training_values")
        cls.development = CUSTODY.dev_projection_payload(cls.source, "development_inputs")

    def test_exact_protocol_record_counts_and_scalar_replay(self):
        train = CUSTODY.dev_verify_payload(self.train, self.source, "training_values")
        dev = CUSTODY.dev_verify_payload(self.development, self.source, "development_inputs")
        self.assertEqual(train["source_cell_records"], 17292)
        self.assertEqual(train["source_time_records"], 0)
        self.assertEqual(train["source_scalar_checks"], 51876)
        self.assertEqual(dev["source_cell_records"], 7040)
        self.assertEqual(dev["source_time_records"], 8448)
        self.assertEqual(dev["source_scalar_checks"], 29568)
        self.assertTrue(all(row["window"] == "prefix" for row in self.development["source_cells"]))
        self.assertTrue(all(set(row) == {"group_id", "cell_index", "frame_1based", "relative_time_min", "time_pointer"} for row in self.development["source_times"]))

    def test_embedded_projection_contracts(self):
        for phase, payload in (("training_values", self.train), ("development_inputs", self.development)):
            fixture = {"document_type": "source_bound_development_projection", "protocol_id": CUSTODY.DEV_PROTOCOL_ID,
                       "phase": phase, "approval_sha256": "0" * 64, "source_asset_id": CUSTODY.DEV_SOURCE_ID,
                       "source_sha256": CUSTODY.DEV_SOURCE_SHA, **payload,
                       "source_binding_audit": {"path": "synthetic-fixture-only.json", "sha256": "0" * 64, "bytes": 1, "role": "source_binding_audit"}}
            CUSTODY.dev_contract(self.protocol, "Projection").validate(fixture)

    def test_null_is_preserved_and_zero_substitution_refused(self):
        self.assertIsNone(self.train["source_cells"][0]["max5"])
        bad = copy.deepcopy(self.train)
        bad["source_cells"][0]["max5"] = 0.0
        with self.assertRaises(CUSTODY.CustodyError):
            CUSTODY.dev_verify_payload(bad, self.source, "training_values")

    def test_tampered_values_and_shuffled_rows_refused(self):
        bad = copy.deepcopy(self.train)
        bad["source_cells"][1]["max5"] += 1.0
        with self.assertRaises(CUSTODY.CustodyError):
            CUSTODY.dev_verify_payload(bad, self.source, "training_values")
        bad = copy.deepcopy(self.train)
        bad["source_cells"][1], bad["source_cells"][2] = bad["source_cells"][2], bad["source_cells"][1]
        with self.assertRaises(CUSTODY.CustodyError):
            CUSTODY.dev_verify_payload(bad, self.source, "training_values")

    def test_disguised_development_response_refused(self):
        bad = copy.deepcopy(self.development)
        bad["source_times"][0]["max5"] = 7.0
        with self.assertRaises(CUSTODY.CustodyError):
            CUSTODY.dev_verify_payload(bad, self.source, "development_inputs")
        bad = copy.deepcopy(self.development)
        bad["source_times"][0]["relative_time_min"] = 7.0
        with self.assertRaises(CUSTODY.CustodyError):
            CUSTODY.dev_verify_payload(bad, self.source, "development_inputs")
        with self.assertRaises(CUSTODY.CustodyError):
            CUSTODY.dev_projection_payload(self.source, "development_responses")

    def test_unapproved_groups_fields_windows_and_pointers_refused(self):
        mutations = (
            ("group_id", CUSTODY.DEV_SOURCE_ID + ":/rep6"),
            ("frame_1based", 49),
            ("max5_pointer", "/rep1/GFP/nucLoc/0/27"),
            ("cell_index", 1),
        )
        for key, value in mutations:
            bad = copy.deepcopy(self.train)
            bad["source_cells"][0][key] = value
            with self.assertRaises(CUSTODY.CustodyError):
                CUSTODY.dev_verify_payload(bad, self.source, "training_values")

    def test_incomplete_or_duplicate_inventory_refused(self):
        bad = copy.deepcopy(self.train)
        bad["source_cells"].pop()
        with self.assertRaises(CUSTODY.CustodyError):
            CUSTODY.dev_verify_payload(bad, self.source, "training_values")
        bad = copy.deepcopy(self.train)
        bad["source_cells"][1] = copy.deepcopy(bad["source_cells"][0])
        with self.assertRaises(CUSTODY.CustodyError):
            CUSTODY.dev_verify_payload(bad, self.source, "training_values")

    def test_strict_json(self):
        for content in ('{"x": 1, "x": 2}', '{"x": NaN}', '{"x": Infinity}', '{"x": 1e999}'):
            with self.assertRaises(CUSTODY.CustodyError):
                CUSTODY.dev_strict_json(content)
        self.assertFalse(CUSTODY.dev_same_scalar(None, 0))
        self.assertFalse(CUSTODY.dev_same_scalar(True, 1))

    def test_mean_of_ratios_prefix_boundary_and_attrition(self):
        group = copy.deepcopy(self.source["rep1"])
        for row in range(CUSTODY.DEV_ROWS["rep1"]):
            group["max5"][row] = [2.0 if row == 0 else 9.0] * 97
            group["median"][row] = [1.0 if row == 0 else 3.0] * 97
            if row >= 2:
                group["max5"][row][27] = None
        group["times"][0][49], group["times"][1][49] = 2.0, 4.0
        group["median"][1][49] = 0.0
        group["max5"][0][50] = group["max5"][1][50] = None
        result = CUSTODY.dev_aggregate_group("rep1", group)
        self.assertEqual(result["eligible_cell_indices"], [0, 1])
        self.assertEqual(result["mu"], 2.5)
        self.assertAlmostEqual(result["s"], 0.5 ** 0.5)
        self.assertEqual(result["frames"][0]["time_min"], 3.0)
        self.assertEqual(result["frames"][0]["observed_mean"], 2.0)
        self.assertEqual(result["frames"][0]["invalid_cells"], 1)
        self.assertIsNone(result["frames"][0]["cell_sample_sd"])
        self.assertIsNone(result["frames"][1]["observed_mean"])
        self.assertEqual(result["frames"][1]["source_missing_cells"], 2)
        self.assertEqual(result["frames"][2]["observed_mean"], 2.5)
        self.assertEqual(len(result["frames"]), 24)
        scaled = copy.deepcopy(group)
        for field in ("max5", "median"):
            scaled[field] = [[None if value is None else 2.0 * value for value in row] for row in scaled[field]]
        scaled_result = CUSTODY.dev_aggregate_group("rep1", scaled)
        self.assertEqual(result["eligible_cell_indices"], scaled_result["eligible_cell_indices"])
        self.assertEqual(result["mu"], scaled_result["mu"])
        self.assertEqual([frame["observed_mean"] for frame in result["frames"]], [frame["observed_mean"] for frame in scaled_result["frames"]])
        poison = copy.deepcopy(group)
        for frame in CUSTODY.DEV_RESPONSE:
            for row in range(CUSTODY.DEV_ROWS["rep1"]):
                poison["max5"][row][frame - 1] = 10000.0
                poison["median"][row][frame - 1] = 1.0
        poisoned = CUSTODY.dev_aggregate_group("rep1", poison)
        self.assertEqual((result["eligible_cell_indices"], result["mu"], result["s"]),
                         (poisoned["eligible_cell_indices"], poisoned["mu"], poisoned["s"]))

    def test_writer_resume_and_path_escape(self):
        with tempfile.TemporaryDirectory(dir=CUSTODY.DEV_ROOT, prefix="synthetic-custody-test-") as directory:
            path = Path(directory) / "fixture.json"
            first = CUSTODY.dev_save(path, {"synthetic_fixture": True})
            self.assertEqual(first, CUSTODY.dev_save(path, {"synthetic_fixture": True}, resume_identical=True))
            with self.assertRaises(CUSTODY.CustodyError):
                CUSTODY.dev_save(path, {"synthetic_fixture": False}, resume_identical=True)
        with self.assertRaises(CUSTODY.CustodyError):
            CUSTODY.dev_save(CUSTODY.DEV_ROOT / "../outside-owner-fixture.json", {"synthetic_fixture": True})


if __name__ == "__main__":
    unittest.main()
