"""Synthetic-only scope/window tests; no native data, fits or scores."""
import copy
import json
from pathlib import Path
import runpy
import unittest

MODULE = runpy.run_path(str(Path(__file__).with_name("release_phase_c.py")), run_name="phase_c_parser_test")
Reader = MODULE["SelectedResponseReader"]
Error = MODULE["P"]["RefreshError"]


def specimen():
    result = {"rep6": {"GFP": {"max5": [[float("nan")]]}, "ignored_string": "braces [ { } ] and escaped quote \""}}
    for rep in ("rep4", "rep5"):
        matrix = [[float(frame) for frame in range(1, 98)] for _ in range(2)]
        # Unapproved columns/fields contain a poison token to demonstrate they
        # are never numerically decoded by this selective reader. Production
        # separately requires the exact already-validated whole-file checksum.
        matrix[0][0] = float("nan")
        result[rep] = {"GFP": {"max5": copy.deepcopy(matrix), "median": copy.deepcopy(matrix), "nucLoc": [[float("nan")]]},
                       "general": {"times": copy.deepcopy(matrix), "origin": 49}}
    return result


class PhaseCParserTests(unittest.TestCase):
    def test_only_authorized_groups_fields_and_frames_are_decoded(self):
        selected = Reader(json.dumps(specimen()), {"rep4": 2, "rep5": 2}).read()
        self.assertEqual(set(selected), {"rep4", "rep5"})
        for rep in selected:
            self.assertEqual(set(selected[rep]["GFP"]), {"max5", "median"})
            for matrix in (selected[rep]["GFP"]["max5"], selected[rep]["GFP"]["median"], selected[rep]["general"]["times"]):
                self.assertEqual(matrix, [list(map(float, range(50, 74)))] * 2)

    def test_null_is_preserved(self):
        data = specimen()
        data["rep4"]["GFP"]["max5"][0][49] = None
        selected = Reader(json.dumps(data), {"rep4": 2, "rep5": 2}).read()
        self.assertIsNone(selected["rep4"]["GFP"]["max5"][0][0])

    def test_nonfinite_selected_scalar_is_rejected(self):
        data = specimen()
        data["rep4"]["GFP"]["max5"][0][49] = float("nan")
        with self.assertRaises(Error):
            Reader(json.dumps(data), {"rep4": 2, "rep5": 2}).read()

    def test_changed_shape_or_origin_is_rejected(self):
        data = specimen()
        data["rep4"]["GFP"]["median"][0].pop()
        with self.assertRaises(Error):
            Reader(json.dumps(data), {"rep4": 2, "rep5": 2}).read()
        data = specimen()
        data["rep5"]["general"]["origin"] = 48
        with self.assertRaises(Error):
            Reader(json.dumps(data), {"rep4": 2, "rep5": 2}).read()

    def test_duplicate_selected_key_is_rejected(self):
        text = json.dumps(specimen()).replace('"origin": 49', '"origin": 49, "origin": 49', 1)
        with self.assertRaises(Error):
            Reader(text, {"rep4": 2, "rep5": 2}).read()


if __name__ == "__main__":
    unittest.main()
