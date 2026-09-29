"""Small synthetic QC checks; real numerical comparison is stage01 --verify."""
import importlib.util
from pathlib import Path
import unittest
import contextlib
import io
from types import SimpleNamespace
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("stage01", Path(__file__).resolve().parents[1] / "scripts" / "01_Data_Quality_and_Preparation.py")
stage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage)


class QCChecks(unittest.TestCase):
    def test_spatial_flags_do_not_modify_observations(self):
        frame = pd.DataFrame(dict(borehole_id=list("abcd"), xcoord=[0, 1, 2, 3],
                                  ycoord=[0, 0, 0, 0], z=[500., 500., 500., 550.]))
        before = frame.copy(deep=True)
        report = stage.Report()
        result = stage.local_plausibility(frame, "z", 10, 3, 20, report, "fixture", "raw")
        pd.testing.assert_frame_equal(frame, before)
        self.assertEqual(result.iloc[-1].local_reference, 500)
        self.assertEqual(result.iloc[-1].local_residual, 50)
        self.assertEqual(len(report.rows), 1)
        self.assertEqual(report.rows[0]["action"], "flag_for_review")

    def test_neighbour_support_excludes_same_borehole(self):
        frame = pd.DataFrame(dict(borehole_id=["a", "a", "b"], xcoord=[0, 0, 1],
                                  ycoord=[0, 0, 0], z=[500., 500., 501.]))
        result = stage.local_plausibility(frame, "z", 10, 2, 10, stage.Report(), "fixture", "raw")
        self.assertTrue(result.local_reference.isna().all())
        self.assertTrue(result.neighbour_count.eq(1).all())

    def test_interval_defects_are_reported_without_sorting_input(self):
        layers = pd.DataFrame(dict(ObjektID=["a"]*4, Obergrenz=[2., 0., 1., 4.],
            Untergrenz=[3., .5, 2.5, 3.5], DIN=["S", None, "?", "TU"], PetBez=["sand", None, "unknown", "clay"]))
        before = layers.copy(deep=True)
        collars = pd.DataFrame(dict(borehole_id=["a"], Endteufe=[5.]))
        report = stage.Report()
        stage.audit_intervals(layers, collars, stage.defaults(), report, "fixture")
        codes = {r["issue_code"] for r in report.rows}
        self.assertTrue({"gap", "overlap", "negative_or_zero_thickness_or_depth", "end_depth_mismatch",
                         "original_interval_order", "missing_DIN", "unknown_DIN_legacy_default_10"} <= codes)
        pd.testing.assert_frame_equal(layers, before)

    def test_numerical_comparison_catches_nan_and_value_changes(self):
        self.assertTrue(stage.compare_array("a", [1., np.nan], [1., np.nan])["passed"])
        self.assertFalse(stage.compare_array("a", [1., np.nan], [1., 0.])["passed"])
        self.assertFalse(stage.compare_array("a", [1.], [1.0000000001])["passed"])

    def test_legacy_swallowed_exceptions_reach_report(self):
        def qnd(*args, **kwargs):
            print("\n VTK - Error in:\n \n['hole-a', 'hole-b']")
            return ("unchanged",)
        u = SimpleNamespace(qnd_compositing=qnd, len_wei_compositing=lambda *a, **kw: "unchanged-tertiary")
        report = stage.Report()
        with contextlib.redirect_stdout(io.StringIO()):
            result = stage.geological_preprocessing(None, None, stage.defaults(), u, report)
        self.assertEqual(result, (("unchanged",), "unchanged-tertiary"))
        self.assertEqual([r["borehole_id"] for r in report.rows], ["hole-a", "hole-b"])
        self.assertTrue(all(r["action"] == "legacy_compositing_skipped" for r in report.rows))


if __name__ == "__main__":
    unittest.main()
