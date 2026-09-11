"""The late window is a choice, so the folds computed inside it have to be swept.

``int(len(tw) * 0.75)`` decides where a trace counts as settled, and every reported fold
and growth rate is computed after that index. A fold near 1.0 at 0.75 and nowhere else is
not the same finding as one near 1.0 across the window, and only the sweep separates them.

These tests exercise the sweep's own logic on constructed frames. They cannot run the
real plates: the four NewProtocol exports are not in the repository, and the script
refuses rather than writing an empty table when they are absent.
"""

import importlib.util
import pathlib

import numpy as np
import pandas as pd
import pytest
from ystwin.readings import CorrectedOD

_SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "run_sensor_characterisation.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("run_sensor_characterisation", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestTheWindowIsAParameterNotALiteral:
    def test_collect_takes_the_fraction(self, script):
        """It was a literal inside a loop, which is why it was never swept."""
        import inspect

        assert "late_fraction" in inspect.signature(script.collect).parameters

    def test_the_default_matches_the_reported_window(self, script):
        import inspect

        default = inspect.signature(script.collect).parameters["late_fraction"].default

        assert default == script.LATE_WINDOW_FRACTION == 0.75

    def test_the_sweep_brackets_the_default(self, script):
        """A sweep that only moves one way would report a trend, not a sensitivity."""
        assert min(script.LATE_WINDOW_SWEEP) < script.LATE_WINDOW_FRACTION
        assert max(script.LATE_WINDOW_SWEEP) > script.LATE_WINDOW_FRACTION
        assert script.LATE_WINDOW_FRACTION in script.LATE_WINDOW_SWEEP

    def test_every_fraction_leaves_points_in_the_window(self, script):
        """At the widest end the window still has to be a window."""
        for fraction in script.LATE_WINDOW_SWEEP:
            assert 0.0 < fraction < 1.0
            # The narrowest real trace that survives growth_window carries 10 points.
            assert int(10 * fraction) < 10


class TestTheFoldIsReadFromTheRightAttribute:
    """The sweep crashed the first time it met real data: `FoldChange` is a dataclass, not
    a float subclass like `PowerEstimate`, so `float(fold)` raises. Cheap to assert, and
    the assumption is the kind that only fails once the plates are present."""

    def test_foldchange_is_not_a_float_subclass(self):
        from ystwin.analysis.uncertainty import FoldChange

        assert not issubclass(FoldChange, float)

    def test_it_carries_the_point_estimate_and_the_interval(self):
        import dataclasses

        from ystwin.analysis.uncertainty import FoldChange

        names = {f.name for f in dataclasses.fields(FoldChange)}

        assert {"point", "low", "high", "n_plates", "estimable"} <= names


class TestTheSweepReportsSpanNotJustAValue:
    @staticmethod
    def _sweep_frame(folds_by_fraction):
        return pd.DataFrame([
            {"late_fraction": f, "construct": "UPRE1", "stressor": "DTT",
             "dose_mM": 2.0, "fold": v, "estimable": True}
            for f, v in folds_by_fraction.items()
        ])

    def test_a_fold_that_crosses_one_is_reported_as_straddling(self, script):
        """The call that matters: 0.97 at the default window is not "no induction" if the
        same fold reads 1.30 at a wider one."""
        frame = self._sweep_frame({0.5: 1.30, 0.6: 1.12, 0.7: 1.02,
                                   0.75: 0.97, 0.8: 0.93, 0.9: 0.88})

        lo, hi, verdict = script.window_verdict(frame.fold)

        assert (lo, hi) == pytest.approx((0.88, 1.30))
        assert verdict == "straddles 1.0 across the sweep"

    def test_a_fold_that_stays_below_one_keeps_its_sign(self, script):
        frame = self._sweep_frame({0.5: 0.71, 0.6: 0.74, 0.7: 0.77,
                                   0.75: 0.79, 0.8: 0.81, 0.9: 0.85})

        lo, hi, verdict = script.window_verdict(frame.fold)

        assert verdict == "keeps its sign"
        assert hi < 1.0

    def test_the_verdict_ignores_a_missing_fold_rather_than_returning_nan(self, script):
        """One unestimable window must not erase the span the others establish."""
        frame = self._sweep_frame({0.5: 1.30, 0.6: float("nan"), 0.75: 0.97})

        lo, hi, _ = script.window_verdict(frame.fold)

        assert (lo, hi) == pytest.approx((0.97, 1.30))

    def test_an_empty_sweep_is_reported_not_crashed(self, script, capsys):
        """No fold estimable at any window is the honest outcome on two plates, and it
        must not look like a clean pass."""
        out = script.report_late_window_sensitivity([])

        assert out.empty
        assert "nothing to sweep" in capsys.readouterr().out


class TestTheWindowedQuantitiesActuallyMove:
    def test_a_later_window_sees_a_different_growth_rate(self):
        """Guards the premise. If the window made no difference there would be nothing to
        sweep, and a trace with curvature is exactly where it does."""
        from ystwin.growth import specific_growth_rate

        t = np.linspace(0.0, 24.0, 60)
        # Growth that slows: mu falls across the trace, so window choice changes mu_late.
        od = 0.05 * np.exp(0.35 * t / (1.0 + t / 12.0))
        mu = specific_growth_rate(t, CorrectedOD(od))
        early = float(np.nanmean(mu[slice(int(len(t) * 0.5), None)]))
        late = float(np.nanmean(mu[slice(int(len(t) * 0.9), None)]))

        assert early > late
