"""The one number carrying the whole stress-to-metabolism coupling.

`latent_bridge` turns an inferred stress state into a single non-growth ATP maintenance
demand, so one coefficient converts promoter activity into ATP. It used to be
`_MAINTENANCE_PER_ACTIVITY = 300.0`, described as "mmol ATP/gDW/h per unit promoter
activity" -- and it multiplied an activity in RFU/OD/h, which is not a unit of anything
per ATP. On the original plate snapshot it demanded 316,000-391,000 mmol ATP/gDW/h
against a resting 0.7. The current measured excess reaches ~1457 RFU/OD/h. Yeast9 is
infeasible at 19.19 and growth is already zero there.

So the branch was not merely unvalidated, it was unrunnable on real data, which is why
feeding it a real activity had never been tried. These tests keep it runnable.
"""

import pytest

from ystwin.bridge.latent_bridge import (
    _ACTIVITY_FULL_SCALE,
    _MAX_STRESS_MAINTENANCE,
    _RESTING_MAINTENANCE,
    LatentState,
    latent_constraints,
)
from ystwin.bridge.physiology_bridge import PhysiologicalConstraints

# Where Yeast9 stops growing at all, at its default glucose bound. Bisected, not asserted.
YEAST9_INFEASIBLE_NGAM = 19.19


@pytest.fixture(scope="module")
def measured_max_excess():
    """Measured input from the existing independent derivation, not another scale literal."""
    return TestTheFullScaleIsStillWhatThePlatesSay._derive()


@pytest.fixture
def physical():
    return PhysiologicalConstraints(
        growth_lower_bound=0.0, glucose_uptake=1.5, oxygen_uptake=2.0, time_h=8.0,
    )


def _ngam(activity, physical, **kw):
    state = LatentState(promoter_activity=activity, construct="UPRE1", time_h=8.0)
    return latent_constraints(
        state, physical, acknowledge_unvalidated=True, **kw
    ).atp_maintenance


class TestTheDemandStaysInsideWhatTheModelCanCarry:
    def test_no_activity_leaves_resting_maintenance_alone(self, physical):
        assert _ngam(0.0, physical) == pytest.approx(_RESTING_MAINTENANCE)

    def test_the_largest_measured_activity_is_still_feasible(self, physical, measured_max_excess):
        """The regression that matters: a measured input, not the production scale as input."""
        assert _ngam(measured_max_excess, physical) < YEAST9_INFEASIBLE_NGAM

    def test_an_absurd_activity_is_clipped_not_extrapolated(self, physical):
        """Extrapolating past the measured range is exactly where the old constant did its
        damage, so the fraction is capped at full scale."""
        assert _ngam(1e9, physical) == pytest.approx(
            _RESTING_MAINTENANCE + _MAX_STRESS_MAINTENANCE)

    def test_the_ceiling_is_below_where_yeast9_stops_growing(self):
        assert _RESTING_MAINTENANCE + _MAX_STRESS_MAINTENANCE < YEAST9_INFEASIBLE_NGAM

    def test_the_old_coefficient_would_not_be(self, measured_max_excess):
        """States the defect as arithmetic so the fix cannot silently revert."""
        assert 0.7 + 300.0 * measured_max_excess > 1e5


class TestItIsDimensionallyCoherent:
    def test_the_coefficient_is_in_atp_units_not_atp_per_rfu(self):
        """Doubling the activity scale halves the demand at a given activity, which is
        what makes the coefficient an ATP quantity rather than a conversion factor."""
        half = _ngam(500.0, physical=PhysiologicalConstraints(0.0, 1.5, 2.0, 8.0),
                     activity_full_scale=_ACTIVITY_FULL_SCALE * 2)
        full = _ngam(500.0, physical=PhysiologicalConstraints(0.0, 1.5, 2.0, 8.0))

        assert (half - _RESTING_MAINTENANCE) == pytest.approx(
            (full - _RESTING_MAINTENANCE) / 2)

    def test_a_zero_scale_is_refused_rather_than_dividing(self, physical):
        with pytest.raises(ValueError, match="divisor"):
            _ngam(500.0, physical, activity_full_scale=0.0)

    def test_activity_below_basal_costs_nothing(self, physical):
        assert _ngam(1e-6, physical) == pytest.approx(_RESTING_MAINTENANCE)

    def test_more_stress_never_costs_less(self, physical, measured_max_excess):
        demands = [_ngam(a, physical) for a in (0.0, 100.0, 500.0, 1000.0, measured_max_excess)]

        assert demands == sorted(demands)


class TestTheRefusalSurvivedTheFix:
    def test_it_still_refuses_without_acknowledgement(self, physical):
        """G4 returned INCONCLUSIVE for all four constructs. Fixing the arithmetic does
        not validate the mapping, and the gate is the only thing saying so."""
        state = LatentState(promoter_activity=500.0, construct="UPRE1", time_h=8.0)

        with pytest.raises(ValueError, match="not validated"):
            latent_constraints(state, physical)

    def test_and_the_refusal_travels_on_the_result(self, physical):
        state = LatentState(promoter_activity=500.0, construct="UPRE1", time_h=8.0)
        got = latent_constraints(state, physical, acknowledge_unvalidated=True)

        assert "REFUSED by G4" in got.validation


@pytest.mark.integration
class TestTheFullScaleIsStillWhatThePlatesSay:
    """`_ACTIVITY_FULL_SCALE` is a DERIVED number living as a literal, and that is the
    shape that goes stale without anyone noticing.

    It sat at 1305.0 from 2026-08-26 to 2026-09-04 while the artifact it was read off had
    moved to 1465.0 -- a 12.3% error in a denominator, reached by nothing that would fail.
    `src/` must not import from `outputs/`, so the constant cannot simply be computed at
    import time; what it can have is a test that recomputes it and refuses to agree.

    If this fails, verify the current `scripts/run_sensor_characterisation.py` output,
    then set `_ACTIVITY_FULL_SCALE` from this independent derivation. The measured-input
    fixture also uses `_derive()`; there is no second full-scale literal to update.
    """

    @staticmethod
    def _derive():
        import pandas as pd

        from ystwin import paths

        table = paths.outputs_dir() / "sensor_characterisation.csv"
        if not table.exists():
            pytest.skip(f"{table} not generated; run scripts/run_sensor_characterisation.py")
        frame = pd.read_csv(table)
        # Excess is per CONSTRUCT, against that construct's own unstressed wells: the four
        # constructs have different basal brightness, so a single pooled baseline would
        # subtract one construct's floor from another's ceiling.
        baseline = frame[frame["dose_mM"] == 0].groupby("construct")["activity_late"].median()
        excess = frame.join(baseline.rename("baseline"), on="construct")
        return float((excess["activity_late"] - excess["baseline"]).max())

    def test_the_constant_matches_the_artifact_it_was_read_off(self):
        derived = self._derive()

        assert _ACTIVITY_FULL_SCALE == pytest.approx(derived, abs=0.05), (
            f"_ACTIVITY_FULL_SCALE = {_ACTIVITY_FULL_SCALE} but "
            f"outputs/sensor_characterisation.csv now gives {derived:.4g}. It is a snapshot "
            f"of that artifact, so re-run scripts/run_sensor_characterisation.py and set the "
            f"constant to the new value -- do not edit this test's expectation")

    def test_the_measured_full_scale_reaches_the_maintenance_ceiling(
            self, physical, measured_max_excess):
        """The independent measured excess must exercise full scale on the default mapping."""
        assert _ngam(measured_max_excess, physical, basal_activity=0.0) == pytest.approx(
            _RESTING_MAINTENANCE + _MAX_STRESS_MAINTENANCE)

    def test_the_full_scale_is_an_excess_and_not_a_raw_reading(self):
        """A guard on the derivation itself. The raw `activity_late` maximum is roughly
        twice the excess, and using it would silently halve every stress fraction."""
        import pandas as pd

        from ystwin import paths

        table = paths.outputs_dir() / "sensor_characterisation.csv"
        if not table.exists():
            pytest.skip("sensor_characterisation.csv not generated")
        raw = float(pd.read_csv(table)["activity_late"].max())

        assert raw > 1.5 * _ACTIVITY_FULL_SCALE
