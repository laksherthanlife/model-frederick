"""The ablation protocol, tested on the failure it exists to prevent.

Every real-data assertion below runs on `data/plates/`, which is tracked, so nothing here
skips on a machine that has never seen a plate reader.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import pytest

from ystwin.generator.panel_experiment import MEASURED_GROWTH_RATE_SE, OBSERVED_ACTIVITY_CV
from ystwin.mech import ablation as ab
from ystwin.mech.params import FreeScalarGateFailed, Target

# The three committed NewProtocol plates that carry a per-construct dose ladder, and the
# four built strains on each. Twelve blocks, 175 rows apiece.
EXPORTS = (
    "20260722_ER&OxidativeStress_NewProtocol_ANALYSED.xlsx",
    "20260803_ER&oxidativestress_Replicate3.xlsx",
    "20260804_ER&OxidativeStress_Replicate4.xlsx",
)
CONSTRUCTS = ("UPRE1", "UPRE2", "NativeYap1", "AlteredYap1")
BLOCKS = tuple((e, c) for e in EXPORTS for c in CONSTRUCTS)


@pytest.fixture(scope="module")
def upre1():
    return ab.reporter_block(EXPORTS[0], "UPRE1")


@pytest.fixture(scope="module")
def sweep():
    """Every block, ablated and compared out of sample once. The evidence for the numbers
    quoted in `docs/ABLATION_PROTOCOL.md` and in this module's own docstring."""
    rows = []
    for export, construct in BLOCKS:
        data = ab.reporter_block(export, construct)
        comparison = ab.nested_out_of_sample(
            ab.MaturingReporter(), ab.SaturatingReporter(), data, ab.REPORTER_ACTIVITY_FLOOR)
        result = ab.ablate(
            ab.MaturingReporter(), ab.SaturatingReporter(), ab.endpoint_fold_induction(data))
        rows.append((export, construct, comparison, result))
    return rows


# --------------------------------------------------------------------------------------
# A floor is a measurement with an assay attached, or it is refused
# --------------------------------------------------------------------------------------

def _observable(name="reporter_activity", units="fold over the zero-dose control",
                assay=None, scored="relative", data=None):
    assay = ab.REPORTER_ACTIVITY_FLOOR.assay if assay is None else assay
    data = pd.DataFrame({"time_h": [0.0, 1.0], "dose_mM": [0.0, 1.0], "signal": [1.0, 2.0]}) \
        if data is None else data
    return ab.Observable(name=name, units=units, assay=assay, scored=scored, data=data,
                         summarise=lambda model, fit: 1.0)


class TestTheFloorMustBeNamed:
    def test_the_two_measured_floors_are_the_repositorys_own_numbers(self):
        assert ab.REPORTER_ACTIVITY_FLOOR.noise_floor == OBSERVED_ACTIVITY_CV == 0.146
        assert ab.GROWTH_RATE_FLOOR.noise_floor == MEASURED_GROWTH_RATE_SE == 0.0117
        assert ab.GROWTH_RATE_FLOOR.units == "1/h"

    def test_an_unregistered_observable_is_refused_and_the_message_names_what_exists(self):
        with pytest.raises(ab.FloorNotNamed) as caught:
            ab.floor_for(_observable(name="cell_wall_thickness"))
        assert "cell_wall_thickness" in str(caught.value)
        assert "reporter_activity" in str(caught.value) and "growth_rate" in str(caught.value)

    @pytest.mark.parametrize("name", sorted(ab.REFUSED_FLOORS))
    def test_an_observable_with_no_assay_here_is_refused_with_its_reason(self, name):
        with pytest.raises(ab.FloorNotNamed) as caught:
            ab.floor_for(_observable(name=name))
        assert ab.REFUSED_FLOORS[name][:40] in str(caught.value)

    def test_and_the_product_ceilings_refusal_names_the_error_it_prevents(self):
        with pytest.raises(ab.FloorNotNamed) as caught:
            ab.floor_for(_observable(name="product_ceiling"))
        message = str(caught.value)
        assert "C2" in message and "plate-reader activity CV" in message

    def test_a_bare_number_is_not_a_floor(self, upre1):
        with pytest.raises(ab.FloorNotNamed) as caught:
            ab.ablate(ab.MaturingReporter(), ab.SaturatingReporter(),
                      ab.endpoint_fold_induction(upre1), floor=0.146)
        assert "names no assay" in str(caught.value)

    def test_a_floor_from_another_instrument_is_refused_even_when_it_is_a_real_measurement(self):
        """The specific defect: scoring one observable against another assay's noise."""
        with pytest.raises(ab.FloorNotNamed) as caught:
            ab.ablate(ab.MaturingReporter(), ab.SaturatingReporter(),
                      _observable(), floor=ab.GROWTH_RATE_FLOOR)
        message = str(caught.value)
        assert ab.GROWTH_RATE_FLOOR.assay in message
        assert ab.REPORTER_ACTIVITY_FLOOR.assay in message

    def test_a_dimensional_floor_cannot_score_a_relative_effect(self):
        pretend = Target(
            name="reporter_activity", observable="anything", units="1/h", noise_floor=0.0117,
            assay=ab.REPORTER_ACTIVITY_FLOOR.assay, source="a floor in the wrong units")
        with pytest.raises(ValueError, match="not one of"):
            ab.ablate(ab.MaturingReporter(), ab.SaturatingReporter(), _observable(), floor=pretend)

    def test_an_absolute_effect_needs_a_floor_in_its_own_units(self):
        observable = _observable(units="mg/gDCW", scored="absolute")
        with pytest.raises(ValueError, match="same units"):
            ab.ablate(ab.MaturingReporter(), ab.SaturatingReporter(), observable,
                      floor=ab.REPORTER_ACTIVITY_FLOOR)


# --------------------------------------------------------------------------------------
# A frozen incumbent is caught by measurement, not by trust
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class _FrozenModel(ab.SaturatingReporter):
    """Declares five free scalars and returns the same five numbers whatever it is shown.

    This is the shape of the incumbent four of the five 2026-09-05 rescues were scored
    against, written down so the check can be proved to fire on it.
    """

    name: str = "FrozenModel"

    def refit(self, data):
        parameters = {"R0": 1.0, "k_basal": 1.0, "k_max": 1.0, "K_dose": 1.0, "lambda_": 1.0}
        return ab.Fit(self.name, parameters, len(self.parameter_names), len(data), 0.0, True)


class _ParameterFreeModel:
    name = "ParameterFree"
    parameter_names: tuple = ()
    response_column = "signal"

    def refit(self, data):
        return ab.Fit(self.name, {}, 0, len(data), 0.0, True)

    def predict(self, data, parameters):
        return np.zeros(len(data))


class TestAFrozenModelIsRefused:
    def test_a_model_that_ignores_its_data_raises_however_it_describes_itself(self, upre1):
        with pytest.raises(ab.FrozenIncumbent) as caught:
            ab.ablate(_FrozenModel(), ab.SaturatingReporter(),
                      ab.endpoint_fold_induction(upre1))
        assert "FrozenModel" in str(caught.value)
        assert "identical parameters" in str(caught.value)

    def test_it_fires_on_whichever_side_is_frozen(self, upre1):
        with pytest.raises(ab.FrozenIncumbent, match="FrozenModel"):
            ab.ablate(ab.MaturingReporter(), _FrozenModel(),
                      ab.endpoint_fold_induction(upre1))

    def test_a_genuinely_parameter_free_model_is_allowed_and_reported_as_such(self, upre1):
        observable = ab.Observable(
            name="reporter_activity", units="fold over the zero-dose control",
            assay=ab.REPORTER_ACTIVITY_FLOOR.assay, scored="relative", data=upre1,
            summarise=lambda model, fit: 1.0 + fit.n_free)
        result = ab.ablate(ab.SaturatingReporter(), _ParameterFreeModel(), observable)
        assert result.reduced_parameter_free is True
        assert result.reduced_n_free == 0
        assert "parameter-free" in result.report()

    def test_removing_a_piece_cannot_add_a_degree_of_freedom(self, upre1):
        with pytest.raises(ValueError, match="not an ablation"):
            ab.ablate(ab.SaturatingReporter(), ab.MaturingReporter(),
                      ab.endpoint_fold_induction(upre1))


# --------------------------------------------------------------------------------------
# Both sides are refit, in the call, on the same rows
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class _RecordingModel(ab.SaturatingReporter):
    name: str = "RecordingModel"

    def __post_init__(self):
        object.__setattr__(self, "seen", [])

    def refit(self, data):
        self.seen.append(data)
        return super().refit(data)


class TestBothModelsAreRefitHere:
    def test_the_call_fits_each_side_rather_than_accepting_a_prepared_fit(self, upre1):
        full, reduced = _RecordingModel(), _RecordingModel()
        observable = ab.endpoint_fold_induction(upre1)
        ab.ablate(full, reduced, observable)
        # One fit on the data, one on the perturbed copy that proves the model responds.
        assert len(full.seen) == 2 and len(reduced.seen) == 2
        assert full.seen[0] is upre1 and reduced.seen[0] is upre1

    def test_and_the_two_sides_see_the_same_rows(self, upre1):
        full, reduced = _RecordingModel(), _RecordingModel()
        ab.ablate(full, reduced, ab.endpoint_fold_induction(upre1))
        assert full.seen[0].equals(reduced.seen[0])

    def test_the_report_carries_both_counts_the_assay_and_the_floors_source(self, upre1):
        result = ab.ablate(ab.MaturingReporter(), ab.SaturatingReporter(),
                           ab.endpoint_fold_induction(upre1))
        report = result.report()
        assert "6 free scalars" in report and "5 free scalars" in report
        assert ab.REPORTER_ACTIVITY_FLOOR.assay in report
        assert "OBSERVED_ACTIVITY_CV" in report
        assert "0.146" in report


# --------------------------------------------------------------------------------------
# Criterion (e), callable
# --------------------------------------------------------------------------------------

def _targets(n, fitted=False):
    return [Target(name=f"t{i}", observable="anything", assay="an assay", noise_floor=0.1,
                   units="relative CV", source="a source", fitted=fitted) for i in range(n)]


class TestTheFreeScalarGate:
    def test_six_scalars_against_seven_targets_passes(self):
        gate = ab.free_scalar_gate(ab.MaturingReporter(), _targets(7))
        assert gate.passes and len(gate.free) == 6 and len(gate.targets) == 7

    def test_six_against_five_refuses_and_names_the_scalars(self):
        with pytest.raises(FreeScalarGateFailed) as caught:
            ab.require_free_scalar_gate(ab.MaturingReporter(), _targets(5))
        assert "tau_mat" in str(caught.value)
        assert "without adding refutability" in str(caught.value)

    def test_a_target_the_model_was_fitted_to_does_not_count(self):
        gate = ab.free_scalar_gate(ab.SaturatingReporter(), _targets(6, fitted=True))
        assert gate.targets == () and gate.fitted_targets != () and not gate.passes

    def test_the_gate_runs_inside_the_out_of_sample_comparison(self, upre1):
        comparison = ab.nested_out_of_sample(
            ab.MaturingReporter(), ab.SaturatingReporter(), upre1, ab.REPORTER_ACTIVITY_FLOOR)
        assert comparison.gate.passes
        assert len(comparison.gate.targets) == len(comparison.groups) == 7

    def test_and_refuses_a_model_with_more_scalars_than_folds(self, upre1):
        three = upre1[upre1.dose_mM.isin([0.0, 1.0, 5.0])]
        with pytest.raises(FreeScalarGateFailed):
            ab.nested_out_of_sample(ab.MaturingReporter(), ab.SaturatingReporter(), three,
                                    ab.REPORTER_ACTIVITY_FLOOR)


# --------------------------------------------------------------------------------------
# Nesting is checked, not assumed
# --------------------------------------------------------------------------------------

class TestNesting:
    def test_two_models_with_the_same_parameters_are_not_a_nested_pair(self, upre1):
        with pytest.raises(ab.NotNested, match="strict restriction"):
            ab.nested_out_of_sample(ab.SaturatingReporter(), ab.SaturatingReporter(), upre1,
                                    ab.REPORTER_ACTIVITY_FLOOR)

    def test_the_reduced_model_is_the_full_one_at_a_point_of_its_box(self, upre1):
        """`tau_mat = 0` is not a limit at infinity but a corner of the fitted box, which is
        what makes the subset check on parameter names mean something."""
        full, reduced = ab.MaturingReporter(), ab.SaturatingReporter()
        parameters = {"R0": 900.0, "k_basal": 400.0, "k_max": 5000.0, "K_dose": 1.4,
                      "lambda_": 0.55}
        assert np.allclose(full.predict(upre1, dict(parameters, tau_mat=0.0)),
                           reduced.predict(upre1, parameters))

    def test_and_a_short_maturation_time_converges_on_it(self, upre1):
        full, reduced = ab.MaturingReporter(), ab.SaturatingReporter()
        parameters = {"R0": 900.0, "k_basal": 400.0, "k_max": 5000.0, "K_dose": 1.4,
                      "lambda_": 0.55}
        near = full.predict(upre1, dict(parameters, tau_mat=1e-6))
        assert np.allclose(near, reduced.predict(upre1, parameters), rtol=1e-4)


# --------------------------------------------------------------------------------------
# The absolute branch: a growth rate is scored in 1/h, never as a CV
# --------------------------------------------------------------------------------------

class _ExponentialGrowth:
    name = "ExponentialGrowth"
    parameter_names = ("od0", "mu")
    response_column = "od"

    def predict(self, data, parameters):
        t = np.asarray(data.time_h, dtype=float)
        return float(parameters["od0"]) * np.exp(float(parameters["mu"]) * t)

    def refit(self, data):
        from scipy.optimize import least_squares
        y = np.asarray(data.od, dtype=float)
        solution = least_squares(
            lambda x: self.predict(data, dict(zip(self.parameter_names, x))) - y,
            [y[0], 0.3], bounds=([1e-6, 1e-4], [10.0, 3.0]))
        parameters = dict(zip(self.parameter_names, (float(v) for v in solution.x)))
        rms = float(np.sqrt(np.mean((self.predict(data, parameters) - y) ** 2)))
        return ab.Fit(self.name, parameters, 2, len(data), rms, bool(solution.status > 0))


class _LogisticGrowth(_ExponentialGrowth):
    name = "LogisticGrowth"
    parameter_names = ("od0", "mu", "capacity")

    def predict(self, data, parameters):
        t = np.asarray(data.time_h, dtype=float)
        od0, mu, capacity = (float(parameters[n]) for n in self.parameter_names)
        return capacity / (1.0 + (capacity - od0) / od0 * np.exp(-mu * t))

    def refit(self, data):
        from scipy.optimize import least_squares
        y = np.asarray(data.od, dtype=float)
        solution = least_squares(
            lambda x: self.predict(data, dict(zip(self.parameter_names, x))) - y,
            [y[0], 0.3, float(y.max()) * 1.5],
            bounds=([1e-6, 1e-4, 1e-3], [10.0, 3.0, 100.0]))
        parameters = dict(zip(self.parameter_names, (float(v) for v in solution.x)))
        rms = float(np.sqrt(np.mean((self.predict(data, parameters) - y) ** 2)))
        return ab.Fit(self.name, parameters, 3, len(data), rms, bool(solution.status > 0))


@pytest.fixture(scope="module")
def logistic_od():
    """A synthetic run from a known logistic truth. Synthetic on purpose: the point is to
    exercise the absolute branch against a mu the test knows, not to measure yeast."""
    t = np.linspace(0.0, 12.0, 49)
    od0, mu, capacity = 0.05, 0.40, 1.20
    od = capacity / (1.0 + (capacity - od0) / od0 * np.exp(-mu * t))
    return pd.DataFrame({"time_h": t, "od": od})


class TestAGrowthRateIsScoredInItsOwnUnits:
    def test_deleting_the_capacity_term_moves_mu_by_more_than_the_growth_assays_se(
            self, logistic_od):
        observable = ab.Observable(
            name="growth_rate", units="1/h", assay=ab.GROWTH_RATE_FLOOR.assay,
            scored="absolute", data=logistic_od,
            summarise=lambda model, fit: fit.parameters["mu"],
            description="fitted specific growth rate")
        result = ab.ablate(_LogisticGrowth(), _ExponentialGrowth(), observable)
        assert result.floor is ab.GROWTH_RATE_FLOOR
        assert result.scored == "absolute"
        assert result.effect > MEASURED_GROWTH_RATE_SE
        assert result.clears_floor
        assert "1/h" in result.report()

    def test_and_the_full_models_mu_is_the_truth_it_was_generated_from(self, logistic_od):
        fitted = _LogisticGrowth().refit(logistic_od)
        assert fitted.parameters["mu"] == pytest.approx(0.40, rel=1e-3)


# --------------------------------------------------------------------------------------
# Tier 0.3, run: the maturation state against the incumbent's algebra, on data in hand
# --------------------------------------------------------------------------------------

class TestTheMaturationStateOnTheCommittedPlates:
    def test_every_block_is_the_shape_the_protocol_expects(self):
        for export, construct in BLOCKS:
            data = ab.reporter_block(export, construct)
            assert len(data) == 175
            assert sorted(pd.unique(data.dose_mM)) == pytest.approx(
                [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0 if construct.startswith("UPRE") else 4.0])

    def test_the_ablation_clears_the_floor_on_none_of_the_twelve(self, sweep):
        assert [r.clears_floor for _, _, _, r in sweep] == [False] * 12

    def test_and_its_largest_effect_is_a_fifth_of_the_floor(self, sweep):
        worst = max(r.ratio for _, _, _, r in sweep)
        assert worst == pytest.approx(0.2038, abs=5e-3)
        assert worst < 1.0

    def test_out_of_sample_the_extra_state_never_earns_its_parameter(self, sweep):
        ratios = [c.ratio(c.interior_groups) for _, _, c, _ in sweep]
        assert max(ratios) == pytest.approx(0.1350, abs=5e-3)
        assert max(ratios) < 1.0

    def test_and_on_two_of_the_twelve_it_predicts_worse_than_the_algebra(self, sweep):
        losses = [(e, c) for e, c, comparison, _ in sweep
                  if comparison.delta(comparison.interior_groups) <= 0.0]
        assert len(losses) == 2
        assert all(construct.endswith("Yap1") for _, construct in losses)

    def test_the_h2o2_blocks_are_a_comparison_between_two_models_that_both_miss(self, sweep):
        """Not a finding about the maturation state. The 2 and 4 mM peroxide wells go
        negative after the dilution correction -- `panel_experiment.PanelDataset.n_unusable`
        records exactly this -- and neither model describes them."""
        by_construct = {}
        for _, construct, comparison, _ in sweep:
            by_construct.setdefault(construct, []).append(comparison.mean_error()[0])
        dtt = [e for c in ("UPRE1", "UPRE2") for e in by_construct[c]]
        peroxide = [e for c in ("NativeYap1", "AlteredYap1") for e in by_construct[c]]
        assert max(dtt) < 0.25
        assert min(peroxide) > 0.60

    def test_the_fitted_maturation_time_is_not_identified_by_these_traces(self):
        """The measurement behind experiment 1.2. A 38x spread across replicates of one
        construct is why `tau_mat` is reported as fitted and never as a maturation time."""
        fitted = {}
        for export, construct in BLOCKS:
            data = ab.reporter_block(export, construct)
            fitted[(export, construct)] = ab.MaturingReporter().refit(data).parameters["tau_mat"]
        spread = max(fitted.values()) / min(fitted.values())
        assert spread > 30.0
        upre1 = [v for (_, c), v in fitted.items() if c == "UPRE1"]
        assert max(upre1) / min(upre1) > 5.0

    def test_the_report_says_which_folds_are_interpolation(self, sweep):
        _, _, comparison, _ = sweep[0]
        assert len(comparison.interior_groups) == 5
        assert "interpolation rather than extrapolation" in comparison.report()


class TestTheDataIsReadFromTheCommittedText:
    def test_a_construct_that_is_not_on_the_plate_is_a_key_error_naming_the_ones_that_are(self):
        with pytest.raises(KeyError, match="UPRE1"):
            ab.reporter_block(EXPORTS[0], "HSE-heat")

    def test_a_dose_that_is_not_on_the_ladder_is_refused(self, upre1):
        with pytest.raises(ValueError, match="not on this ladder"):
            ab.endpoint_fold_induction(upre1, dose=3.0)

    def test_a_frame_without_the_three_columns_says_which_are_missing(self):
        with pytest.raises(KeyError, match="time_h"):
            ab.SaturatingReporter().refit(pd.DataFrame({"t": [0.0], "y": [1.0]}))
