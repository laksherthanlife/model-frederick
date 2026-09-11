"""What the chain can actually answer: which products, and through which channel.

The project's stated goal is environment conditions in, titre out, for any product. These
tests measure how far that currently goes, and they exist because nothing had ever run the
sweep end to end and reported what came back.

They are deliberately written to FAIL when the answer improves. A second calibrated product,
or an environment channel that reaches the product without passing through mu, breaks them --
and that is the point: each one is a claim about a limit, so it should stop being true the
moment the limit moves. `scripts/product_environment_sweep.py` regenerates the tables.

The one channel is mu, and it is THRESHOLD-GATED. A stressor reaches the product only once
it has taken mu_max below the held setpoint; under that dose the content is bit-identical and
over it the incumbent refuses outright. So the tests here come in pairs and every name says
which side of the crossing it is on -- "no effect" is a statement about a dose, not an axis.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import pathlib

import pandas as pd
import pytest

from ystwin.generator.context import CultureContext, context_growth_rate
from ystwin.generator.stress_panel import STRESSORS, viability
from ystwin.pathway import calibrations
from ystwin.pathway.solve import NodeKinetics, solve_pathway
from ystwin.pathway.spec import Fate, RateLaw, available_pathways, load_pathway
from ystwin.predict import (
    Environment,
    Genotype,
    SetpointUnreachable,
    predict_product,
)

KINETICS = calibrations.BETA_CAROTENE_KINETICS
HELD_RATE = 0.18

#: The growth-rate window the carotenoid kinetics were fitted over. A held rate outside it
#: is refused by `pathway/solve.py`, so it also bounds how far the seam can move content.
FIT_FLOOR, FIT_TOP = KINETICS["lycopene"].growth_rate_range

#: What the reference context alone allows: glucose, 30 C, air, exponential. 0.40 /h.
MU_ENV = context_growth_rate(CultureContext())


def _bisect_dose(still_answers, low: float, high: float, tolerance: float) -> float:
    """The lowest dose in ``(low, high]`` at which ``still_answers`` stops being true.

    The crossing is measured rather than typed, so no test here can pin a threshold the
    code does not have. Both ends of the bracket are checked, so a bracket that already
    sits on one side fails loudly instead of returning its own end.
    """
    assert still_answers(low), "the low end of the bracket is already past the crossing"
    assert not still_answers(high), "the high end of the bracket has not reached the crossing"
    while high - low > tolerance:
        middle = 0.5 * (low + high)
        if still_answers(middle):
            low = middle
        else:
            high = middle
    return high


@pytest.fixture(scope="module")
def carotene():
    return load_pathway("beta_carotene")


@pytest.fixture(scope="module")
def script():
    filename = pathlib.Path(__file__).resolve().parents[1] / "scripts" / \
        "product_environment_sweep.py"
    spec = importlib.util.spec_from_file_location("product_environment_sweep", filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _predict(spec, environment):
    return predict_product(spec, Genotype(entry_expression=1.0), environment,
                           calibrations.BETA_CAROTENE_FLUX, KINETICS, mode="empirical")


class TestOnlyOneProductCanBeAnsweredFor:
    def test_four_pathways_are_declared(self):
        assert set(available_pathways()) == {"beta_carotene", "gadusol", "glycogen", "phb"}

    def test_exactly_one_of_them_has_a_flux_calibration(self):
        calibrated = [name for name in available_pathways()
                      if any(getattr(calibrations, attribute).entry_enzyme
                             == load_pathway(name).entry_enzyme
                             for attribute in dir(calibrations)
                             if attribute.isupper()
                             and hasattr(getattr(calibrations, attribute), "entry_enzyme"))]

        assert calibrated == ["beta_carotene"]

    @pytest.mark.parametrize("product", ["gadusol", "glycogen", "phb"])
    def test_the_others_refuse_rather_than_borrowing_a_scalar(self, product):
        """Applying CrtE's fitted alpha to another gene's expression is a units error, and
        the refusal is what stops a number appearing that has no calibration behind it."""
        spec = load_pathway(product)

        with pytest.raises(ValueError, match="was fitted on"):
            predict_product(
                spec, Genotype(entry_expression=1.0),
                Environment(growth_rate_setpoint_per_h=HELD_RATE),
                calibrations.BETA_CAROTENE_FLUX, {}, mode="empirical")

    @pytest.mark.parametrize("product", ["gadusol", "glycogen", "phb"])
    def test_foreign_kinetics_are_refused_independently_of_the_scalar(self, product):
        with pytest.raises(ValueError, match="kinetics supplied for absent nodes"):
            _predict(load_pathway(product), Environment(growth_rate_setpoint_per_h=HELD_RATE))

    @pytest.mark.parametrize("product", ["gadusol", "glycogen", "phb"])
    def test_unused_parameters_on_a_declared_node_are_also_refused(self, product):
        spec = load_pathway(product)
        with pytest.raises(ValueError, match="not used by this node's declared rate law"):
            predict_product(
                spec, Genotype(entry_expression=1.0),
                Environment(growth_rate_setpoint_per_h=HELD_RATE),
                calibrations.BETA_CAROTENE_FLUX,
                {spec.nodes[0].name: NodeKinetics(vmax_per_growth=1e-3, km=1e-4)},
                mode="empirical")

    @pytest.mark.parametrize("product", ["gadusol", "phb"])
    def test_and_they_are_structurally_uncalibratable_not_merely_unfitted(self, product):
        """These diluted chains have no measured intermediate pool to constrain branch
        kinetics. An entry-expression series is a separate dataset requirement; glycogen
        also needs a terminal loss measurement and is not a dilution-only identity."""
        spec = load_pathway(product)

        assert all(node.rate_law == RateLaw.PASSTHROUGH for node in spec.nodes[:-1])
        assert spec.nodes[-1].fate == Fate.DILUTED
        assert spec.calibratability.flux_is_recoverable
        assert not spec.calibratable


class TestTheEnvironmentHasExactlyOneChannel:
    """The fixed empirical comparison is mu-only where a reachable, in-range answer exists."""

    @pytest.fixture(scope="class")
    def held_results(self, carotene):
        answered, refused = [], []
        for carbon in ("glucose", "galactose", "ethanol"):
            for temperature in (25.0, 30.0, 37.0):
                for oxygen in (0.21, 0.05, 0.01):
                    for stressor, dose in ((None, 0.0), ("DTT", 1.0), ("NaCl", 0.4)):
                        environment = Environment(
                            context=CultureContext(carbon_source=carbon,
                                                   temperature_c=temperature,
                                                   oxygen=oxygen),
                            stressor=stressor, dose=dose,
                            growth_rate_setpoint_per_h=HELD_RATE)
                        try:
                            got = _predict(carotene, environment)
                        except SetpointUnreachable as refusal:
                            # A context or dose at or below the setpoint has no held state.
                            # Retain every refused environment and its reason rather than
                            # counting only the conditions that returned a content value.
                            refused.append((environment, refusal))
                        else:
                            answered.append((environment, got))
        return answered, refused

    def test_every_environment_below_the_crossing_at_a_held_rate_gives_the_same_content(
            self, held_results):
        """Only the 35 reachable combinations share a content; 46 return no held state.

        Correct AT THESE DOSES, and the name says so. The seam is threshold-gated: a
        stressor reaches the product only once it has taken mu_max below the held rate, and
        at 0.18 /h that needs 1.6254 M NaCl and 1.6795 mM DTT, against the 0.4 M and 1.0 mM
        dosed here. Every temperature swept is inside the 16.94-42.48 C band as well.
        :class:`TestTheHeldSeamIsThresholdGated` pins the other side of all three.
        """
        answered, refused = held_results
        assert len(answered) == 35
        assert len(refused) == 46
        assert len(answered) + len(refused) == 81
        assert {got.growth_rate_per_h for _, got in answered} == {HELD_RATE}
        assert len({got.content_mg_per_gdcw for _, got in answered}) == 1

    def test_the_refused_combinations_remain_visible(self, held_results):
        answered, refused = held_results
        assert all(environment.context.carbon_source != "ethanol"
                   for environment, _ in answered)
        assert sum(environment.context.carbon_source == "ethanol"
                   for environment, _ in refused) == 27
        assert all("at or below the setpoint" in str(refusal) for _, refusal in refused)

    @pytest.mark.parametrize("stressor,dose", [(None, 0.0), ("DTT", 1.0)])
    def test_a_dose_that_outruns_the_setpoint_is_refused_not_answered(
            self, carotene, stressor, dose):
        """Ethanol cannot hold 0.18 /h even without a stressor: its default capacity is
        0.140 /h, and DTT at 1.0 mM lowers it to about 0.105 /h. Neither returns content."""
        environment = Environment(
            context=CultureContext(carbon_source="ethanol", temperature_c=30.0,
                                   oxygen=0.21),
            stressor=stressor, dose=dose, growth_rate_setpoint_per_h=HELD_RATE)

        with pytest.raises(SetpointUnreachable, match="at or below the setpoint"):
            _predict(carotene, environment)

    @pytest.mark.parametrize("carbon,stressor,dose", [
        ("galactose", "DTT", 1.0), ("galactose", "H2O2", 0.5),
        ("ethanol", None, 0.0), ("ethanol", "DTT", 1.0), ("ethanol", "H2O2", 0.5),
    ])
    def test_content_in_batch_is_reproduced_by_the_growth_rate_alone(
            self, carotene, carbon, stressor, dose):
        """Compare the in-range batch arithmetic with the solver, not a fictitious held
        state at the batch maximum: equality to capacity has no held fixed point."""
        batch = _predict(carotene, Environment(
            context=CultureContext(carbon_source=carbon), stressor=stressor, dose=dose))
        solved = solve_pathway(carotene, batch.flux.flux_mmol_per_gdcw_h,
                               batch.growth_rate_per_h, KINETICS)

        assert abs(solved.terminal.content_mg_per_gdcw - batch.content_mg_per_gdcw) < 1e-9

    @pytest.mark.parametrize("carbon,stressor,dose", [
        ("glucose", None, 0.0), ("glucose", "DTT", 1.0), ("glucose", "H2O2", 0.5),
        ("galactose", None, 0.0),
    ])
    def test_batch_rates_outside_the_fit_are_refused(self, carotene, carbon, stressor, dose):
        with pytest.raises(ValueError, match="outside the range its kinetics were fitted"):
            _predict(carotene, Environment(context=CultureContext(carbon_source=carbon),
                                          stressor=stressor, dose=dose))

    def test_a_batch_maximum_is_not_a_reachable_held_setpoint(self, carotene):
        environment = Environment(context=CultureContext(carbon_source="ethanol"))
        batch = _predict(carotene, environment)

        with pytest.raises(SetpointUnreachable, match="at or below the setpoint"):
            _predict(carotene, dataclasses.replace(
                environment, growth_rate_setpoint_per_h=batch.growth_rate_per_h))

    def test_the_chain_says_so_rather_than_letting_the_silence_speak(self, carotene):
        """Two identical numbers with no note read as 'stress does not affect production',
        which is a claim the chain has not earned."""
        got = _predict(carotene, Environment(growth_rate_setpoint_per_h=HELD_RATE,
                                             stressor="DTT", dose=1.0))

        assert any("no route modelled" in note for note in got.notes)


#: ``(stressor, dose at which the incumbent stops answering at HELD_RATE)``. Bisected, not
#: chosen; the two named agents are the ones the sweep above doses, plus two more.
CROSSINGS_AT_THE_HELD_RATE = [("NaCl", 1.62536659), ("DTT", 1.67954548),
                              ("H2O2", 1.08357773), ("sorbitol", 3.25073319)]


class TestTheHeldSeamIsThresholdGated:
    """The other side of the doses the class above sweeps, and it is not a moved number.

    The seam is real and it is a threshold: a stressor lowers mu_max, and once mu_max drops
    below the held setpoint there is no held state left. What the incumbent does there is
    REFUSE. So the pair is "bit-identical below, no answer above", and neither half is
    evidence for the other -- which is why both are pinned rather than one being inferred.

    Every crossing below is bisected from the incumbent's own behaviour, and then checked
    against `stress_panel.viability`, which is the single function every panel agent reaches
    mu through. That is why one number, the setpoint over mu_env, fixes all of them.
    """

    @staticmethod
    def _answers(carotene, stressor: str, dose: float, setpoint: float) -> bool:
        try:
            _predict(carotene, Environment(stressor=stressor, dose=dose,
                                           growth_rate_setpoint_per_h=setpoint))
        except SetpointUnreachable:
            return False
        return True

    @pytest.mark.parametrize("stressor,crossing", CROSSINGS_AT_THE_HELD_RATE)
    def test_the_crossing_is_bisected_and_lands_where_viability_meets_the_setpoint(
            self, carotene, stressor, crossing):
        """Measured from the incumbent, then explained by the one function behind it."""
        lethal = STRESSORS[stressor].lethal_dose
        measured = _bisect_dose(
            lambda dose: self._answers(carotene, stressor, dose, HELD_RATE),
            0.0, 4.0 * lethal, tolerance=1e-7 * lethal)

        assert measured == pytest.approx(crossing, abs=1e-6)
        assert viability(stressor, measured) == pytest.approx(HELD_RATE / MU_ENV, rel=1e-7)

    @pytest.mark.parametrize("stressor,crossing", CROSSINGS_AT_THE_HELD_RATE)
    def test_below_the_crossing_the_content_is_bit_identical_and_above_it_there_is_none(
            self, carotene, stressor, crossing):
        """Not "close" below and not "smaller" above: the same float, and then no number."""
        unstressed = _predict(carotene, Environment(growth_rate_setpoint_per_h=HELD_RATE))

        for dose in (0.5 * crossing, 0.99 * crossing, crossing - 1e-6):
            got = _predict(carotene, Environment(stressor=stressor, dose=dose,
                                                 growth_rate_setpoint_per_h=HELD_RATE))
            assert got.growth_rate_per_h == HELD_RATE
            assert got.content_mg_per_gdcw == unstressed.content_mg_per_gdcw

        for dose in (crossing + 1e-6, 1.01 * crossing, 2.0 * crossing):
            with pytest.raises(SetpointUnreachable, match="at or below the setpoint"):
                _predict(carotene, Environment(stressor=stressor, dose=dose,
                                               growth_rate_setpoint_per_h=HELD_RATE))

    def test_at_this_held_rate_every_crossing_sits_above_its_agents_lethal_dose(self):
        """The negative, and it is exact rather than a coincidence of the four agents above.

        `stress_panel.py` defines a lethal dose as where growth falls to HALF of control, so
        `viability` is exactly 0.5 there for every agent and mu_max at the lethal dose is
        exactly mu_env/2 = 0.2 /h. A crossing therefore sits at or below the lethal dose if
        and only if the held setpoint is at or above 0.2 /h. HELD_RATE is 0.18, so on this
        sweep every crossing is past the dose the panel already calls lethal -- 1.0836 x it
        for all 25 agents -- and the whole above-crossing branch is arithmetic about the Hill
        rather than a statement about a culture.
        """
        assert HELD_RATE < MU_ENV / 2.0

        for stressor in sorted(STRESSORS):
            lethal = STRESSORS[stressor].lethal_dose
            assert viability(stressor, lethal) == 0.5, stressor
            assert MU_ENV * viability(stressor, lethal) > HELD_RATE, stressor


class TestTheStressRouteWouldBeBelowTheNoiseFloor:
    """The planning question, and the reason this is not simply a to-do item.

    The one candidate channel is stress -> maintenance -> growth -> content. Its size is
    measured in `outputs/maintenance_scale.csv`. Carried through to content it is far smaller
    than the error of the flux scalar it would have to pass through, so wiring it would not
    produce a number anyone could check -- whichever way the unresolved sign goes.
    """

    def test_the_flux_calibration_error_is_the_floor(self):
        assert calibrations.BETA_CAROTENE_FLUX.typical_fold_error > 1.2

    @pytest.mark.parametrize("growth_penalty_percent", [0.85, 8.55])
    def test_a_measured_stress_penalty_below_the_crossing_moves_content_far_less_than_that(
            self, carotene, growth_penalty_percent):
        """0.85% and 8.55% are the measured full-stress growth penalties at glucose -10.0
        and -1.0 (`outputs/maintenance_scale.csv`).

        Below the crossing the penalty does not reach the product AT ALL -- it lands on
        mu_max, which the feed is still holding above. Pricing it as if it reached mu is the
        generous case, and even that is a tenth of the flux scalar's own error.
        """
        base = _predict(carotene, Environment(growth_rate_setpoint_per_h=HELD_RATE))
        stressed = _predict(carotene, Environment(
            growth_rate_setpoint_per_h=HELD_RATE * (1 - growth_penalty_percent / 100.0)))

        moved = abs(stressed.content_mg_per_gdcw - base.content_mg_per_gdcw) \
            / base.content_mg_per_gdcw
        floor = calibrations.BETA_CAROTENE_FLUX.typical_fold_error - 1.0

        assert moved < floor / 10.0

    def test_and_the_largest_move_anywhere_above_the_crossing_is_still_below_that_floor(
            self, carotene):
        """The bound, so the pair does not read as "small penalty, therefore small effect".

        Above the crossing mu IS free to fall, but only to 0.100987987 /h: below that the
        carotenoid kinetics are outside the window they were fitted over and
        `pathway/solve.py` refuses rather than extrapolating. Priced at that floor the
        content moves 8.29% from a 0.18 /h hold and 18.77% from the 0.254320182 /h top -- larger than a
        maintenance penalty by an order of magnitude, and still inside the 1.2216-fold
        typical error of the flux scalar every one of these numbers passes through.
        """
        floor_content = _predict(
            carotene, Environment(growth_rate_setpoint_per_h=FIT_FLOOR)).content_mg_per_gdcw
        flux_floor = calibrations.BETA_CAROTENE_FLUX.typical_fold_error - 1.0

        moves = {}
        for setpoint in (HELD_RATE, FIT_TOP):
            held = _predict(carotene, Environment(growth_rate_setpoint_per_h=setpoint))
            moves[setpoint] = floor_content / held.content_mg_per_gdcw - 1.0

        assert moves[HELD_RATE] == pytest.approx(0.08292027, abs=5e-8)
        assert moves[FIT_TOP] == pytest.approx(0.18766029, abs=5e-8)
        for setpoint, moved in moves.items():
            assert moved > flux_floor / 10.0, setpoint
            assert moved < flux_floor, setpoint

        with pytest.raises(ValueError, match="outside the range its kinetics"):
            _predict(carotene, Environment(growth_rate_setpoint_per_h=FIT_FLOOR * 0.99))


class TestTheSweepReportsEveryAttempt:
    def test_the_shipped_kinetic_window_is_not_removed(self, script):
        assert script.KINETICS is calibrations.BETA_CAROTENE_KINETICS
        assert script.KINETICS["lycopene"].growth_rate_range is not None

    def test_refusals_and_comparison_scope_are_written_and_counted(
            self, script, monkeypatch, tmp_path, capsys):
        monkeypatch.setattr(script.paths, "outputs_dir", lambda: tmp_path)

        assert script.main() == 0
        output = capsys.readouterr().out
        table = pd.read_csv(tmp_path / "product_environment_sweep.csv")
        attempted_per_mode = (len(script.CARBON) * len(script.TEMPERATURE_C)
                              * len(script.OXYGEN) * len(script.STRESSORS))
        assert len(table) == 2 * attempted_per_mode
        assert set(table.prediction_mode) == {"empirical"}
        assert table.prediction_scope.str.contains("not environment validation").all()
        refused = table[table.refused.notna()]
        answered = table[table.refused.isna()]
        assert not refused.empty and not answered.empty
        assert refused.refusal_reason.notna().all()
        assert refused[["growth_rate_per_h", "content_mg_per_gdcw"]].isna().all().all()
        assert answered[["growth_rate_per_h", "content_mg_per_gdcw"]].notna().all().all()
        low, high = script.KINETICS["lycopene"].growth_rate_range
        assert answered.growth_rate_per_h.between(low, high).all()
        held_ethanol = table[(table["mode"] == "held") & (table.carbon_source == "ethanol")]
        assert (held_ethanol.refused == "SetpointUnreachable").all()
        assert held_ethanol.refusal_reason.str.contains("at or below the setpoint").all()
        assert output.count(f"refused out of {attempted_per_mode} attempted") == 2
        assert "SetpointUnreachable:" in output
        assert "not environment validation" in output

    def test_unexpected_errors_are_not_mislabelled_as_domain_refusals(
            self, script, monkeypatch, tmp_path):
        monkeypatch.setattr(script.paths, "outputs_dir", lambda: tmp_path)

        def broken_prediction(*args, **kwargs):
            raise RuntimeError("unexpected solver failure")

        monkeypatch.setattr(script, "predict_product", broken_prediction)
        with pytest.raises(RuntimeError, match="unexpected solver failure"):
            script.main()
