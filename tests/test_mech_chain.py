"""The join, and the four claims the join is scored on.

Everything asserted here is computed by the test rather than repeated from a report:

1. **END TO END.** The incumbent returns ONE content for eight environments at a held
   setpoint -- 0.97656355 mg/gDCW at mu_set = 0.2543 /h, identical to eight decimal places,
   and refuses two more. The chain returns six distinct numbers spanning 1.99x over the same
   ten. The route is not a new mechanism: it is that a feed cannot hold a growth rate the
   strain cannot reach, which `predict.py::_growth_rate` checks for a stressor and not for a
   context. That makes the whole seam THRESHOLD-GATED, so section 1b pins both sides of the
   crossing on every axis it reaches, and every name there says which side it is on.

2. **MECHANISM.** Five differential states the environment drives, against zero in the
   incumbent stress path, where a stressor enters as an instantaneous Hill and the only
   integrated states belong to the culture and the reporter.

3. **TIMESCALES.** The tau/T audit is run on the ASSEMBLED state, in both windows, and the
   sets differ: 4 states integrated on the 4.14 h plate read, 2 on the 5 d fed-batch.

4. **GEOMETRY.** On the one channel the mechanism reaches, the dose x time surface's ray
   share falls from a median 0.999516 to 0.975763 and the off-ray energy rises ~50x. On the
   other 24 stressors it does not move at all, and the panel-wide median is unchanged to six
   decimal places. That negative is asserted here as firmly as the positive.

And the gate: the assembled layer has MORE free scalars than independent targets and
:func:`require_assembled_gate` raises. That refusal is pinned, because a later change that
makes it pass by adding a target nobody measured would otherwise pass silently.
"""

from __future__ import annotations


import numpy as np
import pytest

from ystwin.generator.context import (
    CARDINAL_TEMPERATURES_C,
    CultureContext,
    context_growth_rate,
)
from ystwin.generator.panel_experiment import (
    MEASURED_GROWTH_RATE_SE,
    OBSERVED_ACTIVITY_CV,
)
from ystwin.generator.stress_panel import (
    STRESSORS,
    healthy_ladder,
    transcriptional_reporters,
    viability,
)
from ystwin.mech.chain import (
    CHAIN_PARAMS,
    IN_CHAIN_REGISTRIES,
    LAYER_ORDER,
    SweepPoint,
    assembled_gate,
    assembled_states,
    chain_band,
    dose_time_surface,
    geometry_audit,
    ode_state_census,
    provenance,
    ray_share,
    reaction_depths,
    reachable_reactions,
    reduction_audit,
    require_assembled_gate,
    run_chain,
)
from ystwin.mech.params import FreeScalarGateFailed
from ystwin.mech.state import FEDBATCH_5D, PLATE_READ_4H
from ystwin.pathway import calibrations
from ystwin.pathway.spec import load_pathway
from ystwin.predict import (
    Environment,
    Genotype,
    LayerState,
    SetpointUnreachable,
    predict_product,
)

#: The upper Elizondo state, and the top of the range the carotenoid kinetics were fitted
#: over. Chosen because it is where the incumbent's blindness is widest.
SETPOINT = 0.2543

MID = SweepPoint.midpoint()

#: The growth-rate window `pathway/solve.py` will answer inside. Outside it the product arm
#: refuses, which is what bounds how far past a crossing any test here can go.
FIT_FLOOR, FIT_TOP = calibrations.BETA_CAROTENE_KINETICS["lycopene"].growth_rate_range

#: What the reference context alone allows: glucose, 30 C, air, exponential. 0.40 /h.
MU_ENV = context_growth_rate(CultureContext())


def _bisect(still_held, low: float, high: float, tolerance: float) -> float:
    """The lowest input in ``(low, high]`` at which ``still_held`` stops being true.

    Every threshold in this file is MEASURED with this rather than typed, so a test can
    never assert a crossing the code does not have. The caller supplies the bracket and
    asserts the value that comes back; the two ends are checked here so a bracket that
    already sits on one side of the crossing fails loudly instead of returning its own end.
    """
    assert still_held(low), "the low end of the bracket is already past the crossing"
    assert not still_held(high), "the high end of the bracket has not reached the crossing"
    while high - low > tolerance:
        middle = 0.5 * (low + high)
        if still_held(middle):
            low = middle
        else:
            high = middle
    return high


@pytest.fixture(scope="module")
def spec():
    return load_pathway("beta_carotene")


@pytest.fixture(scope="module")
def product_kwargs(spec):
    return dict(window=FEDBATCH_5D, genotype=Genotype(1.0, "ref"), spec=spec,
                flux_calibration=calibrations.BETA_CAROTENE_FLUX,
                kinetics=calibrations.BETA_CAROTENE_KINETICS)


def _environments() -> dict[str, Environment]:
    """The eight the build brief names, plus two the acid arm needs to be visible."""
    return {
        "reference": Environment(growth_rate_setpoint_per_h=SETPOINT),
        "ethanol": Environment(context=CultureContext(carbon_source="ethanol"),
                               growth_rate_setpoint_per_h=SETPOINT),
        "galactose": Environment(context=CultureContext(carbon_source="galactose"),
                                 growth_rate_setpoint_per_h=SETPOINT),
        "37 C": Environment(context=CultureContext(temperature_c=37.0),
                            growth_rate_setpoint_per_h=SETPOINT),
        "20 C": Environment(context=CultureContext(temperature_c=20.0),
                            growth_rate_setpoint_per_h=SETPOINT),
        "pH 4 medium": Environment(context=CultureContext(ph_medium=4.0),
                                   growth_rate_setpoint_per_h=SETPOINT),
        "anoxia": Environment(
            context=CultureContext(oxygen=0.0, anaerobic_supplements=True),
            growth_rate_setpoint_per_h=SETPOINT),
        "H2O2": Environment(stressor="H2O2", dose=1.0,
                            growth_rate_setpoint_per_h=SETPOINT),
        "acetic 175 pH4": Environment(context=CultureContext(ph_medium=4.0),
                                      stressor="acetic_acid", dose=175.0,
                                      growth_rate_setpoint_per_h=SETPOINT),
        "ethanol 37 C": Environment(
            context=CultureContext(carbon_source="ethanol", temperature_c=37.0),
            growth_rate_setpoint_per_h=SETPOINT),
    }


# --------------------------------------------------------------------------------------
# 1. END TO END
# --------------------------------------------------------------------------------------

class TestTheEnvironmentReachesTheProduct:
    """Criterion 1, measured against the incumbent on the same ten conditions."""

    def test_default_empirical_comparison_refuses_unreachable_or_unmodelled_conditions(
            self, spec):
        """The corrected empirical baseline: supported held rates agree; unsupported inputs refuse."""
        refusals = {
            "ethanol": (SetpointUnreachable, "maximum growth rate"),
            "20 C": (SetpointUnreachable, "maximum growth rate"),
            "pH 4 medium": (ValueError, "ph_medium"),
            "H2O2": (SetpointUnreachable, "maximum growth rate"),
            "acetic 175 pH4": (ValueError, "ph_medium"),
            "ethanol 37 C": (SetpointUnreachable, "maximum growth rate"),
        }
        contents = {}
        for name, environment in _environments().items():
            if name in refusals:
                exception, message = refusals[name]
                with pytest.raises(exception, match=message):
                    predict_product(spec, Genotype(1.0), environment,
                                    calibrations.BETA_CAROTENE_FLUX,
                                    calibrations.BETA_CAROTENE_KINETICS)
                continue
            got = predict_product(spec, Genotype(1.0), environment,
                                  calibrations.BETA_CAROTENE_FLUX,
                                  calibrations.BETA_CAROTENE_KINETICS)
            assert got.mode == "empirical"
            contents[name] = round(got.content_mg_per_gdcw, 8)
        assert set(contents) == {"reference", "galactose", "37 C", "anoxia"}
        assert len(set(contents.values())) == 1
        assert contents["reference"] == pytest.approx(0.97656355, abs=5e-9)

    def test_the_chain_returns_six_distinct_numbers_over_the_same_ten(
            self, product_kwargs):
        results = {name: run_chain(env, sweep=MID, **product_kwargs)
                   for name, env in _environments().items()}
        contents = {round(r.product.population_content_mg_per_gdcw, 8)
                    for r in results.values()}
        assert len(contents) == 6
        assert max(contents) / min(contents) == pytest.approx(1.9867, abs=1e-3)
        # Nothing is refused: the two the incumbent raises on are the two where the feed
        # stops setting mu, and the chain reports the realised rate instead.
        assert not results["H2O2"].feed_sets_mu
        assert not results["acetic 175 pH4"].feed_sets_mu

    def test_the_reference_condition_reproduces_the_incumbent_bit_for_bit(
            self, spec, product_kwargs):
        """The chain must not move the number where no mechanism applies."""
        environment = Environment(growth_rate_setpoint_per_h=0.18)
        incumbent = predict_product(spec, Genotype(1.0), environment,
                                    calibrations.BETA_CAROTENE_FLUX,
                                    calibrations.BETA_CAROTENE_KINETICS)
        chain = run_chain(environment, sweep=MID, **product_kwargs)
        assert chain.product.content_mg_per_gdcw == pytest.approx(
            incumbent.content_mg_per_gdcw, rel=0, abs=1e-12)
        assert incumbent.content_mg_per_gdcw == pytest.approx(1.07098788, abs=5e-9)

    def test_empirical_prediction_refuses_where_the_legacy_chain_reports_a_clipped_rate(
            self, product_kwargs):
        """A reported legacy rate is not evidence that the vessel can hold its requested setpoint."""
        environment = _environments()["ethanol"]
        chain = run_chain(environment, sweep=MID, **product_kwargs)
        assert chain.mu_max_per_h == pytest.approx(0.14, abs=1e-9)
        assert chain.growth_rate_per_h == pytest.approx(0.14, abs=1e-9)
        assert not chain.feed_sets_mu
        # predict.py accepts the same setpoint without checking it.
        with pytest.raises(SetpointUnreachable, match="maximum growth rate"):
            predict_product(
                load_pathway("beta_carotene"), Genotype(1.0), environment,
                calibrations.BETA_CAROTENE_FLUX, calibrations.BETA_CAROTENE_KINETICS)

    def test_the_growth_ablation_clears_its_own_named_floor(self, product_kwargs):
        """OBSERVABLE growth rate, FLOOR MEASURED_GROWTH_RATE_SE. Not the plate CV."""
        chain = run_chain(_environments()["ethanol 37 C"], sweep=MID, **product_kwargs)
        moved = abs(chain.growth_rate_per_h - SETPOINT)
        assert moved / MEASURED_GROWTH_RATE_SE > 10.0
        assert moved == pytest.approx(0.1261, abs=1e-3)

    def test_the_acid_channel_opens_at_a_dose_inside_the_panels_own_ladder(
            self, product_kwargs):
        """Between 160 and 175 mM at pH 4, against a declared lethal dose of 180."""
        def feed_sets_mu(dose: float) -> bool:
            environment = Environment(context=CultureContext(ph_medium=4.0),
                                      stressor="acetic_acid", dose=dose,
                                      growth_rate_setpoint_per_h=SETPOINT)
            return run_chain(environment, sweep=MID, **product_kwargs).feed_sets_mu

        assert feed_sets_mu(160.0)
        assert not feed_sets_mu(175.0)
        assert 175.0 < STRESSORS["acetic_acid"].lethal_dose

    def test_the_acid_growth_arm_is_below_its_floor_wherever_it_is_physical(
            self, product_kwargs):
        """The sharpest negative here, and it is pinned rather than reported.

        The pump-holds bill traps 210x the outside undissociated acid at the resting pH, so
        it passes mech/ph.py's own 500 mM osmotic flag above ~3.7 mM acetate at pH 4.5. At
        the largest dose where the state is still one a cell could hold, the growth effect
        is a third of the growth assay's own noise.
        """
        physical = Environment(context=CultureContext(ph_medium=4.5),
                               stressor="acetic_acid", dose=3.7,
                               growth_rate_setpoint_per_h=SETPOINT)
        chain = run_chain(physical, sweep=MID, **product_kwargs)
        assert chain.bill_is_physical
        moved = abs(chain.mu_env_per_h - chain.mu_max_per_h)
        assert moved / MEASURED_GROWTH_RATE_SE < 1.0

        beyond = run_chain(_environments()["acetic 175 pH4"], sweep=MID, **product_kwargs)
        assert not beyond.bill_is_physical
        assert any("osmotic scale" in note for note in beyond.notes)

    def test_medium_ph_moves_the_growth_rate_only_through_the_acid(self, product_kwargs):
        """pH alone reaches nothing, and that is a refusal rather than an omission."""
        bare = run_chain(_environments()["pH 4 medium"], sweep=MID, **product_kwargs)
        reference = run_chain(_environments()["reference"], sweep=MID, **product_kwargs)
        assert bare.mu_max_per_h == pytest.approx(reference.mu_max_per_h)
        assert bare.layer("weak acid").state == LayerState.NOT_RUN

    def test_a_growth_rate_outside_the_fitted_window_is_still_refused(
            self, product_kwargs):
        environment = Environment(context=CultureContext(growth_phase="stationary"),
                                  growth_rate_setpoint_per_h=SETPOINT)
        with pytest.raises(ValueError, match="outside the range its kinetics"):
            run_chain(environment, sweep=MID, **product_kwargs)


# --------------------------------------------------------------------------------------
# 1b. THE SEAM IS THRESHOLD-GATED
# --------------------------------------------------------------------------------------

#: ``(medium pH, setpoint /h, crossing mM total acetate, content below it, a dose above it,
#: content there)``. The crossings are bisected, not chosen; see the first test below.
ACID_CROSSINGS = (
    (4.0, 0.2543, 161.43215, 0.97656355, 200.0, 1.00301159),
    (4.5, 0.2543, 213.11122, 0.97656355, 300.0, 1.01913211),
    (4.5, 0.18, 454.61452, 1.07098788, 500.0, 1.08230172),
)


def _acid(ph_medium: float, dose_mM: float, setpoint: float, **kwargs):
    """One chain run on the acid arm at one medium pH, dose and held setpoint."""
    kwargs.setdefault("window", FEDBATCH_5D)
    return run_chain(
        Environment(context=CultureContext(ph_medium=ph_medium), stressor="acetic_acid",
                    dose=dose_mM, growth_rate_setpoint_per_h=setpoint),
        sweep=MID, **kwargs)


class TestTheAcidSeamIsThresholdGated:
    """The pair, and why one half of it is not evidence for the other.

    The vessel layer sets mu = min(setpoint, mu_max), so an environment reaches the product
    only once it has pushed mu_max BELOW the setpoint. Under that dose the feed still holds
    mu and the content is bit-identical; over it the content moves. So "the acid does not
    reach the product" and "the acid reaches the product" are both true statements about
    different doses, and both are pinned here rather than one being inferred.

    Every crossing is bisected from `feed_sets_mu`, never typed, and each is quoted with the
    two facts that decide whether it is a prediction or arithmetic: where it sits against
    acetic acid's declared lethal dose, and whether mech/ph.py still calls the bill physical.
    """

    @pytest.mark.parametrize("ph,setpoint,crossing,below,above,moved", ACID_CROSSINGS)
    def test_the_crossing_is_bisected_from_the_feed_rather_than_typed(
            self, ph, setpoint, crossing, below, above, moved):
        """And it is not a dose threshold at all: it is where mu_max meets the setpoint."""
        measured = _bisect(lambda dose: _acid(ph, dose, setpoint).feed_sets_mu,
                           0.0, 4.0 * crossing, tolerance=1e-5)

        assert measured == pytest.approx(crossing, abs=1e-4)
        assert _acid(ph, measured, setpoint).mu_max_per_h == pytest.approx(setpoint,
                                                                          abs=1e-6)

    @pytest.mark.parametrize("ph,setpoint,crossing,below,above,moved", ACID_CROSSINGS)
    def test_below_the_crossing_the_acid_leaves_the_content_bit_identical(
            self, product_kwargs, ph, setpoint, crossing, below, above, moved):
        """Not "close" -- the same float, because mu is the same float."""
        reference = _acid(ph, 0.0, setpoint, **product_kwargs)
        assert reference.product.content_mg_per_gdcw == pytest.approx(below, abs=5e-9)

        for dose in (0.5 * crossing, 0.99 * crossing, crossing - 1e-3):
            got = _acid(ph, dose, setpoint, **product_kwargs)
            assert got.feed_sets_mu
            assert got.growth_rate_per_h == setpoint
            assert (got.product.content_mg_per_gdcw
                    == reference.product.content_mg_per_gdcw)

    @pytest.mark.parametrize("ph,setpoint,crossing,below,above,moved", ACID_CROSSINGS)
    def test_above_the_crossing_the_acid_moves_the_content(
            self, product_kwargs, ph, setpoint, crossing, below, above, moved):
        just_above = _acid(ph, crossing + 0.05, setpoint, **product_kwargs)
        far_above = _acid(ph, above, setpoint, **product_kwargs)

        for got in (just_above, far_above):
            assert not got.feed_sets_mu
            assert got.growth_rate_per_h == got.mu_max_per_h
            assert got.growth_rate_per_h < setpoint
            assert got.product.content_mg_per_gdcw > below
        assert far_above.product.content_mg_per_gdcw == pytest.approx(moved, abs=5e-9)

    def test_which_of_these_doses_are_mathematical_rather_than_physiological(self):
        """Said plainly, because most of them are.

        `stress_panel.py` declares acetic acid lethal at 180 mM (ec50 60.0 x
        _LETHAL_MULTIPLE 3.0). Only the pH 4.0 crossing at 161.4 mM is under it. The pH 4.5
        crossings at 213.1 and 454.6 mM, and both 300 and 500 mM above them, are over it and
        are arithmetic about the rate law rather than statements about a culture. The next
        test takes the 161.4 mM row away too, on the separate ground of the anion pool.
        """
        lethal = STRESSORS["acetic_acid"].lethal_dose
        assert lethal == 180.0

        physiological = [row[2] for row in ACID_CROSSINGS if row[2] < lethal]
        beyond = [row[2] for row in ACID_CROSSINGS if row[2] > lethal]
        assert physiological == [161.43215]
        assert beyond == [213.11122, 454.61452]
        assert all(row[4] > lethal for row in ACID_CROSSINGS)

    @pytest.mark.parametrize("ph,ceiling", [(4.0, 2.809032), (4.5, 3.708284)])
    def test_and_every_acid_crossing_is_far_above_where_the_bill_stays_physical(
            self, ph, ceiling):
        """The harder half of the same honesty, and it costs the acid arm the pH 4.0 row too.

        `mech/ph.py`'s own osmotic flag goes out in single-digit mM, so the anion pool a
        crossing dose implies is one no cell holds. Every crossing here is 57x that ceiling
        or more -- below the lethal dose is not the same as inside the physical range.
        """
        measured = _bisect(lambda dose: _acid(ph, dose, SETPOINT).proton_bill.physical,
                           1e-3, 500.0, tolerance=1e-6)
        assert measured == pytest.approx(ceiling, abs=1e-4)

        for _, setpoint, crossing, *_rest in (r for r in ACID_CROSSINGS if r[0] == ph):
            assert crossing > 50.0 * measured
            assert not _acid(ph, crossing, setpoint).proton_bill.physical

    def test_far_enough_above_a_crossing_the_product_arm_refuses_rather_than_extrapolating(
            self, product_kwargs):
        """The seam does not run out of gate, it runs out of calibration.

        At pH 4.0 and a 0.2543 /h setpoint the content keeps moving up to 834.25 mM, where
        mu_max reaches 0.100987987 /h, the bottom of the window the carotenoid kinetics were
        fitted over. Past that `pathway/solve.py` refuses, which is the right answer.
        """
        refusal = _bisect(lambda dose: _acid(4.0, dose, SETPOINT).mu_max_per_h >= FIT_FLOOR,
                          161.5, 2000.0, tolerance=1e-4)
        assert refusal == pytest.approx(834.249, abs=1e-3)

        inside = _acid(4.0, 833.0, SETPOINT, **product_kwargs)
        assert inside.product.content_mg_per_gdcw == pytest.approx(1.15967912, abs=5e-9)
        with pytest.raises(ValueError, match="outside the range its kinetics"):
            _acid(4.0, 836.0, SETPOINT, **product_kwargs)


class TestTheOtherStressAxesCrossToo:
    """Osmotic, oxidative and temperature-as-a-dose, asked for by name and measured.

    Every agent on the panel reaches mu through ONE function, `stress_panel.viability`, so
    the crossing sits at the same viability for all of them -- the setpoint over mu_env --
    and therefore at the same fraction of each agent's own declared lethal dose. That is why
    these tests select agents by the module a stressor DECLARES rather than by name, and why
    the assertion is against `viability` rather than against a dose anybody wrote down.
    """

    @staticmethod
    def _agents(module: str) -> list[str]:
        return sorted(name for name, spec in STRESSORS.items()
                      if spec.targets.get(module, 0.0) > 0.0)

    @staticmethod
    def _dosed(name: str, dose: float, setpoint: float, **kwargs):
        kwargs.setdefault("window", FEDBATCH_5D)
        return run_chain(Environment(stressor=name, dose=dose,
                                     growth_rate_setpoint_per_h=setpoint),
                         sweep=MID, **kwargs)

    @pytest.mark.parametrize("module", ["osmotic", "oxidative", "heat"])
    @pytest.mark.parametrize("setpoint,fraction", [(0.2543, 0.80028611),
                                                   (0.18, 1.08357773)])
    def test_a_crossing_exists_on_each_axis_at_one_fraction_of_every_lethal_dose(
            self, module, setpoint, fraction):
        """Yes, a crossing exists on all three -- and it is the same crossing wearing
        different units, which is a fact about `viability` rather than about the agents."""
        agents = self._agents(module)
        assert agents

        for name in agents:
            lethal = STRESSORS[name].lethal_dose
            crossing = _bisect(
                lambda dose: self._dosed(name, dose, setpoint).feed_sets_mu,
                0.0, 4.0 * lethal, tolerance=1e-7 * lethal)

            assert crossing / lethal == pytest.approx(fraction, rel=1e-6)
            assert viability(name, crossing) == pytest.approx(setpoint / MU_ENV, rel=1e-6)

    # One content above the crossing for all five, because one viability sets them all.
    @pytest.mark.parametrize("module,name,crossing,above_content", [
        ("osmotic", "NaCl", 1.20042917, 0.97897541),
        ("osmotic", "sorbitol", 2.40085833, 0.97897541),
        ("oxidative", "H2O2", 0.80028611, 0.97897541),
        ("oxidative", "menadione", 240.08583, 0.97897541),
        ("heat", "heat", 14.40515, 0.97897541),
    ])
    def test_below_the_crossing_the_content_is_bit_identical_and_above_it_it_moves(
            self, product_kwargs, module, name, crossing, above_content):
        """The pair on the panel axes, 0.8% either side of each measured crossing."""
        assert name in self._agents(module)
        reference = run_chain(Environment(growth_rate_setpoint_per_h=SETPOINT), sweep=MID,
                              **product_kwargs).product.content_mg_per_gdcw

        below = self._dosed(name, 0.99 * crossing, SETPOINT, **product_kwargs)
        assert below.feed_sets_mu
        assert below.product.content_mg_per_gdcw == reference

        above = self._dosed(name, 1.008 * crossing, SETPOINT, **product_kwargs)
        assert not above.feed_sets_mu
        assert above.product.content_mg_per_gdcw > reference
        assert above.product.content_mg_per_gdcw == pytest.approx(above_content, abs=5e-8)

    def test_no_dose_on_the_panels_own_healthy_ladder_crosses_at_any_admissible_setpoint(
            self):
        """The negative, with the range stated, and it holds for all 25 agents at once.

        `healthy_ladder` tops out at 0.55 x lethal, where `viability` is 0.8168. A crossing
        needs viability down at setpoint/mu_env, and the highest setpoint the carotenoid
        kinetics accept is 0.254320182 /h, i.e. 0.6358. So over the whole admissible window
        [0.100987987, 0.254320182] /h the feed still sets mu at every rung of every ladder,
        and no crossing in this file is reachable by a dose the panel would actually spend.

        acetic_acid is in the loop on its own terms: at the default medium pH its ladder top
        runs through the mechanism's acid arm rather than `viability`, and leaves mu_max at
        0.376 /h. It does not cross either, for a different reason.
        """
        ladder_top_viability = viability("NaCl", 0.55 * STRESSORS["NaCl"].lethal_dose)
        assert ladder_top_viability == pytest.approx(0.81676658, abs=5e-9)
        assert ladder_top_viability > FIT_TOP / MU_ENV

        for name in sorted(STRESSORS):
            top = float(healthy_ladder(name, n_doses=5)[-1])
            assert top == pytest.approx(0.55 * STRESSORS[name].lethal_dose, rel=1e-9)
            assert self._dosed(name, top, FIT_TOP).feed_sets_mu, name


class TestTheTemperatureAxisCrossesOnBothSidesOfTheOptimum:
    """Temperature is a CONTEXT axis rather than a dose, so its admissible range is the
    cardinal one: `generator/context.py` returns a zero growth factor outside (2.8, 45.4) C,
    which is what T_min and T_max mean. Both crossings sit strictly inside it."""

    @staticmethod
    def _at(temperature_c: float, setpoint: float, **kwargs):
        kwargs.setdefault("window", FEDBATCH_5D)
        return run_chain(
            Environment(context=CultureContext(temperature_c=temperature_c),
                        growth_rate_setpoint_per_h=setpoint),
            sweep=MID, **kwargs)

    @pytest.mark.parametrize("setpoint,cold,hot", [(0.2543, 20.529929, 40.852158),
                                                   (0.18, 16.936271, 42.484865)])
    def test_both_crossings_are_bisected_and_sit_inside_the_cardinal_range(
            self, setpoint, cold, hot):
        low, high = CARDINAL_TEMPERATURES_C["min"], CARDINAL_TEMPERATURES_C["max"]
        measured_cold = _bisect(lambda t: not self._at(t, setpoint).feed_sets_mu,
                                low + 0.2, CARDINAL_TEMPERATURES_C["opt"], tolerance=1e-6)
        measured_hot = _bisect(lambda t: self._at(t, setpoint).feed_sets_mu,
                               CARDINAL_TEMPERATURES_C["opt"], high - 0.1, tolerance=1e-6)

        assert measured_cold == pytest.approx(cold, abs=1e-5)
        assert measured_hot == pytest.approx(hot, abs=1e-5)
        assert low < measured_cold < measured_hot < high

    def test_between_the_crossings_the_content_is_bit_identical_and_outside_them_it_moves(
            self, product_kwargs):
        """25 and 37 C are the sweep's own incubator temperatures and both sit between the
        0.2543 /h crossings at 20.53 and 40.85 C; 20 and 41 C sit just outside them."""
        reference = self._at(30.0, SETPOINT, **product_kwargs)
        assert reference.product.content_mg_per_gdcw == pytest.approx(0.97656355, abs=5e-9)

        for temperature in (25.0, 30.0, 37.0):
            held = self._at(temperature, SETPOINT, **product_kwargs)
            assert held.feed_sets_mu
            assert (held.product.content_mg_per_gdcw
                    == reference.product.content_mg_per_gdcw)

        for temperature, content in ((20.0, 0.99072116), (41.0, 0.98439771)):
            moved = self._at(temperature, SETPOINT, **product_kwargs)
            assert not moved.feed_sets_mu
            assert moved.product.content_mg_per_gdcw > reference.product.content_mg_per_gdcw
            assert moved.product.content_mg_per_gdcw == pytest.approx(content, abs=5e-9)

    def test_and_the_outer_edge_of_the_moving_band_is_a_refusal_not_an_extrapolation(
            self, product_kwargs):
        """mu_env leaves the fitted window at 12.85 and 43.89 C, well inside the cardinal
        range, so the band this axis can answer over is 12.85-20.53 C and 40.85-43.89 C."""
        for temperature in (12.9, 43.8):
            assert self._at(temperature, SETPOINT,
                            **product_kwargs).product.content_mg_per_gdcw > 1.15

        for temperature in (12.8, 43.9):
            assert 0.0 < context_growth_rate(CultureContext(temperature_c=temperature))
            with pytest.raises(ValueError, match="outside the range its kinetics"):
                self._at(temperature, SETPOINT, **product_kwargs)


# --------------------------------------------------------------------------------------
# 2. MECHANISM
# --------------------------------------------------------------------------------------

class TestTheStateCount:

    def test_five_mechanistic_states_where_the_old_path_had_none(self):
        states = assembled_states(PLATE_READ_4H, ah_out_mM=38.72, sweep=MID)
        assert states.names == ("A_i", "pH_c", "product_fraction",
                                "plasmid_bearing", "generations")
        assert len(states) == 5

    def test_the_census_separates_the_mechanism_from_the_incumbent(self):
        census = ode_state_census(PLATE_READ_4H, ah_out_mM=38.72, sweep=MID)
        assert census["mechanism"] == 5
        assert census["incumbent generator"] == 3
        assert census["total"] == 8

    def test_every_mechanistic_state_is_integrated_not_asserted(self, product_kwargs):
        environment = _environments()["acetic 175 pH4"]
        chain = run_chain(environment, sweep=MID, **product_kwargs)
        assert chain.ph_states is not None
        assert chain.ph_states.system.names == ("A_i", "pH_c")
        assert chain.ph_states.n_rhs_evals > 0
        assert chain.population.system.names == ("plasmid_bearing", "generations")
        assert chain.population.n_rhs_evals > 0

    def test_the_ph_trace_is_a_transient_and_not_a_constant(self, product_kwargs):
        chain = run_chain(_environments()["acetic 175 pH4"], sweep=MID, **product_kwargs)
        trace = chain.ph_states.of("pH_c")
        assert trace[0] > trace[-1]
        assert trace.std() > 0.01


# --------------------------------------------------------------------------------------
# 3. TIMESCALES
# --------------------------------------------------------------------------------------

class TestTheReductionAudit:

    def test_the_plate_and_the_fedbatch_keep_different_states(self):
        plate = reduction_audit(PLATE_READ_4H, ah_out_mM=38.72, sweep=MID)
        vessel = reduction_audit(FEDBATCH_5D, ah_out_mM=38.72, sweep=MID)
        assert len(plate.integrated) == 3
        assert len(vessel.integrated) == 2
        assert set(vessel.integrated.names) < set(plate.integrated.names)

    def test_both_acid_states_reduce_in_both_windows(self):
        """Post-repair. pH_c used to be KEPT on the plate; that was the charge-balance error.

        `mech/ph.py`'s repair conserves charge, which slaves pH_c to A_i exactly, so its
        tau/T falls from 0.691-1.540 to 0.000402-0.000696 and it eliminates in both windows.
        """
        plate = reduction_audit(PLATE_READ_4H, ah_out_mM=38.72, sweep=MID)
        vessel = reduction_audit(FEDBATCH_5D, ah_out_mM=38.72, sweep=MID)
        assert plate.verdict("A_i").value == "eliminate"
        assert vessel.verdict("A_i").value == "eliminate"
        assert plate.verdict("pH_c").value == "eliminate"
        assert vessel.verdict("pH_c").value == "eliminate"

    def test_only_the_plasmid_states_survive_the_vessel(self):
        vessel = reduction_audit(FEDBATCH_5D, ah_out_mM=38.72, sweep=MID)
        assert set(vessel.integrated.names) == {"plasmid_bearing", "generations"}

    def test_the_audit_table_names_the_criterion_and_its_source(self):
        table = reduction_audit(PLATE_READ_4H, ah_out_mM=38.72, sweep=MID).audit_table()
        assert f"{OBSERVED_ACTIVITY_CV:g}" in table
        assert "panel_experiment.OBSERVED_ACTIVITY_CV" in table


# --------------------------------------------------------------------------------------
# 4. GEOMETRY
# --------------------------------------------------------------------------------------

class TestTheGeometry:

    def test_ray_share_is_one_on_a_rank_one_surface_by_construction(self):
        profile = np.array([1.0, 0.7, 0.4, 0.2])
        surface = np.outer([1.0, 2.0, 5.0], profile)
        assert ray_share(surface) == pytest.approx(1.0, abs=1e-12)

    def test_ray_share_refuses_a_degenerate_or_empty_surface(self):
        with pytest.raises(ValueError, match="at least two doses"):
            ray_share(np.ones((1, 5)))
        with pytest.raises(ValueError, match="identically zero"):
            ray_share(np.zeros((3, 5)))

    def test_the_incumbent_surface_is_a_ray_to_three_decimal_places(self):
        doses = [5.0, 10.0, 20.0, 40.0, 60.0, 90.0]
        shares = [ray_share(dose_time_surface("acetic_acid", name, doses,
                                              mechanism=False))
                  for name in transcriptional_reporters()]
        assert min(shares) > 0.995
        assert np.median(shares) == pytest.approx(0.999516, abs=1e-5)

    def test_the_mechanism_moves_the_geometry_on_the_channel_it_reaches(self):
        doses = [5.0, 10.0, 20.0, 40.0, 60.0, 90.0]
        results = geometry_audit(doses, sweep=MID)
        assert len(results) == len(transcriptional_reporters())
        incumbent = np.array([r.incumbent for r in results])
        mechanistic = np.array([r.mechanistic for r in results])
        assert np.median(mechanistic) < np.median(incumbent)
        assert np.median(mechanistic) == pytest.approx(0.983917, abs=1e-4)
        # The number that matters is the complement: energy a latent state could use.
        gain = np.median(1.0 - mechanistic) / np.median(1.0 - incumbent)
        assert gain > 30.0

    def test_the_quench_does_all_of_it_and_the_cell_dynamics_do_none(self):
        """The decomposition, and it is the sharpest negative in this file.

        The geometry moves because a pH-dependent BRIGHTNESS multiplies the trace, not
        because the cell's own dynamics changed shape. Turning the quench off and leaving
        the growth coupling on returns the incumbent's ray share to five decimal places.
        """
        doses = [5.0, 10.0, 20.0, 40.0, 60.0, 90.0]
        kwargs = dict(sweep=MID, mechanism=True)
        both = ray_share(dose_time_surface("acetic_acid", "UPRE-ER", doses, **kwargs))
        quench_only = ray_share(dose_time_surface(
            "acetic_acid", "UPRE-ER", doses, growth_coupling=False, **kwargs))
        growth_only = ray_share(dose_time_surface(
            "acetic_acid", "UPRE-ER", doses, quench=False, **kwargs))
        incumbent = ray_share(dose_time_surface("acetic_acid", "UPRE-ER", doses,
                                                mechanism=False))
        assert quench_only == pytest.approx(both, abs=1e-4)
        assert growth_only == pytest.approx(incumbent, abs=3e-5)
        # Growth alone adds 2.1% of off-ray energy; the quench adds 33x. Was 1.02 before
        # the pH repair made pH_c algebraic, which shifted growth_only by 1.01e-5.
        assert (1.0 - growth_only) / (1.0 - incumbent) < 1.03

    def test_and_it_is_still_a_ray_which_is_the_honest_half(self):
        """0.9758 is not 0.5. The mechanism did not make the surface two-dimensional."""
        results = geometry_audit([5.0, 10.0, 20.0, 40.0, 60.0, 90.0], sweep=MID)
        assert min(r.mechanistic for r in results) > 0.97

    @pytest.mark.slow
    def test_the_panel_wide_median_does_not_move_because_24_of_25_have_no_mechanism(self):
        old, new, moved = [], [], 0
        for stressor in sorted(STRESSORS):
            ladder = [float(d) for d in healthy_ladder(stressor, n_doses=6)]
            for name in transcriptional_reporters():
                incumbent = dose_time_surface(stressor, name, ladder, mechanism=False)
                if np.linalg.norm(incumbent) <= 0.0:
                    continue
                mechanistic = dose_time_surface(stressor, name, ladder, sweep=MID,
                                                mechanism=True)
                old.append(ray_share(incumbent))
                new.append(ray_share(mechanistic))
                moved += int(not np.allclose(incumbent, mechanistic))
        assert moved == len(transcriptional_reporters())
        assert moved / len(old) < 0.05
        assert np.median(new) == pytest.approx(np.median(old), abs=1e-9)


# --------------------------------------------------------------------------------------
# The free-scalar gate, and the API that enforces it
# --------------------------------------------------------------------------------------

class TestTheAssembledGate:

    def test_every_block_passes_alone(self):
        for registry in IN_CHAIN_REGISTRIES:
            if registry is CHAIN_PARAMS:
                continue
            assert registry.gate().passes, registry.piece

    def test_and_the_union_does_not_because_the_floors_are_shared(self):
        gate = assembled_gate()
        assert not gate.passes
        assert len(gate.free) == 7
        assert set(gate.targets) == {"half_yield_transfer", "growth_rate",
                                     "reporter_activity"}

    def test_the_targets_are_de_duplicated_by_name(self):
        """Summing the per-block gates would count the growth assay four times."""
        naive = sum(len(r.independent_targets()) for r in IN_CHAIN_REGISTRIES)
        assert naive == 8
        assert len(assembled_gate().targets) == 3

    def test_asking_for_a_fitted_point_is_refused_with_the_count_in_the_message(self):
        with pytest.raises(FreeScalarGateFailed, match="7 free scalars against 3"):
            require_assembled_gate()

    def test_run_chain_has_no_default_sweep(self, product_kwargs):
        with pytest.raises(TypeError):
            run_chain(_environments()["reference"], **product_kwargs)

    def test_a_sweep_point_off_a_declared_axis_is_refused(self):
        with pytest.raises(ValueError, match="SWEPT over"):
            SweepPoint(beta_mM_per_ph=1.0, cytosolic_volume_ml_per_gdcw=2.0,
                       p_per_division=0.05, burden=0.1)

    def test_the_band_spans_the_corners_and_is_wider_than_a_point(self, product_kwargs):
        band = chain_band(_environments()["acetic 175 pH4"], **product_kwargs)
        assert len(band) == 16
        rates = [r.mu_max_per_h for r in band]
        assert max(rates) > min(rates)


# --------------------------------------------------------------------------------------
# The layer record
# --------------------------------------------------------------------------------------

class TestTheLayerRecord:

    def test_every_declared_layer_reports_exactly_once(self, product_kwargs):
        chain = run_chain(_environments()["acetic 175 pH4"], sweep=MID,
                          construct="UPRE1", **product_kwargs)
        names = [layer.name for layer in chain.layers]
        assert names == list(LAYER_ORDER)

    def test_a_layer_that_could_not_have_moved_the_number_says_inert(self, product_kwargs):
        chain = run_chain(_environments()["reference"], sweep=MID, construct="UPRE1",
                          **product_kwargs)
        assert chain.layer("photophysics").state == LayerState.INERT
        assert not chain.layer("photophysics").sets_the_number

    def test_the_refused_arms_are_reported_and_never_in_chain(self, product_kwargs):
        chain = run_chain(_environments()["reference"], sweep=MID, **product_kwargs)
        for name in ("allocation", "ethanol"):
            assert chain.layer(name).state == LayerState.REPORTED
            assert not chain.layer(name).sets_the_number

    def test_the_not_run_rows_say_what_would_run_them(self, product_kwargs):
        chain = run_chain(_environments()["reference"], sweep=MID, **product_kwargs)
        for layer in chain.layers:
            if layer.state == LayerState.NOT_RUN:
                assert "pass " in layer.detail

    def test_the_layer_report_prints_one_line_per_layer(self, product_kwargs):
        chain = run_chain(_environments()["reference"], sweep=MID, **product_kwargs)
        assert len(chain.layer_report().splitlines()) == len(LAYER_ORDER)

    def test_the_two_ph_arms_are_reported_as_a_bracket_not_as_one_answer(
            self, product_kwargs):
        chain = run_chain(_environments()["acetic 175 pH4"], sweep=MID, **product_kwargs)
        assert any("bracket" in note for note in chain.notes)
        assert chain.proton_bill.trapped_anion_mM > 0.0
        assert chain.ph_c_final < 7.0

    def test_the_mechanism_replaces_the_panel_curve_on_the_acid_and_keeps_it_elsewhere(
            self, product_kwargs):
        acid = run_chain(_environments()["acetic 175 pH4"], sweep=MID, **product_kwargs)
        peroxide = run_chain(_environments()["H2O2"], sweep=MID, **product_kwargs)
        assert acid.layer("stress panel").state == LayerState.INERT
        assert peroxide.layer("stress panel").state == LayerState.IN_CHAIN
        assert peroxide.layer("stress panel").sets_the_number


# --------------------------------------------------------------------------------------
# The reporter arm and the quench
# --------------------------------------------------------------------------------------

class TestTheReporterArm:

    def test_the_incumbent_generator_does_the_integration(self):
        environment = Environment(context=CultureContext(ph_medium=4.5),
                                  stressor="acetic_acid", dose=60.0)
        chain = run_chain(environment, sweep=MID, window=PLATE_READ_4H,
                          construct="UPRE1")
        assert chain.reporter.construct == "UPRE1"
        assert chain.reporter.times_h[-1] == pytest.approx(PLATE_READ_4H.duration_h)
        assert np.all(chain.reporter.activity > 0)

    def test_the_quench_clears_the_named_reporter_floor_at_every_corner(self):
        environment = Environment(context=CultureContext(ph_medium=4.5),
                                  stressor="acetic_acid", dose=60.0)
        effects = []
        for corner in SweepPoint.corners():
            chain = run_chain(environment, sweep=corner, window=PLATE_READ_4H,
                              construct="UPRE1")
            effects.append(abs(chain.reporter.observed[-1]
                               / chain.reporter.activity[-1] - 1.0))
        assert min(effects) / OBSERVED_ACTIVITY_CV > 3.0

    def test_with_no_acid_the_quench_is_exactly_one(self):
        chain = run_chain(Environment(), sweep=MID, window=PLATE_READ_4H,
                          construct="UPRE1")
        assert np.all(chain.reporter.quench == 1.0)
        assert np.array_equal(chain.reporter.observed, chain.reporter.activity)

    def test_the_nutrient_factor_is_the_only_seam_into_the_generator(self):
        """The mechanism enters simulate_culture through its own declared knob."""
        chain = run_chain(Environment(), sweep=MID, window=PLATE_READ_4H,
                          construct="UPRE1")
        assert chain.reporter.nutrient_factor == pytest.approx(1.0)
        taxed = run_chain(Environment(), sweep=MID, window=PLATE_READ_4H,
                          construct="UPRE1", n_copies=20.0)
        assert 0.0 < taxed.reporter.nutrient_factor < 1.0

    def test_an_unknown_construct_is_refused_by_name(self):
        with pytest.raises(KeyError, match="DEFAULT_PANEL"):
            run_chain(Environment(), sweep=MID, window=PLATE_READ_4H,
                      construct="not-a-strain")


# --------------------------------------------------------------------------------------
# Guards against drift
# --------------------------------------------------------------------------------------

class TestNothingDrifts:

    def test_the_restated_basal_matches_the_generators_private_constant(self):
        from ystwin.generator import stress_panel

        from ystwin.mech.chain import _BASAL_ACTIVITY

        assert _BASAL_ACTIVITY == stress_panel._DEFAULT_BASAL

    def test_the_product_arm_refuses_a_partial_call(self, spec):
        with pytest.raises(ValueError, match="together"):
            run_chain(Environment(growth_rate_setpoint_per_h=SETPOINT), sweep=MID,
                      spec=spec)

    def test_the_product_arm_refuses_a_spec_with_no_genotype(self, spec):
        with pytest.raises(ValueError, match="genotype"):
            run_chain(Environment(growth_rate_setpoint_per_h=SETPOINT), sweep=MID,
                      spec=spec, flux_calibration=calibrations.BETA_CAROTENE_FLUX,
                      kinetics=calibrations.BETA_CAROTENE_KINETICS)

    def test_chain_band_refuses_a_sweep_argument(self):
        with pytest.raises(TypeError, match="runs the corners itself"):
            chain_band(Environment(), sweep=MID)

    def test_the_burden_layer_never_claims_criterion_one(self, product_kwargs):
        chain = run_chain(_environments()["reference"], sweep=MID, **product_kwargs)
        assert "does not claim criterion 1" in chain.layer("burden").detail

    def test_provenance_names_the_gate_and_every_in_chain_registry(self):
        text = provenance()
        assert "7 free scalars against 3 independent targets" in text
        for registry in IN_CHAIN_REGISTRIES:
            assert registry.piece in text

    def test_the_result_is_reproducible(self, product_kwargs):
        first = run_chain(_environments()["acetic 175 pH4"], sweep=MID, **product_kwargs)
        second = run_chain(_environments()["acetic 175 pH4"], sweep=MID, **product_kwargs)
        assert (first.product.population_content_mg_per_gdcw
                == second.product.population_content_mg_per_gdcw)
        assert first.ph_c_final == second.ph_c_final


# --------------------------------------------------------------------------------------
# 5. DEPTH -- the GEM half, which needs the shipped GSMM
# --------------------------------------------------------------------------------------

@pytest.mark.integration
class TestOneInputReachesManyReactions:
    """Criterion 4. The baseline is ONE reaction, at depth 1, in one script."""

    @pytest.fixture(scope="class")
    def _acid_model_template(self, yeast_gem_factory):
        from ystwin.fba.stress_ph import add_weak_acid_uncoupling

        return add_weak_acid_uncoupling(yeast_gem_factory(), "acetic")

    @pytest.fixture
    def acid_model(self, _acid_model_template, model_copy):
        return model_copy(_acid_model_template)

    def test_an_acid_load_moves_hundreds_of_reactions(self, acid_model):
        from ystwin.fba.stress_ph import constrain_proton_budget

        with acid_model as model:
            constrain_proton_budget(model)
            moved = reachable_reactions(model, "acetic", 3.0)
        assert len(moved) > 400

    def test_they_sit_at_several_depths_from_the_entry_step(self, acid_model):
        from ystwin.fba.stress_ph import constrain_proton_budget

        with acid_model as model:
            constrain_proton_budget(model)
            depths = reaction_depths(model, "acetic", 3.0)
        assert max(depths) >= 4
        assert sum(depths.values()) > 400
        assert depths[0] >= 1

    def test_the_moved_set_is_degenerate_and_the_count_must_not_be_pinned(self, acid_model):
        """pFBA has alternate optima here, so the exact count is not a fact about the cell.

        Measured over repeats in one process: 484 on the first call and 487 after, union 494,
        intersection 477. This test exists so that the two inequalities above are never
        "tightened" into an equality the solver does not owe anyone. What IS stable is
        asserted: the depth ceiling, the compartment count, and the size of the core.
        """
        from ystwin.fba.stress_ph import constrain_proton_budget

        runs = []
        for _ in range(3):
            with acid_model as model:
                constrain_proton_budget(model)
                runs.append(set(reachable_reactions(model, "acetic", 3.0)))
        union, core = set.union(*runs), set.intersection(*runs)
        # The stable core is most of it, and the instability is real but small.
        assert len(core) > 400
        assert len(union - core) < 0.15 * len(union)
        compartments = {met.compartment for rid in core
                        for met in acid_model.reactions.get_by_id(rid).metabolites}
        assert len(compartments) >= 6

    def test_the_depth_ceiling_is_stable_even_though_the_count_is_not(self, acid_model):
        """The histogram's SHAPE survives the degeneracy: mass at depths 1-2, ceiling 6."""
        from ystwin.fba.stress_ph import constrain_proton_budget

        for _ in range(2):
            with acid_model as model:
                constrain_proton_budget(model)
                depths = reaction_depths(model, "acetic", 3.0)
            assert max(depths) >= 5
            assert depths[1] + depths[2] > 0.7 * sum(depths.values())

    def test_the_layer_reports_it_when_a_model_is_supplied(self, acid_model,
                                                           product_kwargs):
        from ystwin.fba.stress_ph import constrain_proton_budget

        with acid_model as model:
            constrain_proton_budget(model)
            chain = run_chain(_environments()["acetic 175 pH4"], sweep=MID,
                              gem_model=model, **product_kwargs)
        layer = chain.layer("GEM depth")
        assert layer.state == LayerState.AUDITS
        assert not layer.sets_the_number
