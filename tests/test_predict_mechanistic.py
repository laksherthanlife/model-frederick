"""`predict_product(mech=...)`: what the mechanistic layer buys, measured on the product.

Two claims are pinned here and the second is a negative.

1. **`mech=None` is the incumbent, bit-for-bit.** Every leave-one-strain-out figure this
   package quotes was measured with `mech/` absent, so the seam has to be provably silent
   when nobody opens it -- the content, the layer report and the notes.

2. **With the seam open, the mechanism moves the product through exactly one channel:**
   ``mu = min(mu_set, mu_max)``. `predict.py::_growth_rate` checked reachability on the
   stressor branch only, so a fed-batch setpoint above what the context can support was
   returned unchecked -- which is why seven environments returned one content to eight
   decimal places. Fix that and one of seven moves. Six still return 1.07098788 mg/gDCW to
   eight decimals, and this file says so with numbers rather than leaving the reader to
   assume the layer did more.

The sign bug is here too, because nothing in the suite covered `stress_modules` and the
field dropped every module the panel scored negative -- acetic acid's whole pH arm among
them.
"""

from __future__ import annotations

import pytest

from ystwin.generator.context import CultureContext
from ystwin.generator.stress_panel import STRESSORS, module_response
from ystwin.mech.chain import SweepPoint, run_chain
from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS
from ystwin.pathway.spec import load_pathway
from ystwin.predict import (
    Environment,
    Genotype,
    LayerState,
    SetpointUnreachable,
    predict_product,
)

SETPOINT = 0.18
"""Inside the kinetics' fitted window [0.101, 0.2543] and above the ethanol context's own
maximum of 0.140 /h, which is the only reason any of this is visible."""

# What each condition returns with mech ON at mu_set = 0.18 /h, mg/gDCW. Six of the eight
# are the incumbent's own number to eight decimals; the run that produced them is in the
# report for this task.
EXPECTED_MG_PER_GDCW = {
    "reference": 1.07098788,
    "ethanol": 1.11791951,
    "37 C": 1.07098788,
    "pH 6": 1.07098788,
    "near-anoxia": 1.07098788,
    "H2O2": 1.07098788,
    "acetic acid": 1.07098788,
    "DTT": 1.11977419,
}

HELD_SETPOINT_CONTENT = 1.0709878838715265
"""The whole double behind the 1.07098788 above.

Bit-identity is the point: the acid rows claimed to set a number that does not move in its
last digit across the entire dose axis.
"""

CONDITIONS = {
    "reference": (CultureContext(), None, 0.0),
    "ethanol": (CultureContext(carbon_source="ethanol"), None, 0.0),
    "37 C": (CultureContext(temperature_c=37.0), None, 0.0),
    "pH 6": (CultureContext(ph_medium=6.0), None, 0.0),
    "near-anoxia": (CultureContext(oxygen=0.005), None, 0.0),
    "H2O2": (CultureContext(), "H2O2", 1.0),
    "acetic acid": (CultureContext(ph_medium=4.5), "acetic_acid", 60.0),
    "DTT": (CultureContext(), "DTT", 2.0),
}


@pytest.fixture(scope="module")
def spec():
    return load_pathway("beta_carotene")


@pytest.fixture(scope="module")
def sweep():
    """The centre of the four declared axes. A smoke corner, and every number below is
    checked to be invariant to the choice before it is quoted."""
    return SweepPoint.midpoint()


def _environment(name: str) -> Environment:
    context, stressor, dose = CONDITIONS[name]
    return Environment(context, stressor=stressor, dose=dose,
                       growth_rate_setpoint_per_h=SETPOINT)


def _predict(spec, name, **kwargs):
    result = predict_product(spec, Genotype(1.0, "b-car4"), _environment(name),
                             BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS, **kwargs)
    if kwargs.get("mode") == "legacy":
        assert result.mode == "legacy"
        assert not result.supported
        assert not hasattr(result, "content_mg_per_gdcw")
        assert result.calculation.mode == "legacy"
        return result.calculation
    return result


TOP_OF_FITTED_WINDOW = 0.2543
"""Just inside `BETA_CAROTENE_KINETICS["lycopene"].growth_rate_range` = [0.101, 0.25432].

The highest setpoint `pathway/solve.py` admits, so the first setpoint any acid dose can
clip -- which is why the acid arm's reach is measured against this and not only against
`SETPOINT`.
"""

LETHAL_DOSE_MM = STRESSORS["acetic_acid"].lethal_dose
"""180 mM, read off the panel rather than typed, because every dose below distinguishes a
measured statement from an arithmetic one."""


def _acid(spec, sweep, ph_medium, dose, setpoint=SETPOINT):
    """The acid condition off the `CONDITIONS` table, so dose and medium pH can be swept.

    `dose=0` drops the stressor rather than passing a zero one: `mech/chain.py` routes a
    zero dose down its no-acid branch anyway, and this keeps the [AH]_o = 0 end honest.
    """
    return predict_product(
        spec, Genotype(1.0, "b-car4"),
        Environment(CultureContext(ph_medium=ph_medium),
                    stressor="acetic_acid" if dose > 0.0 else None, dose=dose,
                    growth_rate_setpoint_per_h=setpoint),
        BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS, mode="legacy", mech=sweep).calculation


def _crossing_dose(sweep, ph_medium, target_per_h, lo=1.0, hi=3000.0):
    """The acetate dose at which mu_max falls through `target_per_h`, by bisection.

    On the growth arm alone: past the crossing `pathway/solve.py` refuses the rate outright
    as outside its fitted window, so the product arm cannot be in the loop.
    """
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        chain = run_chain(
            Environment(CultureContext(ph_medium=ph_medium), stressor="acetic_acid",
                        dose=mid, growth_rate_setpoint_per_h=SETPOINT), sweep=sweep)
        lo, hi = (mid, hi) if chain.mu_max_per_h > target_per_h else (lo, mid)
    return 0.5 * (lo + hi)


class TestTheNegativeModulesAreNoLongerDropped:
    """`stress_modules` filtered on ``v > 0.01``, which is one-sided.

    The panel scores a module's activity as a signed deviation, and five of its agents carry
    their most carefully sourced arm on the negative side: acetic acid's cytosolic pH, DTT's
    redox pool, glucose starvation's and antimycin A's ATP, menadione's NADH. A one-sided
    filter deleted all five from the field a caller reads, so a user asking "what does 60 mM
    acetic acid do" got the ESR and the cell wall and no mention of pH at all.
    """

    @pytest.mark.parametrize("stressor,dose,module", [
        ("acetic_acid", 60.0, "ph"),
        ("DTT", 1.0, "redox"),
        ("antimycin_A", 5.0, "atp"),
        ("glucose_starvation", 0.6, "atp"),
        ("menadione", 100.0, "nadh"),
    ])
    def test_a_negative_arm_survives_into_the_field(self, spec, stressor, dose, module):
        got = predict_product(
            spec, Genotype(1.0, "b-car4"),
            Environment(stressor=stressor, dose=dose, growth_rate_setpoint_per_h=SETPOINT),
            BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS)

        assert got.stress_modules[module] < -0.01

    def test_the_acid_arm_is_the_largest_module_it_carries(self, spec):
        """And it was the missing one: 0.696 against the ESR's 0.161."""
        got = _predict(spec, "acetic acid", mode="legacy")

        assert max(got.stress_modules, key=lambda k: abs(got.stress_modules[k])) == "ph"
        assert got.stress_modules["ph"] == pytest.approx(-0.6961, abs=5e-4)

    def test_the_threshold_is_still_a_threshold_and_not_removed(self, spec):
        """`abs` widened the filter; it must not have deleted it."""
        got = _predict(spec, "acetic acid", mode="legacy")
        full = module_response("acetic_acid", 60.0)

        assert {k for k, v in full.items() if abs(v) > 0.01} == set(got.stress_modules)
        assert all(abs(v) > 0.01 for v in got.stress_modules.values())
        assert len(got.stress_modules) < len(full)


class TestTheSeamIsSilentWhenNobodyOpensIt:
    """Default off, and "off" means indistinguishable rather than merely similar."""

    @pytest.mark.parametrize("name", [n for n in CONDITIONS if n != "DTT"])
    def test_default_keeps_supported_content_and_refuses_unsupported_conditions(self, spec, name):
        if name == "ethanol":
            with pytest.raises(SetpointUnreachable):
                _predict(spec, name)
        elif name in ("pH 6", "acetic acid"):
            with pytest.raises(ValueError, match="ph_medium"):
                _predict(spec, name)
        else:
            assert _predict(spec, name).content_mmol_per_gdcw == (
                _predict(spec, "reference").content_mmol_per_gdcw)

    def test_and_the_layer_report_gains_no_row(self, spec):
        got = _predict(spec, "reference")

        assert not [row for row in got.layers if row.name.startswith("mech ")]
        assert got.mech_chain is None

    def test_the_refusal_channel_is_untouched_with_the_seam_shut(self, spec):
        """DTT at 2 mM leaves mu_max = 0.138 /h under a 0.18 /h setpoint."""
        with pytest.raises(SetpointUnreachable):
            _predict(spec, "DTT")


class TestWhatTheMechanismMovesAndWhatItDoesNot:
    """The measured answer, and most of it is "nothing"."""

    @pytest.mark.parametrize("name", sorted(CONDITIONS))
    def test_every_condition_returns_its_measured_content(self, spec, sweep, name):
        got = _predict(spec, name, mode="legacy", mech=sweep)

        assert got.content_mg_per_gdcw == pytest.approx(
            EXPECTED_MG_PER_GDCW[name], abs=5e-9)

    def test_six_of_the_eight_are_still_one_number_to_eight_decimals(self, spec, sweep):
        """THE HONEST RESULT. Wiring 11 mechanistic modules into the product channel left
        six of eight environments returning the identical content, because their mu_max
        stays above the setpoint and mu is the only route the mechanism has to a pool."""
        contents = {name: round(_predict(spec, name, mode="legacy", mech=sweep).content_mg_per_gdcw, 8)
                    for name in CONDITIONS}
        unchanged = [n for n, v in contents.items() if v == 1.07098788]

        assert sorted(unchanged) == ["37 C", "H2O2", "acetic acid", "near-anoxia",
                                     "pH 6", "reference"]

    def test_the_one_environment_that_moves_moves_by_1_0438x(self, spec, sweep):
        """Ethanol, and only because the context caps that culture at 0.140 /h while the
        setpoint asks for 0.180 -- the check `_growth_rate` skipped on this branch."""
        moved = _predict(spec, "ethanol", mode="legacy", mech=sweep)
        base = _predict(spec, "reference", mode="legacy", mech=sweep)

        assert moved.growth_rate_per_h == pytest.approx(0.140, abs=1e-9)
        assert (moved.content_mg_per_gdcw / base.content_mg_per_gdcw
                == pytest.approx(1.0438, abs=5e-5))

    def test_medium_ph_reaches_the_growth_rate_and_never_the_content(self, spec, sweep):
        """pH moves [AH]_o by 11.9x and mu_max by 0.050 /h, and the content not at all.

        `generator/context.py` carries no pH term at all, so pH's only route is the acid
        arm -- and at THIS dose that arm leaves mu_max above every setpoint the kinetics
        admit. The next test says why that is a dose-scoped fact and not a structural one.
        """
        acid = Genotype(1.0, "b-car4")
        low = predict_product(
            spec, acid, Environment(CultureContext(ph_medium=4.5), stressor="acetic_acid",
                                    dose=60.0, growth_rate_setpoint_per_h=SETPOINT),
            BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS, mode="legacy", mech=sweep).calculation
        high = predict_product(
            spec, acid, Environment(CultureContext(ph_medium=6.0), stressor="acetic_acid",
                                    dose=60.0, growth_rate_setpoint_per_h=SETPOINT),
            BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS, mode="legacy", mech=sweep).calculation

        assert low.mech_chain.ah_out_mM / high.mech_chain.ah_out_mM == pytest.approx(
            11.86, abs=0.02)
        assert low.mech_chain.mu_max_per_h == pytest.approx(0.3444, abs=5e-5)
        assert high.mech_chain.mu_max_per_h == pytest.approx(0.3946, abs=5e-5)
        assert low.content_mmol_per_gdcw == high.content_mmol_per_gdcw

    def test_the_acid_arm_cannot_bite_at_this_dose_and_this_medium_ph(self, spec, sweep):
        """A DOSE-SCOPED statement. It read as a universal until 2026-09-09, and as a
        universal it is false.

        At this file's 60 mM and medium pH 4.5 the acid leaves mu_max at 0.344 /h, above
        the top of the window the kinetics were fitted over (0.2543 /h). A clip needs a
        setpoint above mu_max, `pathway/solve.py` admits none that high, so the pH block --
        the two ODE states, the bill, the speciation -- cannot reach this content HERE.

        mu_max falls monotonically with dose, though, so "whatever the dose" was never
        true: the crossings are bisected in the next test but one, and past them the arm
        does move the number. Read as a universal this retires the acid arm, which is the
        only mechanistic route to this content there is; the scoped version is what
        docs/MECHANISTIC_LAYER.md:129-131 already states.
        """
        got = _predict(spec, "acetic acid", mode="legacy", mech=sweep)
        _, highest_fitted = BETA_CAROTENE_KINETICS["lycopene"].growth_rate_range

        assert got.mech_chain.mu_max_per_h > highest_fitted
        assert got.mech_chain.feed_sets_mu

    def test_and_past_the_crossing_dose_it_does_bite(self, spec, sweep):
        """The counterexample, with the caveat that stops it being a titre prediction.

        At medium pH 4.5 / 500 mM the feed lets go and the content moves to 1.0823 mg/gDCW
        against the 1.0710 pinned above. What that settles is only that the pH block's
        route to the content is OPEN -- not that 1.0823 is right. `mech/ph.py` flags the
        trapped anion pool as non-physical at both doses (8090 mM at 60, 67417 mM at 500),
        and 500 mM is past the panel's own declared lethal dose, so the arithmetic runs
        well outside anything measured.
        """
        got = _acid(spec, sweep, 4.5, 500.0)

        assert not got.mech_chain.feed_sets_mu
        assert got.content_mg_per_gdcw == pytest.approx(1.082302, abs=5e-6)
        assert got.content_mg_per_gdcw != pytest.approx(HELD_SETPOINT_CONTENT, abs=1e-6)
        assert 500.0 > LETHAL_DOSE_MM
        assert not got.mech_chain.proton_bill.physical

    @pytest.mark.parametrize("ph_medium,target,expected_mM,sub_lethal", [
        (4.5, TOP_OF_FITTED_WINDOW, 213.11, False),
        (4.5, SETPOINT, 454.61, False),
        (4.0, TOP_OF_FITTED_WINDOW, 161.43, True),
        (4.0, SETPOINT, 344.37, False),
    ])
    def test_the_crossing_doses_are_bisected_and_not_asserted(
            self, sweep, ph_medium, target, expected_mM, sub_lethal):
        """Pinning them is what keeps the docstring above scoped rather than universal.

        One of the four sits below the panel's declared lethal dose: the medium pH 4
        crossing of the fitted ceiling, 161.43 mM, which is the figure
        docs/MECHANISTIC_LAYER.md quotes and the reason the acid arm must not be retired on
        the strength of the old docstring. The other three are arithmetic past a dose the
        panel already calls lethal.
        """
        assert _crossing_dose(sweep, ph_medium, target) == pytest.approx(
            expected_mM, abs=0.01)
        assert (expected_mM < LETHAL_DOSE_MM) is sub_lethal

    def test_the_content_is_invariant_to_all_sixteen_declared_corners(self, spec):
        """What licenses quoting a number out of a REFUSED free-scalar gate.

        `mech/chain.py::assembled_gate` counts 7 free scalars against 3 targets and refuses.
        Buffering capacity and cytosolic volume enter only the acid arm; loss-per-division
        and burden only the bearing fraction. None reaches the content, and this is the
        assertion the `mech free scalars` layer's detail rests on.
        """
        for name in ("reference", "ethanol", "DTT", "acetic acid"):
            contents = {round(_predict(spec, name, mode="legacy", mech=corner).content_mmol_per_gdcw, 12)
                        for corner in SweepPoint.corners()}
            assert len(contents) == 1, name

    def test_the_bearing_fraction_is_where_those_scalars_do_land(self, spec):
        """They are not inert, they are aimed elsewhere: 1.026x across the corners."""
        titres = [_predict(spec, "reference", mode="legacy", mech=corner
                           ).mech_chain.product.population_content_mg_per_gdcw
                  for corner in SweepPoint.corners()]

        assert max(titres) / min(titres) == pytest.approx(1.0263, abs=5e-4)


class TestTheVesselRefusalBecomesAReportedRate:
    """With `mech` on, a setpoint is a ceiling on both branches rather than a refusal on one."""

    def test_dtt_returns_a_number_instead_of_raising(self, spec, sweep):
        got = _predict(spec, "DTT", mode="legacy", mech=sweep)

        assert got.growth_rate_per_h == pytest.approx(0.138349, abs=1e-6)
        assert not got.mech_chain.feed_sets_mu

    def test_and_says_so_where_a_reader_will_see_it(self, spec, sweep):
        got = _predict(spec, "DTT", mode="legacy", mech=sweep)

        assert any("MECHANISM CLIPPED THE SETPOINT" in note for note in got.notes)

    def test_the_stress_panel_goes_back_into_the_chain_when_the_feed_lets_go(
            self, spec, sweep):
        """It is AUDITS under a held setpoint because mu is the feed's; when the dose takes
        mu_max below the setpoint the panel sets the rate the pools are solved at."""
        assert _predict(spec, "DTT", mode="legacy", mech=sweep).layer(
            "stress panel").state == LayerState.IN_CHAIN
        assert _predict(spec, "reference", mode="legacy", mech=sweep).layer(
            "stress panel").state == LayerState.AUDITS

    def test_the_by_construction_note_is_withdrawn_when_it_stops_being_true(
            self, spec, sweep):
        """That note says the stressed number is identical to the unstressed one BY
        CONSTRUCTION. Under a clip it is not, and a stale note here is the exact failure
        this module's layer record exists to prevent."""
        held = _predict(spec, "acetic acid", mode="legacy", mech=sweep)
        clipped = _predict(spec, "DTT", mode="legacy", mech=sweep)

        assert any("BY CONSTRUCTION" in note for note in held.notes)
        assert not any("BY CONSTRUCTION" in note for note in clipped.notes)


class TestTheLayerRecordStaysHonest:
    """The mech rows go through the existing mechanism, and mean what it means."""

    def test_every_mech_row_is_prefixed_and_no_name_is_claimed_twice(self, spec, sweep):
        got = _predict(spec, "acetic acid", mode="legacy", mech=sweep)
        names = [row.name for row in got.layers]

        assert len(names) == len(set(names))
        assert len([n for n in names if n.startswith("mech ")]) == 15

    def test_the_incumbent_rows_still_come_first_and_in_their_own_order(self, spec, sweep):
        with_mech = [row.name for row in _predict(spec, "reference", mode="legacy", mech=sweep).layers]
        without = [row.name for row in _predict(spec, "reference").layers]

        assert with_mech[:len(without)] == without

    def test_in_chain_still_means_the_returned_number_depends_on_it(self, spec, sweep):
        """`mech/chain.py` calls its population and pH-state rows IN_CHAIN, and for ITS
        headline -- population titre -- they are. This class returns content per bearing
        cell, so those rows are demoted rather than copied as a false claim."""
        got = _predict(spec, "acetic acid", mode="legacy", mech=sweep)

        for row in got.layers:
            if row.state == LayerState.IN_CHAIN:
                assert row.sets_the_number, row.name
        assert got.layer("mech population").state == LayerState.REPORTED
        assert "calls this in-chain" in got.layer("mech population").detail

    def test_the_refused_free_scalar_gate_is_attached_rather_than_swallowed(
            self, spec, sweep):
        got = _predict(spec, "reference", mode="legacy", mech=sweep)

        assert "REFUSED" in got.layer("mech free scalars").detail
        assert "7 free scalars against 3 independent targets" in got.layer(
            "mech free scalars").detail

    def test_the_report_still_has_one_line_per_layer(self, spec, sweep):
        got = _predict(spec, "reference", mode="legacy", mech=sweep)

        assert len(got.layer_report().splitlines()) == len(got.layers)


class TestSetsTheNumberIsDecidedPerRunAndNotByAName:
    """`sets_the_number` is a per-run claim everywhere else in this package, and was not here.

    `mech/chain.py` decides it against the run: the stress panel's row is True only when
    ``panel_retained < 1.0`` and the burden row only when ``n_copies != REFERENCE_COPIES``.
    `predict.py` recomputed the copies from a static name allowlist instead, which asks a
    question about the layer rather than about the prediction. Under a held setpoint that
    printed ``mech stress panel ... sets=False`` and, twelve rows later,
    ``mech weak acid ... sets=True`` -- for the same mu, in the same report, with the
    content bit-identical across the entire dose axis.

    The authoritative per-run signal was already in hand ten lines below the bug:
    ``ChainResult.feed_sets_mu``. When the feed holds mu at the setpoint, mu_max is not what
    the pools are solved at, so nothing upstream of mu_max reaches the content.
    """

    MU_UPSTREAM = ("mech environment", "mech stress panel", "mech weak acid",
                   "mech proton bill", "mech burden")
    """The rows whose only route to the content is mu_max, so the feed can cut all five."""

    @pytest.mark.parametrize("dose", [0.0, 30.0, 60.0, 120.0, 200.0, 300.0])
    def test_a_held_setpoint_makes_the_content_deaf_to_the_whole_dose_axis(
            self, spec, sweep, dose):
        """The measurement the claim answers to: bit-identical, not merely close."""
        got = _acid(spec, sweep, 4.5, dose)

        assert got.mech_chain.feed_sets_mu
        assert got.content_mg_per_gdcw == HELD_SETPOINT_CONTENT

    def test_and_the_dose_axis_it_is_deaf_to_is_a_real_one(self, spec, sweep):
        """[AH]_o runs 0 -> 193.6 mM and mu_max 0.400 -> 0.221 /h over that same sweep."""
        ends = [_acid(spec, sweep, 4.5, dose).mech_chain for dose in (0.0, 300.0)]

        assert [round(c.ah_out_mM, 3) for c in ends] == [0.0, 193.606]
        assert [round(c.mu_max_per_h, 4) for c in ends] == [0.4, 0.2214]

    @pytest.mark.parametrize("name", MU_UPSTREAM)
    def test_no_row_upstream_of_mu_claims_the_content_while_the_feed_holds_it(
            self, spec, sweep, name):
        """Was False, False, True, True, False -- three of them false claims."""
        got = _acid(spec, sweep, 4.5, 60.0)
        row = got.layer(name)

        assert got.mech_chain.feed_sets_mu
        assert not row.sets_the_number, row.detail
        assert row.state != LayerState.IN_CHAIN, row.detail

    def test_the_report_no_longer_contradicts_itself_row_to_row(self, spec, sweep):
        """The reported symptom was one report carrying both verdicts about one mu."""
        got = _acid(spec, sweep, 4.5, 60.0)
        claimed = [row.name for row in got.layers
                   if row.name.startswith("mech ") and row.sets_the_number]

        assert claimed == ["mech vessel", "mech metabolism"]

    @pytest.mark.parametrize("name", ["37 C", "near-anoxia"])
    def test_the_context_row_carried_the_same_false_true(self, spec, sweep, name):
        """`mech environment` is hardcoded True in the chain, and the chain is right about
        ITS headline. It is not right about this content while the feed sets mu."""
        got = _predict(spec, name, mode="legacy", mech=sweep)

        assert got.mech_chain.feed_sets_mu
        assert got.content_mg_per_gdcw == pytest.approx(
            EXPECTED_MG_PER_GDCW["reference"], abs=5e-9)
        assert not got.layer("mech environment").sets_the_number

    def test_a_demoted_row_still_says_what_the_chain_called_it(self, spec, sweep):
        """Demotion, not deletion: the chain's own verdict stays readable in the detail."""
        got = _acid(spec, sweep, 4.5, 60.0)

        assert got.layer("mech weak acid").state == LayerState.REPORTED
        assert "calls this in-chain" in got.layer("mech weak acid").detail

    def test_the_acid_rows_claim_it_again_once_the_dose_takes_mu_max_below_the_setpoint(
            self, spec, sweep):
        """THE REASON THE ALLOWLIST MUST NOT SIMPLY LOSE THE ACID NAMES.

        At 300 mM and medium pH 4.5 the feed lets go, mu becomes mu_max, and the acid arm
        does set this content: 1.0191 mg/gDCW against 0.9766 at the same held setpoint.
        A fix that dropped `weak acid` and `proton bill` from the allowlist would print
        False here and be wrong in the other direction.
        """
        bound = _acid(spec, sweep, 4.5, 300.0, setpoint=TOP_OF_FITTED_WINDOW)
        held = _acid(spec, sweep, 4.5, 0.0, setpoint=TOP_OF_FITTED_WINDOW)

        assert not bound.mech_chain.feed_sets_mu
        assert bound.mech_chain.mu_max_per_h == pytest.approx(0.221417, abs=5e-6)
        assert bound.content_mg_per_gdcw == pytest.approx(1.019132, abs=5e-6)
        assert held.content_mg_per_gdcw == pytest.approx(0.976564, abs=5e-6)
        for name in ("mech environment", "mech weak acid", "mech proton bill"):
            assert bound.layer(name).sets_the_number, name
            assert bound.layer(name).state == LayerState.IN_CHAIN, name

    @pytest.mark.parametrize("name", ["mech vessel", "mech metabolism"])
    def test_the_vessel_and_the_solve_claim_it_on_both_branches(self, spec, sweep, name):
        """They are downstream of the min, so no branch of it can cut them."""
        held = _acid(spec, sweep, 4.5, 60.0)
        bound = _acid(spec, sweep, 4.5, 300.0, setpoint=TOP_OF_FITTED_WINDOW)

        assert held.layer(name).sets_the_number
        assert bound.layer(name).sets_the_number

    def test_the_panel_row_still_answers_to_the_chain_when_the_feed_lets_go(
            self, spec, sweep):
        """The gate narrows the chain's verdict and never widens it: DTT clips, so the
        panel's own ``panel_retained < 1.0`` decides, exactly as `mech/chain.py` computed."""
        got = _predict(spec, "DTT", mode="legacy", mech=sweep)

        assert not got.mech_chain.feed_sets_mu
        assert got.layer("mech stress panel").sets_the_number
        assert got.layer("mech stress panel").state == LayerState.IN_CHAIN
