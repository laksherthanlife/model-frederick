"""Move 4: product forecasting comes off FBA and onto kinetics.

D1 showed the constraint-based model gives [0, ceiling] and never a point value.
A kinetic branch does give one, and it also predicts something FBA cannot: the
distribution of carbon across phytoene, lycopene and beta-carotene. Those three
come off one HPLC run, so the intermediates are a nearly free independent anchor
on which enzyme is limiting.

Two structural facts the model has to carry:
  * crtYB is bifunctional. Phytoene synthase and lycopene cyclase draw on one
    enzyme pool, so pushing one activity starves the other.
  * intracellular pools are diluted by growth, exactly as the reporter is.

The second half of this file is new and its job is different. The branch is no longer
parked -- Elizondo & Saa 2025 supplies six chemostat steady states -- and a calibration set
does two things to a model, only one of which is comfortable. It fits the part that is
identifiable, and it *refutes* the part that is not. Both are tested here, and the
refutation tests are deliberately written so that reverting the module to a
growth-rate-independent cyclase fails them.
"""

import numpy as np
import pytest

from ystwin.kinetic.carotenoid import (
    ELIZONDO2025,
    UNIDENTIFIABLE_STEPS,
    CarotenoidKinetics,
    beta_carotene_content,
    branch_rates,
    calibrated_cyclase,
    calibrated_kinetics,
    calibration_states,
    limiting_step,
    predict_beta_carotene_rate,
    simulate_branch,
    solve_branch_from_flux,
    steady_state_pools,
)

BALANCED = CarotenoidKinetics(
    vmax_psy=0.05, km_psy=0.01,
    vmax_crti=0.05, km_crti=0.01,
    vmax_lcy=0.05, km_lcy=0.01,
)


def test_zero_precursor_supply_gives_no_product():
    pools = steady_state_pools(ggpp_supply=0.0, kinetics=BALANCED, growth_rate=0.2)

    assert pools.beta_carotene_flux == pytest.approx(0.0, abs=1e-12)


def test_at_steady_state_product_flux_matches_the_supply_it_can_carry():
    """With no bottleneck and no dilution, carbon in equals carbon out."""
    supply = 0.004  # mmol GGPP / gDCW / h; two GGPP make one phytoene
    pools = steady_state_pools(ggpp_supply=supply, kinetics=BALANCED, growth_rate=0.0)

    assert pools.beta_carotene_flux == pytest.approx(supply / 2, rel=0.02)


def test_growth_dilutes_the_intermediate_pools():
    fast = steady_state_pools(ggpp_supply=0.004, kinetics=BALANCED, growth_rate=0.35)
    slow = steady_state_pools(ggpp_supply=0.004, kinetics=BALANCED, growth_rate=0.05)

    assert fast.beta_carotene < slow.beta_carotene


def test_a_desaturase_bottleneck_accumulates_phytoene():
    """The signature that says crtI is limiting -- visible on the HPLC trace."""
    slow_crti = CarotenoidKinetics(
        vmax_psy=0.05, km_psy=0.01, vmax_crti=0.002, km_crti=0.01, vmax_lcy=0.05, km_lcy=0.01
    )
    pools = steady_state_pools(ggpp_supply=0.004, kinetics=slow_crti, growth_rate=0.1)

    assert pools.phytoene > pools.lycopene * 5


def test_a_cyclase_bottleneck_accumulates_lycopene_instead():
    slow_lcy = CarotenoidKinetics(
        vmax_psy=0.05, km_psy=0.01, vmax_crti=0.05, km_crti=0.01, vmax_lcy=0.002, km_lcy=0.01
    )
    pools = steady_state_pools(ggpp_supply=0.004, kinetics=slow_lcy, growth_rate=0.1)

    assert pools.lycopene > pools.phytoene * 5


def test_intermediate_ratios_identify_which_enzyme_is_limiting():
    """The anchor claim, stated as a test."""
    cases = {
        "crtI": CarotenoidKinetics(0.05, 0.01, 0.002, 0.01, 0.05, 0.01),
        "crtYB_LCY": CarotenoidKinetics(0.05, 0.01, 0.05, 0.01, 0.002, 0.01),
        "crtYB_PSY": CarotenoidKinetics(0.002, 0.01, 0.05, 0.01, 0.05, 0.01),
    }
    for expected, kinetics in cases.items():
        pools = steady_state_pools(ggpp_supply=0.004, kinetics=kinetics, growth_rate=0.1)
        assert limiting_step(pools) == expected, f"misidentified {expected}"


class TestBifunctionalEnzyme:
    def test_the_two_crtyb_activities_share_one_enzyme_pool(self):
        shared = CarotenoidKinetics(0.05, 0.01, 0.05, 0.01, 0.05, 0.01, crtyb_pool=0.05)
        rates = branch_rates(
            {"ggpp": 0.05, "phytoene": 0.05, "lycopene": 0.05}, shared
        )

        assert rates["psy"] + rates["lcy"] <= 0.05 * (1 + 1e-9)

    def test_without_a_shared_pool_the_activities_do_not_compete(self):
        independent = CarotenoidKinetics(0.05, 0.01, 0.05, 0.01, 0.05, 0.01, crtyb_pool=None)
        rates = branch_rates(
            {"ggpp": 0.05, "phytoene": 0.05, "lycopene": 0.05}, independent
        )

        assert rates["psy"] + rates["lcy"] > 0.05

    def test_saturating_one_activity_starves_the_other(self):
        shared = CarotenoidKinetics(0.05, 0.01, 0.05, 0.01, 0.05, 0.01, crtyb_pool=0.05)
        balanced = branch_rates({"ggpp": 0.02, "phytoene": 0.02, "lycopene": 0.02}, shared)
        ggpp_flood = branch_rates({"ggpp": 5.0, "phytoene": 0.02, "lycopene": 0.02}, shared)

        assert ggpp_flood["lcy"] < balanced["lcy"]


def test_product_formation_saturates_as_precursor_supply_rises():
    low = steady_state_pools(ggpp_supply=0.002, kinetics=BALANCED, growth_rate=0.1)
    high = steady_state_pools(ggpp_supply=0.2, kinetics=BALANCED, growth_rate=0.1)

    assert high.beta_carotene_flux < low.beta_carotene_flux * 100


def test_the_dynamic_simulation_settles_on_the_steady_state():
    t = np.linspace(0, 400, 4000)
    traj = simulate_branch(t, ggpp_supply=0.004, kinetics=BALANCED, growth_rate=0.1)
    expected = steady_state_pools(ggpp_supply=0.004, kinetics=BALANCED, growth_rate=0.1)

    assert traj["beta_carotene"][-1] == pytest.approx(expected.beta_carotene, rel=0.02)


def test_a_time_varying_supply_is_followed():
    t = np.linspace(0, 48, 1000)
    supply = np.where(t < 24, 0.001, 0.008)

    traj = simulate_branch(t, ggpp_supply=supply, kinetics=BALANCED, growth_rate=0.1)

    assert traj["beta_carotene"][-1] > traj["beta_carotene"][len(t) // 2] * 2


def test_negative_kinetic_parameters_are_refused():
    with pytest.raises(ValueError, match="positive"):
        CarotenoidKinetics(-0.05, 0.01, 0.05, 0.01, 0.05, 0.01)


# --------------------------------------------------------------------------- #
# what the calibration set measures
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def states():
    return calibration_states()


def test_the_measured_lycopene_rate_is_washout_and_not_desaturase_flux(states):
    """The unit error that would silently halve every fitted cyclase constant.

    Carotenoids are not secreted, so the only way lycopene leaves a chemostat is inside the
    cells. The reported ``q_lycopene`` is therefore growth times the intracellular pool, and
    the desaturase flux is that *plus* everything that went on to beta-carotene. Reading the
    reported rate as the crtI flux understates the pathway by up to 5x here: in b-car2 the
    cyclase carries 62 % of the carbon, so lycopene washout is only 38 % of what crtI made.
    """
    assert np.allclose(states.q_lycopene, states.mu_per_h * states.lycopene_content)
    assert np.allclose(states.desaturase_flux, states.q_lycopene + states.q_betacarotene)
    assert (states.cyclase_fraction > 0.19).all() and (states.cyclase_fraction < 0.71).all()


def test_a_steady_state_reports_washout_rates_not_internal_fluxes():
    """The same distinction, on the model side rather than the data side.

    ``crti_flux`` and ``lycopene_accumulation`` are different numbers whenever any carbon
    reaches the product, and only the second is comparable to an HPLC measurement. Asserting
    the identity here stops a future refactor from aliasing one onto the other, which would
    make the model appear to fit while comparing two different quantities.
    """
    state = steady_state_pools(ggpp_supply=0.004, kinetics=BALANCED, growth_rate=0.1)

    assert state.lycopene_accumulation == pytest.approx(state.crti_flux - state.lcy_flux)
    assert state.beta_carotene_accumulation == pytest.approx(state.lcy_flux, rel=1e-6)
    assert state.lycopene_accumulation < state.crti_flux


def test_the_vendored_states_match_the_published_product_rates(states):
    """Guards the vendored file against a silent unit or ordering change.

    The values are Table 2 of PMID 40891387 in nmol/gDCW/h, re-extracted from the PMC JATS
    independently of the authors' repository the TSV was built from. Two sources agreeing is
    the only reason this file is trusted enough to fit against; if a refresh of the TSV ever
    disagrees with the printed table, the fit is invalid and this says so first.
    """
    published = {
        "2D01": (94.4, 151.3), "2D025": (84.3, 184.2),
        "3D01": (472.0, 227.1), "3D025": (126.5, 294.8),
        "4D01": (760.7, 185.3), "4D025": (392.9, 453.6),
    }
    for _, row in states.iterrows():
        lycopene, betacarotene = published[row.condition]
        assert row.q_lycopene * 1e6 == pytest.approx(lycopene, rel=1e-3)
        assert row.q_betacarotene * 1e6 == pytest.approx(betacarotene, rel=1e-3)


# --------------------------------------------------------------------------- #
# the refutation
# --------------------------------------------------------------------------- #


def test_a_growth_rate_independent_cyclase_is_refuted_by_every_strain(states):
    """The result that forced the model to change, stated so it cannot be quietly reverted.

    Michaelis-Menten says the rate is non-decreasing in substrate concentration, for any
    vmax and any km. Each strain was run at two dilution rates, so each strain is a paired
    test of exactly that with the enzyme complement as fixed as an experiment can hold it.
    All three strains go the wrong way -- 2.8x to 9.4x more lycopene buys 0.41x to 0.82x the
    beta-carotene -- and no parameter choice can rescue that, which is why the module's
    calibrated cyclase carries a growth-rate term at all.
    """
    for strain, group in states.groupby("strain"):
        ordered = group.sort_values("lycopene_content")
        low, high = ordered.iloc[0], ordered.iloc[-1]
        assert high.lycopene_content > low.lycopene_content * 2
        assert high.q_betacarotene < low.q_betacarotene, (
            f"{strain} is consistent with a fixed-vmax cyclase; the refutation has gone away"
        )


def test_the_desaturase_step_is_not_identified_by_this_dataset(states):
    """Why four of the six kinetic parameters stay unfitted.

    Phytoene was never measured, so the desaturase has a pinned flux and an unobserved
    substrate. Two parameter sets three orders of magnitude apart both reproduce every state
    exactly, differing only in a phytoene pool nobody measured. A fit reporting one of them
    would be reporting the optimiser's starting point, not the data.
    """
    flux = states.desaturase_flux.to_numpy()
    for vmax in (1.5 * flux.max(), 1000.0 * flux.max()):
        for km in (1e-5, 1e-3):
            pool = km * flux / (vmax - flux)  # the pool that delivers exactly this flux
            assert np.allclose(vmax * pool / (km + pool), flux)


def test_the_full_parameter_set_is_refused_rather_than_guessed():
    with pytest.raises(NotImplementedError, match="does not identify the full branch"):
        calibrated_kinetics(growth_rate=0.15)


def test_the_refusal_names_the_steps_it_cannot_fit():
    """A refusal that does not say what is missing sends the caller to invent a number.

    The repo has been bitten by plausible substituted constants before, so the message has
    to carry the two step names and the one measurement -- a phytoene peak, already resolved
    at 285 nm in the same injection -- that would lift the refusal.
    """
    with pytest.raises(NotImplementedError) as raised:
        calibrated_kinetics(growth_rate=0.15)

    message = str(raised.value)
    assert all(step in message for step in UNIDENTIFIABLE_STEPS)
    assert "phytoene" in message and "40891387" in message


# --------------------------------------------------------------------------- #
# the fitted cyclase
# --------------------------------------------------------------------------- #


def test_the_fitted_law_reproduces_all_six_steady_states(states):
    """Two free parameters against six observations, on real measurements.

    17 % worst-case is not a triumph; it is what a two-parameter law achieves on data
    spanning three strains and a 2.5-fold range of dilution rate, and the point of the bound
    is that a refit which drifts past it has changed something. The comparison that matters
    is the held-out one below, not this.
    """
    predicted = np.array([
        predict_beta_carotene_rate(row.lycopene_content, row.mu_per_h)
        for _, row in states.iterrows()
    ])
    relative = np.abs(predicted - states.q_betacarotene) / states.q_betacarotene

    assert relative.max() < 0.20
    assert np.sqrt(np.mean((np.log(predicted) - np.log(states.q_betacarotene)) ** 2)) < 0.12


def test_the_held_out_score_beats_predicting_the_training_mean():
    """A fit with no held-out point is not evidence, and this repo says so in CLAIM_BOUNDARY.

    Leave-one-out over all six states, scored against the honest baseline for six numbers:
    the mean of the five the fit did see. Skill is 1 minus the ratio of squared log errors,
    so zero means the law is worth exactly as much as the mean and negative means it is
    worse. Several forms that look reasonable on the full-data fit -- vmax per strain with
    NO growth term, vmax tracking measured CrtYB mRNA -- land negative, which is why this
    bound exists. The growth term is what separates them: the same per-strain capacities WITH
    it land at +0.78.
    """
    assert ELIZONDO2025.skill_vs_training_mean > 0.6
    assert ELIZONDO2025.loo_rmse_log < ELIZONDO2025.baseline_rmse_log


def test_the_calibration_declares_more_observations_than_parameters():
    """Six states, two parameters. The count is carried so a reader need not reconstruct it.

    Recorded as an assertion because the tempting next step is a per-strain capacity, which
    is four parameters against six observations. This docstring said until 2026-09-04 that
    such a form "scores worse than the baseline out of sample". That was the score of a
    candidate which also dropped the GROWTH TERM, so it was reporting
    growth-rate-independence and not gene dosage. With the growth term, a per-strain capacity
    scores 0.2054 held out and BEATS the baseline. The reason to keep two parameters is the
    identifiability margin itself, plus a nested F of 0.5666 at p = 0.6383 --
    never the training residual, and no longer a score that was measuring another factor.
    `test_the_two_per_strain_candidates_are_not_interchangeable` pins the distinction.
    """
    assert ELIZONDO2025.n_states == 6
    assert ELIZONDO2025.n_parameters == 2
    assert ELIZONDO2025.n_states > 2 * ELIZONDO2025.n_parameters


def test_the_fitted_constants_carry_intervals_that_bracket_them():
    """A point estimate dressed as a measurement is the failure mode this repo audits for.

    Both constants are fitted -- no turnover number exists for any of these enzymes -- so
    each one ships with the interval its own residual scatter implies. The km interval spans
    a factor of four, and anything reading the point value alone would not know that.
    """
    for value, (low, high) in (
        (ELIZONDO2025.capacity_mmol_per_gdcw, ELIZONDO2025.capacity_ci95),
        (ELIZONDO2025.km_mmol_per_gdcw, ELIZONDO2025.km_ci95),
    ):
        assert low < value < high
    assert ELIZONDO2025.km_ci95[1] / ELIZONDO2025.km_ci95[0] > 3


def test_the_cyclase_runs_partly_saturated_across_the_measured_states(states):
    """Positive evidence that this step holds control, which is the whole reason to fit it.

    A step running far below saturation carries no control -- it passes on whatever arrives.
    The fitted km sits inside the measured range of lycopene contents, so the cyclase runs
    between 36 % and 93 % saturated across the six states. That is the model-side version of
    Elizondo's MCA result putting the highest flux control on CrtYB.
    """
    saturation = states.lycopene_content / (
        ELIZONDO2025.km_mmol_per_gdcw + states.lycopene_content
    )

    assert 0.3 < saturation.min() < 0.5
    assert 0.9 < saturation.max() < 1.0
    assert states.lycopene_content.min() < ELIZONDO2025.km_mmol_per_gdcw < states.lycopene_content.max()


def test_lycopene_is_the_intermediate_that_accumulates_in_these_strains(states):
    """Our limiting-step diagnostic and Elizondo's control analysis agree on crtYB.

    Elizondo's MCA puts the highest control over beta-carotene on CrtYB, downstream of
    lycopene, in the best producer b-car4. The model's own diagnostic reads the same thing
    off the pellet: b-car4 at the low dilution rate carries four times more lycopene than
    beta-carotene. The agreement is corroboration and not confirmation -- phytoene and GGPP
    were not measured, so nothing here ranks the cyclase against the two upstream steps.
    """
    worst = states.set_index("condition").loc["4D01"]

    assert worst.lycopene_content > 4 * worst.betacarotene_content
    assert limiting_step(
        steady_state_pools(0.004, CarotenoidKinetics(0.05, 0.01, 0.05, 0.01, 0.002, 0.01), 0.1)
    ) == "crtYB_LCY"


def test_the_content_law_saturates_and_passes_through_the_origin():
    zero = beta_carotene_content(0.0)
    huge = beta_carotene_content(1.0)

    assert zero == 0.0
    assert huge == pytest.approx(ELIZONDO2025.capacity_mmol_per_gdcw, rel=1e-3)
    assert beta_carotene_content(ELIZONDO2025.km_mmol_per_gdcw) == pytest.approx(
        ELIZONDO2025.capacity_mmol_per_gdcw / 2
    )


def test_a_negative_lycopene_content_is_refused():
    with pytest.raises(ValueError, match="negative"):
        beta_carotene_content(-1e-6)


def test_the_kinetic_reading_of_the_fit_scales_vmax_with_growth_rate():
    """The uncomfortable half of the result, pinned so it stays visible.

    A growth-rate-independent relation between the two contents is arithmetically the same
    statement as a cyclase vmax proportional to growth rate. That is not what a fixed
    complement of Michaelis-Menten enzyme does, and Elizondo's own RT-qPCR has CrtYB mRNA
    *falling* by a third from the low dilution rate to the high one. The fit is good and its
    mechanism is unestablished; this test asserts the shape rather than endorsing the
    mechanism, so that a later mechanistic form has something concrete to disagree with.
    """
    low, high = ELIZONDO2025.growth_rate_range
    vmax_low, km_low = calibrated_cyclase(low)
    vmax_high, km_high = calibrated_cyclase(high)

    assert km_low == km_high
    assert vmax_high / vmax_low == pytest.approx(high / low, rel=1e-9)


def test_the_calibrated_cyclase_plugs_into_the_rate_law_unchanged(states):
    """The fitted pair has to be usable by ``branch_rates``, not just by the helper.

    Two entry points to the same law drift apart silently. Running the measured lycopene
    content through the general Michaelis-Menten machinery must give what
    ``predict_beta_carotene_rate`` gives, or the twin's dynamic mode and its steady-state
    mode would disagree about the same strain.
    """
    for _, row in states.iterrows():
        vmax_lcy, km_lcy = calibrated_cyclase(row.mu_per_h)
        kinetics = CarotenoidKinetics(0.05, 0.01, 0.05, 0.01, vmax_lcy, km_lcy)
        through_rate_law = branch_rates({"lycopene": row.lycopene_content}, kinetics)["lcy"]

        assert through_rate_law == pytest.approx(
            predict_beta_carotene_rate(row.lycopene_content, row.mu_per_h), rel=1e-12
        )


def test_extrapolating_beyond_the_two_measured_dilution_rates_is_refused():
    """Two points cannot establish a functional form, only interpolate between themselves.

    ``vmax ~ mu`` and any other monotone curve through the same pair are indistinguishable
    on this dataset, so the proportionality has no support outside [0.101, 0.254] /h. A
    batch culture at mu_max ~ 0.4 /h is exactly where a caller would reach for this and
    exactly where it is unsupported, so the refusal has to fire rather than the number
    quietly continuing to scale.
    """
    low, high = ELIZONDO2025.growth_rate_range
    for outside in (low * 0.5, high * 1.5, 0.4):
        with pytest.raises(ValueError, match="outside the calibrated range"):
            calibrated_cyclase(outside)
        with pytest.raises(ValueError, match="outside the calibrated range"):
            predict_beta_carotene_rate(1e-3, outside)

    assert calibrated_cyclase(low)[0] > 0
    assert calibrated_cyclase(high)[0] > 0


def test_a_non_positive_growth_rate_is_refused_before_the_range_check():
    """At mu = 0 the measured "rate" is zero over zero and the content is undefined.

    A chemostat at zero dilution has no washout, so ``q / mu`` is not a content and the
    whole derivation the fit rests on collapses. Refusing here rather than returning zero
    keeps that from reading as a real prediction of no product.
    """
    with pytest.raises(ValueError, match="must be positive"):
        calibrated_cyclase(0.0)
    with pytest.raises(ValueError, match="must be positive"):
        calibrated_cyclase(-0.1)


# --------------------------------------------------------------------------- #
# the two per-strain candidates, which are not interchangeable
# --------------------------------------------------------------------------- #

def _fit_script():
    """`scripts/fit_carotenoid_kinetics.py` as a module. Not importable by name -- scripts/
    is not a package -- so it is loaded by path, the way `conftest` does elsewhere."""
    import importlib.util
    import pathlib as _pathlib

    path = _pathlib.Path(__file__).resolve().parents[1] / "scripts" / "fit_carotenoid_kinetics.py"
    spec = importlib.util.spec_from_file_location("_fit_carotenoid_kinetics", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def fit_frame():
    """The six states with the CrtYB transcript attached, as the script builds them."""
    import pandas as pd

    from ystwin import paths

    frame = calibration_states()
    mrna = pd.read_csv(
        paths.data_dir() / "carotenoid" / "elizondo2025_relative_mrna.tsv", sep="\t")
    crtyb = mrna[mrna.gene == "CrtYB"].set_index("condition").rel_expression
    frame["crtyb_mrna"] = frame.condition.map(crtyb)
    return frame


def test_conditional_candidate_predictions_do_not_read_heldout_beta_outcomes(fit_frame):
    script = _fit_script()
    changed = fit_frame.copy()
    changed["q_betacarotene"] *= 10.0
    changed["desaturase_flux"] = changed.q_lycopene + changed.q_betacarotene
    for name, (_, start, _) in script.CANDIDATES.items():
        before = script.predict(name, np.asarray(start), fit_frame)
        after = script.predict(name, np.asarray(start), changed)
        np.testing.assert_array_equal(before, after, err_msg=name)


def test_no_saturation_candidate_is_first_order_in_lycopene_not_outcome_capped(fit_frame):
    script = _fit_script()
    inputs = fit_frame.drop(columns=["q_betacarotene", "desaturase_flux"]).copy()
    theta = np.array([np.log(2.5)])
    name = "vmax ~ mu, no saturation"
    expected = 2.5 * inputs.mu_per_h.to_numpy() * inputs.lycopene_content.to_numpy()
    np.testing.assert_allclose(script.predict(name, theta, inputs), expected)
    inputs["lycopene_content"] *= 3.0
    np.testing.assert_allclose(script.predict(name, theta, inputs), 3.0 * expected)
    inputs["lycopene_content"] = 0.0
    np.testing.assert_array_equal(script.predict(name, theta, inputs), np.zeros(len(inputs)))


def test_zero_km_limit_does_not_select_the_unsaturated_candidate(fit_frame):
    script = _fit_script()
    inputs = fit_frame.drop(columns=["q_betacarotene", "desaturase_flux"])
    expected = 2.5 * inputs.mu_per_h.to_numpy()
    with np.errstate(under="ignore"):
        predicted = script.predict("vmax ~ mu", np.array([np.log(2.5), -1000.0]), inputs)
    np.testing.assert_allclose(predicted, expected)


def test_the_two_per_strain_candidates_are_not_interchangeable(fit_frame):
    """The one assertion that stops -0.58 being quoted as a fact about gene dosage again.

    Two candidates carry per-strain capacities. One also drops the growth term and loses to
    the training-mean baseline; the other keeps it and beats the baseline comfortably. Until
    2026-09-04 only the first existed, so its score was the only per-strain number available
    and four documents quoted it as evidence about dosage. This test fails the day the two
    are confused again -- and it fails LOUDLY if the growth-carrying row is deleted, because
    the lookup raises.
    """
    script = _fit_script()
    observed = np.log(fit_frame.q_betacarotene.to_numpy())
    baseline = script.leave_one_out("training mean content", fit_frame)
    baseline_rmse = np.sqrt(np.mean((np.log(baseline) - observed) ** 2))

    def rmse(name):
        predicted = script.leave_one_out(name, fit_frame)
        return np.sqrt(np.mean((np.log(predicted) - observed) ** 2))

    with_growth = rmse("capacity per strain x mu")
    without_growth = rmse("vmax per strain, growth-independent")

    assert with_growth < baseline_rmse, (
        "a per-strain capacity WITH the growth term beats the baseline; a claim that it "
        "scores worse out of sample is about the other candidate")
    assert without_growth > baseline_rmse
    assert without_growth > 2 * with_growth


def test_pooling_the_capacity_rests_on_the_f_test_and_not_on_the_loo_score(fit_frame):
    """Why the shipped fit keeps two parameters, stated as the test that could refute it.

    The per-strain form is nested in the pooled one, so the comparison is an F test. It is a
    weak test -- two denominator degrees of freedom -- and the fitted capacities are close
    enough together that it does not fire. Both halves are asserted, because a spread this
    small is the reason to pool and the F alone would not say so.
    """
    script = _fit_script()
    nested = script.nested_f_test("vmax ~ mu", "capacity per strain x mu", fit_frame)
    strains = script.per_strain_capacities(fit_frame)
    spread = (strains.capacity_mmol_per_gdcw.max()
              / strains.capacity_mmol_per_gdcw.min())

    assert nested["p"] > 0.05, "the extra capacities would then be worth their parameters"
    assert spread < 1.5, f"capacities spread {spread:.4g}x across a 3x crtYB copy span"
    assert ELIZONDO2025.n_parameters == 2


def test_the_stored_uncertainty_is_the_fit_that_produced_it(fit_frame):
    """One re-derivation pinning THREE numbers, so none can go stale on its own.

    `capacity_ci95`, `km_ci95` and `log_covariance` all came out of a single
    `uncertainty('vmax ~ mu', ...)` call, and until 2026-09-04 only the first two were kept:
    the covariance was computed, its correlation was printed to a terminal, and the dataclass
    stored the diagonal. Retyping +0.852 as a fourth literal beside three others would
    reproduce the failure `pathway/calibrations.py` already exists to prevent, so this refits
    instead.
    """
    script = _fit_script()
    refit = script.uncertainty("vmax ~ mu", fit_frame, draws=1)
    (capacity, cap_lo, cap_hi), (km, km_lo, km_hi) = refit["asymptotic"]

    assert ELIZONDO2025.capacity_mmol_per_gdcw == pytest.approx(capacity, rel=1e-3)
    assert ELIZONDO2025.km_mmol_per_gdcw == pytest.approx(km, rel=1e-3)
    assert ELIZONDO2025.capacity_ci95 == pytest.approx((cap_lo, cap_hi), rel=1e-3)
    assert ELIZONDO2025.km_ci95 == pytest.approx((km_lo, km_hi), rel=1e-3)
    assert ELIZONDO2025.log_parameter_correlation == pytest.approx(
        refit["log_correlation"], rel=1e-3)


def test_the_two_parameters_are_correlated_enough_that_it_matters():
    """The reason `log_covariance` is carried at all, stated as a threshold."""
    assert ELIZONDO2025.log_parameter_correlation > 0.5

    (vaa, vab), (vba, vbb) = ELIZONDO2025.log_covariance
    assert vab == vba
    assert vaa > 0 and vbb > 0


def test_independent_propagation_overstates_the_error_on_the_shipped_solver():
    """The size of the mistake `joint_log_sd` exists to prevent, measured not asserted.

    The gradient is taken by finite differences through the SHIPPED
    `solve_branch_from_flux`, so this is the propagation a real consumer would do. The two
    log parameters enter with opposite signs, so their errors partly cancel and adding the
    variances overstates the answer at every growth rate the calibration covers.
    """
    import math as _math
    from dataclasses import replace

    states = calibration_states()
    flux = float(states.desaturase_flux.median())

    def gradient(growth_rate):
        out = []
        for field in ("capacity_mmol_per_gdcw", "km_mmol_per_gdcw"):
            base, step = getattr(ELIZONDO2025, field), 1e-4
            high = solve_branch_from_flux(
                flux, growth_rate,
                replace(ELIZONDO2025, **{field: base * _math.exp(step)}))
            low = solve_branch_from_flux(
                flux, growth_rate,
                replace(ELIZONDO2025, **{field: base * _math.exp(-step)}))
            out.append((_math.log(high.beta_carotene_content)
                        - _math.log(low.beta_carotene_content)) / (2 * step))
        return tuple(out)

    for growth_rate in ELIZONDO2025.growth_rate_range:
        g = gradient(growth_rate)

        assert g[0] > 0 > g[1], f"the signs are what makes them cancel: {g}"
        joint = ELIZONDO2025.joint_log_sd(g)
        independent = ELIZONDO2025.independent_log_sd(g)

        assert joint < independent
        assert independent / joint > 1.4, (
            f"at mu = {growth_rate:g} independent propagation is only "
            f"{independent / joint:.3g}x the joint one")


def test_shared_nested_branch_fit_and_analytic_branch_use_the_same_mass_balance(states):
    from dataclasses import replace
    from ystwin.pathway.flux import fit_saturating_branch
    from ystwin.pathway.solve import solve_pathway
    from ystwin.pathway.spec import load_pathway

    spec = load_pathway("beta_carotene")
    train = states[states.strain != "b-car2"]
    kinetics = fit_saturating_branch(spec, train)
    fitted = kinetics["lycopene"]
    calibration = replace(ELIZONDO2025, capacity_mmol_per_gdcw=fitted.vmax_per_growth,
                          km_mmol_per_gdcw=fitted.km, growth_rate_range=fitted.growth_rate_range)
    assert fitted.growth_rate_range == (train.mu_per_h.min(), train.mu_per_h.max())
    for row in states.itertuples():
        flux = row.q_betacarotene + row.q_lycopene
        analytic = solve_branch_from_flux(flux, row.mu_per_h, calibration)
        generic = solve_pathway(spec, flux, row.mu_per_h, kinetics)
        assert analytic.beta_carotene_content == pytest.approx(generic.terminal.content_mmol_per_gdcw)
        assert analytic.lycopene_content == pytest.approx(generic.node("lycopene").content_mmol_per_gdcw)


@pytest.mark.parametrize("flux,growth", [(np.nan, 0.15), (np.inf, 0.15), (1e-4, np.nan)])
def test_nonfinite_branch_inputs_are_refused(flux, growth):
    with pytest.raises(ValueError, match="finite"):
        solve_branch_from_flux(flux, growth)

