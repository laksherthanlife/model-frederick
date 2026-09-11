"""The oxidative block, checked by running it -- including where it fails.

Every number asserted here came out of `ystwin.mech.oxidative` on this checkout. Three of
these tests pin a FAILURE: the bolus/perfusion gap overshoots, the leave-one-dose-out win
sits below the assay floor, and the decay rate the plates prefer collapses to zero if the
growth-halving rung is kept. They are pinned deliberately -- a later change that quietly
turns one of them into a pass has to come through here and say so.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import solve_ivp

from ystwin.generator.panel_experiment import OBSERVED_ACTIVITY_CV
from ystwin.generator.stress_panel import STRESSORS
from ystwin.mech import oxidative as ox
from ystwin.mech.integrate import (
    REDUCTION_RATIO_CEILING,
    Verdict,
    freeze_bias,
    qss_bias,
)
from ystwin.mech.params import RefusedValue, SweptValue, Tag
from ystwin.mech.state import (
    FEDBATCH_5D,
    PLATE_GROWTH_BAND_PER_H,
    PLATE_READ_4H,
    Encoding,
)

MU = 0.3227
"""Median chord growth rate on the three committed plates, from PLATE_GROWTH_BAND_PER_H's
own derivation. Passed explicitly everywhere; nothing here asserts a growth rate."""

class _Fit:
    """A stand-in carrying only what an Observable's summarise() reads: `parameters`."""

    def __init__(self, parameters):
        self.parameters = dict(parameters)


EXPORTS = ox.committed_yap1_blocks()
BAND = ox.CELLS_PER_ML_PER_OD600.bounds


@pytest.fixture(scope="module")
def bionumbers():
    from ystwin import paths

    path = paths.data_dir() / "kaggle" / "BioNumbers_Nov2024.csv"
    if not path.exists():
        pytest.skip("BioNumbers export absent")
    frame = pd.read_csv(path)
    return frame[frame.Organism.astype(str).str.contains("cerevisiae", case=False, na=False)]


def _bnid(yeast, bnid: int) -> str:
    row = yeast[yeast.BNID == bnid]
    assert not row.empty, f"BNID {bnid} missing from the export"
    return str(row.iloc[0]["Value/Range"]).strip("'~ ")


# --------------------------------------------------------------------------------------
# The constants, and the gate that counts them
# --------------------------------------------------------------------------------------

class TestTheConstants:
    def test_the_gate_is_two_free_scalars_against_three_targets(self):
        gate = ox.OXIDATIVE_GATE
        assert gate.passes
        assert set(gate.free) == {"cells_per_ml_per_od600", "K_ex"}
        assert set(gate.targets) == {
            "reporter_activity", "bolus_perfusion_lethal_ratio",
            "cross_protection_generations"}

    def test_every_registered_parameter_carries_a_tag_and_a_source(self):
        for registry in (ox.OXIDATIVE_PARAMS, ox.NOT_BUILT):
            for param in registry.params:
                assert param.tag in Tag.ALL
                assert param.source.strip()
                assert param.units.strip()

    def test_k_ref_is_selvaggios_number_and_the_target_doc_tabulates_it_60x_wrong(self):
        """1.33e-3 /s is 0.0798 per MINUTE. ARCHITECTURE_TARGET.md L1 tabulates it as /h."""
        per_second = 1.33e-3
        assert float(ox.K_REF_PER_H) == pytest.approx(per_second * 3600.0)
        assert per_second * 60.0 == pytest.approx(0.0798, abs=5e-5)
        assert float(ox.K_REF_PER_H) == pytest.approx(4.788, abs=1e-3)
        assert float(ox.K_REF_PER_H) / 60.0 == pytest.approx(0.0798, abs=5e-5)

    def test_the_tabulation_error_did_not_propagate_into_that_documents_time_constants(self):
        """36 min at OD 0.05 and 2.5 min at OD 0.80 are per-HOUR arithmetic, done right."""
        implied = [float(ox.X_REF_CELLS_PER_ML) / (float(ox.K_REF_PER_H) * od * (m / 60.0))
                   for od, m in ((0.05, 36.0), (0.80, 2.5))]
        assert implied == pytest.approx([4.87e7, 4.39e7], rel=0.01)
        assert min(implied) > BAND[1]

    def test_the_documents_taus_are_faster_than_this_modules_swept_band_allows(self):
        ours = [1.0 / ox.dose_decay_rate_per_h(od, BAND[1]) * 60.0 for od in (0.05, 0.80)]
        assert ours == pytest.approx([58.5, 3.7], abs=0.1)
        theirs = [36.0, 2.5]
        assert [o / t for o, t in zip(ours, theirs)] == pytest.approx([1.62, 1.48], abs=0.02)
        assert ours[0] / 36.0 == pytest.approx(1.62, abs=0.02)

    def test_k_ref_carries_selvaggios_own_uncertainty(self):
        low, high = ox.K_REF_PER_H.ci95
        assert (low, high) == pytest.approx((1.25e-3 * 3600.0, 1.41e-3 * 3600.0))

    def test_the_cells_per_od_band_is_two_rows_of_the_vendored_bionumbers_export(
            self, bionumbers):
        assert float(_bnid(bionumbers, 100986)) == pytest.approx(BAND[1])
        assert float(_bnid(bionumbers, 106301)) == pytest.approx(BAND[0])

    def test_guans_growth_band_is_bnid_108255(self, bionumbers):
        text = _bnid(bionumbers, 108255)
        assert "90" in text and "140" in text
        assert ox.GUAN_GROWTH_BAND_PER_H == pytest.approx(
            (math.log(2) / (140 / 60), math.log(2) / (90 / 60)))

    def test_k_ex_is_the_ec50_the_repository_already_carries_for_this_agent(self):
        assert float(ox.K_EX_MM) == pytest.approx(STRESSORS["H2O2"].ec50)
        assert ox.K_EX_MM.tag == Tag.BOUNDED
        assert ox.K_EX_MM.bounds == (0.02, 0.6)
        assert ox.K_EX_MM.missing.strip()

    def test_the_growth_halving_doses_are_quoted_verbatim_from_the_stressor_record(self):
        source = STRESSORS["H2O2"].source
        assert "growth halves at 0.79 and 1.24 mM" in source
        assert ox.GROWTH_HALVING_MM["NativeYap1"] == 0.79
        assert ox.GROWTH_HALVING_MM["AlteredYap1"] == 1.24
        assert float(ox.PLATE_LETHAL_BOLUS_MM) == STRESSORS["H2O2"].lethal_dose

    def test_the_two_protein_turnovers_are_christianos_half_lives(self):
        assert float(ox.PROTEIN_TURNOVER_TSA1_PER_H) == pytest.approx(math.log(2) / 10.9)
        assert float(ox.PROTEIN_TURNOVER_CTT1_PER_H) == pytest.approx(math.log(2) / 108.8)
        assert float(ox.PROTEIN_TURNOVER_CTT1_PER_H) < float(ox.PROTEIN_TURNOVER_TSA1_PER_H)

    def test_the_mrna_decay_is_goulevs_forty_minutes(self):
        assert float(ox.MRNA_DECAY_PER_H) == pytest.approx(math.log(2) / (40.0 / 60.0))

    def test_the_swept_conversion_refuses_to_be_a_float(self):
        with pytest.raises(SweptValue):
            float(ox.CELLS_PER_ML_PER_OD600)
        assert ox.CELLS_PER_ML_PER_OD600.at(1e7) == 1e7
        with pytest.raises(ValueError):
            ox.CELLS_PER_ML_PER_OD600.at(1e9)

    def test_every_refusal_raises_and_names_what_would_close_it(self):
        assert len(ox.NOT_BUILT) >= 15
        for param in ox.NOT_BUILT.params:
            assert param.tag == Tag.REFUSED
            assert param.reason.strip() and param.missing.strip()
            with pytest.raises(RefusedValue):
                float(param)

    def test_the_whole_selvaggio_peroxiredoxin_cycle_is_refused_by_name(self):
        for name in ("kOx", "kCond", "kSulf", "kSrx", "kRed", "PrxT", "TrxT", "VMaxApp",
                     "KM_trxr"):
            assert f"prx_cycle_{name}" in ox.NOT_BUILT

    def test_the_permeability_refusal_carries_all_three_disagreeing_values(self):
        source = ox.NOT_BUILT["k_inf_membrane_permeability"].source
        for value in ("1 /min", "0.331", "9.5"):
            assert value in source

    def test_no_free_scalar_is_left_uncited_except_the_ones_that_cite_a_plate(self):
        plate_sourced = {"plate_inoculum_od600", "bolus_lethal_mM",
                         "growth_halving_nativeyap1_mM", "growth_halving_alteredyap1_mM"}
        assert {p.name for p in ox.OXIDATIVE_PARAMS.uncited()} <= plate_sourced
        assert not ox.NOT_BUILT.uncited()

    def test_the_ratio_floor_is_the_measured_cv_propagated_and_is_not_registered(self):
        assert ox.SPREAD_RATIO_FLOOR.noise_floor == pytest.approx(
            OBSERVED_ACTIVITY_CV * math.sqrt(2.0))
        assert ox.SPREAD_RATIO_FLOOR.assay == ox.REPORTER_TARGET.assay
        assert ox.SPREAD_RATIO_FLOOR.name not in {
            t.name for t in ox.OXIDATIVE_PARAMS.targets}

    def test_the_two_floors_this_block_names_are_measured_inside_their_own_experiments(self):
        assert ox.LETHAL_RATIO_FLOOR.noise_floor == pytest.approx((1.24 - 0.79) / 0.79)
        assert ox.CROSS_PROTECTION_FLOOR.noise_floor == pytest.approx(0.5 / 4.5)
        assert ox.REPORTER_TARGET.noise_floor == OBSERVED_ACTIVITY_CV


# --------------------------------------------------------------------------------------
# The two states criterion (d) will not let this block integrate
# --------------------------------------------------------------------------------------

class TestTheAlgebraThatReplacedTwoStates:
    def test_the_sensor_input_function_is_a_michaelis_form_on_the_external_pool(self):
        k = float(ox.K_EX_MM)
        assert ox.yap1_nuclear_fraction(0.0) == 0.0
        assert float(ox.yap1_nuclear_fraction(k)) == pytest.approx(0.5)
        assert float(ox.yap1_nuclear_fraction(1e6)) > 0.999
        grid = np.linspace(0.0, 2.0, 50)
        assert np.all(np.diff(ox.yap1_nuclear_fraction(grid)) > 0)

    def test_the_membrane_reduces_at_the_slowest_of_the_three_disputed_speeds(self):
        """Goulev 1/min, Tomalin 0.331/s, Selvaggio 9.5/s -- the slowest still reduces."""
        slowest_tau_h = (1.0 / (1.0 / 60.0)) / 3600.0
        assert slowest_tau_h / PLATE_READ_4H.duration_h < REDUCTION_RATIO_CEILING
        assert slowest_tau_h / PLATE_READ_4H.duration_h == pytest.approx(4.03e-3, abs=1e-5)

    def test_the_cytosolic_gradient_refuses_and_names_the_570x(self):
        with pytest.raises(ox.GradientRefused) as caught:
            ox.cytosolic_gradient()
        assert "570x" in str(caught.value)

    def test_yap1_relocation_is_reducible_and_21_to_23x_faster_than_its_transcript(self):
        """The multiple is against the transcript STATE's tau, so it carries the growth band."""
        ratio = ox.YAP1_RELOCATION_TAU_H / PLATE_READ_4H.duration_h
        assert ratio < REDUCTION_RATIO_CEILING
        assert ratio == pytest.approx(0.00805, abs=1e-4)
        multiples = [1.0 / (float(ox.MRNA_DECAY_PER_H) + mu) / ox.YAP1_RELOCATION_TAU_H
                     for mu in PLATE_GROWTH_BAND_PER_H]
        assert sorted(multiples) == pytest.approx([21.3722546966438, 23.3287477679751], abs=0.02)
        assert (1.0 / float(ox.MRNA_DECAY_PER_H)) / ox.YAP1_RELOCATION_TAU_H == pytest.approx(
            28.85, abs=0.02)


# --------------------------------------------------------------------------------------
# L1: the well is not a reservoir
# --------------------------------------------------------------------------------------

class TestTheDepletingWell:
    def test_the_closed_form_solves_the_ode_it_claims_to(self):
        k0, dose = 3.2832, 0.5
        solution = solve_ivp(
            lambda t, y: [-k0 * math.exp(MU * t) * y[0]], (0.0, 4.14), [dose],
            t_eval=np.linspace(0, 4.14, 41), rtol=1e-10, atol=1e-14, method="Radau")
        closed = ox.extracellular_h2o2(solution.t, dose, k0, MU)
        assert np.allclose(closed, solution.y[0], atol=1e-9)

    def test_at_zero_decay_the_dose_is_the_incumbent_exactly(self):
        t = np.linspace(0, 8.0, 17)
        assert np.all(ox.extracellular_h2o2(t, 0.5, 0.0, MU) == 0.5)
        assert ox.dose_half_life_h(0.0, MU) == math.inf

    def test_the_rate_is_linear_in_density_and_in_the_swept_conversion(self):
        assert ox.dose_decay_rate_per_h(0.32, BAND[1]) == pytest.approx(
            2.0 * ox.dose_decay_rate_per_h(0.16, BAND[1]))
        low, high = ox.dose_decay_band_per_h(0.16)
        assert high / low == pytest.approx(BAND[1] / BAND[0])
        assert (low, high) == pytest.approx((0.87552, 3.2832))

    def test_the_measured_half_life_at_the_plates_inoculum(self):
        low, high = ox.dose_decay_band_per_h(float(ox.PLATE_INOCULUM_OD600))
        assert ox.dose_half_life_h(high, MU) * 60.0 == pytest.approx(12.25, abs=0.05)
        assert ox.dose_half_life_h(low, MU) * 60.0 == pytest.approx(42.30, abs=0.05)

    def test_the_note_on_how_fast_the_dose_goes_is_the_number_the_model_gives(self):
        remaining = [float(ox.extracellular_h2o2(
            1.0, 1.0, ox.dose_decay_rate_per_h(0.16, c), MU)) for c in BAND]
        assert remaining == pytest.approx([0.35583, 0.02076], abs=1e-4)
        end = [float(ox.extracellular_h2o2(
            PLATE_READ_4H.duration_h, 1.0, ox.dose_decay_rate_per_h(0.16, c), MU))
            for c in BAND]
        assert max(end) < 6e-4

    def test_dilution_is_a_fifth_of_the_transcripts_loss_and_all_of_the_effectors(self):
        for mu, share in ((0.245, 0.1907), (0.363, 0.2588)):
            assert mu / (float(ox.MRNA_DECAY_PER_H) + mu) == pytest.approx(share, abs=1e-4)
        for mu, share in ((0.245, 0.9747), (0.363, 0.9828)):
            assert mu / (mu + float(ox.PROTEIN_TURNOVER_CTT1_PER_H)) == pytest.approx(
                share, abs=1e-4)

    def test_a_growing_culture_eats_the_dose_sooner_than_a_static_one(self):
        k0 = 1.0
        assert ox.dose_half_life_h(k0, MU) < math.log(2) / k0

    def test_the_exposure_falls_short_of_the_held_dose_and_recovers_it_at_zero_decay(self):
        window, dose = 4.14, 0.5
        held = ox.time_integrated_dose(dose, 0.0, MU, window)
        assert held == pytest.approx(dose * window)
        for k0 in (0.1, 1.0, 3.28):
            assert ox.time_integrated_dose(dose, k0, MU, window) < held

    def test_a_negative_dose_rate_or_a_dead_culture_is_refused(self):
        with pytest.raises(ValueError):
            ox.extracellular_h2o2(1.0, 0.5, -1.0, MU)
        with pytest.raises(ValueError):
            ox.extracellular_h2o2(1.0, 0.5, 1.0, 0.0)
        with pytest.raises(ValueError):
            ox.dose_decay_rate_per_h(-0.1, BAND[1])


# --------------------------------------------------------------------------------------
# Criterion (d)
# --------------------------------------------------------------------------------------

class TestTheReductionAudit:
    def test_there_are_exactly_three_states_and_they_are_named(self):
        system = ox.oxidative_states(PLATE_READ_4H, od600=0.16)
        assert system.names == (
            "h2o2_extracellular", "oxidative_mrna", "oxidative_effector")
        assert len(system) == 3

    def test_every_state_names_an_assay_or_declares_an_axis(self):
        for state in ox.oxidative_states(PLATE_READ_4H, od600=0.16):
            assert state.constrained_by.strip() or state.sweep_axis
            assert state.units.strip() and state.source.strip()

    def test_all_three_states_are_kept_on_the_plate_read(self):
        audit = ox.reduction_audit(PLATE_READ_4H, od600=0.16)
        assert audit.integrated.names == (
            "h2o2_extracellular", "oxidative_mrna", "oxidative_effector")
        assert not audit.eliminated.names and not audit.frozen.names
        for name in audit.integrated.names:
            assert audit.verdict(name) is Verdict.KEEP

    def test_the_fed_batch_eliminates_both_pools_so_the_sets_are_not_nested(self):
        audit = ox.reduction_audit(FEDBATCH_5D, od600=0.16)
        assert audit.verdict("oxidative_mrna") is Verdict.ELIMINATE
        assert audit.verdict("oxidative_effector") is Verdict.ELIMINATE
        assert audit.integrated.names == ("h2o2_extracellular",)

    def test_the_dose_is_integrating_because_it_has_no_quasi_steady_value(self):
        system = ox.oxidative_states(PLATE_READ_4H, od600=0.16)
        assert system["h2o2_extracellular"].encoding is Encoding.INTEGRATING
        assert system["oxidative_mrna"].encoding is Encoding.LEVEL
        assert system["oxidative_effector"].encoding is Encoding.LEVEL

    def test_freezing_the_dose_biases_the_integral_far_above_the_assay_floor(self):
        """The incumbent IS the freeze, and this is what it costs."""
        row = ox.reduction_audit(PLATE_READ_4H, od600=0.16).row("h2o2_extracellular")
        worst = freeze_bias(row.tau_low_h, PLATE_READ_4H.duration_h)
        assert worst > 0.9
        assert worst > 6.0 * REDUCTION_RATIO_CEILING

    def test_the_transcript_pool_is_a_marginal_keep_and_the_margin_is_reported(self):
        row = ox.reduction_audit(PLATE_READ_4H, od600=0.16).row("oxidative_mrna")
        assert REDUCTION_RATIO_CEILING < row.ratio_low < 0.2
        assert row.bias_low == pytest.approx(qss_bias(row.tau_low_h, PLATE_READ_4H.duration_h))

    def test_the_global_mrna_row_would_have_reduced_it_and_the_pathway_number_does_not(self):
        """Chan's 3.6 min median puts tau/T at 0.0126; TSA1's own 40 min puts it at 0.17."""
        fast_tau = (3.6 / 60.0) / math.log(2.0)
        assert fast_tau / PLATE_READ_4H.duration_h < REDUCTION_RATIO_CEILING
        pathway_tau = 1.0 / (float(ox.MRNA_DECAY_PER_H) + MU)
        assert pathway_tau / PLATE_READ_4H.duration_h > REDUCTION_RATIO_CEILING

    def test_the_audit_table_prints_the_ratio_the_bias_and_the_verdict(self):
        table = ox.reduction_audit(PLATE_READ_4H, od600=0.16).audit_table()
        for token in ("h2o2_extracellular", "oxidative_mrna", "oxidative_effector",
                      "keep", "0.146"):
            assert token in table


# --------------------------------------------------------------------------------------
# The integrator and the closed form have to agree
# --------------------------------------------------------------------------------------

def _closed_form_error(*, dose, od, conv, mu, window_h, n_t=61):
    """Worst relative error of the closed form against Radau at rtol 1e-11, one corner."""
    k0 = ox.dose_decay_rate_per_h(od, conv)
    t = np.linspace(0.0, window_h, n_t)[1:]
    reference = solve_ivp(
        lambda tt, y: ox.oxidative_rhs(tt, y, k0_per_h=k0, growth_rate_per_h=mu),
        (0.0, window_h), [dose, 0.0, 0.0], t_eval=t, method="Radau",
        rtol=1e-11, atol=1e-14)
    closed = ox.effector_trace(t, dose_mM=dose, k0_per_h=k0, growth_rate_per_h=mu)
    return float(np.max(np.abs(closed - reference.y[2])) / np.max(reference.y[2]))


def _operating_box(window):
    """Every corner of the swept band and the dose ladder, at the WINDOW's OWN measured mu.

    Sweeping the plate's growth band into a fed-batch would be scoring the quadrature at a
    growth rate that window never runs at.
    """
    mus = sorted({window.growth_rate_low_per_h, window.growth_rate_high_per_h})
    return [{"conv": conv, "od": od, "mu": mu, "dose": dose}
            for conv in BAND for od in (0.05, 0.16, 0.80)
            for mu in mus for dose in (0.1, 0.5)]


class TestTheQuadratureGridTracksTheDriver:
    """The defect this grid was rewritten to fix, pinned so it cannot come back.

    A grid spread uniformly over the REQUESTED window is wrong at exactly this block's
    operating point -- a dose consumed in minutes read out over days -- and on the 120 h
    fed-batch it was 1.10e-1 relative: only 1.3x BELOW the 0.146 reporter floor rather than
    four orders below it, which is no margin at all.
    """

    def test_the_horizon_inverts_the_closed_form_exactly(self):
        for k0, mu in ((3.2832, MU), (16.416, 0.101), (0.876, 0.363)):
            horizon = ox.driver_horizon_h(k0, mu)
            remaining = ox.extracellular_h2o2(horizon, 1.0, k0, mu)
            assert remaining == pytest.approx(ox.DRIVER_FLOOR, rel=1e-9)

    def test_the_neglected_driver_is_orders_below_the_error_it_replaces(self):
        """u <= E/K_ex, so the widest ladder rung over the narrowest K_ex bounds the neglect."""
        widest_dose = 4.0
        neglected = ox.DRIVER_FLOOR * (widest_dose + ox.K_EX_MM.bounds[0]) / ox.K_EX_MM.bounds[0]
        assert neglected == pytest.approx(2.0e-10, rel=0.05)
        assert neglected < 1e-5 / 1000.0

    def test_a_dose_that_never_decays_has_no_horizon(self):
        assert ox.driver_horizon_h(0.0, MU) == math.inf
        with pytest.raises(ValueError):
            ox.driver_horizon_h(-1.0, MU)
        with pytest.raises(ValueError):
            ox.driver_horizon_h(1.0, MU, floor=0.0)

    @pytest.mark.parametrize("window", [PLATE_READ_4H, FEDBATCH_5D])
    def test_the_closed_form_stays_far_below_the_floor_over_the_whole_box(self, window):
        errors = [_closed_form_error(window_h=window.duration_h, **corner)
                  for corner in _operating_box(window)]
        assert max(errors) < 1e-4
        assert OBSERVED_ACTIVITY_CV / max(errors) > 1000.0

    def test_the_fed_batch_corner_that_used_to_fail_is_the_one_that_now_passes(self):
        """OD600 0.80 at the top of the swept band: 0.042 h half-life, 120 h window."""
        error = _closed_form_error(dose=0.1, od=0.80, conv=BAND[1], mu=0.101,
                                   window_h=120.0)
        assert error < 1e-4
        stale = 0.1095
        assert stale > OBSERVED_ACTIVITY_CV / 1.5
        assert error < stale / 1000.0

    def test_the_tail_is_the_free_two_exponential_relaxation(self):
        """Past the horizon the driver is off, so P must fall with no further input."""
        k0 = ox.dose_decay_rate_per_h(0.80, BAND[1])
        horizon = ox.driver_horizon_h(k0, MU)
        t = horizon + np.linspace(0.5, 60.0, 40)
        trace = ox.effector_trace(t, dose_mM=0.5, k0_per_h=k0, growth_rate_per_h=MU)
        assert np.all(np.diff(trace) < 0.0)
        assert trace[-1] > 0.0

    def test_the_plate_window_headline_is_unmoved_by_the_grid_change(self):
        """The read is shorter than the horizon at the plates' own inoculum, so it must be."""
        k0 = ox.dose_decay_rate_per_h(float(ox.PLATE_INOCULUM_OD600), BAND[0])
        assert ox.driver_horizon_h(k0, MU) > PLATE_READ_4H.duration_h
        t = np.linspace(0.0, PLATE_READ_4H.duration_h, 33)
        fine = ox.effector_trace(t, dose_mM=0.5, k0_per_h=k0, growth_rate_per_h=MU,
                                 n_grid=24001)
        default = ox.effector_trace(t, dose_mM=0.5, k0_per_h=k0, growth_rate_per_h=MU)
        assert np.max(np.abs(default - fine)) < 1e-5 * np.max(fine)


class TestTheChain:
    def test_the_convolution_and_the_stiff_solver_agree(self):
        trajectory = ox.simulate_oxidative(
            PLATE_READ_4H, dose_mM=0.5, od600=0.16, cells_per_ml_per_od600=BAND[1],
            growth_rate_per_h=MU, n_points=201)
        k0 = ox.dose_decay_rate_per_h(0.16, BAND[1])
        closed = ox.effector_trace(trajectory.t, dose_mM=0.5, k0_per_h=k0,
                                   growth_rate_per_h=MU, n_grid=4001)
        solved = trajectory.of("oxidative_effector")
        assert np.max(np.abs(closed - solved)) < 0.01 * np.max(solved)

    def test_the_solver_used_the_pinned_implicit_method(self):
        trajectory = ox.simulate_oxidative(
            PLATE_READ_4H, dose_mM=0.5, od600=0.16, cells_per_ml_per_od600=BAND[1],
            growth_rate_per_h=MU)
        assert trajectory.method in ("BDF", "LSODA")

    def test_the_dose_is_gone_and_the_effector_is_not_by_the_end_of_a_plate_read(self):
        trajectory = ox.simulate_oxidative(
            PLATE_READ_4H, dose_mM=0.5, od600=0.16, cells_per_ml_per_od600=BAND[1],
            growth_rate_per_h=MU)
        assert trajectory.of("h2o2_extracellular")[-1] < 1e-6 * 0.5
        assert trajectory.of("oxidative_effector")[-1] > 0.0

    def test_the_effector_peaks_after_the_dose_has_gone(self):
        trajectory = ox.simulate_oxidative(
            PLATE_READ_4H, dose_mM=0.5, od600=0.16, cells_per_ml_per_od600=BAND[1],
            growth_rate_per_h=MU, n_points=801)
        dose = trajectory.of("h2o2_extracellular")
        effector = trajectory.of("oxidative_effector")
        half_dose = trajectory.t[np.argmax(dose < 0.25)]
        assert trajectory.t[int(np.argmax(effector))] > half_dose

    def test_a_pre_treated_culture_starts_the_run_with_a_carrier(self):
        trajectory = ox.simulate_oxidative(
            PLATE_READ_4H, dose_mM=0.5, od600=0.16, cells_per_ml_per_od600=BAND[1],
            growth_rate_per_h=MU, initial_effector=0.02)
        assert trajectory.of("oxidative_effector")[0] == pytest.approx(0.02)

    def test_the_effector_is_zero_without_a_dose(self):
        assert ox.effector_trace(np.linspace(0, 4.14, 9), dose_mM=0.0, k0_per_h=1.0,
                                 growth_rate_per_h=MU).max() == 0.0


# --------------------------------------------------------------------------------------
# Criterion (a)
# --------------------------------------------------------------------------------------

@pytest.mark.integration
class TestCriterionA:
    def test_three_committed_ladders_carry_this_construct(self):
        assert len(EXPORTS) == 3
        assert all(export.endswith(".xlsx") for export in EXPORTS)

    def test_the_ladder_is_truncated_by_a_measured_dose_not_a_chosen_one(self):
        frame = ox.yap1_frame(EXPORTS[0])
        assert frame.dose_mM.max() <= ox.GROWTH_HALVING_MM["NativeYap1"]
        assert sorted(frame.dose_mM.unique()) == [0.0, 0.1, 0.2, 0.5]
        assert set(ox._REQUIRED_COLUMNS) <= set(frame.columns)

    def test_an_unknown_construct_is_refused_rather_than_given_a_default_cutoff(self):
        with pytest.raises(KeyError):
            ox.yap1_frame(EXPORTS[0], "UPRE1")
        assert len(ox.yap1_frame(EXPORTS[0], "UPRE1", max_dose_mM=0.5)) > 0

    def test_the_inoculum_spread_ablation_clears_the_floor_on_every_plate_and_corner(self):
        for export in EXPORTS:
            for conversion in BAND:
                result = ox.inoculum_spread_ablation(
                    export, cells_per_ml_per_od600=conversion, growth_rate_per_h=MU)
                assert result.clears_floor
                assert result.floor is ox.SPREAD_RATIO_FLOOR
                assert 4.5 < result.ratio < 7.5
                assert result.reduced_value == pytest.approx(1.0)
                assert 1.9 < result.full_value < 2.6

    def test_neither_side_of_the_ablation_wins_by_spending_a_parameter(self):
        result = ox.inoculum_spread_ablation(EXPORTS[0], growth_rate_per_h=MU)
        assert result.full_n_free == result.reduced_n_free == 3
        assert not result.full_parameter_free

    def test_the_single_inoculum_ablation_also_clears_at_both_ends(self):
        ratios = [
            ox.inoculum_ablation(export, od600=od, cells_per_ml_per_od600=conversion,
                                 growth_rate_per_h=MU).ratio
            for export in EXPORTS for conversion in BAND for od in (0.05, 0.80)]
        assert min(ratios) > 1.0
        assert 1.7 < min(ratios) < 1.9 and 3.0 < max(ratios) < 3.2

    def test_neither_side_carries_the_panels_viability_term(self):
        """Both models rise with dose; the panel's own response turns down past its peak."""
        from ystwin.generator.stress_panel import module_response

        doses = [0.1, 0.2, 0.5, 1.0, 2.0, 4.0]
        panel = [module_response("H2O2", d)["oxidative"] for d in doses]
        assert min(np.diff(panel)) < 0.0
        for model in (ox.YapReporter(decay="none", growth_rate_per_h=MU),
                      ox.YapReporter(decay="declared", growth_rate_per_h=MU)):
            query = pd.DataFrame({"time_h": [4.14] * len(doses), "dose_mM": doses,
                                  "signal": [0.0] * len(doses),
                                  "od600": [0.16] * len(doses)})
            predicted = model.predict(
                query, {"R0": 0.0, "k_basal": 0.0, "gain": 1.0})
            assert min(np.diff(predicted)) > 0.0

    def test_the_incumbent_cannot_see_the_inoculum_at_all(self):
        model = ox.YapReporter(decay="none", growth_rate_per_h=MU)
        fit = model.refit(ox.yap1_frame(EXPORTS[0]))
        query = pd.DataFrame({"time_h": [4.14, 4.14], "dose_mM": [0.5, 0.5],
                              "signal": [0.0, 0.0], "od600": [0.02, 2.0]})
        predicted = model.predict(query, fit.parameters)
        assert predicted[0] == pytest.approx(predicted[1])

    def test_the_gain_divides_out_only_once_the_two_basal_terms_are_zero(self):
        """Not a claim that the gain cancels: it does not. See TestWhatTheSpreadRatioActuallyCancels."""
        model = ox.YapReporter(decay="declared", cells_per_ml_per_od600=BAND[1],
                               growth_rate_per_h=MU)
        data = ox.yap1_frame(EXPORTS[0])
        fit = model.refit(data)
        observable = ox.inoculum_spread_observable(data)
        doubled = dict(fit.parameters)
        doubled["gain"] = 2.0 * doubled["gain"]
        doubled["R0"] = doubled["k_basal"] = 0.0
        base = dict(fit.parameters)
        base["R0"] = base["k_basal"] = 0.0
        assert observable.summarise(model, ox.Fit(
            "x", doubled, 3, len(data), 0.0, True)) == pytest.approx(
            observable.summarise(model, ox.Fit("x", base, 3, len(data), 0.0, True)))

    def test_the_shape_comparison_wins_out_of_sample_and_stays_below_the_floor(self):
        """PINNED FAILURE. A consistent win, 0.25-0.48x the floor. Not a pass."""
        for export in EXPORTS:
            comparison = ox.shape_comparison(export, growth_rate_per_h=MU)
            assert comparison.full_model_wins
            assert not comparison.clears_floor
            assert 0.2 < comparison.ratio() < 0.5

    def test_the_nesting_is_a_point_of_the_parameter_space_and_not_a_limit(self):
        data = ox.yap1_frame(EXPORTS[0])
        held = ox.YapReporter(decay="none", growth_rate_per_h=MU)
        fitted = ox.YapReporter(decay="fitted", growth_rate_per_h=MU)
        assert set(held.parameter_names) < set(fitted.parameter_names)
        parameters = dict(held.refit(data).parameters)
        parameters["k0_per_h"] = 0.0
        assert np.allclose(fitted.predict(data, parameters),
                           held.predict(data, parameters))

    def test_the_fitted_decay_rate_does_not_depend_on_where_the_search_starts(self):
        """A search that lands somewhere different from each start is a scalar nobody counted."""
        from scipy.optimize import least_squares

        model = ox.YapReporter(decay="fitted", growth_rate_per_h=MU)
        data = ox.yap1_frame(EXPORTS[0])
        landed = [
            float(least_squares(
                lambda x: model._linear_fit(data, max(float(x[0]), 0.0))[1], [start],
                bounds=([0.0], [200.0]), max_nfev=2000).x[0])
            for start in (0.01, 1.0, 100.0)]
        assert landed == pytest.approx([landed[0]] * 3, rel=1e-4)
        spread = max(landed) - min(landed)
        assert spread < 1e-4 * (BAND[1] / BAND[0])

    def test_the_fit_is_deterministic(self):
        data = ox.yap1_frame(EXPORTS[0])
        model = ox.YapReporter(decay="declared", growth_rate_per_h=MU)
        first, second = model.refit(data), model.refit(data)
        assert first.parameters == second.parameters
        assert first.converged


class TestWhatTheSpreadRatioActuallyCancels:
    """The verifier's correction: the plate SCALE cancels in the spread ratio; the gain does not."""

    def test_the_plate_scale_cancels_exactly(self):
        data = ox.yap1_frame(EXPORTS[0])
        model = ox.YapReporter(decay="declared", cells_per_ml_per_od600=BAND[1],
                               growth_rate_per_h=MU)
        fit = model.refit(data)
        observable = ox.inoculum_spread_observable(data)
        base = observable.summarise(model, fit)
        for scale in (0.5, 2.0, 100.0):
            scaled = _Fit({k: v * scale for k, v in fit.parameters.items()})
            assert observable.summarise(model, scaled) == pytest.approx(base, rel=1e-9)

    def test_the_gain_alone_does_not_cancel(self):
        """The two basal terms are common to both wells, so the gain does not divide out."""
        data = ox.yap1_frame(EXPORTS[0])
        model = ox.YapReporter(decay="declared", cells_per_ml_per_od600=BAND[1],
                               growth_rate_per_h=MU)
        fit = model.refit(data)
        observable = ox.inoculum_spread_observable(data)
        moved = []
        for scale in (0.25, 1.0, 4.0):
            parameters = dict(fit.parameters)
            parameters["gain"] *= scale
            moved.append(observable.summarise(model, _Fit(parameters)))
        assert moved[0] < moved[1] < moved[2]
        assert moved == pytest.approx([1.393, 2.453, 5.462], abs=0.01)

    def test_the_incumbent_answers_one_at_any_parameters(self):
        data = ox.yap1_frame(EXPORTS[0])
        model = ox.YapReporter(decay="none", growth_rate_per_h=MU)
        observable = ox.inoculum_spread_observable(data)
        for parameters in ({"R0": 1.0, "k_basal": 0.0, "gain": 5.0},
                           {"R0": 0.0, "k_basal": 3.0, "gain": 0.1}):
            assert observable.summarise(model, _Fit(parameters)) == pytest.approx(1.0)


class TestTheMrnaStateIsTheNarrowestCall:
    def test_the_half_life_at_which_the_transcript_state_would_eliminate(self):
        """29.5-32.2 min. Below it the block is a two-state block; 40 min clears by 1.18-1.29x."""
        window = PLATE_READ_4H
        for mu, expected in ((window.growth_rate_low_per_h, 29.5),
                             (window.growth_rate_high_per_h, 32.2)):
            k_deg = 1.0 / (REDUCTION_RATIO_CEILING * window.duration_h) - mu
            assert math.log(2.0) / k_deg * 60.0 == pytest.approx(expected, abs=0.1)
        margins = [1.0 / (float(ox.MRNA_DECAY_PER_H) + mu) / window.duration_h
                   / REDUCTION_RATIO_CEILING
                   for mu in (window.growth_rate_low_per_h, window.growth_rate_high_per_h)]
        assert min(margins) == pytest.approx(1.18, abs=0.01)
        assert max(margins) == pytest.approx(1.29, abs=0.01)


class TestTheUnfittedRateOutPredictsTheFittedOneAcrossPlates:
    """Leave-one-PLATE-out. The sharpest emergence statement the three plates support."""

    @staticmethod
    def _held_out_rms(model, data, k0):
        _, residual = model._linear_fit(data, k0)
        return float(np.sqrt(np.mean(residual ** 2)))

    def test_branco_beats_a_rate_fitted_on_the_other_two_plates(self):
        frames = {export: ox.yap1_frame(export) for export in EXPORTS}
        for held in EXPORTS:
            train = pd.concat([frames[e] for e in EXPORTS if e != held], ignore_index=True)
            k0 = ox.YapReporter(decay="fitted",
                                growth_rate_per_h=MU).refit(train).parameters["k0_per_h"]
            test = frames[held]
            fitted_elsewhere = self._held_out_rms(
                ox.YapReporter(decay="fitted", growth_rate_per_h=MU), test, k0)
            measured = self._held_out_rms(
                ox.YapReporter(decay="declared", cells_per_ml_per_od600=BAND[0],
                               growth_rate_per_h=MU), test, None)
            incumbent = self._held_out_rms(
                ox.YapReporter(decay="none", growth_rate_per_h=MU), test, None)
            assert measured < fitted_elsewhere < incumbent
            assert incumbent / measured > 1.6

    def test_the_other_end_of_the_swept_band_wins_on_two_plates_of_three(self):
        """Reported because the clean win is at ONE declared point, not across the band."""
        wins = 0
        for export in EXPORTS:
            data = ox.yap1_frame(export)
            high = self._held_out_rms(
                ox.YapReporter(decay="declared", cells_per_ml_per_od600=BAND[1],
                               growth_rate_per_h=MU), data, None)
            incumbent = self._held_out_rms(
                ox.YapReporter(decay="none", growth_rate_per_h=MU), data, None)
            wins += high < incumbent
        assert wins == 2


class TestTheOmittedViabilityTermIsNotCarryingTheResult:
    """The one confound that flattens a dose-response the way a depleting pool does."""

    @staticmethod
    def _viable(**kwargs):
        from ystwin.generator.stress_panel import viability

        class _Viable(ox.YapReporter):
            def _basis(self, data, k0):
                basis = super()._basis(data, k0).copy()
                basis[:, 2] *= np.array(
                    [viability("H2O2", float(d)) for d in data["dose_mM"]])
                return basis

        return _Viable(**kwargs)

    def test_restoring_it_makes_the_incumbent_worse_not_better(self):
        for export in EXPORTS:
            data = ox.yap1_frame(export)
            plain = ox.YapReporter(decay="none", growth_rate_per_h=MU).refit(data)
            viable = self._viable(decay="none", growth_rate_per_h=MU).refit(data)
            assert viable.residual_rms > plain.residual_rms

    def test_restoring_it_barely_moves_the_fitted_decay_rate(self):
        for export in EXPORTS:
            data = ox.yap1_frame(export)
            plain = ox.YapReporter(decay="fitted", growth_rate_per_h=MU).refit(data)
            viable = self._viable(decay="fitted", growth_rate_per_h=MU).refit(data)
            moved = (viable.parameters["k0_per_h"] / plain.parameters["k0_per_h"]) - 1.0
            assert 0.0 < moved < 0.10

    def test_the_omission_is_not_small_on_the_rung_the_headline_is_quoted_at(self):
        """0.850 at 0.5 mM -- a 15% suppression, above the reporter CV. Bounded, not harmless."""
        from ystwin.generator.stress_panel import viability

        assert viability("H2O2", 0.5) == pytest.approx(0.850, abs=0.001)
        assert 1.0 - viability("H2O2", 0.5) > OBSERVED_ACTIVITY_CV


class TestHowManyFittedRatesLandInsideDependsOnMu:
    def test_the_count_moves_across_the_windows_own_growth_band(self):
        """3 of 3 at the band's low end, 2 at the median, 1 at the high end. The middle is quoted."""
        low, high = ox.dose_decay_band_per_h(float(ox.PLATE_INOCULUM_OD600))
        inside = {}
        for mu in (PLATE_READ_4H.growth_rate_low_per_h, MU,
                   PLATE_READ_4H.growth_rate_high_per_h):
            rates = [ox.YapReporter(decay="fitted", growth_rate_per_h=mu)
                     .refit(ox.yap1_frame(export)).parameters["k0_per_h"]
                     for export in EXPORTS]
            inside[mu] = sum(1 for k in rates if low <= k <= high)
        assert inside[PLATE_READ_4H.growth_rate_low_per_h] == 3
        assert inside[MU] == 2
        assert inside[PLATE_READ_4H.growth_rate_high_per_h] == 1

    def test_no_decay_is_rejected_at_every_mu_in_the_band(self):
        """The band-free part of the emergence claim, and the only part that is."""
        for mu in (PLATE_READ_4H.growth_rate_low_per_h, MU,
                   PLATE_READ_4H.growth_rate_high_per_h):
            for export in EXPORTS:
                data = ox.yap1_frame(export)
                fitted = ox.YapReporter(decay="fitted", growth_rate_per_h=mu).refit(data)
                held = ox.YapReporter(decay="none", growth_rate_per_h=mu).refit(data)
                assert fitted.parameters["k0_per_h"] > 0.3
                assert held.residual_rms / fitted.residual_rms > 1.5


class TestTheGapsFailureIsProportionalToTheWindow:
    def test_the_predicted_ratio_grows_with_the_declared_window(self):
        """Reported because the failure's SIZE is a window choice; only its sign is not."""
        ratios = {t: ox.bolus_perfusion_gap(window_h=t, cells_per_ml_per_od600=BAND[1],
                                            growth_rate_per_h=MU).predicted_ratio
                  for t in (0.5, 1.0, 4.14, 24.0)}
        assert ratios[0.5] == pytest.approx(2.09, abs=0.02)
        assert ratios[1.0] == pytest.approx(3.64, abs=0.02)
        assert ratios[4.14] == pytest.approx(14.82, abs=0.02)
        assert ratios[24.0] == pytest.approx(85.94, abs=0.05)
        assert all(v > 1.0 for v in ratios.values())

    def test_a_short_enough_window_would_clear_the_floor_it_fails_on(self):
        short = ox.bolus_perfusion_gap(window_h=0.5, cells_per_ml_per_od600=BAND[1],
                                       growth_rate_per_h=MU)
        plate = ox.bolus_perfusion_gap(cells_per_ml_per_od600=BAND[1], growth_rate_per_h=MU)
        assert short.overshoot - 1.0 < ox.LETHAL_RATIO_FLOOR.noise_floor
        assert plate.overshoot - 1.0 > ox.LETHAL_RATIO_FLOOR.noise_floor


class TestTheDepletingDoseBeatsAnEqualParameterRival:
    def test_a_free_effector_loss_on_a_held_dose_does_not_explain_the_traces(self):
        """The rival is the same 4 scalars: an undepleted dose plus one free protein loss."""
        for export in EXPORTS:
            data = ox.yap1_frame(export)
            depleting = ox.YapReporter(decay="fitted", growth_rate_per_h=MU).refit(data)
            rival = min(
                ox.YapReporter(decay="none", growth_rate_per_h=MU,
                               effector_turnover_per_h=float(delta)).refit(data).residual_rms
                for delta in np.exp(np.linspace(math.log(1e-3), math.log(50.0), 60)))
            assert depleting.residual_rms < rival


class TestTheModelsOwnGuards:
    def test_a_bad_decay_mode_or_a_dead_culture_is_refused_at_construction(self):
        with pytest.raises(ValueError):
            ox.YapReporter(decay="sometimes")
        with pytest.raises(ValueError):
            ox.YapReporter(growth_rate_per_h=0.0)
        with pytest.raises(ValueError):
            ox.YapReporter(decay="declared", cells_per_ml_per_od600=1e9)

    def test_a_frame_without_an_inoculum_column_is_refused(self):
        frame = pd.DataFrame({"time_h": [0.0, 1.0], "dose_mM": [0.5, 0.5],
                              "signal": [1.0, 2.0]})
        with pytest.raises(KeyError):
            ox.YapReporter(growth_rate_per_h=MU).refit(frame)

    def test_the_spread_observable_refuses_an_inverted_inoculum_range(self):
        frame = pd.DataFrame({"time_h": [0.0], "dose_mM": [0.5], "signal": [1.0],
                              "od600": [0.16]})
        with pytest.raises(ValueError):
            ox.inoculum_spread_observable(frame, od_low=0.8, od_high=0.05)


# --------------------------------------------------------------------------------------
# Criterion (c): the plates choose Branco's number for themselves
# --------------------------------------------------------------------------------------

@pytest.mark.integration
class TestCriterionC:
    def test_the_plates_prefer_a_decay_rate_that_agrees_with_branco(self):
        frame = ox.fitted_decay_rates(growth_rate_per_h=MU)
        assert len(frame) == 3
        assert frame.fitted_k0_per_h.tolist() == pytest.approx(
            [1.9656, 0.8818, 0.6695], abs=5e-3)
        assert frame.inside_measured_band.sum() == 2
        outside = frame[~frame.inside_measured_band].iloc[0]
        assert outside.measured_k0_low_per_h / outside.fitted_k0_per_h < 1.4

    def test_letting_the_dose_decay_improves_the_fit_on_every_plate(self):
        frame = ox.fitted_decay_rates(growth_rate_per_h=MU)
        assert frame.rms_improvement.min() > 1.6

    def test_keeping_the_growth_halving_rung_destroys_it(self):
        """PINNED FAILURE, and the reason the truncation is measured rather than chosen."""
        frame = ox.fitted_decay_rates(growth_rate_per_h=MU, max_dose_mM=1.0)
        assert not frame.inside_measured_band.any()
        assert (frame.fitted_k0_per_h < 0.7).all()
        assert (frame.fitted_k0_per_h < 0.05).sum() == 2

    def test_nothing_in_the_fit_was_told_about_branco(self):
        model = ox.YapReporter(decay="fitted", growth_rate_per_h=MU)
        assert "k_ref" not in model.parameter_names
        data = ox.yap1_frame(EXPORTS[0])
        assert model.refit(data).parameters["k0_per_h"] != pytest.approx(
            float(ox.K_REF_PER_H))


class TestCrossProtection:
    def test_the_persistence_ratio_is_threshold_free(self):
        first = ox.memory_persistence_generations(MU, 0.05)["ratio"]
        second = ox.memory_persistence_generations(MU, 0.20)["ratio"]
        assert first == pytest.approx(second)

    def test_the_carrier_outlives_its_transcript_by_three_to_four_fold(self):
        low = ox.memory_persistence_generations(ox.GUAN_GROWTH_BAND_PER_H[1])["ratio"]
        high = ox.memory_persistence_generations(ox.GUAN_GROWTH_BAND_PER_H[0])["ratio"]
        assert (low, high) == pytest.approx((3.2058, 4.4055), abs=1e-3)

    def test_the_memory_is_a_generation_counter_inside_guans_own_resolution(self):
        for mu in ox.GUAN_GROWTH_BAND_PER_H:
            departure = ox.memory_persistence_generations(
                mu)["departure_from_generation_counter"]
            assert departure < ox.CROSS_PROTECTION_FLOOR.noise_floor
            assert 0.013 < departure < 0.022

    def test_a_pre_treated_culture_meets_the_challenge_with_a_carrier_and_a_naive_one_does_not(self):
        result = ox.cross_protection(growth_rate_per_h=MU)
        assert result.effector_at_challenge > 0.0
        assert result.effector_naive == 0.0
        assert 0.03 < result.fraction_retained < 0.06

    def test_the_carrier_choice_is_what_makes_the_memory_long(self):
        ctt1 = ox.cross_protection(growth_rate_per_h=MU)
        assert ctt1.persistence_ratio > 4.0
        tsa1_ratio = (1.0 + float(ox.PROTEIN_TURNOVER_TSA1_PER_H) / MU)
        assert tsa1_ratio > (1.0 + float(ox.PROTEIN_TURNOVER_CTT1_PER_H) / MU)

    def test_a_dead_culture_or_an_impossible_threshold_is_refused(self):
        with pytest.raises(ValueError):
            ox.memory_persistence_generations(0.0)
        with pytest.raises(ValueError):
            ox.memory_persistence_generations(MU, 1.0)

    def test_the_report_says_what_is_not_predicted(self):
        text = ox.cross_protection(growth_rate_per_h=MU).report()
        assert "NOT PREDICTED" in text
        assert "PMID 22851651" in text


# --------------------------------------------------------------------------------------
# The bolus / perfusion gap, which does not fall out
# --------------------------------------------------------------------------------------

class TestTheBolusPerfusionGap:
    def test_the_prediction_overshoots_the_measurement(self):
        """PINNED FAILURE. Right sign, 7.2x too large at the plates' own density."""
        gap = ox.bolus_perfusion_gap(growth_rate_per_h=MU, cells_per_ml_per_od600=BAND[1])
        assert gap.predicted_ratio == pytest.approx(14.82, abs=0.05)
        assert (gap.measured_ratio_low, gap.measured_ratio_high) == pytest.approx(
            (0.79 / 0.6, 1.24 / 0.6))
        assert gap.overshoot == pytest.approx(7.17, abs=0.05)
        assert gap.overshoot > 1.0 + ox.LETHAL_RATIO_FLOOR.noise_floor

    def test_it_overshoots_at_both_ends_of_the_swept_conversion(self):
        ratios = [ox.bolus_perfusion_gap(cells_per_ml_per_od600=c,
                                         growth_rate_per_h=MU).predicted_ratio
                  for c in BAND]
        assert min(ratios) > 4.4 and max(ratios) < 15.0
        assert all(r > 1.24 / 0.6 for r in ratios)

    def test_the_overshoot_survives_the_windows_whole_measured_growth_band(self):
        """The band the docstring quotes: mu is an INPUT here, so it is swept too."""
        ratios = [ox.bolus_perfusion_gap(cells_per_ml_per_od600=c, growth_rate_per_h=mu)
                  .predicted_ratio
                  for c in BAND
                  for mu in (PLATE_READ_4H.growth_rate_low_per_h,
                             PLATE_READ_4H.growth_rate_high_per_h)]
        assert min(ratios) == pytest.approx(4.475, abs=0.01)
        assert max(ratios) == pytest.approx(14.97, abs=0.02)
        assert min(ratios) > 1.24 / 0.6

    def test_a_peak_criterion_is_closer_than_an_exposure_criterion(self):
        gap = ox.bolus_perfusion_gap(growth_rate_per_h=MU, cells_per_ml_per_od600=BAND[1])
        assert gap.peak_criterion_error < gap.overshoot
        assert gap.overshoot / gap.peak_criterion_error == pytest.approx(5.45, abs=0.05)

    def test_tomalins_sink_would_make_it_worse_not_better(self):
        """An extra consumption term only shortens the bolus, so the ratio only grows."""
        base = ox.bolus_perfusion_gap(cells_per_ml_per_od600=BAND[1], growth_rate_per_h=MU)
        faster = ox.bolus_perfusion_gap(od600=2.0 * float(ox.PLATE_INOCULUM_OD600),
                                        cells_per_ml_per_od600=BAND[1],
                                        growth_rate_per_h=MU)
        assert faster.predicted_ratio > base.predicted_ratio
        assert faster.overshoot > base.overshoot

    def test_the_report_names_both_published_explanations_by_their_numbers(self):
        text = ox.bolus_perfusion_gap(growth_rate_per_h=MU).report()
        assert "PERFUSION" in text and "BOLUS" in text and "peak-concentration" in text


# --------------------------------------------------------------------------------------
# The GEM link. stress_ros.py prices the dose; this only converts into its units.
# --------------------------------------------------------------------------------------

class TestTheGemLink:
    def test_the_flux_is_the_consumption_rate_per_gram_of_biomass(self):
        flux = ox.gem_peroxide_uptake_mmol_per_gdcw_h(
            dose_mM=0.5, od600=0.16, cells_per_ml_per_od600=BAND[1])
        k0 = ox.dose_decay_rate_per_h(0.16, BAND[1])
        assert flux == pytest.approx(k0 * 0.5 / (0.42 * 0.16))
        assert flux > 0.0

    def test_the_dose_is_a_large_flux_for_a_short_time(self):
        band = [ox.gem_peroxide_uptake_mmol_per_gdcw_h(
            dose_mM=0.5, od600=0.16, cells_per_ml_per_od600=c) for c in BAND]
        assert band == pytest.approx([6.514, 24.429], abs=0.01)
        assert max(band) > 2.0 * 10.0

    def test_the_flux_per_gram_does_not_depend_on_the_inoculum(self):
        """k scales with OD and so does the biomass it is divided by, so they cancel."""
        first = ox.gem_peroxide_uptake_mmol_per_gdcw_h(
            dose_mM=0.5, od600=0.05, cells_per_ml_per_od600=BAND[1])
        second = ox.gem_peroxide_uptake_mmol_per_gdcw_h(
            dose_mM=0.5, od600=0.80, cells_per_ml_per_od600=BAND[1])
        assert first == pytest.approx(second)

    def test_a_dead_culture_cannot_be_given_a_per_gram_flux(self):
        with pytest.raises(ValueError):
            ox.gem_peroxide_uptake_mmol_per_gdcw_h(
                dose_mM=0.5, od600=0.0, cells_per_ml_per_od600=BAND[1])


@pytest.mark.integration
@pytest.mark.integration
class TestTheOneStructuralAssumptionIsDeclaredAndUntestable:
    """k proportional to X to the first power. Pinned as UNTESTED, not as confirmed."""

    @pytest.fixture(scope="class")
    def identifiability(self):
        return ox.density_exponent_identifiability(growth_rate_per_h=MU)

    def test_no_plate_prefers_a_free_exponent_above_the_reporter_floor(self, identifiability):
        assert not bool(identifiability.clears_floor.any())
        assert list(identifiability.rms_improvement) == pytest.approx(
            [1.124, 1.000, 1.021], abs=2e-3)

    def test_the_three_replicates_disagree_by_six_fold_and_one_sits_on_the_bound(
            self, identifiability):
        exponents = list(identifiability.fitted_exponent)
        assert exponents == pytest.approx([6.000, 1.052, 2.916], abs=2e-3)
        assert max(exponents) / min(exponents) > 5.0
        assert int(identifiability.at_search_bound.sum()) == 1

    def test_the_exponent_is_refused_rather_than_carried_as_a_parameter(self):
        refusal = ox.NOT_BUILT["density_exponent"]
        assert refusal.tag is Tag.REFUSED
        assert "four inocula" in refusal.missing
        assert "density_exponent" not in {p.name for p in ox.OXIDATIVE_PARAMS.free_scalars()}

    def test_the_default_exponent_changes_nothing_anywhere(self):
        t = np.linspace(0.0, 4.14, 25)
        k0 = ox.dose_decay_rate_per_h(0.16, BAND[1])
        assert ox.extracellular_h2o2(t, 0.5, k0, MU, 1.0) == pytest.approx(
            ox.extracellular_h2o2(t, 0.5, k0, MU))
        assert ox.driver_horizon_h(k0, MU, density_exponent=1.0) == pytest.approx(
            ox.driver_horizon_h(k0, MU))

    def test_a_zero_exponent_is_the_ungrowing_well_and_solves_its_own_ode(self):
        """n = 0 means the consuming density never changes, so E is a plain exponential."""
        t = np.linspace(0.0, 4.14, 25)
        assert ox.extracellular_h2o2(t, 0.5, 1.3, MU, 0.0) == pytest.approx(
            0.5 * np.exp(-1.3 * t))
        assert ox.driver_horizon_h(1.3, MU, density_exponent=0.0) == pytest.approx(
            math.log(1.0 / ox.DRIVER_FLOOR) / 1.3)

    def test_a_negative_exponent_can_leave_the_well_a_floor_of_its_own(self):
        """A dying consumer leaves a residue; there is a horizon only if the residue is below it."""
        assert math.exp(1.0 / -MU) > ox.DRIVER_FLOOR
        assert ox.driver_horizon_h(1.0, MU, density_exponent=-1.0) == math.inf
        assert math.exp(10.0 / -MU) < ox.DRIVER_FLOOR
        horizon = ox.driver_horizon_h(10.0, MU, density_exponent=-1.0)
        assert math.isfinite(horizon)
        assert ox.extracellular_h2o2(horizon, 1.0, 10.0, MU, -1.0) == pytest.approx(
            ox.DRIVER_FLOOR, rel=1e-9)


@pytest.mark.integration
class TestTheUnfittedRateBeatsAFittedOne:
    """Criterion (c) at its sharpest: Branco 2004 out-predicts a rate learned on the plates.

    Pinned to the digit because the whole claim is that a number nothing here was fitted to
    beats one that was. A change that quietly moves it has to come through this test.
    """

    @pytest.fixture(scope="class")
    def transfer(self):
        return ox.transfer_comparison(growth_rate_per_h=MU)

    def test_the_transferred_rates_and_their_held_out_residuals(self, transfer):
        assert list(transfer.k0_transferred_per_h) == pytest.approx(
            [0.765, 1.179, 1.351], abs=2e-3)
        assert list(transfer.rms_transferred) == pytest.approx(
            [285.6, 175.4, 184.5], abs=0.2)

    def test_brancos_rate_beats_the_transferred_fit_on_every_plate(self, transfer):
        assert list(transfer.rms_branco_low_conversion) == pytest.approx(
            [269.4, 169.3, 152.0], abs=0.2)
        assert bool(transfer.branco_beats_transferred.all())

    def test_and_beats_the_incumbent_by_between_1_63_and_1_75x(self, transfer):
        gain = transfer.rms_incumbent / transfer.rms_branco_low_conversion
        assert list(transfer.rms_incumbent) == pytest.approx([472.2, 280.2, 247.7], abs=0.2)
        assert gain.min() == pytest.approx(1.63, abs=0.01)
        assert gain.max() == pytest.approx(1.75, abs=0.01)

    def test_the_other_end_of_the_swept_band_is_not_the_interesting_one(self, transfer):
        """It loses on one plate of three, which is what a SWEPT axis is allowed to do."""
        beats = transfer.rms_branco_high_conversion < transfer.rms_incumbent
        assert int(beats.sum()) == 2
        assert list(transfer.rms_branco_high_conversion) == pytest.approx(
            [233.4, 243.6, 265.0], abs=0.2)


class TestTheGemPricesIt:
    @pytest.fixture(scope="class")
    def _ros_template(self, yeast_gem_factory):
        from ystwin.fba.solver import configure
        from ystwin.fba.stress_ros import add_ros_module

        yeast_gem = yeast_gem_factory()
        configure(yeast_gem)
        return add_ros_module(yeast_gem)

    @pytest.fixture
    def ros(self, _ros_template, model_copy):
        return model_copy(_ros_template)

    def test_the_plates_dose_is_a_flux_the_gem_can_carry(self, ros):
        from ystwin.fba.stress_ros import peroxide_dose_response

        flux = ox.gem_peroxide_uptake_mmol_per_gdcw_h(
            dose_mM=0.5, od600=0.16, cells_per_ml_per_od600=BAND[1])
        frame = peroxide_dose_response(ros, [0.0, flux])
        assert bool(frame.feasible.all())

    def test_the_growth_cost_is_zero_with_catalase_and_measurable_without_it(self, ros):
        """Both ends of the swept conversion, priced by stress_ros rather than re-derived."""
        from ystwin.fba.stress_ros import peroxide_dose_response

        fluxes = [ox.gem_peroxide_uptake_mmol_per_gdcw_h(
            dose_mM=0.5, od600=0.16, cells_per_ml_per_od600=c) for c in BAND]
        wild = peroxide_dose_response(ros, [0.0] + fluxes)
        assert wild.growth_cost_fraction.abs().max() < 1e-9
        null = peroxide_dose_response(ros, [0.0] + fluxes, knockouts=("YGR088W",))
        assert list(null.growth_cost_fraction[1:]) == pytest.approx([0.0545, 0.2043],
                                                                   abs=5e-4)

    def test_the_gem_prices_a_wild_type_dose_at_nothing_and_a_catalase_null_at_something(
            self, ros):
        from ystwin.fba.stress_ros import clearance_cost

        flux = ox.gem_peroxide_uptake_mmol_per_gdcw_h(
            dose_mM=0.5, od600=0.16, cells_per_ml_per_od600=BAND[1])
        cost = clearance_cost(ros, ("YGR088W", "YDR256C"), dose=flux)
        assert cost.undosed_cost < 1e-5
        assert cost.dosed_cost > 0.15
        assert cost.dosed_cost > 1e4 * cost.undosed_cost


# --------------------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------------------

class TestProvenance:
    def test_the_summary_carries_the_gate_and_both_registries(self):
        text = ox.provenance_summary()
        assert "yap1_three_state" in text and "not-built" in text
        assert "2 free scalars against 3 independent targets -- PASSES" in text
        assert "**REFUSED**" in text and "**MEASURED**" in text

    def test_the_module_docstring_states_the_numbers_the_tests_pin(self):
        doc = ox.__doc__
        for token in ("0.876-3.283", "1.966, 0.882 and 0.669", "7.2x too large",
                      "3.21-4.41x", "4.73-7.04x", "0.2065", "4.48-14.97x"):
            assert token in doc

    def test_every_public_name_exists(self):
        for name in ox.__all__:
            assert hasattr(ox, name), name
