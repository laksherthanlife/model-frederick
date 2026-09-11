"""The pH / weak-acid block, and every claim its report makes.

Four things are pinned here that nothing else in the suite pins.

  * THE COUNT. `mech/ph.py` is BUILD_REDUCED because `ARCHITECTURE_GAPS.md` 0.3 counted the
    full slice at eight free scalars against two usable targets and criterion (e) refuses
    that. The gate is computed at import; these tests fail if it stops passing, if a free
    scalar creeps back in, or if the six that were dropped stop being recorded as refusals.
  * THE VESICLE TRAP. Gabba 2020 measured the same acids twice -- across the yeast plasma
    membrane and across a POPE:POPG:POPC bilayer -- and the two differ by 192-707x.
    `ARCHITECTURE_TARGET.md` quotes the vesicle time constants (10 ms, 9.8 ms) three lines
    below the sentence forbidding the substitution. A test asserts the yeast values.
  * ABBOTT'S OWN ARITHMETIC. Every partitioning number here was read out of PMID 18676708
    and every one is re-derived from the module's pKa constants, so a drifted pKa fails.
  * THE REGIME. ``s_ATP`` is a two-valued step and the repository's own reference operating
    point sits on the 2.32x branch. A test fails if the block ever prices an ATP debit at
    the respiratory slope.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ystwin.bridge.maintenance_calibration import MAINTENANCE_REGIMES, REFERENCE_REGIME
from ystwin.generator import redox
from ystwin.generator.panel_experiment import (
    MEASURED_GROWTH_RATE_SE,
    OBSERVED_ACTIVITY_CV,
)
from ystwin.generator.stress_panel import STRESSORS
from ystwin.mech import ph
from ystwin.mech.integrate import Verdict
from ystwin.mech.params import RefusedValue, SweptValue, Tag
from ystwin.mech.state import FEDBATCH_5D, PLATE_READ_4H
from ystwin.pathway.thermo_gate import CYTOSOLIC_VOLUMES_ML_PER_GDCW
from ystwin.photophysics import citrine_ph_response

PKA_ACETIC = float(ph.PKA_ACETIC)
PKA_LACTIC = float(ph.PKA_LACTIC)
PKA_BENZOIC = float(ph.PKA_BENZOIC)


# --------------------------------------------------------------------------------------
# Criterion (e): the count that decides whether this block may exist at all
# --------------------------------------------------------------------------------------

class TestFreeScalarGate:

    def test_gate_passes_and_is_computed_not_claimed(self):
        gate = ph.PH_PARAMS.gate()
        assert gate.passes
        assert gate is not ph.PH_GATE or gate.summary() == ph.PH_GATE.summary()

    def test_exactly_two_free_scalars_and_they_are_the_two_that_survived(self):
        free = {p.name for p in ph.PH_PARAMS.free_scalars()}
        assert free == {"beta_cytosolic", "cytosolic_volume"}, (
            "ARCHITECTURE_GAPS.md 0.3 counted eight; this block ships two. Anything else "
            "here is a scalar that crept back in without a target to pay for it")

    def test_three_independent_targets_none_of_them_fitted(self):
        targets = ph.PH_PARAMS.independent_targets()
        assert {t.name for t in targets} == {
            "half_yield_transfer", "growth_rate", "reporter_activity"}
        assert all(not t.fitted for t in targets)

    def test_the_six_dropped_scalars_are_recorded_as_refusals_not_forgotten(self):
        # The eight of GAPS 0.3 minus the two that survived, plus the leak's missing V_m.
        dropped = {"jmax_pump", "pka_pma1", "n_pma1", "vmax_anion_export", "k_anion_export",
                   "g_proton_leak"}
        assert dropped <= {p.name for p in ph.NOT_BUILT.params}
        for name in dropped:
            assert ph.NOT_BUILT[name].tag == Tag.REFUSED

    def test_refusals_live_outside_the_gate(self):
        assert not (set(ph.NOT_BUILT._params) & set(ph.PH_PARAMS._params))
        assert len(ph.NOT_BUILT.refusals()) == len(ph.NOT_BUILT)

    def test_a_third_free_scalar_would_fail_the_gate(self):
        from ystwin.mech.params import FreeScalarGate
        free = tuple(p.name for p in ph.PH_PARAMS.free_scalars()) + ("one_more", "and_another")
        targets = tuple(t.name for t in ph.PH_PARAMS.independent_targets())
        assert not FreeScalarGate("hypothetical", free, targets).passes


class TestProvenance:

    def test_every_pinned_constant_carries_a_resolvable_identifier(self):
        assert ph.PH_PARAMS.uncited() == (), (
            "a number attributed to a paper by prose alone cannot be checked, and that is "
            "the failure mode this project has been burned by repeatedly")

    def test_no_constant_is_borrowed_from_another_organism(self):
        assert ph.PH_PARAMS.borrowed() == (), (
            "every Neurospora and squid-axon number in the Kahm parameter table is in "
            "NOT_BUILT, not in the built block")

    def test_swept_parameters_refuse_to_collapse_to_a_point(self):
        for name in ("beta_cytosolic", "cytosolic_volume"):
            with pytest.raises(SweptValue):
                float(ph.PH_PARAMS[name])

    def test_every_refusal_raises_and_names_what_would_close_it(self):
        for param in ph.NOT_BUILT.params:
            with pytest.raises(RefusedValue) as excinfo:
                float(param)
            assert param.missing and param.missing in str(excinfo.value)

    def test_the_swept_volume_band_is_imported_not_retyped(self):
        low, high = ph.CYTOSOLIC_VOLUME_ML_PER_GDCW.bounds
        assert (low, high) == (min(CYTOSOLIC_VOLUMES_ML_PER_GDCW),
                               max(CYTOSOLIC_VOLUMES_ML_PER_GDCW))

    def test_provenance_summary_names_the_gate_and_the_regime(self):
        text = ph.provenance_summary()
        assert "2 free scalars against 3 independent targets" in text
        assert "MODEL-DERIVED / IN_SILICO" in text


# --------------------------------------------------------------------------------------
# The acids, and the three this block refuses by name
# --------------------------------------------------------------------------------------

class TestAcids:

    def test_only_the_three_acids_with_a_measured_yeast_permeability_are_carried(self):
        assert set(ph.ACIDS) == {"acetic", "formic", "lactic"}

    @pytest.mark.parametrize("name", ["propionic", "benzoic", "sorbic"])
    def test_the_other_three_of_abbotts_series_are_refused_by_name(self, name):
        with pytest.raises(ph.AcidUnmeasured) as excinfo:
            ph.acid(name)
        message = str(excinfo.value)
        assert "VESICLE" in message and "192-707x" in message

    def test_lactic_permeability_is_a_measured_zero_not_an_absent_value(self):
        assert ph.ACIDS["lactic"].permeability.tag == Tag.MEASURED
        assert float(ph.ACIDS["lactic"].permeability) == 0.0
        assert not ph.ACIDS["lactic"].permeates
        assert ph.ACIDS["acetic"].permeates

    def test_permeation_time_constants_are_the_yeast_values_not_the_vesicle_ones(self):
        # ARCHITECTURE_TARGET quotes 10 ms and 9.8 ms; those are the bilayer. The
        # geometry-independent statement is the paper's own P_vesicle/P_yeast ratio.
        assert ph.permeation_tau_s("acetic") == pytest.approx(3.0, abs=0.2)
        assert ph.permeation_tau_s("formic") == pytest.approx(3.8, abs=0.2)
        assert ph.permeation_tau_s("acetic") > 1.0
        vesicle_tau = 3600.0 / (float(ph.SURFACE_TO_VOLUME_PER_CM) * 990e-5 * 3600.0)
        assert ph.permeation_tau_s("acetic") / vesicle_tau == pytest.approx(707.0, rel=1e-3)

    def test_an_impermeant_acid_has_no_permeation_time_constant(self):
        with pytest.raises(ph.AcidUnmeasured):
            ph.permeation_tau_s("lactic")

    def test_surface_to_volume_is_area_over_cell_water_not_area_over_cell(self):
        # Gabba 2020's cell is V0 = 81.9 fL and its OSMOTIC volume is V0 - b = 38.6 fL.
        # 38.6 fL was previously used as the cell volume, which shrank the sphere.
        assert float(ph.CELL_VOLUME_FL) == 81.9
        assert float(ph.CELL_WATER_VOLUME_FL) == 38.6
        radius = (3.0 * float(ph.CELL_VOLUME_FL) * 1e-12 / (4.0 * math.pi)) ** (1.0 / 3.0)
        area = 4.0 * math.pi * radius ** 2
        assert float(ph.SURFACE_TO_VOLUME_PER_CM) == pytest.approx(
            area / (float(ph.CELL_WATER_VOLUME_FL) * 1e-12))
        assert float(ph.SURFACE_TO_VOLUME_PER_CM) == pytest.approx(2.363e4, rel=1e-3)

    def test_the_geometry_reproduces_the_papers_own_cell_diameter(self):
        # The unstirred-layer section of PMC6976801 sets the cell's characteristic length
        # (its diameter) to 5.4 um. A sphere of 38.6 fL would be 4.19 um and is not the cell.
        radius = (3.0 * float(ph.CELL_VOLUME_FL) * 1e-12 / (4.0 * math.pi)) ** (1.0 / 3.0)
        assert radius * 2.0 * 1e4 == pytest.approx(5.4, abs=0.05)

    def test_the_cell_water_volume_is_not_used_as_a_cell_volume(self):
        # The regression this locks: area from V0, concentration in V0 - b.
        assert not hasattr(ph, "REFERENCE_CELL_VOLUME_FL")
        water_sphere = 3.0 / ((3.0 * float(ph.CELL_WATER_VOLUME_FL) * 1e-12
                               / (4.0 * math.pi)) ** (1.0 / 3.0))
        assert float(ph.SURFACE_TO_VOLUME_PER_CM) / water_sphere == pytest.approx(
            1.651, rel=1e-3), (
            "the previous value 1.431e4 /cm was 3/r on the cell WATER volume; the correct "
            "ratio is 1.65x larger and every entry time constant computed from it was that "
            "much too slow")


# --------------------------------------------------------------------------------------
# Henderson-Hasselbalch, re-derived against the papers' own printed arithmetic
# --------------------------------------------------------------------------------------

class TestPartitioning:

    def test_abbotts_lactate_conversions_reproduce(self):
        # "900 mM of lactic acid (~61 mM undissociated acid)" at pH 5.
        assert ph.undissociated_outside(900.0, 5.0, PKA_LACTIC) == pytest.approx(61.0, abs=0.5)
        # "85 mM (~60 mM undissociated acid at pH 3.5)".
        assert ph.undissociated_outside(85.0, 3.5, PKA_LACTIC) == pytest.approx(60.0, abs=1.0)

    def test_abbotts_benzoate_conversions_reproduce_and_pin_the_pka(self):
        # "2 mM (pH 5) and 0.3 mM (pH 3.5), corresponding to 0.27 mM and 0.25 mM".
        assert ph.undissociated_outside(2.0, 5.0, PKA_BENZOIC) == pytest.approx(0.27, abs=0.005)
        assert ph.undissociated_outside(0.3, 3.5, PKA_BENZOIC) == pytest.approx(0.25, abs=0.005)

    def test_the_repos_own_acetic_ec50_agrees_with_abbott_in_the_permeating_variable(self):
        ec50 = STRESSORS["acetic_acid"].ec50
        repo = ph.undissociated_outside(ec50, 4.5, PKA_ACETIC)
        abbott = ph.undissociated_outside(105.0, 5.0, PKA_ACETIC)
        assert repo == pytest.approx(abbott, rel=0.01), (
            "the panel's 60 mM 'at pH 4.5' is 38.7 mM undissociated and Abbott's measured "
            "half-yield 105 mM at pH 5.0 is 38.4 mM: the dose is corroborated to 1%. The "
            "defect is the missing citation and the non-separability, not the number")

    def test_a_fixed_ph_steady_state_demands_an_impossible_anion_pool(self):
        ah_out = ph.undissociated_outside(105.0, 5.0, PKA_ACETIC)
        assert ph.trapped_total(ah_out, 7.08, PKA_ACETIC) == pytest.approx(8051.0, rel=1e-3)
        bill = ph.proton_bill(ah_out, 0.10, 2.0, PKA_ACETIC)
        assert not bill.physical, (
            "8 M internal acetate is the argument for the two bounds; if this ever passes "
            "as physical the bracket has quietly become a point estimate")

    def test_dose_and_medium_ph_are_not_separable(self):
        at_low = ph.undissociated_outside(40.0, 4.0, PKA_ACETIC)
        at_high = ph.undissociated_outside(40.0, 5.5, PKA_ACETIC)
        assert at_low / at_high == pytest.approx(5.53, rel=0.01), (
            "one Hill EC50 on total dose cannot express a 5.5x change in the permeating "
            "species at a fixed dose")

    def test_the_trapped_acid_buffers_itself_by_the_standard_monoprotic_formula(self):
        beta, a_i, ph_c = 145.0, 6704.0, 7.0
        f = 1.0 - ph.undissociated_fraction(ph_c, PKA_ACETIC)
        expected = beta + math.log(10.0) * a_i * f * (1.0 - f)
        assert ph.buffer_capacity(beta, a_i, ph_c, PKA_ACETIC) == pytest.approx(expected)
        assert ph.buffer_capacity(beta, a_i, ph_c, PKA_ACETIC) > 1.5 * beta


# --------------------------------------------------------------------------------------
# The two bounds
# --------------------------------------------------------------------------------------

class TestBounds:

    def test_no_pump_bound_solves_its_own_proton_balance(self):
        ah_out = ph.undissociated_outside(40.0, 4.0, PKA_ACETIC)
        bound = ph.no_pump_ph(ah_out, 145.0, PKA_ACETIC)
        assert 145.0 * bound.drop == pytest.approx(bound.trapped_anion_mM, rel=1e-6)
        assert bound.ph_c == pytest.approx(5.568, abs=0.002)

    def test_a_bigger_buffer_gives_a_smaller_drop(self):
        ah_out = ph.undissociated_outside(40.0, 4.0, PKA_ACETIC)
        drops = [ph.no_pump_ph(ah_out, b, PKA_ACETIC).drop
                 for b in ph.CYTOSOLIC_BUFFER_CAPACITY.sweep_points(5)]
        assert drops == sorted(drops, reverse=True)

    def test_zero_dose_moves_nothing(self):
        bound = ph.no_pump_ph(0.0, 145.0, PKA_ACETIC)
        assert bound.drop == 0.0
        assert bound.ph_c == float(ph.RESTING_CYTOSOLIC_PH)

    def test_the_bill_is_the_mass_balance_and_scales_with_growth_and_volume(self):
        ah_out = ph.undissociated_outside(105.0, 5.0, PKA_ACETIC)
        one = ph.proton_bill(ah_out, 0.10, 2.0, PKA_ACETIC)
        faster = ph.proton_bill(ah_out, 0.20, 2.0, PKA_ACETIC)
        bigger = ph.proton_bill(ah_out, 0.10, 4.0, PKA_ACETIC)
        assert faster.mmol_atp_per_gdcw_h == pytest.approx(2 * one.mmol_atp_per_gdcw_h)
        assert bigger.mmol_atp_per_gdcw_h == pytest.approx(2 * one.mmol_atp_per_gdcw_h)
        assert one.mmol_atp_per_gdcw_h == pytest.approx(
            0.10 * one.trapped_anion_mM * 2.0e-3)

    def test_the_bill_at_abbotts_half_yield_is_a_sixth_of_the_catabolic_supply(self):
        ah_out = ph.undissociated_outside(105.0, 5.0, PKA_ACETIC)
        fractions = [ph.proton_bill(ah_out, 0.10, v, PKA_ACETIC).fraction_of_catabolic_atp
                     for v in CYTOSOLIC_VOLUMES_ML_PER_GDCW]
        assert min(fractions) == pytest.approx(0.066, abs=0.002)
        assert max(fractions) == pytest.approx(0.179, abs=0.002)
        assert all(f < 0.60 for f in fractions), (
            "Holyoak 1996's 'up to 60% of ATP under weak acid stress' was inferred from ATP "
            "depletion rather than measured as a flux. This is the model's PREDICTION "
            "against it, not a parameter taken from it")

    def test_the_catabolic_denominator_is_two_atp_per_measured_glucose(self):
        assert float(ph.ANAEROBIC_CATABOLIC_ATP) == pytest.approx(
            2.0 * float(ph.ABBOTT_REFERENCE_GLUCOSE))


# --------------------------------------------------------------------------------------
# B1 / B3 / B4 / B9 -- and the regime correction
# --------------------------------------------------------------------------------------

class TestBranches:

    def test_growth_uses_the_regime_slope_and_not_the_respiratory_scalar(self):
        ah_out = ph.undissociated_outside(40.0, 4.5, PKA_ACETIC)
        bill = ph.proton_bill(ah_out, 0.245, 2.0, PKA_ACETIC)
        fermentative = ph.growth_under_bill(0.245, bill)
        respiratory = ph.growth_under_bill(
            0.245, bill, regime=MAINTENANCE_REGIMES["respiratory_glucose_limited"])
        assert fermentative < respiratory, (
            "the old scalar 0.004642 priced every ATP debit at 43% of its value at the "
            "operating point this repository validates its GEM against")
        assert (REFERENCE_REGIME.ngam_growth_slope_per_h_per_mmol_atp
                / MAINTENANCE_REGIMES["respiratory_glucose_limited"]
                .ngam_growth_slope_per_h_per_mmol_atp) == pytest.approx(2.32, abs=0.01)

    def test_s_atp_is_flagged_as_model_derived_wherever_it_is_taken_from(self):
        assert "IN_SILICO" in REFERENCE_REGIME.provenance
        assert "No culture was read" in REFERENCE_REGIME.provenance

    def test_no_dose_costs_no_growth(self):
        bill = ph.proton_bill(0.0, 0.245, 2.0, PKA_ACETIC)
        assert ph.growth_under_bill(0.245, bill) == pytest.approx(0.245)

    def test_growth_falls_monotonically_with_dose_and_never_goes_negative(self):
        rates = [ph.growth_under_bill(
            0.245, ph.proton_bill(ph.undissociated_outside(d, 4.5, PKA_ACETIC),
                                  0.245, 2.0, PKA_ACETIC))
            for d in (0.0, 5.0, 10.0, 20.0, 40.0)]
        assert rates == sorted(rates, reverse=True)
        assert all(r > 0.0 for r in rates)

    def test_the_atp_bill_is_far_too_weak_at_the_measured_half_effect_dose(self):
        """B1's own half-effect dose against the two independent measured ones.

        Abbott 2008's half-yield is 38.35 mM undissociated acetic acid and this repository's
        panel EC50 is 38.72 mM; the ATP bill needs 4.3-11.6x that to halve anything.
        """
        from scipy.optimize import brentq
        halves = []
        for volume in (min(CYTOSOLIC_VOLUMES_ML_PER_GDCW), max(CYTOSOLIC_VOLUMES_ML_PER_GDCW)):
            model = ph._SeparabilityGrowth(cytosolic_volume_ml_per_gdcw=volume,
                                           reference_growth_per_h=0.245)
            dose = brentq(lambda d: model.mu(d, 4.5) - 0.5 * 0.245, 1.0, 1e6)
            halves.append(ph.undissociated_outside(dose, 4.5, PKA_ACETIC))
        assert min(halves) == pytest.approx(164.5, rel=0.01)
        assert max(halves) == pytest.approx(444.1, rel=0.01)
        measured = ph.undissociated_outside(105.0, 5.0, PKA_ACETIC)
        assert min(halves) / measured == pytest.approx(4.29, abs=0.05)
        assert max(halves) / measured == pytest.approx(11.58, abs=0.05)

    def test_orij_law_is_the_divisor_form_and_carries_its_measured_slope(self):
        assert ph.orij_growth_per_h(float(ph.ORIJ_INTERCEPT)) == pytest.approx(1.0)
        decade = math.log10(ph.orij_growth_per_h(7.0) / ph.orij_growth_per_h(6.0))
        assert decade == pytest.approx(1.0 / float(ph.ORIJ_SLOPE))
        assert decade == pytest.approx(1.370, abs=0.001)
        # The paper's Figure 4a data sit near 0.48 /h at pH_c 7.0; the multiplier reading
        # of the same fit would give 0.67 and is what this check excludes.
        assert ph.orij_growth_per_h(7.0) == pytest.approx(0.48, abs=0.02)

    def test_orij_band_brackets_the_point_and_is_about_1_4x_wide_per_0_1_ph(self):
        low, high = ph.orij_growth_band_per_h(6.5, ph_uncertainty=0.0)
        point = ph.orij_growth_per_h(6.5)
        assert low <= point <= high
        wide_low, wide_high = ph.orij_growth_band_per_h(6.5, ph_uncertainty=0.1)
        assert wide_low < low and wide_high > high

    def test_the_atp_route_opposes_orij_rather_than_double_counting_it(self):
        ah_out = ph.undissociated_outside(105.0, 5.0, PKA_ACETIC)
        slope = ph.atp_route_slope_decades_per_ph(ah_out, 7.08, 2.0, PKA_ACETIC)
        assert slope < 0.0, (
            "ARCHITECTURE_TARGET B9 says the ATP route and Orij's law double-count. With a "
            "weak acid present they oppose: a higher pH_c traps more anion and costs MORE "
            "pumping, while Orij says a higher pH_c grows faster")
        assert abs(slope) * float(ph.ORIJ_SLOPE) == pytest.approx(0.126, abs=0.005)

    def test_the_atp_route_is_identically_zero_without_an_acid(self):
        assert ph.atp_route_slope_decades_per_ph(0.0, 7.08, 2.0, PKA_ACETIC) == 0.0

    def test_e6_rescoped_reports_the_atp_borne_fraction_as_a_band(self):
        ah_out = ph.undissociated_outside(105.0, 5.0, PKA_ACETIC)
        low, high = ph.atp_borne_fraction_of_orij_slope(ah_out)
        assert (low, high) == pytest.approx((0.0630, 0.1702), abs=5e-4), (
            "REVISED_BUILD_LIST re-scopes E6 to 'how much of the pH-growth coupling is "
            "ATP-borne'. This block's answer is 6.3-17.0%, computed, not the 34% quoted")
        assert ph.atp_borne_fraction_of_orij_slope(0.0) == (0.0, 0.0), (
            "Orij titrated pH_c with HCl and no weak acid; weak-acid influx is the only "
            "proton load this module carries, so its ATP route is exactly zero there")

    def test_the_quoted_thirty_four_percent_needs_a_volume_outside_the_declared_band(self):
        ah_out = ph.undissociated_outside(105.0, 5.0, PKA_ACETIC)
        anion = ph.trapped_anion(ah_out, float(ph.RESTING_CYTOSOLIC_PH), PKA_ACETIC)
        slope = REFERENCE_REGIME.ngam_growth_slope_per_h_per_mmol_atp
        needed = 0.473 / (slope * anion / 1000.0)
        assert needed == pytest.approx(5.477, abs=0.01)
        assert needed / max(CYTOSOLIC_VOLUMES_ML_PER_GDCW) == pytest.approx(2.03, abs=0.01), (
            "so the 0.473 decades is not reachable from these equations at any point in the "
            "band pathway/thermo_gate.py declares, and this pass does not adopt it")

    def test_citrine_quench_routes_through_photophysics_and_is_not_reimplemented(self):
        expected = citrine_ph_response(5.6) / citrine_ph_response(float(ph.RESTING_CYTOSOLIC_PH))
        assert ph.citrine_quench(5.6) == pytest.approx(expected)
        assert ph.citrine_quench(float(ph.RESTING_CYTOSOLIC_PH)) == pytest.approx(1.0)

    def test_the_quench_is_a_third_of_the_signal_at_a_survivable_dose(self):
        bound = ph.no_pump_ph(ph.undissociated_outside(40.0, 4.5, PKA_ACETIC), 145.0,
                              PKA_ACETIC)
        assert ph.citrine_quench(bound.ph_c) < 0.7, (
            "measured induction folds here are ~1.5, so a quench below 0.7 is a dose-"
            "structured artefact comparable to the whole signal")

    def test_redox_shift_routes_through_redox_py_and_exceeds_the_real_response(self):
        bound = ph.no_pump_ph(ph.undissociated_outside(40.0, 4.0, PKA_ACETIC), 145.0,
                              PKA_ACETIC)
        shift = ph.apparent_redox_shift_mv(bound.ph_c)
        assert shift == pytest.approx(
            redox.apparent_shift_from_ph(float(ph.RESTING_CYTOSOLIC_PH), bound.ph_c))
        assert shift > 50.0, (
            "a genuine peroxide response is 40-50 mV and carries the same sign, so the two "
            "cannot be told apart without measuring pH")

    def test_both_correcting_functions_are_actually_called_not_merely_reproduced(
            self, monkeypatch):
        # Equal arithmetic is not a call. The point of this block is that these two
        # previously-uncalled functions become reachable, so patch them and require a change.
        monkeypatch.setattr(ph, "apparent_shift_from_ph", lambda a, b: 12345.0)
        monkeypatch.setattr(ph, "citrine_ph_response", lambda p: 2.0)
        assert ph.apparent_redox_shift_mv(6.0) == 12345.0, (
            "mech/ph.py must call generator/redox.py::apparent_shift_from_ph, not re-derive "
            "NERNST_MV_PER_PH * (ref - pH_c) and claim it routes through it")
        assert ph.citrine_quench(6.0) == pytest.approx(1.0), (
            "a constant citrine_ph_response must make the normalised quench exactly 1"
        )


# --------------------------------------------------------------------------------------
# Criterion (d): the tau/T audit, and the fact that the two windows do not nest
# --------------------------------------------------------------------------------------

class TestReduction:

    @pytest.fixture(scope="class")
    def ah_out(self):
        return ph.undissociated_outside(40.0, 4.0, PKA_ACETIC)

    def test_the_fast_state_reduces_in_both_windows(self, ah_out):
        for window in (PLATE_READ_4H, FEDBATCH_5D):
            audit = ph.reduction_audit(window, ah_out_mM=ah_out)
            assert audit.verdict("A_i") is Verdict.ELIMINATE

    def test_ph_eliminates_in_both_windows_now_that_charge_is_conserved(self, ah_out):
        """RETRACTS the block's criterion-4 headline, which was that the verdict FLIPS.

        It read KEEP on the 4.14 h plate and ELIMINATE in the 5 d vessel, and the block
        called that flip its answer to criterion (4). The flip was the missing dilution term
        in `weak_acid_rhs`: the tau that produced it was the charge-balance error's own
        relaxation, which scales as 1/mu because the error does. Repaired, pH_c is algebra in
        both windows -- and the two windows' reduction sets are now identical, so the claim
        that this block's irreducible sets are non-nested goes with it.
        """
        plate = ph.reduction_audit(PLATE_READ_4H, ah_out_mM=ah_out)
        vessel = ph.reduction_audit(FEDBATCH_5D, ah_out_mM=ah_out)
        assert plate.verdict("pH_c") is Verdict.ELIMINATE
        assert vessel.verdict("pH_c") is Verdict.ELIMINATE
        assert {s.name for s in plate.integrated} == {s.name for s in vessel.integrated}
        assert not plate.integrated.variables

    def test_the_eliminated_state_biases_the_time_integral_far_below_the_assay_floor(
            self, ah_out):
        # It used to bias it by 53-74%, i.e. above the floor, which is why it was KEPT.
        row = ph.reduction_audit(PLATE_READ_4H, ah_out_mM=ah_out).row("pH_c")
        assert row.bias_high < OBSERVED_ACTIVITY_CV / 100.0
        assert 4e-4 < row.ratio_low < row.ratio_high < 1e-3

    def test_the_mode_is_entry_driven_and_no_longer_scales_with_growth(self, ah_out):
        """RETRACTS "the slow mode is growth-driven", which was true of the error only.

        The old test asserted tau_vessel/tau_plate == mu_plate/mu_vessel = 2.43. The one
        physical mode is set by entry, so the same dose gives the same band in a vessel
        growing 2.4x slower.
        """
        plate = ph.reduction_audit(PLATE_READ_4H, ah_out_mM=ah_out).row("pH_c")
        vessel = ph.reduction_audit(FEDBATCH_5D, ah_out_mM=ah_out).row("pH_c")
        assert vessel.tau_high_h / plate.tau_high_h == pytest.approx(1.0, abs=0.01)
        assert vessel.tau_low_h / plate.tau_low_h == pytest.approx(1.0, abs=0.01)

    def test_the_tau_band_in_the_provenance_string_is_computed_not_quoted(self):
        # A hardcoded band goes stale at another dose; this one has to move with it.
        working = ph.ph_states(PLATE_READ_4H, ah_out_mM=ph.undissociated_outside(
            40.0, 4.0, PKA_ACETIC))["pH_c"].source
        abbott = ph.ph_states(PLATE_READ_4H, ah_out_mM=ph.undissociated_outside(
            105.0, 5.0, PKA_ACETIC))["pH_c"].source
        assert "0.001665-0.002883 h" in working
        assert "0.001585-0.002665 h" in abbott, (
            "the band for Abbott's half-yield acetate; if this moves, the module docstring "
            "is stale and so is the chain diagram")
        assert "ENTRY-driven" in working and "1/mu" in working

    def test_the_tau_band_is_the_same_at_both_windows_within_one_percent(self):
        # Entry-driven, so it is tau itself and no longer tau*mu that is window-invariant.
        for total, ph_ex in ((40.0, 4.0), (105.0, 5.0), (20.0, 4.5)):
            ah = ph.undissociated_outside(total, ph_ex, PKA_ACETIC)
            bands = [ph.ph_states(w, ah_out_mM=ah)["pH_c"].tau_h
                     for w in (PLATE_READ_4H, FEDBATCH_5D)]
            for plate_end, vessel_end in zip(*bands):
                assert plate_end == pytest.approx(vessel_end, rel=0.01)

    def test_the_quoted_tau_bands_are_the_working_dose_and_are_current(self):
        # RETRACTED: these read (2.9, 6.4) h and (10.3, 15.4) h before the repair.
        ah = ph.undissociated_outside(40.0, 4.0, PKA_ACETIC)
        plate = ph.ph_states(PLATE_READ_4H, ah_out_mM=ah)["pH_c"].tau_h
        vessel = ph.ph_states(FEDBATCH_5D, ah_out_mM=ah)["pH_c"].tau_h
        assert (round(plate[0], 4), round(plate[1], 4)) == (0.0017, 0.0029)
        assert (round(vessel[0], 4), round(vessel[1], 4)) == (0.0017, 0.0029), (
            "the module docstring's chain diagram quotes this band and names this dose; tau "
            "is a property of the operating point, so both must be stated together")

    def test_the_hand_derived_manifold_rate_matches_the_finite_difference_jacobian(self):
        """The one place this module trades a numerical Jacobian for hand algebra.

        `_manifold_mode_tau_h` substitutes dpH/dA_i = -f_d/D into dA_i/dt instead of taking
        an eigenvalue. On the charge-balance manifold that must BE the fast eigenvalue of
        `_jacobian`, and the other one is the -mu decay of an imbalance the physical
        trajectory never carries.
        """
        for beta in (90.0, 145.0, 200.0):
            for mu in (0.101, 0.245, 0.363):
                ah = ph.undissociated_outside(40.0, 4.0, PKA_ACETIC)
                bound = ph.no_pump_ph(ah, beta, PKA_ACETIC)
                a_i = ph.trapped_total(ah, bound.ph_c, PKA_ACETIC)
                k = ph.entry_rate_constant_per_h("acetic")
                eigs = np.linalg.eigvals(
                    ph._jacobian(a_i, bound.ph_c, ah, beta, PKA_ACETIC, mu, k)).real
                fast = float(min(eigs))
                tau = ph._manifold_mode_tau_h(ah, beta, PKA_ACETIC, mu, bound.ph_c, k)
                assert 1.0 / tau == pytest.approx(-fast, rel=1e-3)
                assert float(max(eigs)) == pytest.approx(-mu, rel=0.1)

    def test_the_states_declare_how_they_would_be_constrained(self, ah_out):
        states = ph.ph_states(PLATE_READ_4H, ah_out_mM=ah_out)
        assert "pHluorin" in states["pH_c"].constrained_by
        assert "390/470" in states["pH_c"].constrained_by, (
            "criterion (b) is not satisfied by a table entry: the manifest is single-channel "
            "mCitrine 480/530 and this needs dual excitation, i.e. new optics"
        )
        assert "34477863" in states["A_i"].constrained_by

    def test_time_constants_are_computed_from_the_windows_growth_band(self, ah_out):
        plate = ph.ph_states(PLATE_READ_4H, ah_out_mM=ah_out)["A_i"]
        vessel = ph.ph_states(FEDBATCH_5D, ah_out_mM=ah_out)["A_i"]
        assert plate.tau_h != vessel.tau_h

    def test_an_impermeant_acid_gives_a_stationary_ph_not_a_divide_by_zero(self):
        # P_AH = 0 zeroes the pH_c row of the Jacobian, so its eigenvalue is exactly 0.
        lactic_ah = ph.undissociated_outside(900.0, 5.0, PKA_LACTIC)
        states = ph.ph_states(PLATE_READ_4H, ah_out_mM=lactic_ah, acid_name="lactic")
        assert states["pH_c"].tau_h == (math.inf, math.inf)
        assert "STATIONARY" in states["pH_c"].source

    def test_the_impermeant_acid_freezes_ph_in_both_windows(self):
        lactic_ah = ph.undissociated_outside(900.0, 5.0, PKA_LACTIC)
        for window in (PLATE_READ_4H, FEDBATCH_5D):
            audit = ph.reduction_audit(window, ah_out_mM=lactic_ah, acid_name="lactic")
            assert audit.verdict("pH_c") is Verdict.FREEZE, (
                "the headline in reduction language: 900 mM lactic acid at pH_ex 5.0 leaves "
                "pH_c where it started, in both vessels, with no fitted parameter")

    def test_the_impermeant_acids_pool_relaxes_by_dilution_alone(self):
        lactic_ah = ph.undissociated_outside(900.0, 5.0, PKA_LACTIC)
        state = ph.ph_states(PLATE_READ_4H, ah_out_mM=lactic_ah, acid_name="lactic")["A_i"]
        assert state.tau_h == pytest.approx(
            (1.0 / PLATE_READ_4H.growth_rate_high_per_h,
             1.0 / PLATE_READ_4H.growth_rate_low_per_h))


# --------------------------------------------------------------------------------------
# Criterion (c): the ODE reproduces a measurement nobody fitted it to
# --------------------------------------------------------------------------------------

class TestIntegration:

    @pytest.mark.parametrize("beta", [90.0, 145.0, 200.0])
    def test_the_partition_drop_is_complete_within_two_minutes_at_every_swept_buffer(
            self, beta):
        trace = ph.simulate_weak_acid(PLATE_READ_4H, total_acid_mM=40.0, ph_ex=4.0,
                                      beta_mM_per_ph=beta)
        bound = ph.no_pump_ph(ph.undissociated_outside(40.0, 4.0, PKA_ACETIC), beta,
                              PKA_ACETIC)
        values = trace.of("pH_c")
        target = float(ph.RESTING_CYTOSOLIC_PH) - 0.95 * bound.drop
        reached_h = trace.t[int(np.argmax(values <= target))]
        assert 0.0 < reached_h * 60.0 < 2.0, (
            "Ullah 2012 (PMID 23001666) reads 'the minimum pH_i reached within 2 min', "
            "sampled at 1 s. What arrives by then here is 95% of the partition drop; "
            "nothing is fitted to it")

    @pytest.mark.parametrize("beta,minutes", [(90.0, 0.355), (200.0, 0.628)])
    def test_the_quoted_ninety_five_percent_time_is_the_model_and_not_the_output_grid(
            self, beta, minutes):
        # The docstring band was 0.37-0.75 min, which is 3 and 6 steps of the default
        # 2001-point grid. Resolved, it is 0.36-0.63.
        trace = ph.simulate_weak_acid(PLATE_READ_4H, total_acid_mM=40.0, ph_ex=4.0,
                                      beta_mM_per_ph=beta, n_points=200001)
        bound = ph.no_pump_ph(ph.undissociated_outside(40.0, 4.0, PKA_ACETIC), beta,
                              PKA_ACETIC)
        values = trace.of("pH_c")
        target = float(ph.RESTING_CYTOSOLIC_PH) - 0.95 * bound.drop
        reached = trace.t[int(np.argmax(values <= target))] * 60.0
        assert reached == pytest.approx(minutes, abs=0.01)

    @pytest.mark.parametrize("beta", [90.0, 145.0, 200.0])
    def test_the_trace_plateaus_and_there_is_still_no_recovery_arm(self, beta):
        """RETRACTS "pH_c falls monotonically for the whole window and its minimum is last".

        That was the charge-balance error accumulating, not acidification. Repaired, pH_c
        reaches its plateau within about three minutes and stays there to within 1e-5 pH
        units. What does NOT change is the Ullah 2012 caveat: there is no pump in this bound,
        so the model matches his 2 min READ and cannot reproduce the recovery that follows.
        """
        trace = ph.simulate_weak_acid(PLATE_READ_4H, total_acid_mM=40.0, ph_ex=4.0,
                                      beta_mM_per_ph=beta)
        values = trace.of("pH_c")
        after = values[trace.t >= 3.0 / 60.0]
        assert np.max(after) - np.min(after) < 1e-5
        assert float(ph.RESTING_CYTOSOLIC_PH) - values[-1] > 1.4, (
            "it acidifies and stays acidified: no Pma1 here, so no recovery arm")

    def test_the_integrated_two_minute_ph_matches_the_algebraic_bound(self):
        trace = ph.simulate_weak_acid(PLATE_READ_4H, total_acid_mM=40.0, ph_ex=4.0,
                                      beta_mM_per_ph=145.0)
        bound = ph.no_pump_ph(ph.undissociated_outside(40.0, 4.0, PKA_ACETIC), 145.0,
                              PKA_ACETIC)
        at_two_min = float(np.interp(2.0 / 60.0, trace.t, trace.of("pH_c")))
        assert at_two_min == pytest.approx(bound.ph_c, abs=0.01)

    def test_the_no_pump_system_does_reach_a_steady_state_while_the_cell_grows(self):
        """RETRACTS "with mu > 0 the anion is diluted away and every replacement delivers
        another proton, so there is no steady state". The protons are diluted with the anion
        -- that is the repair -- so the system settles."""
        trace = ph.simulate_weak_acid(PLATE_READ_4H, total_acid_mM=40.0, ph_ex=4.0,
                                      beta_mM_per_ph=145.0)
        values = trace.of("pH_c")
        assert abs(values[-1] - values[len(values) // 2]) < 1e-6

    def test_the_stiff_driver_is_used_and_the_run_is_cheap(self):
        trace = ph.simulate_weak_acid(PLATE_READ_4H, total_acid_mM=40.0, ph_ex=4.0,
                                      beta_mM_per_ph=145.0)
        assert trace.method in ("BDF", "LSODA")
        assert trace.n_rhs_evals < 5000

    def test_an_impermeant_acid_never_enters(self):
        trace = ph.simulate_weak_acid(PLATE_READ_4H, total_acid_mM=900.0, ph_ex=5.0,
                                      beta_mM_per_ph=145.0, acid_name="lactic")
        assert np.allclose(trace.of("A_i"), 0.0)
        assert np.allclose(trace.of("pH_c"), float(ph.RESTING_CYTOSOLIC_PH))

    @pytest.mark.parametrize("beta", [90.0, 200.0])
    @pytest.mark.parametrize("mu", [0.0, 0.101, 0.245, 0.363])
    @pytest.mark.parametrize("window", [PLATE_READ_4H, FEDBATCH_5D])
    def test_the_growing_case_conserves_charge_at_every_growth_rate(self, beta, mu, window):
        """THE REPAIR, PINNED. Replaces a test that asserted the imbalance was 1.49x at
        mu = 0.101 and 2.45x at mu = 0.245 and told the next reader to delete it if the
        ratio ever became 1.0. It is 1.0.

        Only the neutral species crosses, so every A- inside is balanced by one H+ on the
        buffer and ``beta (pH_0 - pH_c) = [A-]_i`` is charge balance rather than an
        equilibrium condition. `weak_acid_rhs` now dilutes the buffered protons with the
        anion, so the identity holds at every mu and in both windows -- it used to fail by
        60-65x by the end of FEDBATCH_5D.
        """
        from scipy.integrate import solve_ivp
        ah_out = ph.undissociated_outside(40.0, 4.0, PKA_ACETIC)
        p0 = float(ph.RESTING_CYTOSOLIC_PH)
        solved = solve_ivp(
            lambda t, y: ph.weak_acid_rhs(
                t, y, ah_out_mM=ah_out, k_entry_per_h=ph.entry_rate_constant_per_h("acetic"),
                pka=PKA_ACETIC, beta_mM_per_ph=beta, growth_rate_per_h=mu),
            (0.0, window.duration_h), [0.0, p0], method="BDF", rtol=1e-10, atol=1e-12)
        a_i, ph_c = solved.y[0, -1], solved.y[1, -1]
        anion = a_i * (1.0 - ph.undissociated_fraction(ph_c, PKA_ACETIC))
        assert beta * (p0 - ph_c) / anion == pytest.approx(1.0, abs=1e-6)
        assert ph.charge_imbalance_mM(a_i, ph_c, beta, PKA_ACETIC) == pytest.approx(
            0.0, abs=1e-5)

    @pytest.mark.parametrize("mu", [0.101, 0.245])
    def test_an_imbalance_the_trajectory_never_carries_decays_at_exactly_mu(self, mu):
        # dE/dt = -mu E identically, which is why E(0) = 0 makes pH_c algebra for good.
        ah_out = ph.undissociated_outside(40.0, 4.0, PKA_ACETIC)
        beta = 200.0
        k = ph.entry_rate_constant_per_h("acetic")
        bound = ph.no_pump_ph(ah_out, beta, PKA_ACETIC)
        a_i = ph.trapped_total(ah_out, bound.ph_c, PKA_ACETIC)
        for offset in (0.3, -0.3):
            y = (a_i, bound.ph_c + offset)
            d_a, d_ph = ph.weak_acid_rhs(0.0, y, ah_out_mM=ah_out, k_entry_per_h=k,
                                         pka=PKA_ACETIC, beta_mM_per_ph=beta,
                                         growth_rate_per_h=mu)
            e = ph.charge_imbalance_mM(y[0], y[1], beta, PKA_ACETIC)
            f_d = 1.0 - ph.undissociated_fraction(y[1], PKA_ACETIC)
            d_e = (-beta * d_ph - f_d * d_a
                   - y[0] * math.log(10.0) * f_d * (1.0 - f_d) * d_ph)
            assert d_e == pytest.approx(-mu * e, rel=1e-9)

    @pytest.mark.parametrize("beta", [90.0, 145.0, 200.0])
    def test_the_bound_is_settled_on_rather_than_crossed(self, beta):
        """RETRACTS "the trajectory crosses the bound it is bracketed by, by 0.31 pH units".

        A bound a model's own trajectory passes through is not a bound, and that was the
        visible face of the charge-balance defect. Repaired, the trace settles ON
        `no_pump_ph` and approaches it FROM ABOVE, because growth holds A_i a hair below
        partition equilibrium: J_in = mu A_i > 0 requires [AH]_i < [AH]_o.
        """
        trace = ph.simulate_weak_acid(PLATE_READ_4H, total_acid_mM=40.0, ph_ex=4.0,
                                      beta_mM_per_ph=beta)
        bound = ph.no_pump_ph(ph.undissociated_outside(40.0, 4.0, PKA_ACETIC), beta,
                              PKA_ACETIC)
        settled = trace.of("pH_c")[trace.t >= 3.0 / 60.0]
        assert np.all(settled >= bound.ph_c)
        assert np.max(settled - bound.ph_c) < 1e-3

    @pytest.mark.parametrize("mu", [0.0, 0.245])
    def test_the_right_hand_side_conserves_protons_against_anion_formed(self, mu):
        # The buffer must absorb exactly the anion that formed. The old test ran mu = 0 only,
        # which is the one growth rate at which the broken equation also passed.
        ah_out = ph.undissociated_outside(40.0, 4.0, PKA_ACETIC)
        k = ph.entry_rate_constant_per_h("acetic")
        a_i, ph_c, dt = 0.0, float(ph.RESTING_CYTOSOLIC_PH), 1e-7
        beta = 145.0
        for _ in range(200000):
            d_a, d_ph = ph.weak_acid_rhs(0.0, (a_i, ph_c), ah_out_mM=ah_out,
                                         k_entry_per_h=k, pka=PKA_ACETIC,
                                         beta_mM_per_ph=beta, growth_rate_per_h=mu)
            a_i += d_a * dt
            ph_c += d_ph * dt
        anion = a_i * (1.0 - ph.undissociated_fraction(ph_c, PKA_ACETIC))
        absorbed = beta * (float(ph.RESTING_CYTOSOLIC_PH) - ph_c)
        assert absorbed == pytest.approx(anion, rel=1e-3)


# --------------------------------------------------------------------------------------
# Criterion (a): three ablations, three floors, three assays
# --------------------------------------------------------------------------------------

class TestAblations:

    def test_the_lactate_floor_is_measured_on_the_positive_control_in_the_same_experiment(self):
        predicted = (ph.undissociated_outside(2.0, 5.0, PKA_BENZOIC)
                     / ph.undissociated_fraction(3.5, PKA_BENZOIC))
        assert ph.LACTATE_TRANSFER_FLOOR.noise_floor == pytest.approx(
            abs(predicted - 0.3) / 0.3)
        assert ph.LACTATE_TRANSFER_FLOOR.noise_floor == pytest.approx(0.0940, abs=1e-3)
        assert ph.LACTATE_TRANSFER_FLOOR.noise_floor != OBSERVED_ACTIVITY_CV, (
            "criterion (a) forbids scoring a chemostat quantity with the plate-reader CV; "
            "REVISED_BUILD_LIST.md records that borrowing as the error that sank Route C2")

    def test_lactate_ablation_clears_the_floor_with_no_parameter_on_either_side(self):
        result = ph.lactate_ablation()
        assert result.full_parameter_free and result.reduced_parameter_free
        assert result.full_n_free == 0 and result.reduced_n_free == 0
        assert result.effect == pytest.approx(9.305, rel=1e-3)
        assert result.ratio == pytest.approx(99.0, rel=0.02)
        assert result.clears_floor

    def test_the_reduced_model_reproduces_abbotts_own_printed_prediction(self):
        reduced = ph._HalfYieldTransfer(name="hh", use_permeability=False)
        assert reduced.transfer("lactic", 5.0, 900.0, 3.5) == pytest.approx(85.0, rel=0.03), (
            "Abbott wrote 85 mM; the reduced model is the paper's own null hypothesis")

    def test_the_full_model_leaves_the_permeant_positive_control_untouched(self):
        full = ph._HalfYieldTransfer(name="full", use_permeability=True)
        reduced = ph._HalfYieldTransfer(name="hh", use_permeability=False)
        assert (full.transfer("benzoic", 5.0, 2.0, 3.5)
                == pytest.approx(reduced.transfer("benzoic", 5.0, 2.0, 3.5))), (
            "a piece that changes the answer where it should not is not an ablation, it is "
            "a different model")

    def test_the_full_model_is_better_but_not_exonerated(self):
        full = ph._HalfYieldTransfer(name="full", use_permeability=True)
        reduced = ph._HalfYieldTransfer(name="hh", use_permeability=False)
        full_error = abs(full.transfer("lactic", 5.0, 900.0, 3.5) - 750.0) / 750.0
        reduced_error = abs(reduced.transfer("lactic", 5.0, 900.0, 3.5) - 750.0) / 750.0
        assert full_error == pytest.approx(0.20, abs=0.01)
        assert reduced_error == pytest.approx(0.884, abs=0.01)
        assert full_error > ph.LACTATE_TRANSFER_FLOOR.noise_floor, (
            "the defensible claim is the NEGATIVE -- Henderson-Hasselbalch must fail for "
            "lactate and hold for benzoate -- not that the fallback predicts 750 mM"
        )

    @pytest.mark.parametrize("volume", list(CYTOSOLIC_VOLUMES_ML_PER_GDCW))
    @pytest.mark.parametrize("mu_ref", [0.245, 0.363])
    def test_growth_separability_clears_at_every_corner_of_both_swept_axes(self, volume,
                                                                          mu_ref):
        result = ph.growth_separability_ablation(cytosolic_volume_ml_per_gdcw=volume,
                                                 reference_growth_per_h=mu_ref)
        assert result.clears_floor
        assert result.reduced_value == 0.0, (
            "CultureContext.ph_medium reaches nothing today, so the incumbent predicts the "
            "same growth rate at every medium pH. That is the null this ablation is against")

    def test_growth_separability_fails_on_the_other_oxygen_regime(self):
        """The pass is conditional on a regime a 96-well plate does not measure, and the
        regime is now a DECLARED argument rather than an axis a monkeypatch had to reach.

        REVISED_BUILD_LIST.md E2: "for a 96-well plate the oxygen regime is UNKNOWN, so
        s_ATP(regime) is a refusal there until a kLa is measured, not a lookup".
        """
        respiratory = MAINTENANCE_REGIMES["respiratory_glucose_limited"]
        result = ph.growth_separability_ablation(regime=respiratory)
        assert not result.clears_floor
        assert result.ratio == pytest.approx(0.548394, abs=0.01)
        assert result.effect == pytest.approx(0.006416207351594316, rel=1e-3)
        old_reference = ph.growth_separability_ablation(regime=respiratory, reference_growth_per_h=0.245)
        assert old_reference.effect == pytest.approx(0.006383735282193798, rel=1e-3)
        assert ph.growth_separability_ablation(
            regime=respiratory,
            cytosolic_volume_ml_per_gdcw=max(CYTOSOLIC_VOLUMES_ML_PER_GDCW)
        ).ratio == pytest.approx(1.383, abs=0.01)

    def test_the_declared_default_regime_is_the_one_every_quoted_ratio_was_computed_at(self):
        assert ph.growth_separability_ablation(regime=REFERENCE_REGIME).ratio == (
            pytest.approx(ph.growth_separability_ablation().ratio))
        assert REFERENCE_REGIME.ngam_growth_slope_per_h_per_mmol_atp == pytest.approx(
            0.010778, rel=1e-3)

    def test_growth_separability_default_is_the_smallest_corner_of_the_sweep(self):
        default = ph.growth_separability_ablation()
        biggest = ph.growth_separability_ablation(
            cytosolic_volume_ml_per_gdcw=max(CYTOSOLIC_VOLUMES_ML_PER_GDCW),
            reference_growth_per_h=0.363)
        assert default.ratio == pytest.approx(1.21, abs=0.02)
        assert biggest.ratio > default.ratio
        assert default.floor.noise_floor == MEASURED_GROWTH_RATE_SE

    @pytest.mark.parametrize("beta", [90.0, 145.0, 200.0])
    def test_reporter_quench_clears_at_every_swept_buffer(self, beta):
        result = ph.reporter_quench_ablation(beta_mM_per_ph=beta)
        assert result.clears_floor
        assert result.floor.noise_floor == OBSERVED_ACTIVITY_CV

    def test_reporter_quench_default_is_the_smallest_corner(self):
        assert ph.reporter_quench_ablation().ratio == pytest.approx(3.38, abs=0.02)

    def test_the_state_ablation_is_retracted_and_fails_at_every_swept_buffer(self):
        """RETRACTED. It reported 2.34x the reporter floor at beta 200 and 2.86x at 90, and
        the block called it criterion (a)'s answer for the ODE. It was scoring
        `weak_acid_rhs`'s charge-balance error. Both sides are parameter-free, so no fit can
        be blamed for the collapse."""
        for beta, expected in ((90.0, 0.004034), (145.0, 0.004676), (200.0, 0.005111)):
            result = ph.reporter_state_ablation(beta_mM_per_ph=beta)
            assert result.full_n_free == 0 and result.reduced_n_free == 0
            assert result.ratio == pytest.approx(expected, rel=0.02)
            assert not result.clears_floor

    def test_no_ablation_in_this_block_removes_a_state_and_clears_its_floor(self):
        # The consequence of the retraction, stated as the assertion the block now makes.
        assert not ph.reporter_state_ablation().clears_floor
        for run in (ph.lactate_ablation, ph.growth_separability_ablation,
                    ph.reporter_quench_ablation):
            assert run().clears_floor

    def test_the_ablations_are_deterministic(self):
        for run in (ph.lactate_ablation, ph.growth_separability_ablation,
                    ph.reporter_quench_ablation):
            assert run().ratio == pytest.approx(run().ratio)


# --------------------------------------------------------------------------------------
# The published rows this block is scored on, read from the paper rather than a summary
# --------------------------------------------------------------------------------------

class TestWhatTheAblationsActuallyAblate:
    """The adversarial pass's own findings, locked so they cannot be quietly undone."""

    def test_the_growth_branch_is_exactly_an_n_equals_one_hill_in_dose(self):
        # Not a slur on the branch, a description of it: what the mechanism buys is that K
        # is computed from a measured pKa, not that the shape is anything but saturating.
        from ystwin.bridge.maintenance_calibration import REFERENCE_REGIME
        model = ph._SeparabilityGrowth(partition=True, cytosolic_volume_ml_per_gdcw=1.0,
                                       reference_growth_per_h=0.245)
        slope = REFERENCE_REGIME.ngam_growth_slope_per_h_per_mmol_atp
        for ph_ex in (4.0, 4.5, 5.0, 5.5):
            k_half = 1.0 / (slope * (1.0 / 1000.0)
                            * 10.0 ** (float(ph.RESTING_CYTOSOLIC_PH) - PKA_ACETIC)
                            * ph.undissociated_fraction(ph_ex, PKA_ACETIC))
            for dose in (0.0, 20.0, 40.0, 160.0, 1000.0):
                assert model.mu(dose, ph_ex) == pytest.approx(0.245 / (1.0 + dose / k_half),
                                                              rel=1e-12)

    def test_the_computed_ec50_moves_with_medium_ph_which_is_the_actual_claim(self):
        from ystwin.bridge.maintenance_calibration import REFERENCE_REGIME
        slope = REFERENCE_REGIME.ngam_growth_slope_per_h_per_mmol_atp

        def k_half(ph_ex):
            return 1.0 / (slope * 1e-3
                          * 10.0 ** (float(ph.RESTING_CYTOSOLIC_PH) - PKA_ACETIC)
                          * ph.undissociated_fraction(ph_ex, PKA_ACETIC))

        assert k_half(4.0) == pytest.approx(521.3, rel=0.01)
        assert k_half(4.5) == pytest.approx(688.1, rel=0.01)
        assert k_half(5.0) == pytest.approx(1215.8, rel=0.01)
        assert k_half(5.5) == pytest.approx(2884.5, rel=0.01)
        assert k_half(5.5) / k_half(4.0) == pytest.approx(5.534, rel=0.01), (
            "a single EC50 on total dose cannot do this; that is the block's claim, and it "
            "is a claim about a constant rather than about a state")

    def test_no_constant_scoring_ablation_integrates_a_state(self, monkeypatch):
        # The three original ablations are instantaneous algebra end to end.
        def refuse(*args, **kwargs):
            raise AssertionError("an ablation reached the ODE")

        monkeypatch.setattr(ph, "simulate_weak_acid", refuse)
        monkeypatch.setattr(ph, "weak_acid_rhs", refuse)
        for run in (ph.lactate_ablation, ph.growth_separability_ablation,
                    ph.reporter_quench_ablation):
            assert run().ratio > 1.0

    @pytest.mark.parametrize("ph_ex,worst", [(4.0, 0.0350), (5.5, 0.0104)])
    def test_the_quench_branch_is_absorbed_by_a_fitted_hill_at_one_medium_ph(self, ph_ex,
                                                                            worst):
        """B3 scored against a FITTED incumbent rather than one frozen at 1.0.

        `reporter_quench_ablation` reports 3.38x the floor against a reduced model that
        predicts no quench at any dose. Every plate in this repository is at ONE medium pH,
        and there the whole branch is a two-parameter Hill in total dose to under 0.04x the
        floor over a 300-fold range -- so the incumbent's own fitted EC50 absorbs it and the
        artefact is not identifiable. It separates only on the medium-pH-crossed plate that
        does not exist yet, where the same comparison gives 1.37x.
        """
        from scipy.optimize import curve_fit
        doses = np.logspace(-1.0, 2.5, 40)
        quench = np.array([
            ph.citrine_quench(ph.no_pump_ph(
                ph.undissociated_outside(d, ph_ex, PKA_ACETIC), 200.0, PKA_ACETIC).ph_c)
            for d in doses])

        def hill(dose, ec50, n):
            return 1.0 / (1.0 + (dose / ec50) ** n)

        fitted, _ = curve_fit(hill, doses, quench, p0=[40.0, 1.0], maxfev=200000)
        residual = float(np.max(np.abs(hill(doses, *fitted) - quench)))
        assert residual / OBSERVED_ACTIVITY_CV == pytest.approx(worst, abs=0.002)

    def test_the_state_ablation_does_integrate_and_now_fails_the_reporter_floor(self):
        """RETRACTED, with the numbers on both sides of the retraction.

        It read full 0.33312 against reduced 0.50642 -- a 17 percentage-point gap on a
        channel whose plate-to-plate CV is 14.6% -- and scored 2.344x the floor. The gap was
        the charge-balance error: the integrated arm was drifting off the manifold while the
        algebraic arm sat on it. Repaired, the two arms agree to 4e-4 of the signal.
        """
        result = ph.reporter_state_ablation()
        assert result.full_n_free == 0 and result.reduced_n_free == 0
        assert result.full_value == pytest.approx(0.50680, abs=1e-4)
        assert result.reduced_value == pytest.approx(0.50642, abs=1e-4)
        assert result.ratio == pytest.approx(0.005111, rel=0.02)
        assert not result.clears_floor

    def test_slaving_ph_to_the_anion_exactly_reproduces_the_integrated_trace(self):
        """What is left of the 2.34x is not a state and not solver noise.

        `reporter_state_ablation`'s reduced arm uses `no_pump_ph`, which puts A_i at
        partition equilibrium; growth holds it just below. Slave pH_c to the A_i the ODE
        actually reaches and the two arms agree to integrator tolerance, at 4.4e-12 on the
        endpoint against the 3.78e-4 the ablation still scores.
        """
        from scipy.integrate import solve_ivp
        from scipy.optimize import brentq
        beta, p0 = 200.0, float(ph.RESTING_CYTOSOLIC_PH)
        ah = ph.undissociated_outside(40.0, 4.0, PKA_ACETIC)
        k = ph.entry_rate_constant_per_h("acetic")
        grid = np.linspace(0.0, PLATE_READ_4H.duration_h, 2001)
        solved = solve_ivp(
            lambda t, y: ph.weak_acid_rhs(t, y, ah_out_mM=ah, k_entry_per_h=k,
                                          pka=PKA_ACETIC, beta_mM_per_ph=beta,
                                          growth_rate_per_h=PLATE_READ_4H.growth_rate_low_per_h),
            (0.0, PLATE_READ_4H.duration_h), [0.0, p0], t_eval=grid, method="BDF",
            rtol=1e-12, atol=1e-14)

        def slaved(a_i):
            if a_i <= 0.0:
                return p0
            return brentq(lambda q: ph.charge_imbalance_mM(a_i, q, beta, PKA_ACETIC),
                          1.0, p0 - 1e-15, xtol=1e-15)

        integrated = np.array([ph.citrine_quench(q) for q in solved.y[1]])
        algebraic = np.array([ph.citrine_quench(slaved(a)) for a in solved.y[0]])
        assert abs(integrated[-1] - algebraic[-1]) < 1e-9
        assert abs(np.trapezoid(integrated, grid)
                   - np.trapezoid(algebraic, grid)) < 1e-8
        by_bound = ph.citrine_quench(ph.no_pump_ph(ah, beta, PKA_ACETIC).ph_c)
        assert abs(integrated[-1] - by_bound) == pytest.approx(3.78e-4, rel=0.02), (
            "the dilution offset the ablation still scores, and it is tolerance-independent")

    def test_the_reference_growth_rate_is_imported_from_the_window_not_retyped(self):
        # It used to be a retyped (0.245, 0.363) in the signature. Same failure shape as the
        # tau bands: a number that is right today and silently wrong at another window.
        assert ph.growth_separability_ablation(window=PLATE_READ_4H).effect == pytest.approx(
            ph.growth_separability_ablation(
                reference_growth_per_h=PLATE_READ_4H.growth_rate_low_per_h).effect)
        assert ph.growth_separability_ablation(window=FEDBATCH_5D).effect != pytest.approx(
            ph.growth_separability_ablation(window=PLATE_READ_4H).effect)

    def test_the_full_lactate_arm_is_a_choice_among_three_not_a_deduction(self):
        # P_AH = 0 rules the neutral species out; it does not rule the total in. The rival
        # "the anion is the agent" is worse on Abbott's own rows, and that is why, not P_AH.
        pka = float(ph.PKA_LACTIC)
        anion_at_5 = 900.0 * (1.0 - ph.undissociated_fraction(5.0, pka))
        anion_model = anion_at_5 / (1.0 - ph.undissociated_fraction(3.5, pka))
        assert anion_model == pytest.approx(2761.7, rel=1e-3)
        assert abs(anion_model - 750.0) / 750.0 > abs(900.0 - 750.0) / 750.0

    def test_the_growth_ablation_is_scored_at_an_unphysical_anion_pool(self):
        # Not a bug: holding pH_c at Orij's resting value under this much acid is exactly
        # what the module says is impossible. It must travel with the 1.21x, not behind it.
        acidic, mild = ph.growth_separability_physicality()
        assert acidic.trapped_anion_mM == pytest.approx(7119.9, rel=1e-3)
        assert mild.trapped_anion_mM == pytest.approx(1286.6, rel=1e-3)
        assert not acidic.physical and not mild.physical
        assert ph.growth_separability_ablation().clears_floor

    def test_the_state_ablation_scores_exactly_the_reduction_the_audit_licenses(self):
        # The audit and the ablation now agree instead of contradicting each other: the
        # verdict is eliminate, and eliminating costs 0.0051x the measured floor.
        ah = ph.undissociated_outside(40.0, 4.0, PKA_ACETIC)
        assert ph.reduction_audit(
            PLATE_READ_4H, ah_out_mM=ah).verdict("pH_c").value == "eliminate"
        assert not ph.reporter_state_ablation().clears_floor

    def test_the_quench_branch_pooled_across_medium_ph_is_the_number_that_travels(self):
        """`reporter_quench_ablation`'s 3.38x is against a reduced model frozen at 1.0, which
        is not the incumbent. Against a fitted two-parameter Hill in total dose it is
        0.0104-0.0350x at any single medium pH and 1.33x pooled across four."""
        worst, best, pooled = ph.quench_identifiability()
        assert worst == pytest.approx(0.0104, abs=0.002)
        assert best == pytest.approx(0.0350, abs=0.002)
        assert pooled == pytest.approx(1.330, abs=0.02)
        assert ph.quench_identifiability(90.0)[2] == pytest.approx(1.334, abs=0.02)


class TestAbbottTable:

    def test_the_table_carries_the_rows_the_paper_prints(self):
        table = ph.abbott_half_yield_table()
        assert list(table.columns) == list(ph.ABBOTT_HALF_YIELD_COLUMNS)
        lactic = table[table["acid"] == "lactic"]
        assert sorted(lactic["total_mM"]) == [500.0, 750.0, 900.0]
        benzoic = table[table["acid"] == "benzoic"]
        assert sorted(benzoic["total_mM"]) == [0.3, 2.0]

    def test_the_nine_fold_failure_is_reproduced_from_the_pka_alone(self):
        table = ph.abbott_half_yield_table().set_index(["acid", "ph_ex"])
        anchor = table.loc[("lactic", 5.0), "undissociated_mM"]
        predicted_total = anchor / ph.undissociated_fraction(3.5, PKA_LACTIC)
        observed_total = table.loc[("lactic", 3.5), "total_mM"]
        assert observed_total / predicted_total == pytest.approx(8.59, abs=0.05), (
            "the paper says 'almost 9-fold higher (750 mM)'")

    def test_benzoate_obeys_henderson_hasselbalch_to_within_ten_percent(self):
        table = ph.abbott_half_yield_table().set_index(["acid", "ph_ex"])
        high = table.loc[("benzoic", 5.0), "undissociated_mM"]
        low = table.loc[("benzoic", 3.5), "undissociated_mM"]
        assert abs(low - high) / high < 0.10


# --------------------------------------------------------------------------------------
# B5: one pH_c reaching eleven reactions, with the signs the measurement actually gives
# --------------------------------------------------------------------------------------

class TestEnzymeBranch:

    def test_the_branch_carries_the_eleven_reactions_luzia_actually_assayed(self):
        # PMC9790636's own denominator is ten: "four out of ten enzymes exhibiting a
        # decrease in activity above 60%", with GAPDH counted once but assayed both ways.
        assert len(ph.LUZIA_ENZYMES) == 11
        assert {"HXK", "GAPDH_fwd", "GAPDH_rev", "PDC"} <= set(ph.LUZIA_ENZYMES)

    def test_the_coupling_enzymes_are_not_carried_as_assayed_ones(self):
        # PGK, ADH and G3PDH appear in the paper's Methods as auxiliary enzymes added in
        # excess, never as enzymes whose own V_max was read against pH.
        for coupling in ("PGK", "ADH", "G3PDH"):
            assert coupling not in ph.LUZIA_ENZYMES, (
                f"{coupling} is a coupling enzyme in Luzia 2022 (PMID 35429225), not one "
                "whose V_max was assayed against pH")
            assert coupling not in ph.ENZYME_SIGN_ABOVE_REFERENCE

    def test_every_sign_and_every_interval_names_an_assayed_reaction(self):
        assert set(ph.ENZYME_SIGN_ABOVE_REFERENCE) <= set(ph.LUZIA_ENZYMES)
        assert {e for e, _, _ in ph.ENZYME_VMAX_CHANGES} <= set(ph.LUZIA_ENZYMES)

    def test_the_two_directions_of_gapdh_move_opposite_ways_over_one_excursion(self):
        forward = ph.enzyme_vmax_change("GAPDH_fwd", 7.1, 6.4)
        reverse = ph.enzyme_vmax_change("GAPDH_rev", 7.1, 6.4)
        assert forward == pytest.approx(-0.83)
        assert reverse == pytest.approx(+0.40)
        assert forward * reverse < 0.0, (
            "one enzyme, one measured pH excursion, opposite signs: no scalar 'acid inhibits "
            "metabolism' multiplier can represent this table, which is the whole reason the "
            "branch is a table and not a number")

    def test_the_enzymes_do_not_share_a_sign_above_the_reference_ph(self):
        signs = set(ph.ENZYME_SIGN_ABOVE_REFERENCE.values())
        assert signs == {+1, -1}
        assert ph.ENZYME_SIGN_ABOVE_REFERENCE["GAPDH_fwd"] == +1
        assert ph.ENZYME_SIGN_ABOVE_REFERENCE["ALD"] == -1

    def test_the_alkaline_interval_is_the_verified_one_and_reproduces_six_fold(self):
        change = ph.enzyme_vmax_change("GAPDH_fwd", 7.8, 6.8)
        assert 1.0 / (1.0 + change) == pytest.approx(6.0), (
            "the harvest recorded this 6-fold as spanning 'the full pH range'; the "
            "verification pass read PMC9790636 and found the paper says 6.8 to 7.8"
        )
        assert ("6.8", "7.8") != ("6.19", "7.9")
        assert ("GAPDH_fwd", 6.19, 7.9) not in ph.ENZYME_VMAX_CHANGES

    def test_an_unmeasured_interval_is_refused_rather_than_interpolated(self):
        with pytest.raises(ph.EnzymeIntervalUnmeasured) as excinfo:
            ph.enzyme_vmax_change("PYK", 6.8, 6.2)
        assert "pH_kinetics" in str(excinfo.value)

    def test_a_ph_outside_the_assayed_range_is_refused_by_name(self):
        with pytest.raises(ph.EnzymeIntervalUnmeasured) as excinfo:
            ph.enzyme_vmax_change("GAPDH_fwd", 7.08, 5.6)
        assert "outside" in str(excinfo.value) and "35429225" in str(excinfo.value)
        assert not ph.enzyme_branch_in_range(5.6)
        assert ph.enzyme_branch_in_range(float(ph.RESTING_CYTOSOLIC_PH))

    def test_the_branch_costs_no_free_scalar(self):
        assert ph.PH_PARAMS.gate().passes
        assert {p.name for p in ph.PH_PARAMS.free_scalars()} == {
            "beta_cytosolic", "cytosolic_volume"}
        for name in ("vmax_gapdh_fwd_7p1_to_6p4", "vmax_gapdh_rev_7p1_to_6p4",
                     "luzia_assay_ph_low", "luzia_assay_ph_high"):
            assert ph.PH_PARAMS[name].tag == Tag.MEASURED

    def test_the_enzyme_assay_floor_is_refused_not_borrowed_from_the_plate(self):
        refusal = ph.NOT_BUILT["enzyme_vmax_assay_floor"]
        assert refusal.tag == Tag.REFUSED
        assert "OBSERVED_ACTIVITY_CV" in refusal.reason, (
            "criterion (a) says NAME ITS ASSAY and its floor; this block names the assay and "
            "refuses the floor rather than reaching for the plate CV"
        )
        with pytest.raises(RefusedValue):
            float(refusal)

    def test_the_continuous_curve_is_refused_and_names_the_dataset_that_closes_it(self):
        refusal = ph.NOT_BUILT["enzyme_vmax_curve"]
        assert refusal.tag == Tag.REFUSED
        assert "DavidLaoM/pH_kinetics" in refusal.missing

    def test_the_branch_is_unusable_at_both_bounds_at_the_working_dose(self):
        verdict = ph.enzyme_branch_verdict(40.0, 4.0)
        assert not verdict.no_pump_in_range, (
            "the no-pump bound drives pH_c to 5.41-5.68, below Luzia's 6.19")
        assert verdict.pump_holds_in_range
        assert verdict.pump_holds_ph_c == float(ph.RESTING_CYTOSOLIC_PH), (
            "at the pumped bound pH_c does not move, so every enzyme factor is exactly 1")

    def test_the_impermeant_acid_has_no_partition_bound_and_no_dose_ceiling(self):
        assert ph.enzyme_branch_dose_ceiling_mM("lactic") == (math.inf, math.inf)
        verdict = ph.enzyme_branch_verdict(900.0, 5.0, acid_name="lactic")
        assert verdict.no_pump_ph_c == (float(ph.RESTING_CYTOSOLIC_PH),) * 2, (
            "P_AH = 0 means the neutral species never equilibrates, so applying partition "
            "equilibrium to lactate would quietly undo the headline")
        assert verdict.no_pump_in_range

    def test_the_dose_ceiling_is_computed_and_sits_far_below_abbotts_half_yield(self):
        low, high = ph.enzyme_branch_dose_ceiling_mM()
        assert (low, high) == pytest.approx((2.976, 6.613), abs=0.005)
        abbott = ph.undissociated_outside(105.0, 5.0, PKA_ACETIC)
        assert abbott / high > 5.0
        # Below the ceiling the branch IS evaluable, so the refusal is dose-scoped.
        assert ph.enzyme_branch_verdict(5.0, 5.0).no_pump_in_range
