"""The burden block: a declared genotype axis, and what it is and is not allowed to claim.

Six things are pinned here and each is a number this file computes rather than repeats.

1. **The two corrections the build list mandates.** The Metzl-Raz ribosome slope is 0.35 per
   GENERATION per hour, so against mu it is 0.35/ln(2) = 0.504943 -- imported from
   `fba/stress_proteostasis.py` and never retyped, which this file checks by scanning the
   source. And ``b`` is SWEPT over 0.43-0.53, computed from the two rows it is a ratio of, so
   ``float(b)`` raises and 0.72 and 0.53-as-a-point are both outside the band.

2. **The identity that makes the sweep a reparameterisation.** ``b * phi_per_copy`` is Kafri's
   measured 0.01 per copy along the PAIRED diagonal of the two bands, so the growth channel
   does not move when you slide along it. Pinned apart the two corners give 0.00826 and
   0.01211, which is why the one function that used both now derives ``b`` from the ``phi``
   it was handed. 0.72 breaks the identity by 37%, which is the whole reason it is refused.

3. **Criterion (a) on real wells, and its condition.** The burden term is deleted and the two
   models are refit on the same 183 gate-1 wells at the same one free scalar each. The effect
   is EXACTLY ZERO at the panel's own single-copy genotype and first clears the measured
   0.0117 /h growth floor at 5 declared copies, reaching 1.85x at 8. A block that only matters
   above 5 copies is a finding, not a failure, and it is stated as a number here. What the
   wells do NOT do is also pinned: the copy column is constant across them, so the two fits
   are one fit reparameterised and no well discriminates the models. The effect is the closed
   form, and a test asserts it to machine precision rather than letting "183 real wells" carry
   an implication the arithmetic does not support.

4. **Criterion (d) is window-dependent and the two verdicts differ.** The one state has
   tau = 1/(mu + k_dH), which is 0.664-0.981 of a 4.14 h plate read -- KEEP -- and 0.0825 of a
   120 h fed-batch -- ELIMINATE. Audited at k_dH = 0 because that is the slowest case.

5. **E7 is a consistency check and E8 does not reproduce.** Growth is invariant to a 20-fold
   destabilisation at zero parameter cost, because the tax is on the synthesis fraction --
   but that structure was chosen from the same figure, so it is not held out and the module
   now says so. And no route to the E8 lag number gives the +14% the architecture bills, the
   five candidates being 17.4, 21.8, 24.2, 30.9 and 17.6 per cent.

6. **The citations resolve.** Every registered parameter names a literature identifier, and
   two misattributions found while sourcing this block are frozen out by name.

The three gate-1 tables and the GSMM are real data; the classes that need them skip when they
are absent.
"""

from __future__ import annotations

import ast
import math
import pathlib

import numpy as np
import pytest

from ystwin.fba.stress_proteostasis import (
    EGUCHI_BURDEN_LIMIT,
    METZL_RAZ_RIBOSOME_RESERVE,
    METZL_RAZ_RIBOSOME_SLOPE,
    METZL_RAZ_RIBOSOME_SLOPE_PER_GENERATION,
    add_heterologous_protein_sink,
)
from ystwin.generator.panel_experiment import (
    MEASURED_GROWTH_RATE_SE,
    OBSERVED_ACTIVITY_CV,
)
from ystwin.generator.stress_panel import STRESSORS
from ystwin.mech import burden
from ystwin.mech.ablation import GROWTH_RATE_FLOOR, REPORTER_ACTIVITY_FLOOR
from ystwin.mech.integrate import Verdict
from ystwin.mech.params import (
    FreeScalarGateFailed,
    RefusedValue,
    SweptValue,
    Tag,
)
from ystwin.mech.state import (
    FEDBATCH_5D,
    FEDBATCH_GROWTH_PER_H,
    PLATE_GROWTH_BAND_PER_H,
    PLATE_READ_4H,
    Encoding,
)

SOURCE = pathlib.Path(burden.__file__).read_text(encoding="utf-8")

CODE = SOURCE.split('"""', 2)[2]
"""The module with its own docstring removed. The docstring RECORDS the two misattributions
by name, so a "this string is gone" check has to read the code below it."""

NUMERIC_LITERALS = tuple(
    float(node.value) for node in ast.walk(ast.parse(SOURCE))
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float))
    and not isinstance(node.value, bool))
"""Every number actually written in the module, so a "never retyped" check reads the
code rather than the prose that explains it."""

PLATE_MEDIAN_MU = 0.32482577175526683
"""Median specific growth rate over the 183 passing gate-1 wells, 1/h.

Recomputed from the committed tables by ``TestTheGrowthWellsAreTheRepositorysOwn`` rather than
trusted; carried here so the classes that do not touch real data can still state a plate mu.
"""


def plate_wells():
    """The real gate-1 rows, or a skip naming what is missing."""
    try:
        return burden.plate_growth_rates()
    except FileNotFoundError as exc:
        pytest.skip(str(exc))


# --------------------------------------------------------------------------------------
# 1. The two corrections the assignment mandates
# --------------------------------------------------------------------------------------

class TestTheRibosomeSlopeIsPerGenerationNotPerMu:
    """Correction 1. The published 0.35 is against generations/h; against mu it is 0.504943."""

    def test_the_module_uses_the_corrected_slope(self):
        assert float(burden.RIBOSOME_SLOPE) == pytest.approx(
            METZL_RAZ_RIBOSOME_SLOPE_PER_GENERATION / math.log(2.0))
        assert float(burden.RIBOSOME_SLOPE) == pytest.approx(0.504943, abs=1e-6)

    def test_the_uncorrected_slope_would_be_a_thirty_percent_error(self):
        ratio = float(burden.RIBOSOME_SLOPE) / METZL_RAZ_RIBOSOME_SLOPE_PER_GENERATION
        assert ratio == pytest.approx(1.0 / math.log(2.0))
        assert ratio == pytest.approx(1.4427, abs=1e-3)

    def test_it_is_imported_and_never_retyped(self):
        """Rule: import from the module that already carries it and its regression test.

        Checked over the parsed numeric literals rather than the raw text, so the prose that
        explains the correction in the docstring does not count as retyping it.
        """
        assert "from ..fba.stress_proteostasis import" in SOURCE
        assert "METZL_RAZ_RIBOSOME_SLOPE," in SOURCE
        for value in NUMERIC_LITERALS:
            assert value != pytest.approx(0.35, abs=1e-9)
            assert value != pytest.approx(METZL_RAZ_RIBOSOME_SLOPE, abs=1e-6)

    def test_the_slope_row_is_derived_and_names_its_import(self):
        assert burden.RIBOSOME_SLOPE.tag == Tag.DERIVED
        assert "stress_proteostasis" in burden.RIBOSOME_SLOPE.source
        assert "GENERATIONS per hour" in burden.RIBOSOME_SLOPE.source


class TestBIsSweptAndNotFixed:
    """Correction 2. 0.43-0.53, computed from its two rows, and 0.72 is a category error."""

    def test_b_has_no_single_float(self):
        with pytest.raises(SweptValue):
            float(burden.BURDEN_SLOPE)

    def test_the_band_is_the_mandated_one_and_is_computed_not_typed(self):
        low, high = burden.burden_slope_band()
        assert (round(low, 2), round(high, 2)) == (0.43, 0.53)
        kafri = float(burden.KAFRI_GROWTH_LOSS_PER_COPY)
        assert low == pytest.approx(kafri / 0.023)
        assert high == pytest.approx(kafri / 0.019)
        assert "bounds=(float(KAFRI_GROWTH_LOSS_PER_COPY) / 0.023," in SOURCE

    @pytest.mark.parametrize("refused", [0.72, 0.68, 0.53, 0.40])
    def test_the_values_outside_the_band_cannot_be_pinned(self, refused):
        """0.53 included on purpose: the band's top is 0.5263 and 0.53 is already outside it."""
        with pytest.raises(ValueError):
            burden.BURDEN_SLOPE.at(refused)

    def test_the_product_b_times_phi_is_kafris_measurement_everywhere_in_the_band(self):
        measured = float(burden.KAFRI_GROWTH_LOSS_PER_COPY)
        for b in burden.BURDEN_SLOPE.sweep_points(9):
            assert b * burden.consistent_phi_per_copy(b) == pytest.approx(measured)

    def test_the_b_band_maps_onto_the_phi_band_rather_than_somewhere_else(self):
        """The previous test is an identity by construction; this one has content."""
        low, high = burden.BURDEN_SLOPE.bounds
        assert burden.consistent_phi_per_copy(high) == pytest.approx(
            burden.PHI_PER_COPY.bounds[0])
        assert burden.consistent_phi_per_copy(low) == pytest.approx(
            burden.PHI_PER_COPY.bounds[1])

    def test_the_identity_holds_only_at_the_paired_corners_and_the_module_says_so(self):
        """b and phi are two spellings of one measurement; pinned apart they disagree."""
        measured = float(burden.KAFRI_GROWTH_LOSS_PER_COPY)
        b_low, b_high = burden.BURDEN_SLOPE.bounds
        phi_low, phi_high = burden.PHI_PER_COPY.bounds
        assert b_low * phi_high == pytest.approx(measured)
        assert b_high * phi_low == pytest.approx(measured)
        assert b_low * phi_low == pytest.approx(0.00826, abs=1e-5)
        assert b_high * phi_high == pytest.approx(0.01211, abs=1e-5)
        assert "only\n    the paired corners satisfy" in SOURCE

    def test_the_refuted_slopes_carry_the_arithmetic_that_refutes_them(self):
        frame = burden.refuted_burden_slopes()
        by_value = {round(float(v), 2): row for v, row in
                    zip(frame["value"], frame.to_dict("records"), strict=True)}
        assert by_value[0.72]["verdict"].startswith("REFUTED")
        assert by_value[0.72]["in_the_band"] is False
        assert by_value[0.68]["verdict"].startswith("EXCLUDED")
        assert bool(frame["in_the_band"].sum()) == 1

    def test_the_seventy_two_arithmetic_reproduces_from_the_two_papers_numbers(self):
        """0.18 at 18 copies over 0.25 at a strain that is 13.2 copies: a 1.37x inflation."""
        copies_at_25_percent = 0.25 / 0.019
        assert copies_at_25_percent == pytest.approx(13.158, abs=1e-3)
        inflation = 18.0 / copies_at_25_percent
        assert inflation == pytest.approx(1.368, abs=1e-3)
        assert 0.18 / 0.25 == pytest.approx(0.72)
        assert (float(burden.KAFRI_GROWTH_LOSS_PER_COPY) / 0.019) * inflation == pytest.approx(
            0.72, abs=5e-3)

    def test_seventy_two_would_inflate_kafris_measured_per_copy_cost_by_thirty_seven_percent(self):
        assert 0.72 * 0.019 == pytest.approx(0.01368)
        assert 0.72 * 0.019 / float(burden.KAFRI_GROWTH_LOSS_PER_COPY) == pytest.approx(
            1.368, abs=1e-3)


# --------------------------------------------------------------------------------------
# 2. The growth arm: one measured number, and no threshold anywhere inside the range
# --------------------------------------------------------------------------------------

class TestTheGrowthArmIsParameterFree:

    def test_the_only_constant_is_kafris_one_percent(self):
        assert float(burden.KAFRI_GROWTH_LOSS_PER_COPY) == 0.01
        assert burden.KAFRI_GROWTH_LOSS_PER_COPY.tag == Tag.MEASURED
        assert "26725116" in burden.KAFRI_GROWTH_LOSS_PER_COPY.source

    @pytest.mark.parametrize("n", [0, 1, 2.5, 8, 18, 20])
    def test_the_law_is_exactly_linear_with_no_threshold(self, n):
        assert burden.growth_loss(n) == pytest.approx(0.01 * n)

    def test_kafris_own_eighteen_copy_strain_cross_checks(self):
        assert burden.growth_loss(18) == pytest.approx(0.18)

    @pytest.mark.parametrize("bad", [-1e-9, -1.0])
    def test_a_negative_genotype_is_refused(self, bad):
        with pytest.raises(ValueError, match="non-negative"):
            burden.growth_loss(bad)

    def test_total_arrest_is_refused_rather_than_extrapolated(self):
        with pytest.raises(ValueError, match="past total"):
            burden.growth_loss(100)
        assert burden.growth_loss(99.9) == pytest.approx(0.999)

    def test_the_burden_is_an_increment_and_vanishes_at_the_reference(self):
        assert burden.burdened_growth_rate(0.32, 1, reference_copies=1) == pytest.approx(0.32)
        assert burden.burdened_growth_rate(0.32, 5, reference_copies=5) == pytest.approx(0.32)

    def test_the_increment_is_the_documented_algebra(self):
        mu_env, query, ref = 0.32, 8.0, 1.0
        loss = float(burden.KAFRI_GROWTH_LOSS_PER_COPY)
        expected = mu_env * (1 - loss * query) / (1 - loss * ref)
        assert burden.burdened_growth_rate(mu_env, query, ref) == pytest.approx(expected)

    @pytest.mark.parametrize("bad", [0.0, -0.1])
    def test_a_non_positive_environment_rate_is_refused(self, bad):
        with pytest.raises(ValueError, match="must be positive"):
            burden.burdened_growth_rate(bad, 5)

    def test_the_stress_multiplier_is_refused_and_says_what_would_close_it(self):
        with pytest.raises(RefusedValue):
            float(burden.BURDEN_STRESS_MULTIPLIER)
        assert "COLONY SIZE AFTER 48 h ON SOLID MEDIUM" in (
            burden.BURDEN_STRESS_MULTIPLIER.reason)
        assert "specific growth rate in liquid" in burden.BURDEN_STRESS_MULTIPLIER.missing

    def test_the_refusal_does_not_rest_on_the_stressor_being_off_the_panel(self):
        """Farkas's heat ladder is 30/37/40 degC and stress_panel carries heat as 'degC above
        30' with a lethal dose of 18, so +7 and +10 are inside its range. A refusal that said
        'not on the panel' would be false, and it is the observable that actually refuses."""
        heat = STRESSORS["heat"]
        assert heat.units == "degC above 30"
        assert 7.0 < heat.lethal_dose and 10.0 < heat.lethal_dose
        assert "Heat IS on this panel" in burden.BURDEN_STRESS_MULTIPLIER.reason

    def test_the_multiplier_records_all_four_farkas_measured_not_the_two_usually_quoted(self):
        source = burden.BURDEN_STRESS_MULTIPLIER.source
        for agent in ("cycloheximide", "MPA", "AZC", "HEAT"):
            assert agent in source
        assert "3.7-fold" in source and "2.8-fold at 40" in source

    def test_the_mpa_dose_is_carried_as_the_papers_own_inconsistency_not_as_a_number(self):
        """Materials and methods say 0.30 ug/ml, the Fig 2B legend says 30 ug/ml. Picking
        either would be choosing a side of a 100x spread on the reader's behalf."""
        source = burden.BURDEN_STRESS_MULTIPLIER.source
        assert "INTERNALLY INCONSISTENT" in source
        assert "0.30 ug/ml" in source and "30 ug/ml" in source
        assert "neither is usable" in source

    def test_kafris_own_slope_is_declared_rich_medium_and_not_condition_free(self):
        """The linearity generalises across media and the SLOPE does not, and Kafri says so in
        the same paragraph. Quoting the linearity without the qualifier widened the claim that
        licenses treating 0.01 as environment-independent."""
        source = burden.KAFRI_GROWTH_LOSS_PER_COPY.source
        assert "although to different extents" in source
        assert "low-phosphate" in source
        assert "possibly reflecting" in source
        assert "constant RELATIVE loss assumes away" in source
        assert "RICH MEDIUM" in source
        assert "NO PER-CONDITION SLOPE IS STATED IN " in source
        assert "Fig 4D plots" in source
        assert "0.01 per copy, MEASURED IN RICH MEDIUM." in SOURCE

    def test_the_refusal_does_not_claim_the_multiplier_is_one(self):
        """Kafri measured a medium x burden interaction on the RIGHT observable -- relative
        growth rate in liquid -- so m != 1 is measured, and only its size is missing."""
        reason = burden.BURDEN_STRESS_MULTIPLIER.reason
        assert "AND THE REFUSAL IS NOT THAT m = 1" in reason
        assert "MEDIUM x burden" in reason
        assert "states no per-condition slope in its text" in reason


# --------------------------------------------------------------------------------------
# 3. Criterion (e): the free-scalar gate, one arm passing and one failing on purpose
# --------------------------------------------------------------------------------------

class TestCriterionEIsComputedForBothArms:

    def test_the_growth_arm_passes_one_against_one(self):
        gate = burden.growth_gate()
        assert gate.passes
        assert len(gate.free) == 1 and len(gate.targets) == 1
        assert gate.targets == (GROWTH_RATE_FLOOR.name,)
        burden.GROWTH_PARAMS.require_gate()

    def test_the_growth_arms_one_free_scalar_is_a_refusal_not_a_fitted_number(self):
        free = burden.GROWTH_PARAMS.free_scalars()
        assert [p.name for p in free] == ["m_stress"]
        assert free[0].tag == Tag.REFUSED

    def test_the_allocation_arm_fails_eight_against_zero(self):
        gate = burden.allocation_gate()
        assert not gate.passes
        assert len(gate.free) == 8
        assert len(gate.targets) == 0
        with pytest.raises(FreeScalarGateFailed):
            burden.ALLOCATION_PARAMS.require_gate()

    def test_the_allocation_arm_is_worse_than_the_ph_slice_the_gate_was_calibrated_on(self):
        """ARCHITECTURE_GAPS 0.3 counted 8 against 2 there; here it is 8 against 0."""
        gate = burden.allocation_gate()
        assert len(gate.free) >= 8 and len(gate.targets) < 2

    def test_the_eight_is_an_upper_bound_and_the_verdict_survives_the_tighter_count(self):
        """b is a function of phi_per_copy and c_r of phi_max_library, so the registry's 8
        double-counts two rows and four more are REFUSED with no value at all. The honest
        irreducible count is two swept axes -- and against zero targets it still REFUSES,
        which is why the loose count was never load-bearing."""
        free = {p.name for p in burden.ALLOCATION_PARAMS.free_scalars()}
        refused = {p.name for p in burden.ALLOCATION_PARAMS.params if p.tag == Tag.REFUSED}
        derived_from_another_row = {"b", "c_r"}
        assert derived_from_another_row < free
        assert len(refused) == 4
        irreducible = free - refused - derived_from_another_row
        assert irreducible == {"phi_per_copy", "phi_max_library"}
        assert len(irreducible) > len(burden.allocation_gate().targets)

    def test_the_two_registries_are_separate_so_the_failure_cannot_hide(self):
        assert burden.GROWTH_PARAMS is not burden.ALLOCATION_PARAMS
        assert set(p.name for p in burden.GROWTH_PARAMS.params).isdisjoint(
            p.name for p in burden.ALLOCATION_PARAMS.params)


# --------------------------------------------------------------------------------------
# 4. Criterion (a): the growth ablation, on real wells, against the measured growth floor
# --------------------------------------------------------------------------------------

class TestTheGrowthWellsAreTheRepositorysOwn:

    def test_the_three_committed_plates_give_one_hundred_and_eighty_three_wells(self):
        wells = plate_wells()
        assert len(wells) == 183
        assert wells["plate"].nunique() == 3
        assert set(wells.columns) >= {"plate", "well", "growth_rate"}

    def test_the_median_matches_the_band_the_state_module_carries(self):
        wells = plate_wells()
        assert wells["growth_rate"].median() == pytest.approx(PLATE_MEDIAN_MU)
        low, high = np.percentile(wells["growth_rate"], [25, 75])
        assert (low, high) == pytest.approx(PLATE_GROWTH_BAND_PER_H, rel=1e-12)

    def test_a_missing_export_is_a_named_refusal_not_a_silent_empty_frame(self):
        assert "FileNotFoundError" in SOURCE
        assert "not scoreable without them" in SOURCE

    def test_the_read_length_is_imported_from_the_window_and_never_retyped(self):
        """4.14 h is mech/state.py's, with its own source. Retyping it here would be the
        same defect as retyping the ribosome slope."""
        assert burden._PLATE_READ_H == PLATE_READ_4H.duration_h
        for value in NUMERIC_LITERALS:
            assert value != pytest.approx(PLATE_READ_4H.duration_h, abs=1e-9)


class TestCriterionAOnTheGrowthArm:

    def test_the_effect_is_exactly_zero_at_the_panels_own_genotype(self):
        result = burden.score_growth_ablation(1.0, data=plate_wells())
        assert result.effect == pytest.approx(0.0, abs=1e-12)
        assert not result.clears_floor

    def test_it_first_clears_the_measured_floor_at_five_declared_copies(self):
        wells = plate_wells()
        assert burden.copies_to_clear_the_floor(data=wells) == 5
        assert not burden.score_growth_ablation(4, data=wells).clears_floor
        assert burden.score_growth_ablation(5, data=wells).clears_floor

    def test_at_eight_copies_it_is_one_point_eight_times_the_floor(self):
        result = burden.score_growth_ablation(8, data=plate_wells())
        assert result.effect == pytest.approx(0.02158843621636264, abs=1e-5)
        assert result.ratio == pytest.approx(1.8451654885780033, abs=5e-3)
        assert result.clears_floor

    def test_the_effect_is_linear_in_the_declared_increment(self):
        wells = plate_wells()
        effects = [burden.score_growth_ablation(n, data=wells).effect for n in (2, 3, 5, 9)]
        increments = [n - 1 for n in (2, 3, 5, 9)]
        per_copy = [e / i for e, i in zip(effects, increments, strict=True)]
        assert per_copy == pytest.approx([per_copy[0]] * 4, rel=1e-9)

    def test_neither_model_wins_by_spending_a_parameter(self):
        result = burden.score_growth_ablation(8, data=plate_wells())
        assert result.full_n_free == result.reduced_n_free == 1
        assert result.n_rows == 183

    def test_no_well_discriminates_the_two_models_and_the_module_says_so(self):
        """The copy column is constant, so the two fits are one fit reparameterised. The wells
        supply mu_bar; they do not weigh the burden term. Claiming otherwise would be the
        frozen-incumbent error's mirror image -- a comparison that cannot lose."""
        wells = plate_wells()
        observable = burden.growth_rate_observable(8, data=wells)
        full, reduced = burden.BurdenedGrowth(), burden.UnburdenedGrowth()
        fit_full = full.refit(observable.data)
        fit_reduced = reduced.refit(observable.data)
        assert fit_full.residual_rms == pytest.approx(fit_reduced.residual_rms, rel=1e-12)
        assert full.predict(observable.data, fit_full.parameters) == pytest.approx(
            reduced.predict(observable.data, fit_reduced.parameters), rel=1e-12)
        assert "no well discriminates the two models" in burden.score_growth_ablation.__doc__

    def test_the_effect_is_exactly_the_closed_form_the_docstring_gives(self):
        wells = plate_wells()
        mu_bar = float(wells["growth_rate"].mean())
        loss = float(burden.KAFRI_GROWTH_LOSS_PER_COPY)
        for n, ref in ((5, 1.0), (8, 1.0), (9, 2.0)):
            closed = mu_bar * loss * (n - ref) / (1.0 - loss * ref)
            got = burden.score_growth_ablation(n, reference_copies=ref, data=wells).effect
            assert got == pytest.approx(closed, rel=1e-12)

    def test_it_is_scored_absolutely_against_the_growth_floor_not_the_plate_cv(self):
        result = burden.score_growth_ablation(8, data=plate_wells())
        assert result.scored == "absolute"
        assert result.floor is GROWTH_RATE_FLOOR
        assert result.floor.noise_floor == MEASURED_GROWTH_RATE_SE == 0.0117
        assert result.floor is not REPORTER_ACTIVITY_FLOOR
        assert result.floor.noise_floor != OBSERVED_ACTIVITY_CV

    def test_the_report_names_the_assay_and_the_verdict(self):
        text = burden.score_growth_ablation(8, data=plate_wells()).report()
        assert "log OD600" in text
        assert "CLEARS THE FLOOR" in text

    def test_no_copy_number_clears_it_from_its_own_reference(self):
        wells = plate_wells()
        with pytest.raises(ValueError, match="no copy number up to"):
            burden.copies_to_clear_the_floor(reference_copies=1, data=wells, max_copies=4)

    def test_the_observable_declares_the_axis_it_is_evaluated_along(self):
        observable = burden.growth_rate_observable(8, data=plate_wells())
        assert observable.units == "1/h"
        assert "declared" in observable.description
        assert set(observable.data["copies"]) == {1.0}


# --------------------------------------------------------------------------------------
# 5. Criterion (b) and (d): one state, its identifiability, and the two windows
# --------------------------------------------------------------------------------------

class TestTheOneStateIsDeclaredNotInvented:

    def test_there_is_exactly_one_state_and_the_growth_arm_has_none(self):
        assert len(burden.BURDEN_STATES) == 1
        assert burden.BURDEN_STATES.names == ("product_fraction",)

    def test_it_is_a_declared_sweep_axis_because_nothing_measures_a_copy_number_here(self):
        assert burden.PRODUCT_FRACTION.sweep_axis is True
        assert burden.PRODUCT_FRACTION.constrained_by == ""
        assert "genotype" in burden.PRODUCT_FRACTION.note

    def test_it_is_a_level_and_not_an_integrating_pool(self):
        assert burden.PRODUCT_FRACTION.encoding is Encoding.LEVEL

    def test_its_time_constant_is_one_over_mu(self):
        assert burden.PRODUCT_FRACTION.tau_growth_multiple == 1.0
        low, high = burden.PRODUCT_FRACTION.taus_in(PLATE_READ_4H)
        assert low == pytest.approx(1.0 / PLATE_GROWTH_BAND_PER_H[1])
        assert high == pytest.approx(1.0 / PLATE_GROWTH_BAND_PER_H[0])


class TestCriterionDReductionIsWindowDependent:

    def test_on_a_four_hour_plate_read_it_must_be_integrated(self):
        reduced = burden.reduction_audit(PLATE_READ_4H)
        assert reduced.verdict("product_fraction") is Verdict.KEEP
        row = reduced.row("product_fraction")
        assert row.ratio_low == pytest.approx(0.664, abs=5e-3)
        assert row.ratio_high == pytest.approx(0.981, abs=5e-3)
        assert len(reduced.integrated) == 1

    def test_on_a_five_day_fed_batch_it_collapses_to_algebra(self):
        reduced = burden.reduction_audit(FEDBATCH_5D)
        assert reduced.verdict("product_fraction") is Verdict.ELIMINATE
        row = reduced.row("product_fraction")
        assert row.ratio_high == pytest.approx(1.0 / FEDBATCH_GROWTH_PER_H / 120.0, rel=1e-9)
        assert row.ratio_high == pytest.approx(0.0825, abs=5e-4)
        assert row.ratio_high < OBSERVED_ACTIVITY_CV
        assert len(reduced.eliminated) == 1

    def test_the_two_verdicts_are_different_which_is_why_the_window_is_an_input(self):
        assert burden.reduction_audit(PLATE_READ_4H).verdict(
            "product_fraction") is not burden.reduction_audit(FEDBATCH_5D).verdict(
            "product_fraction")

    def test_the_audit_is_run_at_the_slowest_case(self):
        """tau = 1/(mu + k_dH). Any product turnover only shortens it, so k_dH = 0 is the
        conservative audit -- a state that must be kept at k_dH = 0 can only become more
        eliminable, never less."""
        mu = PLATE_GROWTH_BAND_PER_H[0]
        assert 1.0 / mu > 1.0 / (mu + 0.1) > 1.0 / (mu + 1.0)
        assert "SLOWEST case" in burden.PRODUCT_FRACTION.source

    def test_the_audit_table_prints_the_bias_the_ratio_bounds(self):
        table = burden.reduction_audit(PLATE_READ_4H).audit_table()
        assert "product_fraction" in table
        assert "0.146" in table


class TestThePoolIntegratesToItsClosedForm:

    @pytest.mark.parametrize("k_dH", [0.0, 0.3, 2.0])
    def test_the_trajectory_matches_the_analytic_solution(self, k_dH):
        f_H, mu = 8 * 0.019, 0.32
        trajectory = burden.integrate_product_fraction(
            PLATE_READ_4H, f_H=f_H, mu=mu, k_dH=k_dH)
        tau = 1.0 / (mu + k_dH)
        steady = burden.steady_state_fraction(f_H, mu, k_dH)
        expected = steady * (1.0 - np.exp(-trajectory.t / tau))
        assert trajectory.of("product_fraction") == pytest.approx(expected, abs=1e-6)

    def test_a_stable_product_relaxes_to_its_synthesis_fraction(self):
        assert burden.steady_state_fraction(0.152, 0.32, 0.0) == pytest.approx(0.152)

    def test_the_fed_batch_run_is_effectively_at_steady_state(self):
        f_H = 8 * 0.019
        trajectory = burden.integrate_product_fraction(
            FEDBATCH_5D, f_H=f_H, mu=FEDBATCH_GROWTH_PER_H)
        assert trajectory.of("product_fraction")[-1] == pytest.approx(f_H, rel=1e-4)

    def test_the_plate_run_is_not(self):
        f_H = 8 * 0.019
        trajectory = burden.integrate_product_fraction(PLATE_READ_4H, f_H=f_H, mu=0.32)
        assert trajectory.of("product_fraction")[-1] == pytest.approx(0.1116, abs=1e-3)
        assert trajectory.of("product_fraction")[-1] < 0.75 * f_H

    @pytest.mark.parametrize("bad", [(0.0, 0.0), (-0.1, 0.0), (0.3, -1.0)])
    def test_the_right_hand_side_refuses_impossible_rates(self, bad):
        mu, k_dH = bad
        with pytest.raises(ValueError):
            burden.product_fraction_rhs(0.1, mu, k_dH)
        with pytest.raises(ValueError):
            burden.steady_state_fraction(0.1, mu, k_dH)


# --------------------------------------------------------------------------------------
# 6. The ceiling: a refusal at the edge, never a Hill term inside the range
# --------------------------------------------------------------------------------------

class TestTheCeilingIsRefusedNotSaturated:

    def test_the_two_ceilings_are_the_measured_ones(self):
        assert float(burden.EGUCHI_CEILING) == EGUCHI_BURDEN_LIMIT == 0.15
        assert float(burden.FUJITA_CEILING) == 0.40

    def test_the_non_toxic_ceiling_is_more_than_two_and_a_half_times_the_egfp_like_one(self):
        assert float(burden.FUJITA_CEILING) / float(burden.EGUCHI_CEILING) > 2.5

    def test_a_load_past_its_class_ceiling_raises_rather_than_bending_the_law(self):
        with pytest.raises(burden.BurdenAboveCeiling, match="above the measured ceiling"):
            burden.require_admissible(0.151, burden.ProductClass.EGFP_LIKE)
        assert burden.require_admissible(0.15, burden.ProductClass.EGFP_LIKE) == 0.15

    def test_the_same_load_is_admissible_for_a_declared_non_toxic_product(self):
        assert burden.require_admissible(0.35, burden.ProductClass.NON_TOXIC) == 0.35
        with pytest.raises(burden.BurdenAboveCeiling):
            burden.require_admissible(0.41, burden.ProductClass.NON_TOXIC)

    def test_the_class_must_be_declared_because_the_gate_is_folding_not_amount(self):
        with pytest.raises(ValueError, match="not one of"):
            burden.require_admissible(0.05, "some-other-protein")
        with pytest.raises(ValueError, match="non-negative"):
            burden.require_admissible(-0.01)

    def test_eight_copies_at_the_mass_spec_end_is_already_past_the_egfp_like_ceiling(self):
        f_H = burden.proteome_fraction(8, 0.019)
        assert f_H == pytest.approx(0.152)
        with pytest.raises(burden.BurdenAboveCeiling):
            burden.require_admissible(f_H)

    def test_no_saturating_term_was_installed_below_the_ceiling(self):
        """Kafri measured linearity with no threshold, so the law must stay a straight line."""
        assert burden.growth_loss(14) / burden.growth_loss(7) == pytest.approx(2.0)
        assert "Hill" not in SOURCE.split("BurdenAboveCeiling")[0]


# --------------------------------------------------------------------------------------
# 7. E7 for free, E8 not reproducing
# --------------------------------------------------------------------------------------

class TestE7GrowthIsInvariantToProductStability:
    """Kafri 2016 Fig 3B: a ~20-fold destabilised GFP cost the same fitness. Zero parameters."""

    def test_growth_does_not_move_when_the_product_is_destabilised(self):
        frame = burden.destabilisation_invariance(8, 0.32, 0.019)
        assert frame["growth_loss"].nunique() == 1
        assert frame["growth_loss"].iloc[0] == pytest.approx(0.08)

    def test_a_pool_taxing_model_would_get_it_wrong_by_the_turnover_ratio(self):
        mu = 0.32
        frame = burden.destabilisation_invariance(8, mu, 0.019, k_dH_multiples=(0.0, 20.0))
        stable, unstable = frame["pool_taxing_model_would_predict"]
        assert stable / unstable == pytest.approx(21.0, rel=1e-9)

    @pytest.mark.parametrize("phi", [0.019, 0.021, 0.023])
    def test_the_counterfactual_agrees_at_zero_turnover_at_every_pinned_phi(self, phi):
        """It disagreed by 17% at phi = 0.019 while b was pinned at an unpaired endpoint."""
        frame = burden.destabilisation_invariance(8, 0.32, phi, k_dH_multiples=(0.0,))
        assert frame["pool_taxing_model_would_predict"].iloc[0] == pytest.approx(
            frame["growth_loss"].iloc[0])

    def test_the_pool_itself_does_move(self):
        frame = burden.destabilisation_invariance(8, 0.32, 0.019)
        assert frame["steady_state_pool"].iloc[1] < 0.05 * frame["steady_state_pool"].iloc[0]

    def test_the_state_is_named_for_the_synthesis_flux_the_tax_is_on(self):
        assert "SYNTHESIS" in burden.proteome_fraction.__doc__

    def test_e7_is_not_claimed_as_held_out_because_the_structure_came_from_the_figure(self):
        """The tax was put on synthesis BECAUSE of Fig 3B. Reproducing Fig 3B afterwards is a
        consistency check, and E8 gets that disclosure two sections down, so E7 gets it too."""
        doc = burden.destabilisation_invariance.__doc__
        assert "IT IS NOT HELD OUT" in doc
        assert "consistency check" in doc
        assert "chosen from" in doc.lower()


class TestE8IsAFittedCheckAndDoesNotReproduce:

    def test_no_route_gives_the_architectures_fourteen_percent(self):
        frame = burden.metzl_raz_lag_check()
        assert not frame["reproduces_the_claim"].any()
        assert (frame["architecture_claimed"] == 0.14).all()

    def test_the_five_routes_are_the_numbers_the_docstring_names(self):
        frame = burden.metzl_raz_lag_check()
        got = sorted(round(float(v), 4) for v in frame["fractional_increase"])
        assert got == pytest.approx(
            [0.15, 0.15, 0.1735, 0.1765, 0.218, 0.2425, 0.3093], abs=5e-4)

    def test_only_the_measurement_is_held_out(self):
        frame = burden.metzl_raz_lag_check()
        held = frame[frame["held_out"]]
        assert len(held) == 1
        assert held["fractional_increase"].iloc[0] == pytest.approx(0.15)
        assert "MEASURED" in held["route"].iloc[0]

    def test_the_measured_target_names_recovery_time_and_not_lag_from_inoculum(self):
        """The paper's sentence is about RECOVERY time. Calling it lag would silently swap
        the observable, which is the substitution LAG_ASSAY_FLOOR exists to refuse."""
        source = SOURCE.split("LAG_MEASURED_FRACTIONAL_INCREASE_AT_8_COPIES = 0.15")[1]
        assert "recovery time of eight-copy burden cells" in source
        assert "RECOVERY TIME, not lag from a fresh inoculum" in source

    def test_the_check_is_not_held_out_because_c_r_comes_off_the_same_library(self):
        frame = burden.metzl_raz_lag_check()
        algebra = frame[frame["route"].str.startswith("algebra")]
        assert len(algebra) == 4
        assert all("same library" in note for note in algebra["note"])

    def test_the_papers_own_fifteen_percent_route_gives_seventeen_point_six(self):
        assert 1.0 / (1.0 - burden.LAG_MEASURED_FRACTIONAL_INCREASE_AT_8_COPIES) - 1.0 == (
            pytest.approx(0.17647, abs=1e-5))

    def test_the_lag_ratio_is_r0_over_the_burdened_reserve(self):
        f_H, c_r = 8 * 0.019, 0.0788
        assert burden.lag_ratio(f_H, c_r) == pytest.approx(
            float(burden.R0_UNBURDENED) / burden.ribosome_reserve(f_H, c_r))

    def test_a_reserve_driven_to_zero_is_refused_not_clipped(self):
        with pytest.raises(ValueError, match="extrapolation"):
            burden.ribosome_reserve(1.0, 0.104)

    def test_it_cannot_be_scored_under_criterion_a_because_its_assay_has_no_floor_here(self):
        with pytest.raises(RefusedValue):
            float(burden.LAG_ASSAY_FLOOR)
        assert "lag" in burden.LAG_ASSAY_FLOOR.reason
        assert "re-inoculated from stationary" in burden.LAG_ASSAY_FLOOR.missing
        assert str(OBSERVED_ACTIVITY_CV) in burden.LAG_ASSAY_FLOOR.source


# --------------------------------------------------------------------------------------
# 8. What is REPORTED, and every refusal that keeps it from feeding a prediction
# --------------------------------------------------------------------------------------

class TestTheAllocationArmFeedsNothing:

    def test_the_module_says_so_in_its_own_docstring(self):
        assert "DOES NOT CLAIM CRITERION 1" in burden.__doc__
        assert "TORC1" in burden.__doc__
        assert "Dot6/Tod6" in burden.__doc__

    def test_the_phi_band_is_the_two_calibrations_of_one_quantity(self):
        assert burden.PHI_PER_COPY.bounds == (0.019, 0.023)
        assert burden.PHI_PER_COPY.tag == Tag.SWEPT
        with pytest.raises(SweptValue):
            float(burden.PHI_PER_COPY)

    def test_the_two_r0_endpoints_are_the_only_measured_excursion(self):
        assert float(burden.R0_UNBURDENED) == 0.081
        assert float(burden.R0_HIGHEST_BURDEN) == 0.055

    def test_the_two_r0_calibrations_do_not_drift_apart(self):
        """8.1% here against the 8% fba/stress_proteostasis.py rounds to: 1.2% apart."""
        relative = abs(float(burden.R0_UNBURDENED) - METZL_RAZ_RIBOSOME_RESERVE) / float(
            burden.R0_UNBURDENED)
        assert relative == pytest.approx(0.0123, abs=1e-3)
        assert relative < 0.02

    def test_the_phi_max_band_is_the_papers_own_inconsistency(self):
        assert burden.PHI_MAX_LIBRARY.bounds == (0.25, 0.33)
        assert "not jointly consistent" in burden.PHI_MAX_LIBRARY.source
        assert "copy number of Metzl-Raz's highest-burden strain" in burden.PHI_MAX_LIBRARY.missing

    def test_the_three_statements_really_are_inconsistent(self):
        """(i)+(ii) put the top strain at 13.2 copies and force a 19.5% r0 fall at 8, against
        the paper's own ~15%; (iii)+(ii) force ~17 copies instead."""
        fall_at_8 = (0.081 - 0.055) / 0.25 * 8 * 0.019
        assert fall_at_8 / 0.081 == pytest.approx(0.195, abs=2e-3)
        top_from_iii = 8 * ((0.081 - 0.055) / 0.081) / 0.15
        assert top_from_iii == pytest.approx(17.1, abs=0.2)
        assert top_from_iii * 0.019 == pytest.approx(0.325, abs=5e-3)

    def test_c_r_is_computed_from_the_rows_above_it(self):
        low, high = burden.RIBOSOME_RESERVE_SLOPE.bounds
        assert low == pytest.approx((0.081 - 0.055) / 0.33)
        assert high == pytest.approx((0.081 - 0.055) / 0.25)
        assert (round(low, 3), round(high, 3)) == (0.079, 0.104)

    def test_the_twenty_five_percent_is_the_heterologous_fraction_not_the_ribosome_one(self):
        """Metzl-Raz states TWO different ~25% quantities and this block depends on using
        the right one. PHI_MAX_LIBRARY's is 'mCherry levels reaching ~25% of the total
        proteome in the highest burden cells' (Results, verbatim). The other -- ~25% of
        RIBOSOMAL proteins inactive, r0 ~8% over r ~30% -- is a different quantity and is
        not what this row means. Recorded because the two were confused during this build."""
        assert "the total proteome in the highest burden cells" in burden.PHI_MAX_LIBRARY.source
        assert "ribosom" not in burden.PHI_MAX_LIBRARY.source.lower()
        doc = " ".join(burden.misattributed_quotes.__doc__.split())
        assert "INACTIVE-RIBOSOME ~25%" in doc
        assert "is verbatim in the Results" in doc

    @pytest.mark.parametrize("name", ["phi_per_copy_from_paxdb", "tai_burden_gate",
                                      "b_secreted", "lag_time_noise_floor"])
    def test_every_refusal_names_a_reason_and_a_measurement(self, name):
        param = burden.ALLOCATION_PARAMS[name]
        assert param.tag == Tag.REFUSED
        with pytest.raises(RefusedValue):
            float(param)
        assert len(param.reason) > 40 and len(param.missing) > 20

    def test_the_paxdb_route_is_refused_with_the_size_of_the_error_attached(self):
        frame = burden.paxdb_route_refused()
        factors = dict(zip(frame["protein"], frame["mass_over_molar_fraction"], strict=True))
        assert factors["crtE"] == pytest.approx(0.84306, abs=1e-5)
        assert factors["crtI"] == pytest.approx(1.30182, abs=1e-5)
        assert factors["crtYB"] == pytest.approx(1.49472, abs=1e-5)
        assert frame["average_protein_kda"].eq(50.0).all()

    def test_the_molar_for_mass_substitution_is_up_to_one_and_a_half_fold(self):
        frame = burden.paxdb_route_refused()
        assert frame["mass_over_molar_fraction"].max() > 1.49
        assert frame["mass_over_molar_fraction"].min() < 0.85

    def test_farkas_is_carried_as_an_upper_bound_on_a_different_observable(self):
        assert float(burden.FARKAS_COLONY_SIZE_SLOPE) == pytest.approx(0.6757, abs=1e-4)
        assert "colony" in burden.FARKAS_COLONY_SIZE_SLOPE.name
        assert float(burden.FARKAS_COLONY_SIZE_SLOPE) > burden.burden_slope_band()[1]


# --------------------------------------------------------------------------------------
# 9. The GEM coupling: through fba/stress_proteostasis.py, and by how much it differs
# --------------------------------------------------------------------------------------

class TestTheGemCouplingIsNotDuplicated:

    def test_the_stoichiometry_is_imported_rather_than_retyped(self):
        for name in ("burden_curve", "grams_per_gdcw_at_protein_fraction",
                     "read_protein_composition", "ribosome_allocation_penalty"):
            assert name in SOURCE
        assert "add_reaction" not in SOURCE
        assert "r_4047" not in SOURCE.replace("``r_4047``", "")

    def test_the_displacement_law_is_steeper_than_the_measured_burden_everywhere(self):
        for mu in (*PLATE_GROWTH_BAND_PER_H, FEDBATCH_GROWTH_PER_H, PLATE_MEDIAN_MU):
            for phi in burden.PHI_PER_COPY.bounds:
                assert burden.displacement_over_measured(mu, phi) > 2.5

    def test_it_is_two_point_eight_to_three_point_four_at_the_plates_median_rate(self):
        got = [burden.displacement_over_measured(PLATE_MEDIAN_MU, phi)
               for phi in burden.PHI_PER_COPY.bounds]
        assert got == pytest.approx([2.8267, 3.4218], abs=1e-3)

    def test_it_is_four_point_nine_to_five_point_nine_at_the_fed_batch_setpoint(self):
        got = [burden.displacement_over_measured(FEDBATCH_GROWTH_PER_H, phi)
               for phi in burden.PHI_PER_COPY.bounds]
        assert got == pytest.approx([4.8804, 5.9079], abs=1e-3)

    def test_over_the_plates_whole_growth_band_it_widens_to_two_point_seven_to_three_point_eight(
            self):
        got = [burden.displacement_over_measured(mu, phi)
               for mu in PLATE_GROWTH_BAND_PER_H for phi in burden.PHI_PER_COPY.bounds]
        assert min(got) == pytest.approx(2.7271, abs=1e-3)
        assert max(got) == pytest.approx(3.7798, abs=1e-3)

    def test_the_ratio_is_independent_of_copy_number_because_both_laws_are_linear(self):
        one = burden.displacement_over_measured(0.32, 0.019)
        assert burden.displacement_over_measured(0.32, 0.023) / one == pytest.approx(
            0.023 / 0.019)

    def test_a_phi_outside_the_band_cannot_be_slipped_in(self):
        with pytest.raises(ValueError, match="outside it"):
            burden.displacement_over_measured(0.32, 0.037)


@pytest.mark.integration
class TestTheGemPricesTheLoadBesideTheMeasurement:
    """The one comparison that needs the real GSMM. Skips when it is absent."""

    @pytest.fixture(scope="class")
    def _sink_model_template(self, yeast_gem_factory):
        return add_heterologous_protein_sink(yeast_gem_factory())

    @pytest.fixture
    def sink_model(self, _sink_model_template, model_copy):
        return model_copy(_sink_model_template)

    POSTMA_MAX_OXYGEN_UPTAKE = -12.0
    """mmol O2/gDW/h, Postma 1989 (PMID 2566299) in CBS 8066, as `docs/FINDINGS.md` carries it.

    The fixture ran at ``r_1992 = -1000`` until this pass, which is the branch that same
    document calls "more than double the highest oxygen uptake ever measured in this
    organism" and where the model shows no Crabtree effect at all -- it grows at 0.888 /h on
    10 mmol glucose. That is the `s_ATP(regime)` defect `REVISED_BUILD_LIST.md` §1 found in the
    pH block, and it was here too.
    """

    @staticmethod
    def carbon_limited(model):
        model.reactions.get_by_id("r_1714").lower_bound = -10.0
        model.reactions.get_by_id("r_1992").lower_bound = (
            TestTheGemPricesTheLoadBesideTheMeasurement.POSTMA_MAX_OXYGEN_UPTAKE)
        return model

    @staticmethod
    def unbounded_oxygen(model):
        model.reactions.get_by_id("r_1714").lower_bound = -10.0
        model.reactions.get_by_id("r_1992").lower_bound = -1000.0
        return model

    def test_the_chemistry_alone_undercharges_against_the_measured_burden(self, sink_model):
        frame = burden.gem_burden_comparison(
            sink_model, 8, 0.32, 0.019, constrain=self.carbon_limited,
            product_class=burden.ProductClass.NON_TOXIC)
        chemistry = frame.set_index("source").loc[
            "this GEM: amino acids + translation ATP", "over_measured"]
        assert chemistry == pytest.approx(0.6553, abs=5e-3)
        assert chemistry < 1.0

    @pytest.mark.parametrize("oxygen_lb", [-12.0, -8.0, -3.7])
    def test_the_number_is_invariant_across_the_whole_admissible_oxygen_band(
            self, sink_model, oxygen_lb):
        """The ratio does not move between 3.7 and 12 mmol O2/gDW/h, so the fixture is not
        choosing a point inside the measured band -- only refusing the branch outside it."""
        def constrain(model):
            model.reactions.get_by_id("r_1714").lower_bound = -10.0
            model.reactions.get_by_id("r_1992").lower_bound = oxygen_lb
            return model

        frame = burden.gem_burden_comparison(
            sink_model, 8, 0.32, 0.019, constrain=constrain,
            product_class=burden.ProductClass.NON_TOXIC)
        chemistry = frame.set_index("source").loc[
            "this GEM: amino acids + translation ATP", "over_measured"]
        assert chemistry == pytest.approx(0.6553, abs=5e-3)

    def test_the_unbounded_oxygen_branch_inflates_it_by_a_third(self, sink_model):
        """Recorded rather than deleted: the old fixture's 0.8665 was 1.32x this one, and the
        branch it came from grows faster than any measured mu_max for this organism."""
        frame = burden.gem_burden_comparison(
            sink_model, 8, 0.32, 0.019, constrain=self.unbounded_oxygen,
            product_class=burden.ProductClass.NON_TOXIC)
        chemistry = frame.set_index("source").loc[
            "this GEM: amino acids + translation ATP", "over_measured"]
        assert chemistry == pytest.approx(0.8665, abs=5e-3)
        assert chemistry / 0.6553 == pytest.approx(1.322, abs=5e-3)

    def test_the_measured_row_is_the_axis_the_others_are_scored_against(self, sink_model):
        frame = burden.gem_burden_comparison(
            sink_model, 5, 0.32, 0.019, constrain=self.carbon_limited)
        measured = frame[frame["source"].str.startswith("Kafri")]
        assert len(measured) == 1
        assert measured["over_measured"].iloc[0] == pytest.approx(1.0)
        assert measured["relative_growth_loss"].iloc[0] == pytest.approx(0.05)

    def test_the_displacement_row_agrees_with_the_standalone_ratio(self, sink_model):
        frame = burden.gem_burden_comparison(
            sink_model, 5, 0.32, 0.019, constrain=self.carbon_limited)
        displacement = frame[frame["source"].str.startswith("proteome displacement")]
        assert displacement["over_measured"].iloc[0] == pytest.approx(
            burden.displacement_over_measured(0.32, 0.019), rel=1e-9)

    def test_a_load_past_the_declared_class_ceiling_is_refused_before_the_lp_runs(
            self, sink_model):
        with pytest.raises(burden.BurdenAboveCeiling):
            burden.gem_burden_comparison(
                sink_model, 8, 0.32, 0.019, constrain=self.carbon_limited)


# --------------------------------------------------------------------------------------
# 10. Provenance: every number resolves, and two misattributions stay dead
# --------------------------------------------------------------------------------------

class TestEveryNumberCarriesAResolvableCitation:

    @pytest.mark.parametrize("registry_name", ["GROWTH_PARAMS", "ALLOCATION_PARAMS"])
    def test_no_registered_parameter_is_uncited(self, registry_name):
        registry = getattr(burden, registry_name)
        assert [p.name for p in registry.uncited()] == []

    def test_nothing_is_borrowed_from_another_organism(self):
        assert burden.GROWTH_PARAMS.borrowed() == ()
        assert burden.ALLOCATION_PARAMS.borrowed() == ()

    def test_nothing_is_merely_asserted(self):
        for registry in (burden.GROWTH_PARAMS, burden.ALLOCATION_PARAMS):
            assert registry.by_tag()[Tag.ASSERTED] == 0

    def test_the_provenance_table_carries_both_arms_and_both_gate_lines(self):
        table = burden.provenance()
        assert "PASSES" in table and "REFUSED" in table
        for registry in (burden.GROWTH_PARAMS, burden.ALLOCATION_PARAMS):
            for param in registry.params:
                assert param.name in table

    def test_farkas_is_cited_as_the_elife_paper_the_reference_list_resolves_to(self):
        """eLife 7:e29845, PMID 29377792. Checked against the reference list of the scraped
        PMC copy in scratchpad/burden/namba2022.md, which is the only source for it here."""
        assert "29377792" in CODE
        assert "e29845" in CODE
        assert "PLoS" not in CODE

    def test_the_non_toxic_ceiling_is_attributed_to_fujita_not_kintaka(self):
        """The harvest calls it 'Kintaka/Moriya 2025'. eLife 13:RP99572 (PMID 40960085) lists
        Fujita, Namba, Kamada and Moriya; Kintaka R is not an author. The >40% sentence is
        verbatim from that abstract."""
        assert "Kintaka" not in burden.FUJITA_CEILING.source
        assert "40960085" in burden.FUJITA_CEILING.source
        assert "Fujita" in burden.FUJITA_CEILING.source
        assert "40% of total protein while maintaining yeast growth" in (
            burden.FUJITA_CEILING.source)
        row = burden.misattributed_quotes().iloc[0]
        assert "Kintaka/Moriya 2025" in row["claim"]
        assert "Kintaka R is not an author" in row["source_really_says"]
        assert "Fujita Y, Namba S, Kamada Y and Moriya H" in row["source_really_says"]

    def test_the_2016_kintaka_citation_that_is_correct_is_left_alone(self):
        assert "27538565" in burden.B_SECRETED.source

    def test_the_module_records_the_misattribution_so_it_cannot_come_back(self):
        assert "MISATTRIBUTION" in burden.__doc__
        assert "Kintaka R is not on it" in " ".join(burden.__doc__.split())

    def test_phi_per_copy_no_longer_asserts_a_basis_neither_paper_states(self):
        """Metzl-Raz's abundance is top-3-peptide intensity normalised to PaxDb and Kafri
        Fig 1E's axis is PaxDb/CYCLOPS ppm, which PHI_FROM_PAXDB itself calls MOLAR. Both
        verified against the EuropePMC full texts for PMID 28857745 and PMID 26725116."""
        source = burden.PHI_PER_COPY.source
        assert "MASS fraction, not molecule fraction" not in source
        assert "BASIS NOT STATED BY EITHER PAPER" in source
        assert "three most" in source and "PaxDb" in source
        row = burden.misattributed_quotes().iloc[1]
        assert "MASS fraction" in row["claim"]
        assert "PaxDb ppm counts MOLECULES" in row["source_really_says"]
        assert "NOT rescaled" in row["consequence"]

    def test_the_molar_exposure_is_on_the_reported_arm_and_not_the_growth_arm(self):
        """The growth channel is b*phi = Kafri's 0.01/copy on the growth axis, so a basis
        correction cannot move it; the ceiling comparison is where it would land."""
        for b in burden.BURDEN_SLOPE.sweep_points(5):
            assert b * burden.consistent_phi_per_copy(b) == pytest.approx(
                float(burden.KAFRI_GROWTH_LOSS_PER_COPY))
        f_H = burden.proteome_fraction(8, 0.019)
        assert f_H == pytest.approx(0.152)
        assert f_H * (27.0 / 50.0) == pytest.approx(0.08208)
        assert f_H * (27.0 / 50.0) < float(burden.EGUCHI_CEILING) < f_H

    def test_geiler_samerotte_is_not_credited_with_the_full_hsf1_regulon(self):
        """PMID 21187411 reports 20 of 25 differentially regulated proteins as Hsf1 targets
        and its abstract says 'in the absence of a wider stress response'."""
        doc = " ".join(burden.ProductClass.__doc__.split())
        assert "NOT 'the full Hsf1 regulon'" in doc
        assert "20 have Hsf1p-bound promoters" in doc
        assert "in the absence of a wider stress response" in doc
        assert "3.2% of growth rate" in doc


    def test_every_public_name_is_exported(self):
        for name in burden.__all__:
            assert hasattr(burden, name), name
