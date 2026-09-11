"""The UPR block, tested on the plates it was built to predict.

Every number in `mech/upr.py`'s docstring is re-derived here from the committed exports, so
a stale claim in the prose fails a test rather than surviving in a paragraph -- which is not
a hypothetical: seven assertions in an earlier version of this file were the only thing that
caught a wrong observation model in the module.

Several of these are deliberately assertions that the block MISSES something -- criterion
(a), the two 20260804 blocks, Pincus's deactivation bracket, nine of the ten rows of
DEACTIVATION_VERSUS_RELAXATION, and the two scalars that sit on their box. A test suite that
only encodes the wins is how a model's known limits get lost, and this block's limits are
its most useful output.
"""
from __future__ import annotations

import inspect
import math
import pathlib

import numpy as np
import pandas as pd
import pytest

from ystwin.generator.panel_experiment import OBSERVED_ACTIVITY_CV
from ystwin.mech import upr
from ystwin.mech.integrate import FALLBACK_METHOD, PINNED_METHOD, Verdict
from ystwin.mech.params import FreeScalarGateFailed, RefusedValue, SweptValue, Tag
from ystwin.mech.state import FEDBATCH_5D, PLATE_READ_4H

FIT_EXPORT = "20260722_ER&OxidativeStress_NewProtocol_ANALYSED.xlsx"
REPLICATE_EXPORT = "20260804_ER&OxidativeStress_Replicate4.xlsx"
ANOMALOUS_EXPORT = "20260803_ER&oxidativestress_Replicate3.xlsx"
NO_SHEET_EXPORT = "20260728_ER&OxidativeStress_NewProtocol_Replicate2.xlsx"
SHEET_EXPORTS = ((FIT_EXPORT, "20260722"), (REPLICATE_EXPORT, "20260804"),
                 (ANOMALOUS_EXPORT, "20260803"))
FLOOR = OBSERVED_ACTIVITY_CV


HELD_OUT = [(FIT_EXPORT, "UPRE2"),
            (NO_SHEET_EXPORT, "UPRE1"), (NO_SHEET_EXPORT, "UPRE2"),
            (ANOMALOUS_EXPORT, "UPRE1"), (ANOMALOUS_EXPORT, "UPRE2"),
            (REPLICATE_EXPORT, "UPRE1"), (REPLICATE_EXPORT, "UPRE2")]


@pytest.fixture(scope="module")
def ladder():
    """One fit on 20260722/UPRE1 and every other block scored against it, with no refit."""
    return upr.score_dtt_ladder(FIT_EXPORT, "UPRE1", held_out=HELD_OUT)


def _module_docstring_of(name: str) -> str:
    """The docstring literal that follows a module-level assignment, which `inspect` cannot
    reach because a data docstring is not stored on the object."""
    source = inspect.getsource(upr)
    return source.split(f"\n{name} = ", 1)[1].split('"""')[1]


@pytest.fixture(scope="module")
def fitted_entry(ladder):
    fit, _ = ladder
    return upr.DoseEntry(top_load=fit.parameters["top_load"], top_dose_mM=1.0,
                         declared=True, note="the fitted axis, still a declared axis")


# --------------------------------------------------------------------------------------
# Rule 3: no invented numbers
# --------------------------------------------------------------------------------------

def test_every_constant_carries_a_tag_and_a_source():
    for registry in (upr.UPR_PARAMS, upr.PLATE_PARAMS):
        assert len(registry) > 0
        for param in registry.params:
            assert param.tag in Tag.ALL
            assert param.source.strip()


def test_the_only_uncited_rows_are_a_bnid_and_this_project_s_own_plates():
    """`citations.scan_text` resolves PMIDs and DOIs, not BNIDs or plate exports."""
    assert {p.name for p in upr.UPR_PARAMS.uncited()} == {
        "haploid_cell_volume", "dtt_lethal_dose"}
    assert "BNID 100427" in upr.CELL_VOLUME_UM3.source
    assert "generator/stress_panel.py" in upr.DTT_LETHAL_MM.source


def test_nothing_is_borrowed_from_another_organism():
    """The mammalian rows in Pincus's chain -- R = 12.5 and the BiP decay -- are REFUSED or
    replaced, not imported with a flag."""
    assert upr.UPR_PARAMS.borrowed() == ()
    assert upr.IRE1_BIP_RATIO.tag == Tag.REFUSED
    assert "Bertolotti" in upr.IRE1_BIP_RATIO.source


def test_the_measured_rows_are_the_ones_the_verifier_confirmed():
    assert float(upr.HAC1_MRNA_DECAY_PER_H) == pytest.approx(3.0)
    assert float(upr.HAC1_SPLICING_PER_H) == pytest.approx(60.0 / 11.0)
    assert float(upr.KAR2_PROTEIN_DECAY_PER_H) == pytest.approx(math.log(2.0) / 27.0)
    assert float(upr.HAC1_UPRE2_KD_NM) == 427.0
    assert upr.HAC1_UPRE2_KD_NM.ci95 == (390.0, 464.0)
    for pmid in ("8898194", "9348528", "25466257", "23054834", "17108329", "29361465"):
        assert any(pmid in p.source for p in upr.UPR_PARAMS.params), pmid


def test_hac1_protein_half_life_bracket_is_pal_2007():
    """1-1.5 min, so the rate is ln2 over the midpoint and the ci95 is the bracket."""
    low, high = upr.HAC1_PROTEIN_DECAY_PER_H.ci95
    assert low == pytest.approx(math.log(2.0) / (1.5 / 60.0))
    assert high == pytest.approx(math.log(2.0) / (1.0 / 60.0))
    assert low < float(upr.HAC1_PROTEIN_DECAY_PER_H) < high


def test_the_occupancy_replaces_two_fitted_scalars_with_measured_arithmetic():
    """Pincus's f(H) = H^2/(a0 + a1*H + H^2) fits a0 = 296.5 and a1 = 5.26. This is the
    same step built from Fordyce's Kd and Ho's abundance instead, and the arithmetic is
    re-derived here rather than trusted."""
    nm = (float(upr.HAC1_ABUNDANCE_MOLECULES) / 6.02214076e23
          / (float(upr.CELL_VOLUME_UM3) * 1e-15) * 1e9)
    assert nm == pytest.approx(354.6, abs=0.5)
    assert float(upr.HAC1_MAX_OCCUPANCY) == pytest.approx(
        nm / (float(upr.HAC1_UPRE2_KD_NM) + nm))
    assert upr.HAC1_MAX_OCCUPANCY.tag == Tag.DERIVED
    assert "296.5" in upr.HAC1_MAX_OCCUPANCY.source


def test_the_declared_compartment_is_free_and_the_ladder_saturates_the_ceiling(ladder):
    """The claim HAC1_MAX_OCCUPANCY makes, MEASURED both halves. Hac1 acts in the nucleus
    and the whole-cell volume is a DECLARED compartment, so the ceiling could be 14x higher
    in concentration. Refit there, nothing moves -- because theta enters only as
    (b + theta)/(b + theta0) and basal_share rescales with it. An earlier draft justified
    the same conclusion by saying the knee was never reached; it is, from 0.5 mM up."""
    fit, scores = ladder
    mu = upr._measured_growth_rates(FIT_EXPORT, "UPRE1")
    theta0 = upr._occupancy(upr._spliced_fraction(upr._ire1_active(
        fit.parameters["basal_load"] / (1.0 - fit.parameters["basal_load"]), 1.0,
        fit.parameters["hill_n"]), mu[0.0]))
    entry = upr.DoseEntry(top_load=fit.parameters["top_load"], top_dose_mM=1.0,
                          declared=True, note="the fitted axis")
    model = upr.PlateUpr()
    at_top = model._promoter_trace(
        min(entry.client_load(1.0, fit.parameters["basal_load"]), 1.0 - 1e-9),
        fit.parameters["hill_n"], fit.parameters["basal_share"], mu[1.0], theta0)[1]
    saturated = float(np.max(at_top)) - fit.parameters["basal_share"]
    assert saturated == pytest.approx(float(upr.HAC1_MAX_OCCUPANCY), rel=1e-6)

    nuclear_nm = (float(upr.HAC1_ABUNDANCE_MOLECULES) / 6.02214076e23
                  / (3.0 * 1e-15) * 1e9)
    ceiling = nuclear_nm / (float(upr.HAC1_UPRE2_KD_NM) + nuclear_nm)
    whole_cell = upr._P_MAX
    try:
        upr._P_MAX = ceiling / (1.0 - ceiling)
        nuclear_fit, nuclear_scores = upr.score_dtt_ladder(
            FIT_EXPORT, "UPRE1", held_out=HELD_OUT)
    finally:
        upr._P_MAX = whole_cell
    moved = max(abs(a.relative_rms - b.relative_rms)
                for a, b in zip(scores, nuclear_scores))
    assert moved < FLOOR / 1000.0
    assert (nuclear_fit.parameters["basal_share"] / fit.parameters["basal_share"]
            == pytest.approx(ceiling / float(upr.HAC1_MAX_OCCUPANCY), rel=0.01))
    assert "it is false" in upr.HAC1_MAX_OCCUPANCY.source


def test_the_langmuir_is_single_site_and_saturates_at_the_derived_occupancy():
    """A single-site Langmuir is a hyperbola, so it is CONCAVE: half the driver gives more
    than half the occupancy, never less. Pinned as the half-max identity, which fixes the
    whole curve from one number, rather than as an inequality that could hold by accident."""
    assert upr._occupancy(0.0) == 0.0
    assert upr._occupancy(1.0) == pytest.approx(float(upr.HAC1_MAX_OCCUPANCY))
    assert upr._occupancy(1.0 / upr._P_MAX) == pytest.approx(0.5)
    assert upr._occupancy(0.5) == pytest.approx(0.2934, abs=1e-4)
    assert 0.5 * upr._occupancy(1.0) < upr._occupancy(0.5) < upr._occupancy(1.0)
    assert upr._occupancy(20.0) < 1.0


def test_pincus_refusals_are_listed_where_a_reader_will_find_them():
    joined = " ".join(upr.REFUSED_FROM_PINCUS)
    assert "Personal Communication" in joined
    assert "Stroberg 2018" in joined
    assert "2.15 um^3" in joined
    assert len(upr.REFUSED_FROM_PINCUS) == 8


# --------------------------------------------------------------------------------------
# Condition (1): the entry function is REFUSED
# --------------------------------------------------------------------------------------

def test_the_dose_entry_refuses_to_be_a_number():
    for param in (upr.DTT_ENTRY, upr.TUNICAMYCIN_ENTRY):
        assert param.tag == Tag.REFUSED
        with pytest.raises(RefusedValue) as raised:
            float(param)
        assert "REFUSED" in str(raised.value)
        assert param.missing.strip()


def test_the_refusal_names_what_would_close_it():
    assert "lumenal" in upr.DTT_ENTRY.missing
    assert "Alg7" in upr.TUNICAMYCIN_ENTRY.missing
    assert "stress_panel" in upr.TUNICAMYCIN_ENTRY.provenance()
    assert "stress_panel" in upr.TUNICAMYCIN_ENTRY.reason


def test_the_missing_deactivation_mechanism_is_refused_not_fitted():
    """The block cannot terminate the response on the timescale the plates show, and the
    refusal names both measurements that would settle why."""
    assert upr.DEACTIVATION_MECHANISM.tag == Tag.REFUSED
    with pytest.raises(RefusedValue):
        float(upr.DEACTIVATION_MECHANISM)
    assert "Ellman" in upr.DEACTIVATION_MECHANISM.missing
    assert "OD600" in upr.DEACTIVATION_MECHANISM.missing
    assert "DEACTIVATION_VERSUS_RELAXATION" in upr.DEACTIVATION_MECHANISM.source


def test_the_deactivation_table_is_re_derived_and_only_one_row_of_ten_clears():
    """The claim the module makes is that the comparison is suggestive on ONE block and
    inconclusive on the other nine. Every entry is recomputed here, and the count is
    asserted, so the table cannot quietly be summarised to its best row."""
    export = {"20260722": FIT_EXPORT, "20260804": REPLICATE_EXPORT,
              "20260803": ANOMALOUS_EXPORT}
    above = 0
    for stamp, construct, dose, decay, relax, conservative in upr.DEACTIVATION_VERSUS_RELAXATION:
        measured = upr.deactivation_rate(export[stamp], construct, dose)
        bound = upr.chaperone_relaxation_rate(export[stamp], construct, dose)
        assert measured == pytest.approx(decay, abs=5e-3)
        assert bound == pytest.approx(relax, abs=5e-3)
        assert measured[0] / bound[1] == pytest.approx(conservative, abs=0.02)
        above += conservative > 1.0
    assert above == 1
    assert len(upr.DEACTIVATION_VERSUS_RELAXATION) == 10


def test_the_two_blocks_with_no_decay_to_fit_raise_rather_than_return_a_number():
    """20260803/UPRE2 never rises 2% above its control after 2 h. That is why it is absent
    from the table, and it is refused rather than fitted to noise."""
    for dose in (0.5, 1.0):
        with pytest.raises(ValueError, match="no decay to fit"):
            upr.deactivation_rate(ANOMALOUS_EXPORT, "UPRE2", dose)


def test_the_growth_rate_argument_will_not_accept_a_rate_that_is_not_on_the_plate():
    with pytest.raises(ValueError, match="MEASURED column"):
        upr.upre_activity_trace(FIT_EXPORT, "UPRE1", rate="mu_guess")


def test_an_undeclared_dose_axis_raises():
    with pytest.raises(upr.EntryRefused) as raised:
        upr.DoseEntry(top_load=0.4, top_dose_mM=1.0)
    assert "DECLARED SENSITIVITY AXIS" in str(raised.value)
    assert isinstance(raised.value, RefusedValue)


def test_a_declared_axis_is_admissible_only_below_saturation_and_below_the_lethal_dose():
    with pytest.raises(ValueError, match="top_load"):
        upr.DoseEntry(top_load=1.0, top_dose_mM=1.0, declared=True)
    with pytest.raises(ValueError, match="lethal"):
        upr.DoseEntry(top_load=0.4, top_dose_mM=2.0, declared=True)
    entry = upr.DoseEntry(top_load=0.4, top_dose_mM=1.0, declared=True)
    assert entry.client_load(0.0, 0.1) == pytest.approx(0.1)
    assert entry.client_load(1.0, 0.1) == pytest.approx(0.4)
    assert entry.client_load(0.5, 0.1) == pytest.approx(0.25)
    assert not hasattr(entry, "load")  # see DoseEntry.client_load for why the name moved


def test_the_block_refuses_the_supra_lethal_rungs_of_its_own_ladder():
    assert upr.DOSE_LADDER_MM[-2:] == (2.0, 5.0)
    assert 2.0 not in upr.SUBLETHAL_DOSE_MM and 5.0 not in upr.SUBLETHAL_DOSE_MM
    data = upr.upre_activity_fold(FIT_EXPORT, "UPRE1")
    assert set(np.unique(data.dose_mM)) == {0.1, 0.2, 0.5, 1.0}
    trace = upr.upre_activity_trace(FIT_EXPORT, "UPRE1")
    assert set(np.unique(trace.dose_mM)) == {0.1, 0.2, 0.5, 1.0}


# --------------------------------------------------------------------------------------
# Condition (2): the Hac1 protein and UPRE occupancy states are back
# --------------------------------------------------------------------------------------

def test_the_hac1_protein_state_exists_and_is_measured():
    state = upr.UPR_STATES["hac1_protein"]
    assert "17108329" in state.source
    assert state.constrained_by
    low, high = state.taus_in(PLATE_READ_4H)
    assert low == pytest.approx((1.0 / 60.0) / math.log(2.0))
    assert high == pytest.approx((1.5 / 60.0) / math.log(2.0))


def test_the_chain_carries_the_occupancy_explicitly(fitted_entry, ladder):
    fit, _ = ladder
    chain = upr.UprChain(basal_load=fit.parameters["basal_load"],
                         hill_n=fit.parameters["hill_n"],
                         basal_share=fit.parameters["basal_share"],
                         folding_rate_per_h=3.0, mu_per_h=0.214)
    trajectory = chain.simulate(fitted_entry, 1.0)
    theta = upr.UprChain.occupancy(trajectory)
    assert theta.shape == trajectory.t.shape
    assert np.all((theta >= 0.0) & (theta <= float(upr.HAC1_MAX_OCCUPANCY)))


def test_the_declared_basal_state_is_actually_a_fixed_point():
    """The test that found the defect the module docstring records. `basal()` returns a
    state and `rhs()` returns the equations that state is meant to sit still in; nothing
    was checking that it did, and it did not -- the Hac1 pool was normalised two different
    ways and the zero-dose control decayed 37% across the read."""
    chain = upr.UprChain(basal_load=0.2424, hill_n=8.0, basal_share=1.8093,
                         folding_rate_per_h=3.0, mu_per_h=0.214)
    start = chain.basal()
    derivative = chain.rhs(chain.basal_load, 0.0)(
        0.0, [start[name] for name in upr.UPR_STATES.names])
    assert np.max(np.abs(derivative)) < 1e-12
    flat = chain.simulate(
        upr.DoseEntry(top_load=0.9, top_dose_mM=1.0, declared=True, note="unused at 0 mM"),
        0.0, n_points=51)
    for name in upr.UPR_STATES.names:
        trace = flat.of(name)
        assert abs(trace[-1] - trace[0]) / trace[0] < 1e-9


def test_the_chain_and_the_shape_agree_on_the_occupancy_they_share():
    """The chain integrates the Hac1 pool and the plate shape solves it algebraically. Both
    are normalised so the occupancy reaches HAC1_MAX_OCCUPANCY at FULL splicing -- that is
    what makes Ho 2018's INDUCED abundance the right number to derive it from -- and if the
    two conventions drift apart, `reduction_bias` measures the drift and calls it the price
    of the reduction."""
    mu = 0.214
    chain = upr.UprChain(basal_load=0.2424, hill_n=8.0, basal_share=0.0,
                         folding_rate_per_h=3.0, mu_per_h=mu)
    start = chain.basal()
    # The chain's own fixed point for the protein pool, at the basal Ire1 activation.
    assert start["hac1_protein"] == pytest.approx(
        upr._P_MAX * start["hac1_spliced"] / start["hac1_spliced_full"], rel=1e-12)
    # And at FULL splicing the same normalisation lands exactly on the measured occupancy.
    assert upr._occupancy(upr._spliced_fraction(1.0, mu)) == pytest.approx(
        float(upr.HAC1_MAX_OCCUPANCY), rel=1e-12)
    assert upr._spliced_fraction(1.0, mu) == pytest.approx(1.0, rel=1e-12)


def test_criterion_b_every_state_names_an_assay_or_a_declared_axis():
    for state in upr.UPR_STATES:
        assert state.constrained_by or state.sweep_axis
        assert not (state.constrained_by and state.sweep_axis)


def test_ire1_is_absent_as_a_state_because_no_rate_for_it_exists():
    assert "ire1" not in " ".join(upr.UPR_STATES.names)
    assert upr.IRE1_DEACTIVATION.tag == Tag.REFUSED
    assert "There are no known measurements" in upr.IRE1_DEACTIVATION.source


# --------------------------------------------------------------------------------------
# Criterion (d): the reduction, audited and then priced
# --------------------------------------------------------------------------------------

def test_on_the_plate_read_exactly_one_state_survives_the_ratio_test():
    reduced = upr.upr_reduction(PLATE_READ_4H)
    assert reduced.integrated.names == ("chaperone",)
    assert len(reduced.eliminated) == 5
    assert reduced.verdict("chaperone") is Verdict.KEEP
    row = reduced.row("chaperone")
    assert 0.146 < row.ratio_low < row.ratio_high < 1.0


def test_the_two_windows_disagree_which_is_why_the_window_is_an_argument():
    plate = upr.upr_reduction(PLATE_READ_4H)
    fedbatch = upr.upr_reduction(FEDBATCH_5D)
    assert len(fedbatch.eliminated) == 6
    assert len(fedbatch.integrated) == 0
    assert plate.integrated.names != fedbatch.integrated.names


def test_the_response_deactivation_window_is_offered_not_applied():
    """Scoring the ratio against Pincus's 2 h asks a different question and gets a
    different answer; the module says so and this pins both halves."""
    strict = upr.upr_reduction(PLATE_READ_4H,
                               driver_persists_h=upr.RESPONSE_DEACTIVATION_H)
    assert len(strict.eliminated) == 2
    assert strict.verdict("hac1_spliced") is Verdict.KEEP
    assert strict.verdict("er_client") is Verdict.STRADDLES
    assert "RESPONSE_DEACTIVATION_H" in upr.upr_reduction.__doc__
    assert "different question" in upr.upr_reduction.__doc__
    assert upr.upr_reduction(PLATE_READ_4H).integrated.names == ("chaperone",)


def test_the_reduction_is_free_to_the_assay_and_the_price_is_measured(fitted_entry, ladder):
    """The finding, asserted as a finding: eliminating five of the six states costs less
    than the plate can see, at every sub-lethal rung."""
    fit, _ = ladder
    bias = upr.reduction_bias(fitted_entry, (0.1, 0.2, 0.5, 1.0), mu_per_h=0.214,
                              hill_n=fit.parameters["hill_n"],
                              basal_load=fit.parameters["basal_load"],
                              basal_share=fit.parameters["basal_share"])
    assert all(v < FLOOR for v in bias.values())
    assert max(bias.values()) / FLOOR == pytest.approx(0.34, abs=0.04)
    assert bias[0.1] == pytest.approx(0.021, abs=0.005)
    assert bias[1.0] == pytest.approx(0.001, abs=0.001)
    # A hump, not a monotone: two earlier drafts asserted the wrong shape here, the second
    # of them because the chain's own control was decaying.
    assert bias[0.2] > bias[0.1] > bias[0.5] > bias[1.0]


def test_the_docstring_quotes_the_measured_bias(fitted_entry, ladder):
    fit, _ = ladder
    bias = upr.reduction_bias(fitted_entry, (0.1, 1.0), mu_per_h=0.214,
                              hill_n=fit.parameters["hill_n"],
                              basal_load=fit.parameters["basal_load"],
                              basal_share=fit.parameters["basal_share"])
    assert f"{bias[0.1]:.3f}" in upr.reduction_bias.__doc__
    assert f"{bias[1.0]:.3f}" in upr.reduction_bias.__doc__


def test_the_whole_chain_is_algebra_as_far_as_this_plate_can_tell(ablation, fitted_entry,
                                                                  ladder):
    """The sharpest thing this block learned, and it needs both halves. Eliminating the five
    states the ratio test drops costs at most 0.34x the floor; keeping the one it holds buys
    0.0029x. Neither side of the audit's cut is visible to the assay."""
    fit, _ = ladder
    reduced = upr.upr_reduction(PLATE_READ_4H)
    assert reduced.integrated.names == ("chaperone",) == upr.PLATE_STATE.names
    bias = upr.reduction_bias(fitted_entry, (0.1, 0.2, 0.5, 1.0), mu_per_h=0.214,
                              hill_n=fit.parameters["hill_n"],
                              basal_load=fit.parameters["basal_load"],
                              basal_share=fit.parameters["basal_share"])
    assert max(bias.values()) < FLOOR
    assert ablation.ratio < max(bias.values()) / FLOOR < 1.0


# --------------------------------------------------------------------------------------
# Criterion (e): the free-scalar gate
# --------------------------------------------------------------------------------------

def test_the_reference_chain_is_refused_by_its_own_gate():
    gate = upr.UPR_PARAMS.gate()
    assert len(gate.free) == 13
    assert len(gate.targets) == 4
    assert not gate.passes
    with pytest.raises(FreeScalarGateFailed):
        upr.UPR_PARAMS.require_gate()


def test_the_shipped_shape_passes_and_the_two_counts_agree():
    registry = upr.PLATE_PARAMS.gate()
    fitted = upr.plate_gate()
    assert len(registry.free) == len(fitted.free) == 4
    assert len(registry.targets) == len(fitted.targets) == 4
    assert registry.passes and fitted.passes
    upr.PLATE_PARAMS.require_gate()
    assert set(fitted.free) == set(upr.PlateUpr().parameter_names)


def test_the_dropped_rows_are_named_with_a_reason():
    dropped = dict(upr.DROPPED_TO_PASS_THE_GATE)
    assert len(dropped) == 9
    reference = {p.name for p in upr.UPR_PARAMS.free_scalars()}
    shipped = {p.name for p in upr.PLATE_PARAMS.free_scalars()}
    assert reference - shipped == set(dropped)
    for reason in dropped.values():
        assert len(reason) > 30


def test_the_prose_counts_match_the_registries_that_produce_them():
    """Both registry docstrings quote their own counts in words, and both had gone stale.
    Pin them to the count rather than to the word."""
    words = {4: "four", 9: "nine", 13: "thirteen"}
    free_chain = len(upr.UPR_PARAMS.free_scalars())
    chain_doc = _module_docstring_of("UPR_PARAMS")
    plate_doc = _module_docstring_of("PLATE_PARAMS")
    dropped_doc = _module_docstring_of("DROPPED_TO_PASS_THE_GATE")
    assert f"{free_chain} free against {len(upr.UPR_PARAMS.independent_targets())}" in chain_doc
    dropped = len(upr.DROPPED_TO_PASS_THE_GATE)
    assert f"{words[dropped]} fewer free ones" in plate_doc
    assert f"The {words[dropped]} rows the plate shape drops" in dropped_doc
    assert dropped == free_chain - len(upr.PLATE_PARAMS.free_scalars())


def test_the_fitted_target_is_not_counted_in_the_denominator():
    assert upr.UPRE1_LADDER_FITTED.fitted
    assert upr.UPRE1_LADDER_FITTED.name not in upr.plate_gate().targets
    assert upr.UPRE1_LADDER_FITTED.name in upr.plate_gate().fitted_targets


def test_the_literature_target_does_not_borrow_the_plate_floor():
    assert upr.HAC1_DEACTIVATION_TIME.noise_floor == pytest.approx(2.0 / 3.0)
    assert upr.HAC1_DEACTIVATION_TIME.noise_floor != FLOOR
    assert "Northern" in upr.HAC1_DEACTIVATION_TIME.assay
    for target in (upr.UPRE2_HELD_OUT, upr.ACTIVITY_TRANSIENT,
                   upr.PLATE_REPLICATE_HELD_OUT):
        assert target.noise_floor == FLOOR


# --------------------------------------------------------------------------------------
# The observable and the fit
# --------------------------------------------------------------------------------------

def test_the_scored_observable_comes_out_of_the_committed_file_untouched():
    """`upre_activity_fold` reads activity_late and divides. If it ever starts re-deriving
    the activity instead, this catches it."""
    from ystwin import paths

    data = upr.upre_activity_fold(FIT_EXPORT, "UPRE1")
    assert set(upr._PLATE_COLUMNS) <= set(data.columns)
    assert len(data) == 12
    table = pd.read_csv(paths.outputs_dir() / "sensor_characterisation.csv")
    rows = table[(table.plate.str.startswith("20260722")) & (table.construct == "UPRE1")
                 & (table.stressor == "DTT")]
    reference = float(rows[rows.dose_mM == 0.0].activity_late.mean())
    for _, row in data.iterrows():
        committed = rows[(rows.well == row.well)].activity_late
        assert float(committed.iloc[0]) / reference == pytest.approx(row.signal)
    assert (data.mu_per_h > 0).all()
    assert "read_h" not in data.columns


def test_the_plate_with_no_dose_response_sheet_is_still_scored():
    """The second defect this pass found. `read_h` was in the scored frame, no equation read
    it, and it came from `plate/replay.py` -- so the whole 20260728 replicate was silently
    excluded although its activity_late is committed alongside every other plate's."""
    rows = upr.upre_activity_fold(NO_SHEET_EXPORT, "UPRE1")
    assert len(rows) == 12
    assert "read_h" not in rows.columns
    with pytest.raises(ValueError, match="no per-construct dose-response sheet"):
        upr.read_duration_h(NO_SHEET_EXPORT)
    assert "No equation read it" in _module_docstring_of("_PLATE_COLUMNS")


def test_the_declared_window_matches_the_reads_that_can_check_it():
    """4.14 h is declared in `mech/state.py`; the three committed sheets end at 4.0 h. The
    3.5% is not hidden -- it is measured on the prediction, and it is 0.84% of the floor."""
    from ystwin.mech.state import Window

    for export, _ in SHEET_EXPORTS:
        assert upr.read_duration_h(export) == pytest.approx(4.0)
    shorter = Window(name="the reads themselves", duration_h=4.0,
                     growth_rate_low_per_h=PLATE_READ_4H.growth_rate_low_per_h,
                     growth_rate_high_per_h=PLATE_READ_4H.growth_rate_high_per_h,
                     source="max time_h of the committed per-construct sheets")
    data = upr.upre_activity_fold(FIT_EXPORT, "UPRE1")
    guess = {"basal_load": 0.2424, "top_load": 0.99, "hill_n": 8.0, "basal_share": 1.8093}
    declared = upr.PlateUpr().predict(data, guess)
    measured = upr.PlateUpr(window=shorter).predict(data, guess)
    moved = float(np.max(np.abs(measured - declared) / declared))
    assert moved < FLOOR / 100.0
    assert f"{moved:.4f}" in _module_docstring_of("_PLATE_COLUMNS")


def test_the_committed_sheet_is_per_cell_which_is_why_the_accumulation_model_went():
    """The root cause, pinned. The late-window mean of replay's `signal` reproduces
    `naive_late` -- the PER-CELL column -- not a total. An observation model that treats it
    as a biomass integral is wrong, and the test says so in the repository rather than in a
    commit message."""
    from ystwin import paths

    table = pd.read_csv(paths.outputs_dir() / "sensor_characterisation.csv")
    worst = 0.0
    for export, stamp in SHEET_EXPORTS:
        block = upr.replay.load_doses(export)
        for construct in ("UPRE1", "UPRE2"):
            wide = block[block.construct == construct].pivot_table(
                index="time_h", columns="dose_mM", values="signal")
            late = slice(int(len(wide) * upr.LATE_WINDOW_FRACTION), None)
            rows = table[(table.plate.str.startswith(stamp))
                         & (table.construct == construct) & (table.stressor == "DTT")]
            committed = rows.groupby("dose_mM").naive_late.mean()
            for dose in wide.columns:
                if dose not in committed.index:
                    continue
                mine = float(wide[dose].iloc[late].mean())
                worst = max(worst, abs(mine - committed[dose]) / committed[dose])
    assert worst < 0.013
    assert "NORMALISED" in upr.upre_activity_fold.__doc__


def test_the_diagnostic_trace_reproduces_the_committed_activity_in_the_late_window():
    """The one approximation the diagnostic makes -- one scalar mu per dose instead of the
    time-resolved rate -- is MEASURED against the committed file rather than argued for."""
    from ystwin.reporter import SpecificFluorescence, promoter_activity
    from ystwin import paths

    table = pd.read_csv(paths.outputs_dir() / "sensor_characterisation.csv")
    worst = 0.0
    for export, stamp in SHEET_EXPORTS:
        for construct in ("UPRE1", "UPRE2"):
            block = upr.replay.load_doses(export)
            wide = block[block.construct == construct].pivot_table(
                index="time_h", columns="dose_mM", values="signal")
            times = wide.index.to_numpy(float)
            rows = table[(table.plate.str.startswith(stamp))
                         & (table.construct == construct) & (table.stressor == "DTT")]
            mu = rows.groupby("dose_mM").mu_late.mean()
            reference = rows.groupby("dose_mM").activity_late_percell.mean()
            late = slice(int(len(times) * upr.LATE_WINDOW_FRACTION), None)
            for dose in upr.SUBLETHAL_DOSE_MM:  # the two supra-lethal rungs reach 2.4%
                activity = promoter_activity(
                    times, SpecificFluorescence(wide[dose].to_numpy(float)),
                    np.full_like(times, float(mu[dose])), upr._PLATE_KINETICS)
                mine = float(np.nanmean(activity[late]))
                worst = max(worst, abs(mine - reference[dose]) / reference[dose])
    assert worst < 0.014
    assert worst < FLOOR / 10.0


def test_growth_rate_is_measured_and_travels_with_the_rows():
    data = upr.upre_activity_fold(FIT_EXPORT, "UPRE1")
    measured = upr._measured_growth_rates(FIT_EXPORT, "UPRE1")
    for dose, mu in zip(data.dose_mM, data.mu_per_h):
        assert mu == pytest.approx(measured[float(dose)])
    assert data.mu_control_per_h.nunique() == 1
    assert "mu" not in " ".join(upr.PlateUpr().parameter_names)


def test_the_fit_is_deterministic():
    data = upr.upre_activity_fold(FIT_EXPORT, "UPRE1")
    model = upr.PlateUpr()
    first, second = model.refit(data), model.refit(data)
    assert first.parameters == second.parameters


def test_the_fit_is_weighted_the_way_the_floor_is_declared(ladder):
    """The floor is a CV, so the residual is relative. The reported fit RMS and the score
    are then the same quantity, which is what makes them comparable."""
    fit, scores = ladder
    fitted = next(s for s in scores if not s.held_out)
    assert fit.residual_rms == pytest.approx(fitted.relative_rms, rel=1e-6)
    assert "multiplicative" in upr.PlateUpr.refit.__doc__


def test_the_fit_lands_inside_every_declared_bracket(ladder):
    fit, _ = ladder
    low, high = upr.IRE1_HILL_N.bounds
    assert low <= fit.parameters["hill_n"] <= high
    assert 0.0 < fit.parameters["basal_load"] < fit.parameters["top_load"] < 1.0
    assert fit.n_free == 4


def test_two_of_the_four_scalars_sit_on_their_declared_box_and_the_module_says_so(ladder):
    """The rail, pinned. Twelve wells at four distinct doses do not identify four scalars,
    and the docstring says which two go to the edge rather than leaving it in the dict."""
    fit, _ = ladder
    low, high = upr.PlateUpr().bounds()
    names = upr.PlateUpr().parameter_names
    at_edge = {n for n, lo, hi in zip(names, low, high)
               if fit.parameters[n] in (pytest.approx(hi, rel=1e-6),
                                        pytest.approx(lo, abs=1e-9))}
    assert at_edge == {"hill_n", "top_load"}
    assert fit.parameters["hill_n"] == pytest.approx(upr.IRE1_HILL_N.bounds[1])
    assert "TWO OF THE FOUR SIT ON THE BOX" in upr.PlateUpr.__doc__
    assert "0.2424" in upr.PlateUpr.__doc__ and "0.9900" in upr.PlateUpr.__doc__


def test_the_in_sample_and_held_out_construct_are_reported_whole(ladder):
    _, scores = ladder
    by_key = {(s.export, s.construct): s for s in scores}
    fitted = by_key[(FIT_EXPORT, "UPRE1")]
    held = by_key[(FIT_EXPORT, "UPRE2")]
    assert not fitted.held_out and held.held_out
    assert fitted.relative_rms == pytest.approx(0.0374, abs=0.005)
    assert held.relative_rms == pytest.approx(0.1417, abs=0.005)
    assert held.clears_floor


def test_criterion_c_five_of_seven_held_out_blocks_clear_and_two_do_not(ladder):
    """UPRE1 and UPRE2 share one Hac1 pool and one upre_basal_share, so UPRE2 costs no
    parameter at all. Fordyce 2012 says the two site geometries differ, which is why this
    is a test and not an assumption -- and the two misses are asserted, not hidden."""
    fit, scores = ladder
    held = [s for s in scores if s.held_out]
    assert len(held) == 7
    assert sum(s.clears_floor for s in held) == 5
    missed = {(s.export, s.construct) for s in held if not s.clears_floor}
    assert missed == {(REPLICATE_EXPORT, "UPRE1"), (REPLICATE_EXPORT, "UPRE2")}
    assert "same Hac1 pool" in upr.UPRE_BASAL_SHARE.source
    assert "two site geometries" in upr.UPRE_BASAL_SHARE.source
    assert fit.n_free == len(upr.PlateUpr().parameter_names)


def test_the_plate_the_old_observable_missed_is_now_the_best_predicted(ladder):
    """20260803 missed by 4.7x on the accumulated increment and clears comfortably on the
    growth-corrected fold. Its anomaly is in the raw control accumulation, which the fold
    does not use, so this is evidence about the OBSERVABLE and not about the mechanism."""
    _, scores = ladder
    anomalous = [s for s in scores if s.export == ANOMALOUS_EXPORT]
    assert len(anomalous) == 2
    assert all(s.clears_floor for s in anomalous)
    increments = {}
    for export in (FIT_EXPORT, REPLICATE_EXPORT, ANOMALOUS_EXPORT):
        block = upr.replay.load_doses(export)
        control = block[(block.construct == "UPRE1") & (block.dose_mM == 0.0)]
        control = control.sort_values("time_h")
        increments[export] = float(control.signal.iloc[-1] - control.signal.iloc[0])
    assert increments[ANOMALOUS_EXPORT] < increments[FIT_EXPORT] / 3.0
    assert increments[ANOMALOUS_EXPORT] < increments[REPLICATE_EXPORT] / 3.0
    assert "565 RFU" in upr.score_dtt_ladder.__doc__


# --------------------------------------------------------------------------------------
# Criterion (a): the ablation, which this block FAILS
# --------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ablation():
    return upr.upr_ablation()


def test_removing_the_chaperone_pool_does_not_move_the_plate_reading(ablation):
    """The block's headline result, asserted as a failure. Sixty times below the floor."""
    assert not ablation.clears_floor
    assert ablation.effect == pytest.approx(0.00042, abs=0.0002)
    assert ablation.ratio == pytest.approx(0.0029, abs=0.0015)
    assert ablation.floor.noise_floor == FLOOR
    assert "FAILS IT" in upr.__doc__


def test_the_ablation_scores_against_the_plate_assay_and_no_other(ablation):
    assert ablation.floor.assay == "96-well plate reader, mCitrine 480/530, OD600-corrected"
    assert ablation.scored == "relative"
    assert ablation.full_n_free == ablation.reduced_n_free == 4
    assert ablation.n_rows == 12


def test_the_frozen_model_predicts_the_held_out_blocks_just_as_well(ladder):
    """The ablation summarises one number; this checks the whole held-out column, which is
    where a real difference would have to show up. It does not."""
    live_fit, live_scores = ladder
    frozen_fit, frozen_scores = upr.score_dtt_ladder(
        FIT_EXPORT, "UPRE1", held_out=HELD_OUT, model=upr.FrozenChaperoneUpr())
    held = [(a, b) for a, b in zip(live_scores, frozen_scores) if a.held_out]
    assert len(held) == 7
    for live, frozen in held:
        assert live.export == frozen.export and live.construct == frozen.construct
        assert abs(live.relative_rms - frozen.relative_rms) < FLOOR / 50.0
    # The frozen model is the BETTER of the two on six of the seven, so the sign does not
    # rescue the state either.
    assert sum(frozen.relative_rms < live.relative_rms for live, frozen in held) == 6


def test_the_measured_transient_is_far_larger_than_the_live_model_produces(ladder):
    """Why the ablation is small: the live model barely makes a transient on this window,
    while the plate's fold at 1 mM falls by more than half over the same hours."""
    fit, _ = ladder
    trace = upr.upre_activity_trace(FIT_EXPORT, "UPRE1")
    top = trace[np.isclose(trace.dose_mM, 1.0)].sort_values("time_h")
    late = top[top.time_h >= top.time_h.max() / 2.0]
    observed_fall = float(late.signal.iloc[0] / late.signal.iloc[-1])
    assert observed_fall > 2.0
    model = upr.PlateUpr()
    read_h = upr.read_duration_h(FIT_EXPORT)
    mu = upr._measured_growth_rates(FIT_EXPORT, "UPRE1")
    theta0 = upr._occupancy(upr._spliced_fraction(upr._ire1_active(
        fit.parameters["basal_load"] / (1.0 - fit.parameters["basal_load"]), 1.0,
        fit.parameters["hill_n"]), mu[0.0]))
    load = fit.parameters["top_load"]
    grid, rate = model._promoter_trace(load, fit.parameters["hill_n"],
                                       fit.parameters["basal_share"], mu[1.0], theta0)
    window = grid >= read_h / 2.0
    predicted_fall = float(rate[window][0] / rate[window][-1])
    assert predicted_fall < 1.5
    assert observed_fall / predicted_fall > 1.8


# --------------------------------------------------------------------------------------
# The chain, and the transient it produces
# --------------------------------------------------------------------------------------

def test_the_chain_produces_a_step_not_a_pulse_at_the_fitted_parameters(fitted_entry,
                                                                       ladder):
    """The chain side of the same shortfall. At the parameters the committed fold selects,
    the occupancy is STILL RISING at the end of the 4.14 h window -- Pincus's measured pulse
    at 1.5 mM DTT is not reproduced, and the module refuses a deactivation rate rather than
    fitting one that would produce it. `DEACTIVATION_MECHANISM` is that refusal."""
    fit, _ = ladder
    chain = upr.UprChain(basal_load=fit.parameters["basal_load"],
                         hill_n=fit.parameters["hill_n"],
                         basal_share=fit.parameters["basal_share"],
                         folding_rate_per_h=3.0, mu_per_h=0.214)
    trajectory = chain.simulate(fitted_entry, 1.0, n_points=201)
    theta = upr.UprChain.occupancy(trajectory)
    assert int(np.argmax(theta)) == len(theta) - 1
    assert theta[-1] > 10.0 * theta[0]
    assert trajectory.of("chaperone")[-1] > 1.0
    assert upr.DEACTIVATION_MECHANISM.tag == Tag.REFUSED


def test_the_chain_does_produce_a_pulse_on_a_gentler_declared_axis(ladder):
    """The chain is structurally capable of the pulse -- give the entry a gentler declared
    axis and the occupancy peaks inside the window and comes back down. What no admissible
    setting buys is a pulse at the top rung of the ladder the plate actually ran."""
    fit, _ = ladder
    gentle = upr.DoseEntry(top_load=0.45, top_dose_mM=1.0, declared=True,
                           note="a declared axis, not a measurement")
    chain = upr.UprChain(basal_load=0.20, hill_n=fit.parameters["hill_n"],
                         basal_share=0.0, folding_rate_per_h=3.0, mu_per_h=0.214)
    trajectory = chain.simulate(gentle, 1.0, n_points=201)
    theta = upr.UprChain.occupancy(trajectory)
    peak = int(np.argmax(theta))
    assert 0 < peak < len(theta) - 1
    assert theta[-1] < theta[peak]
    assert trajectory.of("chaperone")[-1] > 1.0


def test_the_chain_refuses_a_supra_lethal_dose(fitted_entry, ladder):
    fit, _ = ladder
    chain = upr.UprChain(basal_load=fit.parameters["basal_load"],
                         hill_n=fit.parameters["hill_n"],
                         basal_share=fit.parameters["basal_share"],
                         folding_rate_per_h=3.0, mu_per_h=0.214)
    with pytest.raises(ValueError, match="lethal"):
        chain.simulate(fitted_entry, 2.0)


def test_the_chain_refuses_a_folding_rate_off_its_declared_axis():
    with pytest.raises(SweptValue):
        float(upr.FOLDING_RATE_PER_H)
    with pytest.raises(ValueError, match="swept"):
        upr.UprChain(basal_load=0.1, hill_n=4.0, basal_share=0.0,
                     folding_rate_per_h=0.5, mu_per_h=0.2)
    with pytest.raises(ValueError, match="Korennykh"):
        upr.UprChain(basal_load=0.1, hill_n=12.0, basal_share=0.0,
                     folding_rate_per_h=3.0, mu_per_h=0.2)
    with pytest.raises(ValueError, match="MEASURED"):
        upr.UprChain(basal_load=0.1, hill_n=4.0, basal_share=0.0,
                     folding_rate_per_h=3.0, mu_per_h=0.0)


def test_the_chain_misses_pincus_s_deactivation_bracket_at_the_fitted_parameters(
        fitted_entry, ladder):
    """Nothing in the fit sees a deactivation time, so this is an unfitted prediction -- and
    it MISSES. At the parameters the committed fold selects, the chain is still above a
    tenth of its peak at the end of the 4.14 h window, where Pincus MEASURED 2 h at 1.5 mM.
    That is the same shortfall the plate shows, seen from the chain's side, and the function
    reports it rather than clipping to the window's end."""
    fit, _ = ladder
    with pytest.raises(ValueError, match="still at more than"):
        upr.deactivation_time(fitted_entry, 1.5, mu_per_h=0.214,
                              hill_n=fit.parameters["hill_n"],
                              basal_load=fit.parameters["basal_load"],
                              basal_share=fit.parameters["basal_share"])
    assert upr.HAC1_DEACTIVATION_TIME.noise_floor == pytest.approx(2.0 / 3.0)


def test_the_chain_does_deactivate_at_a_gentler_declared_axis(ladder):
    """The refusal is about RATE, not about the chain being structurally incapable: on a
    gentler declared entry the six-state chain deactivates inside the window. What no
    admissible setting buys is a rate fast enough for the top rung of this ladder."""
    fit, _ = ladder
    gentle = upr.DoseEntry(top_load=0.45, top_dose_mM=1.0, declared=True,
                           note="a declared axis, not a measurement")
    times = [upr.deactivation_time(gentle, 1.5, mu_per_h=0.214,
                                   hill_n=fit.parameters["hill_n"], basal_load=0.20,
                                   basal_share=0.0, fraction=f)
             for f in (0.5, 0.2, 0.1)]
    assert times == sorted(times)
    assert all(0.0 < t < PLATE_READ_4H.duration_h for t in times)


def test_the_deactivation_threshold_is_named_as_a_convention():
    """It has to be chosen and said out loud, and the docstring points at the plate's own
    deactivation, which is measured elsewhere and is faster than the chain can go."""
    assert upr.DEACTIVATION_FRACTION == 0.1
    quoted = inspect.getsource(upr).split("DEACTIVATION_FRACTION = ")[1].split('"""')[1]
    assert "ASSERTED as a" in quoted
    assert "deactivation_rate" in quoted


# --------------------------------------------------------------------------------------
# Numerics: the settings that are not model parameters, pinned by measurement
# --------------------------------------------------------------------------------------

def test_the_integrator_choice_does_not_move_the_answer(ladder):
    fit, _ = ladder
    data = upr.upre_activity_fold(FIT_EXPORT, "UPRE1")
    pinned = upr.PlateUpr().predict(data, fit.parameters)
    assert upr.PlateUpr().method == PINNED_METHOD
    fast = upr.PlateUpr(method=FALLBACK_METHOD).predict(data, fit.parameters)
    assert float(np.max(np.abs(fast - pinned) / pinned)) < FLOOR / 1000.0


def test_the_integration_grid_is_converged(ladder, monkeypatch):
    fit, _ = ladder
    data = upr.upre_activity_fold(FIT_EXPORT, "UPRE1")
    coarse = upr.PlateUpr().predict(data, fit.parameters)
    monkeypatch.setattr(upr, "INTERNAL_GRID_POINTS", 4 * upr.INTERNAL_GRID_POINTS)
    fine = upr.PlateUpr().predict(data, fit.parameters)
    assert float(np.max(np.abs(coarse - fine) / fine)) < FLOOR / 100.0


def test_the_late_window_is_a_convention_the_prediction_does_not_hinge_on(ladder,
                                                                          monkeypatch):
    """The one place this module's window can differ from the committed pipeline's, priced
    by sweeping it rather than argued about."""
    fit, _ = ladder
    data = upr.upre_activity_fold(FIT_EXPORT, "UPRE1")
    reference = upr.PlateUpr().predict(data, fit.parameters)
    for fraction in (0.60, 0.70, 0.85):
        monkeypatch.setattr(upr, "LATE_WINDOW_FRACTION", fraction)
        moved = upr.PlateUpr().predict(data, fit.parameters)
        assert float(np.max(np.abs(moved - reference) / reference)) < FLOOR / 20.0


def test_the_late_window_matches_the_script_that_produced_the_column():
    assert upr.LATE_WINDOW_FRACTION == 0.75
    script = (pathlib.Path(__file__).resolve().parents[1]
              / "scripts" / "run_sensor_characterisation.py").read_text()
    assert "LATE_WINDOW_FRACTION = 0.75" in script
    assert upr._PLATE_KINETICS.k_deg == 0.0
    assert not upr._PLATE_KINETICS.has_maturation
    assert "ReporterKinetics(k_deg=0.0)" in script


def test_the_declared_starting_points_are_a_literal_not_a_search():
    assert len(upr.STARTING_POINTS) == 2
    for start in upr.STARTING_POINTS:
        assert len(start) == len(upr.PlateUpr().parameter_names)
        low, high = upr.PlateUpr().bounds()
        assert all(lo <= v <= hi for v, lo, hi in zip(start, low, high))


def test_a_frame_without_measured_growth_is_refused():
    data = upr.upre_activity_fold(FIT_EXPORT, "UPRE1").drop(columns=["mu_per_h"])
    with pytest.raises(KeyError, match="MEASURED"):
        upr.PlateUpr().refit(data)


def test_the_provenance_table_renders_both_registries_and_their_gates():
    text = upr.provenance()
    assert "six-state reference chain" in text and "plate-window shape" in text
    assert "REFUSED" in text and "PASSES" in text
    assert "hac1_upre2_kd" in text


# --------------------------------------------------------------------------------------
# The adversarial pass: what the held-out column does NOT show, and three prose repairs
# --------------------------------------------------------------------------------------

def _hill_of_dose(x, dose):
    """The incumbent shape ARCHITECTURE_GAPS.md 0.5 warns about: a saturating function of
    dose alone, with no state, no mu and no measured constant in it."""
    amplitude, half, order = x
    return 1.0 + amplitude * dose ** order / (half ** order + dose ** order)


def _fit_hill_of_dose(dose, observed):
    from scipy.optimize import least_squares
    solution = least_squares(
        lambda x: (_hill_of_dose(x, dose) - observed) / observed, [1.0, 0.3, 2.0],
        bounds=([0.0, 1e-3, 0.5], [100.0, 100.0, 20.0]), xtol=1e-14, ftol=1e-14)
    return solution.x


def test_a_bare_hill_in_dose_matches_the_whole_held_out_column(ladder):
    """CRITERION (c), SCORED AGAINST THE INCUMBENT IT HAS TO BEAT AND NOT ONLY AGAINST ITSELF.

    Three free scalars against this block's four, no chaperone state, no measured Kd, no mu.
    If it tracks the block across every block, then the held-out column measures that a dose
    shape transfers between plates -- not that this mechanism is the one producing it."""
    fit, scores = ladder
    model = upr.PlateUpr()
    data = upr.upre_activity_fold(FIT_EXPORT, "UPRE1")
    hill = _fit_hill_of_dose(np.asarray(data.dose_mM, dtype=float),
                             np.asarray(data.signal, dtype=float))
    assert len(hill) == 3 < len(model.parameter_names)

    gaps, better = [], 0
    for score in scores:
        rows = upr.upre_activity_fold(score.export, score.construct)
        observed = np.asarray(rows.signal, dtype=float)
        dose = np.asarray(rows.dose_mM, dtype=float)
        mechanistic = upr._relative_rms(model.predict(rows, fit.parameters), observed)
        algebraic = upr._relative_rms(_hill_of_dose(hill, dose), observed)
        gaps.append(abs(mechanistic - algebraic))
        if score.held_out and algebraic < mechanistic:
            better += 1
    assert len(gaps) == 8
    assert max(gaps) < FLOOR / 40.0, max(gaps)
    assert better == 4


def test_the_shipped_shape_is_a_saturating_function_of_dose_to_below_the_floor(ladder):
    """The same statement from the model's side rather than the data's: the mechanism's own
    dose-response curve is a Hill. This is what makes criterion (a)'s failure structural and
    not an artefact of the summary statistic it was scored on."""
    fit, _ = ladder
    model = upr.PlateUpr()
    dose = np.linspace(1e-3, 1.0, 400)
    frame = pd.DataFrame({"dose_mM": dose, "signal": np.ones_like(dose),
                          "mu_per_h": np.full_like(dose, 0.214),
                          "mu_control_per_h": np.full_like(dose, 0.214)})
    curve = model.predict(frame, fit.parameters)
    hill = _fit_hill_of_dose(dose, curve)
    assert upr._relative_rms(_hill_of_dose(hill, dose), curve) < FLOOR / 20.0

    quoted = upr.score_dtt_ladder.__doc__
    assert "0.022x the assay floor" in quoted and "0.038x the floor" in quoted


def test_the_measured_growth_rate_is_not_what_makes_a_plate_held_out(ladder):
    """PLATE_REPLICATE_HELD_OUT's stated reason, measured. mu is MEASURED and travels with
    the rows, but across the whole committed range it moves the prediction by a tenth of the
    floor -- so the second plate is a prediction because it is different wells."""
    fit, scores = ladder
    rates = [score.growth_rate_per_h for score in scores]
    assert min(rates) < 0.21 and max(rates) > 0.40
    model = upr.PlateUpr()
    dose = np.array([0.1, 0.2, 0.5, 1.0])

    def at(mu):
        return model.predict(pd.DataFrame({
            "dose_mM": dose, "signal": np.ones(4), "mu_per_h": np.full(4, mu),
            "mu_control_per_h": np.full(4, mu)}), fit.parameters)

    moved = float(np.max(np.abs(at(min(rates)) - at(max(rates))) / at(max(rates))))
    assert moved < FLOOR / 8.0, moved
    assert f"{moved:.4f}" == "0.0155"
    assert "0.106x the assay floor" in upr.PLATE_REPLICATE_HELD_OUT.source


def test_the_kar2_contradiction_is_quoted_rate_against_rate(ladder):
    """A half-life divided by a mean lifetime is 27/2 = 13.5x and it is the wrong comparison.
    Christiano's 27.0 h is a HALF-life; Pincus's g[B] = 1.39e-4 /s is a rate."""
    christiano = float(upr.KAR2_PROTEIN_DECAY_PER_H)
    pincus = 1.39e-4 * 3600.0
    assert christiano == pytest.approx(math.log(2.0) / 27.0)
    assert pincus == pytest.approx(0.5004, abs=1e-4)
    assert pincus / christiano == pytest.approx(19.5, abs=0.05)
    assert 27.0 / 2.0 == 13.5
    for text in (upr.KAR2_PROTEIN_DECAY_PER_H.provenance(),
                 upr.DEACTIVATION_MECHANISM.provenance()):
        assert "19.5x" in text
    share = christiano / (0.28 + christiano)
    assert f"{share:.1%}" == "8.4%"
    assert "8.4% of the loss" in upr.KAR2_PROTEIN_DECAY_PER_H.provenance()


def test_which_end_of_ho_s_hac1_range_is_taken_does_not_move_the_ladder(ladder):
    """The reason the TOP of Ho 2018's range is admissible. NOT 'it is an induced abundance'
    -- this project's own verification pass says none of the five datasets is -- but the
    degeneracy that also frees the compartment: theta enters only as (b+theta)/(b+theta0)."""
    _, scores = ladder
    baseline = np.array([score.relative_rms for score in scores])
    original = upr._P_MAX
    moved = {}
    try:
        for count in (2069.0, 1778.0):
            nanomolar = count / 6.02214076e23 / (42.0 * 1e-15) * 1e9
            ceiling = nanomolar / (float(upr.HAC1_UPRE2_KD_NM) + nanomolar)
            upr._P_MAX = ceiling / (1.0 - ceiling)
            _, refit = upr.score_dtt_ladder(FIT_EXPORT, "UPRE1", held_out=HELD_OUT)
            moved[count] = float(np.max(np.abs(
                np.array([score.relative_rms for score in refit]) - baseline)))
    finally:
        upr._P_MAX = original
    assert moved[2069.0] == pytest.approx(0.0044, abs=5e-4)
    assert moved[1778.0] == pytest.approx(0.00024, abs=5e-5)
    assert max(moved.values()) < FLOOR / 30.0
    provenance = upr.HAC1_ABUNDANCE_MOLECULES.provenance()
    assert "0.030x the assay floor" in provenance
    assert "nearly absent" in provenance


def test_the_refused_entry_still_leaves_a_michaelis_menten_in_dose():
    """What DTT_ENTRY refuses is the SCALE, not the FORM. The chain's capacity loss is
    identically D*dose/(1+D*dose); the module says so rather than letting a reader assume
    the refusal removed the saturating function ARCHITECTURE_GAPS.md 0.5 warns about."""
    basal, top, top_dose = 0.2424, 0.99, 1.0
    entry = upr.DoseEntry(top_load=top, top_dose_mM=top_dose, declared=True,
                          note="a declared axis, not a measurement")
    strength = (top - basal) / (basal * top_dose)
    for dose in (0.0, 0.1, 0.2, 0.5, 1.0, 1.5):
        load = entry.client_load(dose, basal)
        loss = 0.0 if load <= basal else 1.0 - basal / load
        assert loss == pytest.approx(strength * dose / (1.0 + strength * dose), abs=1e-12)
    assert strength == pytest.approx(3.084, abs=5e-3)
    assert "Michaelis-Menten" in upr.DoseEntry.__doc__
    assert upr.DTT_ENTRY.tag == Tag.REFUSED


def test_criterion_a_fails_at_every_rung_not_only_the_saturated_one():
    """The default ablation dose is where both models sit against HAC1_MAX_OCCUPANCY, which
    is the one place a difference could not show. Run the whole ladder: it fails everywhere,
    and the worst case for the block is the BOTTOM rung."""
    scored = {}
    for dose in (0.1, 0.2, 0.5, 1.0):
        result = upr.upr_ablation(dose_mM=dose)
        scored[dose] = result.effect / result.floor.noise_floor
    assert max(scored.values()) < 1.0 / 100.0
    assert max(scored, key=scored.get) == 0.1
    assert scored[0.1] == pytest.approx(0.0062, abs=5e-4)
    assert scored[1.0] == pytest.approx(0.0029, abs=5e-4)
    assert "0.0062x the floor at 0.1 mM" in upr.upr_ablation.__doc__
