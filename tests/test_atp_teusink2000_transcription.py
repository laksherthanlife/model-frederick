"""Does the Teusink 2000 transcription reproduce the source's own published numbers?

A transcription nobody checked against the source is not evidence. Every number asserted in
``mech/atp_teusink2000.py``'s docstring is recomputed here, including the two that do NOT
match the printed article and the four quantities the source does not contain.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from ystwin.mech import atp_teusink2000 as teusink
from ystwin.mech.params import RefusedValue, Tag

ROOT = Path(__file__).resolve().parents[1]
_ASSET = ROOT / teusink.SOURCE_PATH
requires_asset = pytest.mark.skipif(
    not _ASSET.is_file(),
    reason=f"{teusink.SOURCE_PATH} absent; it is CC0 and refetchable from BIOMD0000000064")


@pytest.fixture(scope="module")
def steady():
    return teusink.steady_state()


@requires_asset
def test_the_pinned_source_is_the_one_that_was_transcribed():
    assert teusink.verify_source(ROOT) == teusink.SOURCE_DIGEST


@requires_asset
def test_the_repository_loader_still_refuses_the_deposit():
    """Condition (a): the file is transcribed BECAUSE the loader refuses it. If this test
    ever fails, the loader was loosened and the transcription is no longer the reason."""
    from ystwin.mech.kinetic_sbml import KineticModel, UnsupportedSBMLError

    with pytest.raises(UnsupportedSBMLError) as raised:
        KineticModel.from_sbml(str(_ASSET))
    assert "must be a nonconstant global parameter" in str(raised.value)


def test_the_adenylate_pool_reproduces_table_4(steady):
    """THE published value this build exists for: ATP 2.51 / ADP 1.29 / AMP 0.30 mM."""
    assert steady.atp == pytest.approx(2.51, abs=0.005)
    assert steady.adp == pytest.approx(1.29, abs=0.005)
    assert steady.amp == pytest.approx(0.30, abs=0.005)
    assert steady.energy_charge == pytest.approx(0.769, abs=0.001)


def test_it_is_a_steady_state_and_not_a_horizon_artefact(steady):
    assert steady.residual < 1e-9


def test_every_table_4_concentration_reproduces_to_the_precision_it_is_printed_at(steady):
    """Table 4 prints two decimals, so the test is at two decimals. Reporting 0.04444
    against a printed 0.04 as an 11% error would be scoring the article's print precision."""
    rows = steady.table4()
    for name, published in teusink.PUBLISHED_TABLE4.items():
        if name == "G6P":
            continue
        assert round(rows[name], 2) == pytest.approx(published, abs=0.005), name


def test_the_one_row_that_disagrees_is_the_deposit_s_own_declared_divergence(steady):
    """The CC0 encoding and the printed article disagree on G6P and the curator says so.
    This transcription must reproduce the ENCODING, and inherit the gap, not close it."""
    g6p = steady.concentrations["G6P"]
    assert g6p == pytest.approx(teusink.DEPOSIT_G6P, abs=0.005)
    assert not math.isclose(g6p, teusink.PUBLISHED_TABLE4["G6P"], abs_tol=0.01)


def test_published_fluxes_reproduce(steady):
    for name, published in teusink.PUBLISHED_FLUXES.items():
        assert steady.fluxes[name] == pytest.approx(published, rel=0.02), name


def test_the_assignment_rules_are_pure_algebra_and_conserve_the_moiety():
    """The three species-targeted rules the loader refuses, checked as algebra."""
    atp, adp, amp = teusink.adenylate(6.31)
    assert atp + adp + amp == pytest.approx(float(teusink.SUM_P), rel=1e-12)
    assert 2.0 * atp + adp == pytest.approx(6.31, rel=1e-12)
    # KeqAK is ATP*AMP/ADP^2, i.e. written for 2 ADP -> ATP + AMP. The engine writes the
    # same reaction the other way round, so the direction is pinned here deliberately.
    assert (atp * amp) / (adp ** 2) == pytest.approx(float(teusink.KEQ_AK), rel=1e-9)


def test_glucose_starvation_moves_the_pool_monotonically():
    """The one panel-dosable arm. Not a claim that the panel's dose is representable."""
    charges = [teusink.steady_state(glco=g).energy_charge for g in (50.0, 20.0, 5.0, 2.0)]
    assert charges == sorted(charges, reverse=True)
    assert charges[0] == pytest.approx(0.769, abs=0.001)
    assert charges[-1] < 0.45


def test_below_the_validated_floor_it_refuses_instead_of_returning_a_number():
    """The source does not settle below ~1.5 mM glucose. An unconverged number would be
    worse than a refusal, so the deep-starvation regime is refused outright."""
    with pytest.raises(ValueError, match="validated floor"):
        teusink.steady_state(glco=1.0)


@pytest.mark.parametrize("name", [
    "teusink.stress_to_atp_demand",
    "teusink.cytosol_l_per_gdw",
    "teusink.respiratory_atp_yield",
    "teusink.growth_dilution",
])
def test_the_four_absent_quantities_raise_rather_than_defaulting(name):
    param = teusink.TEUSINK_PARAMS[name]
    assert param.tag == Tag.REFUSED
    assert param.value is None
    with pytest.raises(RefusedValue):
        float(param)


def test_the_mM_to_gdw_bridge_cannot_be_used_by_accident():
    """Condition (b): the conversion runs through an asserted engine prior, so it refuses."""
    with pytest.raises(RefusedValue):
        teusink.to_mmol_per_gdw(2.51)


def test_no_constant_is_graded_above_asserted():
    """Following ph_ke2013, which grades all 78 of Ke's transcribed constants ASSERTED. The
    article was not retrievable, so no row may claim a measurement it was not checked against."""
    tags = {tag: count for tag, count in teusink.TEUSINK_PARAMS.by_tag().items() if count}
    assert set(tags) == {Tag.ASSERTED, Tag.REFUSED}, tags
    assert tags[Tag.REFUSED] == 4


def test_the_deposit_s_valued_parameter_count_is_all_present():
    """85 valued parameters (70 reaction-local + 15 global), plus SUM_P and F26BP."""
    valued = [p for p in teusink.TEUSINK_PARAMS.params if p.tag == Tag.ASSERTED]
    assert len(valued) == 87
    assert all(p.value is not None for p in valued)
    assert all("BIOMD0000000064" in p.source for p in valued)


def test_the_gate_refuses_and_says_why():
    """Reproducing a source's own model output scores the port, not the biology, so no
    independent target is registered and the free-scalar gate must REFUSE."""
    verdict = teusink.gate()
    assert not verdict.passes
    assert verdict.targets == ()
    assert len(verdict.free) == 91


def test_provenance_never_claims_upstream_availability():
    record = teusink.provenance()
    assert record["availability"] == "transcribed_here"
    assert "CC0" in record["licence"]
    assert record["source_sha256"] == teusink.SOURCE_DIGEST
    assert set(record["transfer_assumptions"]) == {
        "non_growing_to_fed_batch", "fermentative_to_aerobic", "mM_to_mmol_per_gdw",
        "no_camp_reservoir", "adenylate_kinase_equilibrium"}


def test_the_panel_limit_is_stated_rather_than_implied():
    assert "antimycin_A" in teusink.PANEL_LIMIT
    assert "does not close the atp module" in teusink.PANEL_LIMIT
