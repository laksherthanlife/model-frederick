"""The GSMM must reproduce a measured phenotype before any bound it yields means anything.

Reference: S. cerevisiae CEN.PK113-7D, aerobic glucose-excess batch
(van Hoek, van Dijken & Pronk 1998, Appl Environ Microbiol 64:4226):
    mu = 0.40 /h, qGlc = 21.3, qEtOH = 27.4, qO2 = 7.8, qCO2 = 20.4 mmol/gDW/h.

Plain Yeast9 cannot hit this even when handed the measured uptake rates: FBA takes
the optimal-yield route that real yeast does not. The enzyme-constrained model
recovers growth, glucose and ethanol from its protein pool alone -- but still
misses respiration badly, which bounds what any oxygen- or energy-linked claim
can rest on.
"""

import pytest

from ystwin.fba.physiology import (
    REFERENCE_AEROBIC_BATCH,
    aerobic_batch_constraints,
    validate_physiology,
)

pytestmark = pytest.mark.integration

FERMENTATIVE = ("growth_rate", "glucose", "ethanol")
RESPIRATORY = ("oxygen", "co2")


def test_free_oxygen_abolishes_overflow_metabolism_entirely(yeast_gem):
    with yeast_gem as m:
        m.reactions.get_by_id("r_1714").lower_bound = -REFERENCE_AEROBIC_BATCH.glucose_uptake
        sol = m.optimize()

    assert sol.fluxes["r_1761"] == pytest.approx(0.0, abs=1e-6)
    # Roughly 2.8x the oxygen-limited growth below, which is the point: respiration is
    # far more efficient and the model takes it whenever oxygen allows.
    assert sol.objective_value > 0.9


def test_plain_fba_underpredicts_growth_given_the_measured_uptake_rates(yeast_gem):
    """This asserted the opposite, and the reversal is a citation fix, not a model change.

    The reference used to carry glucose 21.3 and ethanol 27.4 attributed to van Hoek 1998.
    Those strings occur zero times in that paper, which ran no batch culture and used
    strain DS28911. With its actual Table 1 values at D = 0.40 -- glucose 11.1, oxygen 3.7
    -- plain Yeast9 grows at 0.346 against a measured 0.40. It undershoots by 13.5%; the
    docs reported it overshooting by 73%.
    """
    with yeast_gem as m:
        aerobic_batch_constraints(m)
        sol = m.optimize()

    assert sol.fluxes["r_1761"] > 5.0, "overflow should reappear once oxygen is limited"
    assert sol.objective_value < 0.40, "and growth should now fall short of the measured rate"
    assert sol.objective_value == pytest.approx(0.346, abs=0.01)


def test_the_enzyme_constrained_model_reaches_the_right_growth_by_the_wrong_route(ec_yeast_gem):
    """Its protein pool alone puts growth within 5.8% of the measured rate, which is the
    reason this model is used here. But left unconstrained it gets there on 61% more
    glucose and 113% more ethanol than the culture it is compared against -- it is finding
    a glucose-excess optimum against a glucose-limited chemostat. The growth agreement is
    real and it is not evidence that the flux distribution is."""
    report = validate_physiology(ec_yeast_gem, apply_constraints=False, keys=FERMENTATIVE)

    assert report.relative_error["growth_rate"] < 0.10
    assert report.relative_error["glucose"] > 0.5
    assert report.relative_error["ethanol"] > 1.0
    assert not report.passed


def test_the_enzyme_constrained_model_needs_no_hand_set_uptake_bounds(ec_yeast_gem):
    """Overflow emerges from the protein pool, not from a fitted oxygen bound."""
    report = validate_physiology(ec_yeast_gem, apply_constraints=False, keys=FERMENTATIVE)

    assert report.observed["ethanol"] > 20.0


def test_the_enzyme_constrained_model_still_misses_respiration(ec_yeast_gem):
    """The documented miss was 65% on oxygen. It is 27% against the corrected reference,
    and CO2 is still far out. Two artefacts inflated the old figure: the reference fluxes
    were unsourced and about twice the paper's own, and cap_uptake was a no-op on this
    model because it is irreversibly split, so oxygen was never actually constrained."""
    report = validate_physiology(ec_yeast_gem, apply_constraints=False, keys=RESPIRATORY)

    assert not report.passed
    assert 0.2 < report.relative_error["oxygen"] < 0.4
    assert report.relative_error["co2"] > 0.5


def test_validation_reports_every_requested_flux_with_its_relative_error(ec_yeast_gem):
    report = validate_physiology(ec_yeast_gem, apply_constraints=False)

    assert set(report.relative_error) == {"growth_rate", "glucose", "oxygen", "ethanol", "co2"}


def test_validation_is_read_only(yeast_gem):
    before = yeast_gem.reactions.get_by_id("r_1992").bounds
    validate_physiology(yeast_gem)

    assert yeast_gem.reactions.get_by_id("r_1992").bounds == before


def test_an_unknown_reference_key_is_rejected(ec_yeast_gem):
    with pytest.raises(KeyError, match="acetate"):
        validate_physiology(ec_yeast_gem, keys=("acetate",))


@pytest.mark.parametrize("fixture_name", ["yeast_gem", "ec_yeast_gem"])
def test_physiology_is_invariant_to_exchange_stoichiometry_rescaling(request, fixture_name):
    model = request.getfixturevalue(fixture_name)
    reference = validate_physiology(model)

    with model:
        for base in ("r_1714", "r_1992", "r_1761", "r_1672"):
            for rid in (base, f"{base}_REV"):
                if rid not in model.reactions:
                    continue
                reaction = model.reactions.get_by_id(rid)
                reaction.add_metabolites({
                    metabolite: coefficient * 2.0
                    for metabolite, coefficient in reaction.metabolites.items()
                }, combine=False)
                reaction.bounds = tuple(bound / 2.0 for bound in reaction.bounds)

        observed = validate_physiology(model)

    assert observed.observed == pytest.approx(reference.observed, rel=1e-6, abs=1e-6)
