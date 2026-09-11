"""Move 5, encoded as a check rather than a note in a document.

A latent module may only be mapped to an FBA constraint if the model reproduces the
fluxes that module's anchor is measured in. Both available models miss respiration
by more than half, so an oxygen- or OUR-linked module has nothing to validate
against and must be refused until that changes.

Written as a guard because a deferral recorded only in prose gets forgotten the
next time someone adds a module.
"""

import pytest

from ystwin.gates.module_admission import (
    ModuleSpec,
    admit_module,
    ModuleRefused,
)
from ystwin.fba.physiology import PhysiologyReport, REFERENCE_AEROBIC_BATCH

FERMENTATIVE_OK = PhysiologyReport(
    reference=REFERENCE_AEROBIC_BATCH,
    observed={"growth_rate": 0.377, "glucose": 17.9, "ethanol": 29.6},
    expected={"growth_rate": 0.40, "glucose": 21.3, "ethanol": 27.4},
    tolerance=0.35,
    relative_error={"growth_rate": 0.058, "glucose": 0.158, "ethanol": 0.079},
)
RESPIRATORY_BAD = PhysiologyReport(
    reference=REFERENCE_AEROBIC_BATCH,
    observed={"oxygen": 2.70, "co2": 32.99},
    expected={"oxygen": 7.80, "co2": 20.40},
    tolerance=0.35,
    relative_error={"oxygen": 0.654, "co2": 0.617},
)

PATHWAY_MODULE = ModuleSpec(
    name="pathway_capacity",
    anchor_fluxes=("growth_rate", "ethanol"),
    constraint="upper bound on the crtYB/crtI pathway",
)
REDOX_MODULE = ModuleSpec(
    name="atp_redox_burden",
    anchor_fluxes=("oxygen", "co2"),
    constraint="bounded respiratory capacity",
)


def test_a_module_whose_anchors_the_model_reproduces_is_admitted():
    decision = admit_module(PATHWAY_MODULE, [FERMENTATIVE_OK])

    assert decision.admitted


def test_an_oxygen_linked_module_is_refused_against_these_models():
    decision = admit_module(REDOX_MODULE, [FERMENTATIVE_OK, RESPIRATORY_BAD])

    assert not decision.admitted
    assert "oxygen" in decision.reason


def test_the_refusal_names_the_flux_and_its_error():
    decision = admit_module(REDOX_MODULE, [RESPIRATORY_BAD])

    assert "65" in decision.reason or "0.65" in decision.reason


def test_a_module_whose_anchor_the_model_never_reports_is_refused():
    unknown = ModuleSpec("mystery", anchor_fluxes=("acetate",), constraint="something")

    decision = admit_module(unknown, [FERMENTATIVE_OK])

    assert not decision.admitted
    assert "acetate" in decision.reason


def test_building_a_constraint_from_a_refused_module_raises():
    with pytest.raises(ModuleRefused, match="atp_redox_burden"):
        admit_module(REDOX_MODULE, [RESPIRATORY_BAD]).require()


def test_an_admitted_module_passes_through_require():
    assert admit_module(PATHWAY_MODULE, [FERMENTATIVE_OK]).require() is PATHWAY_MODULE


def test_a_module_with_no_anchor_flux_at_all_is_refused():
    """A latent state with nothing to validate against is not admissible."""
    decision = admit_module(ModuleSpec("vibes", (), "a multiplier"), [FERMENTATIVE_OK])

    assert not decision.admitted
    assert "no anchor" in decision.reason


def test_the_tolerance_can_be_tightened_for_a_confirmatory_claim():
    """Its anchors sit at 5.8% (growth) and 7.9% (ethanol); 7% refuses the second."""
    lenient = admit_module(PATHWAY_MODULE, [FERMENTATIVE_OK], max_relative_error=0.10)
    strict = admit_module(PATHWAY_MODULE, [FERMENTATIVE_OK], max_relative_error=0.07)

    assert lenient.admitted
    assert not strict.admitted
    assert "ethanol" in strict.reason


def test_a_flux_the_module_does_not_use_cannot_refuse_it():
    """glucose is off by 15.8% but is not an anchor of this module."""
    decision = admit_module(PATHWAY_MODULE, [FERMENTATIVE_OK], max_relative_error=0.12)

    assert decision.admitted
