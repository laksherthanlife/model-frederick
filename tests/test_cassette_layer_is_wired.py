"""The cassette layer reports what `pathway/capacity.py` returned, not what somebody typed.

`tests/test_everything_is_wired.py` pins that every named layer is REACHABLE and that only
the layers which may move the number do. This file pins the thing one layer down, which that
one cannot see: whether a layer's status was COMPUTED or WRITTEN.

The cassette layer was written. Until 2026-09-02 both branches of `_cassette_layer` returned
a hand-built `LayerStatus`, and the string NAMED the module it was speaking for --
"pathway/capacity.py names crtYB as the enzyme that sets the ceiling and refuses to scale by
it" -- while `predict.py` imported `.pathway.flux`, `.pathway.solve`, `.pathway.spec` and
`.pathway.thermo_gate` and NOT `.pathway.capacity`, and `capacity_for` was called from
nowhere in `src/` or `scripts/` at all. The layer was a quotation. A quotation is right until
the quoted module changes, and nothing in the suite would have noticed the day it did --
which is the same failure as a docstring reporting what the code would do, in the one
mechanism this package built to stop exactly that.

What these tests pin, and what breaks if each stops holding:

1. **`predict.py` imports capacity.py and calls it.** If `_cassette_layer` goes back to a
   literal, the detail can no longer carry a number the call returned and these fail.

2. **Both of capacity.py's answers are reachable from a prediction.** A cassette it can
   speak to comes back INERT with the calibration strain's own capacity; one it cannot comes
   back REPORTED carrying `CapacityUnmeasured` verbatim. A refusal nothing can reach is a
   refusal nobody is being told about.

3. **Wiring it in did not let the cassette move the number, and must not.** capacity.py
   ACCEPTS a cassette and REFUSES to scale by it on purpose: the flux law is fitted on crtE,
   the ceiling is set at crtYB, and the coefficient between cyclase dosage and capacity does
   not exist in the literature. A change that makes any of the contents below differ is a
   Tier 0 constant somebody invented, and `==` rather than `approx` is the assertion so that
   no contribution is small enough to slip past.

Nothing here needs a GEM, a thermodynamic table or a plate.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS
from ystwin.pathway.capacity import CapacityUnmeasured, capacity_for
from ystwin.pathway.flux import FluxCalibration
from ystwin.pathway.solve import content_ceiling
from ystwin.pathway.spec import load_pathway
from ystwin.predict import (
    Environment, Genotype, LayerState, LegacyProductDiagnostic, predict_product,
)

#: A cassette naming only the enzyme the FLUX law is fitted on. capacity.py can answer it,
#: and the answer is the anchor unchanged -- raising crtE does not raise the ceiling, which
#: is the whole diagnosis.
FLUX_ENZYME_ONLY = {"crtE": 10.0}

#: The calibration strain itself, stated explicitly. Answerable by definition.
CALIBRATION_STRAIN = {"crtYB": 1.0}

#: The cyclase at a dosage nothing measured. This is the one capacity.py refuses.
UNMEASURED_DOSAGE = {"crtYB": 4.0}

REPO = pathlib.Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def spec():
    return load_pathway("beta_carotene")


@pytest.fixture(scope="module")
def phb_spec():
    return load_pathway("phb")


@pytest.fixture(scope="module")
def phb_calibration():
    """Illustrative, and it does not matter: nothing here reads the flux's magnitude.

    Same fixture as `tests/test_everything_is_wired.py`, and for the same reason -- `phb`
    ships no calibration because Kocharin 2012 measures no expression.
    """
    return FluxCalibration(alpha=9.93e-2, entry_enzyme="PhaA", loso_rmse_log=0.2,
                           loso_skill=0.0, n_states=2, source="illustrative, not fitted")


def _predict(spec, cassette=None, **kwargs):
    return predict_product(spec, Genotype(1.0, "b-car4", cassette=cassette),
                           Environment(growth_rate_setpoint_per_h=0.18),
                           BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS, **kwargs)


def _legacy(spec, cassette=None):
    result = _predict(spec, cassette, mode="legacy")
    assert isinstance(result, LegacyProductDiagnostic)
    assert result.supported is False
    assert result.calculation.mode == "legacy"
    return result


def _layer(spec, cassette=None):
    return _predict(spec, cassette).layer("cassette dosage")


def _legacy_layer(spec, cassette=None):
    return _legacy(spec, cassette).calculation.layer("cassette dosage")


class TestPredictAsksCapacityPyRatherThanQuotingIt:
    def test_predict_imports_the_module_its_status_line_names(self):
        """The defect in one line. `predict.py` printed "pathway/capacity.py names crtYB..."
        out of a literal while importing four `pathway` modules and not that one, so the
        sentence and the code behind it could disagree with nothing to catch it."""
        source = (REPO / "src" / "ystwin" / "predict.py").read_text(encoding="utf-8")
        imported = {node.module for node in ast.walk(ast.parse(source))
                    if isinstance(node, ast.ImportFrom)}

        assert "pathway.capacity" in imported

    def test_the_capacity_in_the_detail_is_the_one_capacity_for_returned(self, spec):
        """Derived, not typed. The anchor is THIS call's own ceiling rather than a constant
        copied into `predict.py`, so a recalibration moves the detail with it."""
        anchor = content_ceiling(spec, BETA_CAROTENE_KINETICS)
        answered = capacity_for(spec.product, anchor, None)

        assert f"{answered:.4g}" in _layer(spec).detail

    def test_a_cassette_capacity_py_can_speak_to_is_inert(self, spec):
        """INERT is the enum's "ran, and provably could not have changed anything". The
        capacity came back equal to the anchor, so nothing could have moved."""
        got = _layer(spec, CALIBRATION_STRAIN)

        assert got.state == LayerState.INERT
        assert "crtYB" in got.detail

    def test_a_gene_it_does_not_name_is_answered_rather_than_refused(self, spec):
        """crtE at ten times dosage is answerable and the answer is the anchor. Refusing it
        would be a guess in the other direction -- capacity.py refuses one enzyme, not every
        cassette."""
        got = _legacy_layer(spec, FLUX_ENZYME_ONLY)

        assert got.state == LayerState.INERT
        assert "crtE 10x" in got.detail

    def test_a_dosage_nothing_measured_comes_back_reported(self, spec):
        """REPORTED is the enum's state for output "attached and labelled but feeds
        nothing... the state that keeps a refused layer visible instead of deleted", which
        is what a `CapacityUnmeasured` is. Same treatment the G4-refused latent branch
        gets."""
        assert _legacy_layer(spec, UNMEASURED_DOSAGE).state == LayerState.REPORTED

    def test_the_two_cassettes_cannot_both_come_from_one_literal(self, spec):
        """**The check can actually fail.** A hardcoded `LayerStatus` cannot do this: the
        same layer must report a different state AND a different string for a cassette
        capacity.py can answer and one it refuses. If someone replaces the call with a
        literal again, one of these two collapses onto the other."""
        answered = _layer(spec, CALIBRATION_STRAIN)
        refused = _legacy_layer(spec, UNMEASURED_DOSAGE)

        assert answered.state != refused.state
        assert answered.detail != refused.detail


class TestTheRefusalIsReachableAndVerbatim:
    def test_the_detail_is_what_capacity_for_raised(self, spec):
        """Verbatim, not paraphrased. A paraphrase is a second copy of the refusal that
        drifts from the first, which is the defect this whole file is about."""
        anchor = content_ceiling(spec, BETA_CAROTENE_KINETICS)
        with pytest.raises(CapacityUnmeasured) as raised:
            capacity_for(spec.product, anchor, UNMEASURED_DOSAGE)
        collapsed = " ".join(str(raised.value).split())

        assert collapsed in _legacy_layer(spec, UNMEASURED_DOSAGE).detail

    def test_and_it_names_the_experiment_that_would_settle_it(self, spec):
        """The reason the refusal is worth carrying at all: it names a measurement rather
        than only declining. A layer that said "unmeasured" and stopped would be worse than
        the literal it replaced."""
        detail = _legacy_layer(spec, UNMEASURED_DOSAGE).detail

        assert "WHAT WOULD SETTLE IT" in detail
        assert "dosage ALONE" in detail
        assert "39215465" in detail

    def test_it_says_the_number_above_it_is_still_the_calibration_strains(self, spec):
        """The honest half of catching the refusal instead of raising it. The content is
        still solved against the anchor, and a caller who reads the number without reading
        that has been misled."""
        assert "still solved against the calibration strain's capacity" in _legacy_layer(
            spec, UNMEASURED_DOSAGE).detail

    def test_the_refusal_reuses_an_existing_state_rather_than_inventing_a_sixth(self, spec):
        """`LayerState.ALL` is consumed elsewhere -- `LayerStatus.__post_init__` validates
        against it and `tests/test_everything_is_wired.py` reads it -- so a sixth member
        would make every consumer wrong."""
        assert len(LayerState.ALL) == 5
        assert _legacy_layer(spec, UNMEASURED_DOSAGE).state in LayerState.ALL

    def test_the_whole_refusal_still_fits_on_one_line(self, spec):
        """`layer_report` is one line per layer and a multi-paragraph detail would silently
        break that count for every caller printing it."""
        got = _legacy(spec, UNMEASURED_DOSAGE).calculation

        assert len(got.layer_report().splitlines()) == len(got.layers)


class TestWiringItInDoesNotMoveTheNumber:
    """`==`, not `approx`. capacity.py accepts a cassette and refuses to SCALE by it, and a
    tolerance here would let an invented coefficient through unnoticed."""

    def test_a_cassette_it_answers_changes_nothing(self, spec):
        assert (_predict(spec, CALIBRATION_STRAIN).content_mmol_per_gdcw
                == _predict(spec).content_mmol_per_gdcw)

    def test_nor_does_raising_the_flux_enzyme_tenfold(self, spec):
        """The diagnosis, as a prediction: the cassette entry for the enzyme the flux law is
        fitted on does not reach the ceiling either, because `entry_expression` is the field
        that feeds the flux and this mapping feeds nothing."""
        assert (_legacy(spec, FLUX_ENZYME_ONLY).calculation.content_mmol_per_gdcw
                == _predict(spec).content_mmol_per_gdcw)

    def test_nor_does_one_it_refuses(self, spec):
        """A refusal must not become a silent scaling in the other direction either."""
        assert (_legacy(spec, UNMEASURED_DOSAGE).calculation.content_mmol_per_gdcw
                == _predict(spec).content_mmol_per_gdcw)

    def test_nor_the_flux_the_content_is_solved_from(self, spec):
        assert (_legacy(spec, UNMEASURED_DOSAGE).calculation.flux.flux_mmol_per_gdcw_h
                == _predict(spec).flux.flux_mmol_per_gdcw_h)

    def test_the_layer_never_claims_to_set_the_number(self, spec):
        for cassette in (None, FLUX_ENZYME_ONLY, CALIBRATION_STRAIN, UNMEASURED_DOSAGE):
            assert not _legacy_layer(spec, cassette).sets_the_number, cassette

    def test_and_the_set_of_layers_that_do_is_unchanged_by_a_cassette(self, spec):
        """The same set `tests/test_everything_is_wired.py` pins, asserted under the input
        that could plausibly have joined it."""
        setting = {layer.name for layer in _legacy(spec, UNMEASURED_DOSAGE).calculation.layers
                   if layer.sets_the_number}

        assert setting == {"metabolism", "flux law"}


class TestAPathwayWithNoCeilingSaysThereIsNothingToAsk:
    """`phb` run with no kinetics has no saturating step at the end of its chain, so
    `content_ceiling` returns `None` and there is no measured anchor to hand capacity.py.
    Passing a stand-in for one would be the invention this layer exists not to make."""

    def test_the_ceiling_this_rests_on_is_really_absent(self, phb_spec):
        assert content_ceiling(phb_spec, {}) is None

    def test_the_layer_says_so_instead_of_naming_a_capacity(self, phb_spec, phb_calibration):
        got = predict_product(phb_spec, Genotype(1.0, "SCKK006"),
                              Environment(growth_rate_setpoint_per_h=0.1), phb_calibration, {})
        layer = got.layer("cassette dosage")

        assert layer.state == LayerState.INERT
        assert "no saturating step" in layer.detail

    def test_even_when_a_cassette_is_stated(self, phb_spec, phb_calibration):
        """`phb` has no entry in `CAPACITY_SETTING_ENZYME` and no anchor either. Two
        independent reasons to answer nothing, and the layer must not fall through to a
        capacity it does not have."""
        got = predict_product(phb_spec, Genotype(1.0, "SCKK006", cassette={"phaA": 4.0}),
                              Environment(growth_rate_setpoint_per_h=0.1), phb_calibration, {},
                              mode="legacy")
        assert isinstance(got, LegacyProductDiagnostic)
        assert got.supported is False
        layer = got.calculation.layer("cassette dosage")

        assert layer.state == LayerState.INERT
        assert "phaA 4x" in layer.detail

    def test_default_refuses_cassette_dosage_without_a_calibrated_route(self, phb_spec, phb_calibration):
        with pytest.raises(ValueError, match="cassette dosage"):
            predict_product(phb_spec, Genotype(1.0, "SCKK006", cassette={"phaA": 4.0}),
                            Environment(growth_rate_setpoint_per_h=0.1), phb_calibration, {})


def test_default_propagates_the_capacity_refusal_verbatim(spec):
    with pytest.raises(CapacityUnmeasured) as direct:
        capacity_for(spec.product, content_ceiling(spec, BETA_CAROTENE_KINETICS), UNMEASURED_DOSAGE)
    with pytest.raises(CapacityUnmeasured) as public:
        _predict(spec, UNMEASURED_DOSAGE)
    assert str(public.value) == str(direct.value)


def test_default_refuses_copy_number_where_measured_expression_is_required(spec):
    with pytest.raises(ValueError, match="measured expression, not gene copy number"):
        _predict(spec, FLUX_ENZYME_ONLY)
