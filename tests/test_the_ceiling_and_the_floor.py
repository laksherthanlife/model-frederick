"""Two hard limits the model imposes on every prediction, neither of which was visible.

Both were found by an adversarial pass that was told to make the model produce a confidently
wrong answer, and both are properties of the architecture rather than of one calibration.

**The floor: no pool had a physical bound.** At an entry expression of 1e4 the lycopene pool
came back as 36,093 mg/gDCW -- 3,609% of the cell's own dry weight -- as an ordinary float,
with no note and no refusal. The steady-state solve has exactly one positive root and no
opinion about whether that root describes anything.

**The ceiling: beta-carotene content can never exceed 1.2483 mg/gDCW.** For a terminal node
fed by a saturating step, ``content = flux_out/mu`` and ``flux_out <= vmax_per_growth * mu``,
so ``content <= vmax_per_growth`` and **mu cancels exactly**. No genotype at any growth rate
can pass it. That is the sharpest falsifiable claim in the package and nobody had noticed the
model was making it.

The ceiling is not an artefact of one state -- leave-one-out moves it between 1.11 and 1.41
mg/gDCW. It is 17x below Lopez 2019's 21 mg/gDCW, which is a red flag rather than a
refutation: different strain, flask rather than chemostat, total carotenoid rather than
beta-carotene. Somebody should close that gap, and until they do the model is asserting
something about every strain that will ever be built.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ystwin import paths
from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS
from ystwin.pathway.solve import (
    DRY_WEIGHT_G_PER_GDCW,
    ImplausibleContent,
    content_ceiling,
)
from ystwin.pathway.spec import load_pathway
from ystwin.predict import Environment, Genotype, predict_product

MOLAR_MASS = 536.87


@pytest.fixture(scope="module")
def spec():
    return load_pathway("beta_carotene")


def run(spec, expression, rate=0.15):
    return predict_product(spec, Genotype(expression), Environment(growth_rate_setpoint_per_h=rate),
                           BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS)


class TestNoPoolMayOutweighTheCellHoldingIt:
    def test_an_absurd_expression_is_refused_rather_than_answered(self, spec):
        with pytest.raises(ImplausibleContent, match="of the cell's own mass"):
            run(spec, 1e4)

    def test_the_refusal_says_the_arithmetic_is_right_and_the_inputs_are_not(self, spec):
        """Which is true and is the useful thing to tell someone. The solve is exact; what
        is wrong is an entry expression four orders of magnitude outside the fit."""
        with pytest.raises(ImplausibleContent, match="inputs are not"):
            run(spec, 1e4)

    def test_the_bound_is_arithmetic_and_needs_no_citation(self):
        """One gram per gram of dry weight is the cell made entirely of one pool. A tighter
        bound would be a claim about how much product a yeast can hold -- a real and
        contested number -- and asserting one in a solver would smuggle in a Tier 0
        constant. The ceiling is deliberately absurd so that it is uncontroversial."""
        assert DRY_WEIGHT_G_PER_GDCW == 1.0

    def test_an_ordinary_prediction_is_nowhere_near_it(self, spec):
        got = run(spec, 1.0)

        assert got.intermediates()["lycopene"] * MOLAR_MASS < 10.0

    def test_a_pool_with_no_molar_mass_cannot_be_checked_and_is_not(self):
        """Honest gap, stated so it is not mistaken for coverage: the bound is in grams, so a
        node that declares no molar mass is unbounded. Every shipped spec declares one."""
        for name in ("beta_carotene",):
            for node in load_pathway(name).nodes:
                assert node.molar_mass_g_per_mol is not None


class TestTheCeilingThisCalibrationImposes:
    """A property of the MODEL, and known since 2026-08-28 to be false of the organism.

    Every assertion below still holds -- they are statements about what this calibration
    can produce, and the algebra is unchanged. What changed is the claim they were once
    read as supporting. Arhar 2024 (PMID 39215465) measured **79 mg/gDCW** of beta-carotene
    by HPLC on gravimetric dry weight with lycopene below detection, about 63x this ceiling,
    by retargeting CrtYB. Capacity here was fitted flat across strains differing only in the
    ENTRY enzyme, and the cyclase dosage that actually sets it is not a field `Genotype`
    carries -- so the repair is structural and no refit expresses it. See docs/HARD_TESTS.md.

    The class was called `TestTheCeilingNoGenotypeCanPass`. That name asserted the refuted
    thing, so it is gone.
    """

    def test_it_exists_and_is_the_fitted_capacity(self, spec):
        ceiling = content_ceiling(spec, BETA_CAROTENE_KINETICS)

        assert ceiling == BETA_CAROTENE_KINETICS["lycopene"].vmax_per_growth

    def test_growth_rate_cancels_out_of_it(self, spec):
        """The reason it is a hard limit rather than a soft one. If mu did not cancel, a
        slow chemostat would reach higher contents and the ceiling would be per-condition."""
        ceiling = content_ceiling(spec, BETA_CAROTENE_KINETICS) * MOLAR_MASS

        for rate in (0.101, 0.15, 0.20, 0.254):
            assert run(spec, 100.0, rate).content_mg_per_gdcw < ceiling

    def test_no_expression_however_large_passes_it(self, spec):
        ceiling = content_ceiling(spec, BETA_CAROTENE_KINETICS) * MOLAR_MASS
        best = max(run(spec, e).content_mg_per_gdcw for e in (1.0, 10.0, 100.0))

        assert best < ceiling

    def test_and_it_is_close_enough_to_matter(self, spec):
        """1.2476 against a ceiling of 1.2483 at an expression of 100. The limit is not
        academic -- a fourfold increase in cassette dosage buys 0.06%."""
        ceiling = content_ceiling(spec, BETA_CAROTENE_KINETICS) * MOLAR_MASS

        assert run(spec, 100.0).content_mg_per_gdcw > 0.999 * ceiling

    def test_a_prediction_near_it_says_so(self, spec):
        """Without the note, two strains differing by 0.4% at the asymptote read as a
        ranking. They are both just reporting the fitted capacity."""
        near = run(spec, 10.0)


        assert any("THIS CALIBRATION imposes" in note for note in near.notes)

    def test_and_an_ordinary_one_does_not(self, spec):
        assert not any("THIS CALIBRATION imposes" in note
                       for note in run(spec, 1.0).notes)

    def test_the_note_says_the_ceiling_is_refuted_in_reality(self, spec):
        """The note is where a caller meets this number, so it is where the caveat belongs."""
        note = next(n for n in run(spec, 10.0).notes if "THIS CALIBRATION imposes" in n)

        assert "39215465" in note and "79" in note, (
            "the near-ceiling note must name Arhar 2024 and its 79 mg/gDCW, or it tells a "
            "caller the bound is unreachable when it has been exceeded 63-fold")

    def test_a_chain_with_nothing_saturating_has_no_ceiling(self):
        """PHB is passthrough end to end, so nothing limits its content and
        `content_ceiling` must say so rather than inventing one."""
        phb = load_pathway("phb")

        assert content_ceiling(phb, {}) is None


class TestTheCeilingIsNotAnArtefactOfOneState:
    """If it rested on a single measurement it would be a fitting accident. It does not."""

    @pytest.fixture(scope="class")
    def measured(self):
        path = paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv"
        if not path.exists():
            pytest.skip(f"{path} not present")
        return pd.read_csv(path, sep="\t")

    @staticmethod
    def _refit(measured, keep):
        from scipy.optimize import least_squares

        mu = measured.mu_per_h.to_numpy()[keep]
        flux = (measured.q_lycopene + measured.q_betacarotene).to_numpy()[keep]
        observed = measured.q_betacarotene.to_numpy()[keep]

        def predicted(p):
            vmax = p[0] * mu
            a, b, c = mu, mu * p[1] + vmax - flux, -flux * p[1]
            lycopene = (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
            return vmax * lycopene / (p[1] + lycopene)

        fitted = least_squares(lambda p: np.log(predicted(p)) - np.log(observed),
                               x0=[2.3e-3, 6e-4], bounds=([1e-8, 1e-8], [1.0, 1.0]))
        return float(fitted.x[0])

    def test_leave_one_out_moves_it_less_than_twenty_percent(self, measured):
        everything = np.arange(len(measured))
        base = self._refit(measured, everything)

        for held_out in everything:
            refitted = self._refit(measured, np.setdiff1d(everything, [held_out]))

            assert abs(refitted / base - 1.0) < 0.20

    def test_the_calibration_states_do_approach_it(self, measured):
        """So it is set by data rather than extrapolated into empty space -- one state sits
        at 97% of it. That is what makes it a claim rather than an accident."""
        ceiling = BETA_CAROTENE_KINETICS["lycopene"].vmax_per_growth
        content = measured.q_betacarotene / measured.mu_per_h

        assert content.max() / ceiling > 0.9
