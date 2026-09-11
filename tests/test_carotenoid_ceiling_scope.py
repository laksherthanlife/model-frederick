"""The refuted ceiling must not be readable as a number without its refutation.

``pathway/solve.py::content_ceiling`` returned a bare ``float``. Everything that made the
number safe to read -- that it bounds THIS CALIBRATION rather than the organism, that the
enzyme setting it is not the enzyme the flux law reads, that Arhar 2024 measured 63x past it
-- lived in docstrings, in `pathway/capacity.py` and in one `predict.py` note that fires only
within 10% of the bound. None of it travelled with the value. A caller holding the float held
a number and no scope, which is exactly how a refuted bound gets reported as a prediction.

That is the same failure shape `pathway/published_cassettes.py` was written to end: a claim
about a number, copied into prose, with nothing that rechecks it. So the refutation here is
carried as the MEASUREMENT and cross-checked against the vendored survey row, and the fold is
computed from the shipped constant rather than written down.

WHY THE REFIT WAS NOT DONE, asserted below rather than described. Making the capacity
proportional to an enzyme's measured expression is the only refit the six Elizondo states
could support, and every version of it fits WORSE than the flat capacity and still lands
within 1.4x of the same refuted ceiling. See `TestNoRefitOnTheseSixStatesReachesTheMeasurement`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import least_squares

from ystwin import paths
from ystwin.pathway.calibrations import BETA_CAROTENE_KINETICS
from ystwin.pathway.solve import NodeKinetics, content_ceiling
from ystwin.pathway.spec import Fate, Node, PathwaySpec, RateLaw, load_pathway

MOLAR_MASS = 536.87
ARHAR_PMID = "39215465"
ARHAR_CONTENT_MG_PER_GDCW = 79.0


@pytest.fixture(scope="module")
def spec():
    return load_pathway("beta_carotene")


@pytest.fixture(scope="module")
def ceiling(spec):
    return content_ceiling(spec, BETA_CAROTENE_KINETICS)


class TestTheCeilingSaysWhatItIsACeilingOf:
    def test_the_number_itself_is_untouched(self, ceiling):
        """Inertness first. This changes what the value CARRIES, never what it IS."""
        assert ceiling == 2.3252e-3
        assert isinstance(ceiling, float)

    def test_it_names_the_pool_it_bounds(self, ceiling):
        assert ceiling.bounded_node == "beta_carotene"

    def test_and_the_step_whose_capacity_sets_it(self, ceiling):
        assert ceiling.capacity_node == "lycopene"
        assert ceiling.capacity_enzyme == "CrtYB"

    def test_which_is_not_the_enzyme_the_flux_law_reads(self, spec, ceiling):
        """The whole diagnosis, computed from the spec rather than asserted in prose: the
        ceiling is set at CrtYB and the only genotype input the model accepts is CrtE."""
        assert ceiling.capacity_enzyme != spec.entry_enzyme

    def test_it_reports_the_unit_the_literature_uses(self, ceiling):
        assert ceiling.mg_per_gdcw == pytest.approx(1.2483, abs=1e-3)


class TestTheCalibrationIsMarkedLegacyOnly:
    def test_the_ceiling_carries_the_measurement_that_refutes_it(self, ceiling):
        assert ceiling.refuted_by is not None
        assert ARHAR_PMID in ceiling.refuted_by.source

    def test_and_therefore_declares_itself_legacy_only(self, ceiling):
        """`predict_product(mode='empirical')` returns a content solved against this
        capacity. A caller must not be able to hold that number without holding this."""
        assert ceiling.legacy_only is True

    def test_the_marking_travels_on_the_shipped_constant_itself(self):
        """Not beside it. A field on `NodeKinetics` would be dropped by every caller that
        reads `vmax_per_growth`, and `predict.py` refuses unrecognised fields outright."""
        capacity = BETA_CAROTENE_KINETICS["lycopene"].vmax_per_growth

        assert capacity == 2.3252e-3
        assert capacity.refuted_by is not None

    def test_the_fold_is_computed_from_the_shipped_constant(self, ceiling):
        """Not written down. A literal 63x goes stale the first time anyone refits."""
        assert ceiling.refuted_fold == pytest.approx(63.3, rel=0.01)

    def test_an_undeclared_calibration_is_not_marked(self):
        """The marking is a declaration about a specific fit, not a default. A pathway that
        has not been refuted must not inherit beta-carotene's refutation."""
        chain = PathwaySpec(
            product="prod", organism="S. cerevisiae", entry_enzyme="E0",
            precursor_metabolite="p",
            nodes=(Node("n0", RateLaw.SATURATING, Fate.DILUTED, "E0"),
                   Node("prod", RateLaw.PASSTHROUGH, Fate.DILUTED, "",
                        molar_mass_g_per_mol=MOLAR_MASS)))
        plain = content_ceiling(chain, {"n0": NodeKinetics(vmax_per_growth=1e-3, km=1e-4)})

        assert plain == 1e-3
        assert plain.refuted_by is None
        assert plain.legacy_only is False
        assert plain.refuted_fold is None


class TestTheMarkingSurvivesTheCapacityLayer:
    """`pathway/capacity.py` is the one place the ceiling changes hands.

    `predict.py`'s cassette layer passes the ceiling to `capacity_for` and reports what comes
    back. A `float()` cast there returns a bare number, so the layer that exists to say the
    capacity was not measured on this strain would have reported it stripped of the fact that
    the capacity is refuted outright.
    """

    def test_an_unstated_cassette_returns_the_ceiling_still_marked(self, ceiling):
        from ystwin.pathway.capacity import capacity_for

        answered = capacity_for("beta_carotene", ceiling, None)

        assert answered == ceiling
        assert answered.refuted_by is ceiling.refuted_by

    def test_and_so_does_the_calibration_strain_s_own_dosage(self, ceiling):
        from ystwin.pathway.capacity import capacity_for

        assert capacity_for("beta_carotene", ceiling, {"crtYB": 1.0}).legacy_only is True

    def test_a_bare_anchor_still_normalises(self):
        """The cast had a job -- an int or a numpy scalar -- and it still does it."""
        from ystwin.pathway.capacity import capacity_for

        assert capacity_for("beta_carotene", 1, None) == 1.0
        assert isinstance(capacity_for("beta_carotene", np.float64(2.0), None), float)


class TestTheRefutationIsTheVendoredMeasurementAndNotACopyOfIt:
    """Cross-checked against `data/carotenoid_batch/published_batch_titres.tsv`.

    Two copies of a number in two files is how the last set of claims went stale -- see
    `pathway/published_cassettes.py`. This makes the survey row the authority.
    """

    @pytest.fixture(scope="class")
    def survey(self):
        from ystwin.pathway.published_cassettes import published_cassettes

        return {row.pmid: row for row in published_cassettes()}

    def test_the_content_matches_the_surveyed_row(self, ceiling, survey):
        row = survey[ARHAR_PMID]

        assert ceiling.refuted_by.measured_mg_per_gdcw == row.content_mg_per_gdcw
        assert row.content_mg_per_gdcw == ARHAR_CONTENT_MG_PER_GDCW

    def test_the_refuting_row_is_a_specific_assay_for_the_same_analyte(self, ceiling, survey):
        """A ceiling on beta-carotene is refuted only by beta-carotene. The survey's
        `measurand` column is the authoritative signal for that, not the notes text."""
        assert survey[ARHAR_PMID].measurand == "beta_carotene"
        assert ceiling.refuted_by.analyte == "beta_carotene"

    def test_a_larger_absorbance_sum_would_not_have_done(self, survey):
        """Lopez 2019 reports 21 mg/gDCW, 17x the ceiling, and refutes nothing: a plate
        reader sums beta-carotene with every other carotenoid in the extract."""
        lopez = survey["31380362"]

        assert lopez.content_mg_per_gdcw > 1.25
        assert lopez.measurand == "total_carotenoid"


class TestNoRefitOnTheseSixStatesReachesTheMeasurement:
    """Why option (a) -- refit the capacity on another enzyme -- was refused.

    The capacity is fitted flat across strains that differ in expression. The only refit the
    vendored data supports is to let it scale with a measured relative mRNA, and
    `data/carotenoid/elizondo2025_relative_mrna.tsv` carries CrtE, CrtI and CrtYB for all six
    states. Every such refit is tried here. All of them fit worse than flat, and all of them
    land within 1.4x of the ceiling Arhar 2024 passes 63-fold.
    """

    @pytest.fixture(scope="class")
    def states(self):
        directory = paths.data_dir() / "carotenoid"
        steady = directory / "elizondo2025_steady_states.tsv"
        mrna = directory / "elizondo2025_relative_mrna.tsv"
        if not (steady.exists() and mrna.exists()):
            pytest.skip("Elizondo calibration data not vendored")
        expression = pd.read_csv(mrna, sep="\t").pivot_table(
            index="condition", columns="gene", values="rel_expression")
        return pd.read_csv(steady, sep="\t").set_index("condition").join(expression)

    @staticmethod
    def _fit(states, scale):
        """Capacity ``c * scale`` and one shared km, on the shipped saturating law."""
        mu = states.mu_per_h.to_numpy()
        flux = (states.q_lycopene + states.q_betacarotene).to_numpy()
        observed = states.q_betacarotene.to_numpy()

        def predicted(p):
            vmax = p[0] * scale * mu
            a, b, c = mu, mu * p[1] + vmax - flux, -flux * p[1]
            lycopene = (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
            return vmax * lycopene / (p[1] + lycopene)

        fit = least_squares(lambda p: np.log(predicted(p)) - np.log(observed),
                            x0=[2.3e-3, 6e-4], bounds=([1e-8, 1e-8], [1.0, 1.0]))
        return float(fit.x[0]), float(np.sqrt(np.mean(fit.fun ** 2)))

    def test_scaling_the_capacity_by_any_measured_expression_fits_worse(self, states):
        flat_capacity, flat_rmse = self._fit(states, np.ones(len(states)))

        for gene in ("CrtE", "CrtYB", "CrtI"):
            _, rmse = self._fit(states, states[gene].to_numpy())

            assert rmse > flat_rmse, (
                f"a capacity proportional to {gene} expression now fits the six states "
                f"better than the flat capacity ({rmse:.4f} against {flat_rmse:.4f}); the "
                f"refusal to refit rests on it fitting worse, so revisit it")
        assert flat_capacity * MOLAR_MASS == pytest.approx(1.2, abs=0.1)

    def test_and_none_of_them_comes_within_forty_fold_of_the_measurement(self, states):
        """The refit is not merely worse, it is nowhere near. The gap is structural: the
        capacity-setting dosage is not a field any of these six states varies."""
        for gene in ("CrtE", "CrtYB", "CrtI"):
            capacity, _ = self._fit(states, states[gene].to_numpy())
            reached = capacity * states[gene].max() * MOLAR_MASS

            assert ARHAR_CONTENT_MG_PER_GDCW / reached > 40.0
