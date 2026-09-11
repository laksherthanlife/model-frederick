"""Arithmetic of a fixed expression law and its exploratory fourteen-gene comparison.

``Environment.pathway_flux`` used to be a caller-supplied constant, and every entry point
obtained it as ``q_lycopene + q_betacarotene`` off the very chemostat state being predicted.
The model measured the product in order to predict the product; ``scripts/predict_product.py``
printed "this is a residual and not a prediction" in its own output, which was honest and
still left the tree with no non-circular source of a flux.

``pathway/flux.py`` supplies one from the genotype instead: relative expression of the first
committed heterologous step times a single fitted scalar. This file does two jobs.

The first is arithmetic bookkeeping -- the shipped ``BETA_CAROTENE_FLUX`` must be what the
vendored data gives, so a hand-edited constant cannot drift away from the file it claims to
come from.

The second is an exploratory comparison in :class:`TestOnlyTheHeterologousGenesCarryIt`,
whose historical name is cited by ``data/pathways/beta_carotene.toml``. ERG9 actually has
slightly positive skill against the training-only baseline. A gap between gene groups
does not prove causal cassette dosage or exclude expression confounding. CrtE was ranked
on these same three strains; the default forward validator therefore selects BOTH gene
and model inside training folds, tested separately in ``test_product_nested.py``.
"""

from __future__ import annotations

import math

import pandas as pd
import pytest

from ystwin import paths
from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX
from ystwin.pathway.flux import (
    FluxCalibration,
    fit_flux_law,
    predict_flux_from_expression,
)

SOURCE = "PMID 40891387 (doi 10.1021/acssynbio.5c00256), Elizondo 2025"

# Native to the host, so they are the control and not the candidates. Every one of these is
# transcribed from the same RNA as CrtE and normalised the same way.
NATIVE = ("ERG10", "ERG12", "ERG13", "ERG20", "ERG8", "ERG9",
          "HMG1", "HMG2", "BTS1", "IDI1", "MVD1")
HETEROLOGOUS = ("CrtE", "CrtI", "CrtYB")


@pytest.fixture(scope="module")
def measurements():
    """The six steady states with all fourteen relative-expression columns joined on.

    Flux is ``q_lycopene + q_betacarotene``: carotenoids are not secreted, so both measured
    rates are growth washing out an intracellular pool, and their sum is everything the
    desaturase made.
    """
    directory = paths.data_dir() / "carotenoid"
    states = pd.read_csv(directory / "elizondo2025_steady_states.tsv", sep="\t")
    states["flux"] = states.q_lycopene + states.q_betacarotene
    mrna = pd.read_csv(directory / "elizondo2025_relative_mrna.tsv", sep="\t")
    wide = mrna.pivot_table(index="condition", columns="gene",
                            values="rel_expression").reset_index()
    return states[["condition", "strain", "mu_per_h", "flux"]].merge(wide, on="condition")


@pytest.fixture(scope="module")
def genes(measurements):
    return [c for c in measurements.columns
            if c not in {"condition", "strain", "mu_per_h", "flux"}]


@pytest.fixture(scope="module")
def scores(measurements, genes):
    """One leave-one-strain-out fit per measured gene, keyed by gene."""
    return {
        gene: fit_flux_law(list(measurements[gene]), list(measurements.flux),
                           list(measurements.strain), gene, SOURCE)
        for gene in genes
    }


@pytest.fixture(scope="module")
def crte(scores):
    return scores["CrtE"]


class TestTheShippedConstantIsWhatTheDataGives:
    def test_the_shipped_alpha_is_a_fresh_fit_of_the_vendored_states(self, crte):
        """The constant in ``calibrations.py`` is the only copy anything imports.

        It was produced by ``scripts/fit_pathway_flux.py`` from the two TSVs in
        ``data/carotenoid/``. Nothing else re-derives it at import time, so if the file is
        refreshed, or the constant is hand-edited, only this comparison notices.
        """
        assert BETA_CAROTENE_FLUX.alpha == pytest.approx(crte.alpha, rel=1e-3)

    def test_the_shipped_leave_one_strain_out_error_is_the_measured_one(self, crte):
        assert BETA_CAROTENE_FLUX.loso_rmse_log == pytest.approx(crte.loso_rmse_log,
                                                                 abs=1e-4)

    def test_the_leave_one_strain_out_error_is_about_a_fifth_of_a_log(self, crte):
        assert crte.loso_rmse_log == pytest.approx(0.2002, abs=1e-3)

    def test_the_shipped_skill_is_the_measured_one(self, crte):
        assert BETA_CAROTENE_FLUX.loso_skill == pytest.approx(crte.loso_skill, abs=1e-3)

    def test_the_held_out_skill_beats_the_geometric_mean_of_the_training_strains(self,
                                                                                crte):
        """Fixed-model skill compares against each fold's training geometric mean.

        +0.745 is the reduction in log RMSE, not squared-error skill or an independently
        validated score for choosing CrtE from fourteen genes on these same strains.
        """
        assert crte.loso_skill == pytest.approx(0.7446797787533761, abs=1e-12)

    def test_the_calibration_is_fitted_on_the_six_published_states(self, crte):
        assert crte.n_states == BETA_CAROTENE_FLUX.n_states == 6

    def test_the_calibration_names_the_gene_its_numbers_belong_to(self):
        """A scalar carries its gene so it cannot be applied to another gene's expression.

        CrtYB fits an alpha within 3% of CrtE's, so the two are numerically
        interchangeable and biologically not -- nothing but this field would catch the
        swap.
        """
        assert BETA_CAROTENE_FLUX.entry_enzyme == "CrtE"


class TestOnlyTheHeterologousGenesCarryIt:
    """Cited by name from ``data/pathways/beta_carotene.toml``; do not rename lightly."""

    def test_all_fourteen_measured_genes_are_scored(self, scores):
        """The control is only a control if it covers everything the paper measured.

        Dropping a native gene from the sweep -- by a column rename, or by a filter that
        skips a gene with a zero somewhere -- would narrow the comparison until the result
        became a selection of the genes that agree with it.
        """
        assert len(scores) == 14

    def test_the_entry_enzyme_ranks_first_of_all_fourteen(self, scores):
        """The historical CrtE choice is data-selected, not prespecified.

        This ranking reproduces an exploratory observation. It does not identify a
        causal flux-control step, and its winner must be selected again inside training
        folds before being scored as a new-strain prediction.
        """
        best = max(scores, key=lambda gene: scores[gene].loso_skill)

        assert best == "CrtE"

    def test_native_erg9_is_slightly_positive_against_the_training_baseline(self, scores):
        """The old all-native-negative claim relied on a baseline containing test outcomes.

        Correcting that baseline leaves ERG9 at +0.003, rather than -0.46. The gene-group
        gap remains descriptive evidence, but neither a zero threshold nor this ranking
        establishes that cassette dosage is the cause.
        """
        positive = {gene: scores[gene].loso_skill
                    for gene in NATIVE if scores[gene].loso_skill >= 0}

        assert set(positive) == {"ERG9"}
        assert positive["ERG9"] == pytest.approx(0.0030432865356496697, abs=1e-12)

    def test_the_best_native_gene_is_beaten_by_the_worst_heterologous_one(self, scores):
        """A gap, not an ordering -- the three cassette genes clear the native block whole.

        CrtI is the weakest of the three at +0.359 and sits above ERG9's +0.003. The
        descriptive separation survives the training-only baseline correction. Selection
        uncertainty still requires nested validation, not a causal interpretation.
        """
        best_native = max(scores[gene].loso_skill for gene in NATIVE)
        worst_heterologous = min(scores[gene].loso_skill for gene in HETEROLOGOUS)

        assert worst_heterologous > best_native

    def test_the_native_genes_are_all_present_in_the_measured_set(self, genes):
        assert set(NATIVE) <= set(genes)

    def test_the_cassette_genes_are_all_present_in_the_measured_set(self, genes):
        assert set(HETEROLOGOUS) <= set(genes)


class TestTheFitRefusesWhatItCannotScore:
    def test_a_single_strain_cannot_be_scored_leave_one_strain_out(self, measurements):
        """With one strain there is nothing to hold out and the score is in sample.

        Leave-one-STATE-out would still run here, and would look respectable, because the
        strain's own flux level stays in the training set -- which is most of the answer to
        "what does a strain I have never measured produce".
        """
        one = measurements[measurements.strain == "b-car2"]

        with pytest.raises(ValueError, match="at least two strains"):
            fit_flux_law(list(one.CrtE), list(one.flux), list(one.strain), "CrtE", SOURCE)

    def test_fewer_than_two_states_are_refused(self):
        with pytest.raises(ValueError, match="at least two states"):
            fit_flux_law([1.0], [1e-3], ["b-car2"], "CrtE", SOURCE)

    def test_mismatched_list_lengths_are_refused(self, measurements):
        """The three lists are positional, so a short one silently re-pairs the data.

        ``zip`` would truncate and the strain labels would drift off their states, giving a
        held-out score computed on the wrong partition.
        """
        with pytest.raises(ValueError, match="same length"):
            fit_flux_law(list(measurements.CrtE), list(measurements.flux)[:-1],
                         list(measurements.strain), "CrtE", SOURCE)

    def test_a_zero_expression_state_is_refused(self, measurements):
        """The fit is a mean of log ratios, and log 0 is not a large negative number here.

        Zero expression is a strain that does not carry the cassette, which is a different
        object from a strain expressing it weakly. Refusing says so.
        """
        expression = list(measurements.CrtE)
        expression[0] = 0.0

        with pytest.raises(ValueError, match="log space"):
            fit_flux_law(expression, list(measurements.flux), list(measurements.strain),
                         "CrtE", SOURCE)

    def test_a_negative_flux_state_is_refused(self, measurements):
        flux = list(measurements.flux)
        flux[0] = -1e-4

        with pytest.raises(ValueError, match="log space"):
            fit_flux_law(list(measurements.CrtE), flux, list(measurements.strain),
                         "CrtE", SOURCE)

    def test_a_non_positive_alpha_cannot_be_declared(self):
        with pytest.raises(ValueError, match="alpha must be positive"):
            FluxCalibration(0.0, "CrtE", 0.2, 0.6, 6, SOURCE)

    def test_a_calibration_fitted_on_one_state_cannot_be_declared(self):
        with pytest.raises(ValueError, match="fewer than two states"):
            FluxCalibration(1e-3, "CrtE", 0.2, 0.6, 1, SOURCE)


class TestPredictionCarriesItsError:
    def test_the_flux_is_the_scalar_times_the_expression(self):
        predicted = predict_flux_from_expression(0.5, BETA_CAROTENE_FLUX)

        assert predicted.flux_mmol_per_gdcw_h == BETA_CAROTENE_FLUX.alpha * 0.5

    def test_zero_expression_is_refused_rather_than_returning_zero_flux(self):
        """Zero is not a small flux; it is a strain without the pathway.

        Returning ``alpha * 0`` would look like a prediction of "makes almost none", which
        is a statement about a producer, from a law fitted in log space that has nothing to
        say about a non-producer.
        """
        with pytest.raises(ValueError, match="log space"):
            predict_flux_from_expression(0.0, BETA_CAROTENE_FLUX)

    def test_a_negative_expression_is_refused(self):
        with pytest.raises(ValueError, match="must be positive"):
            predict_flux_from_expression(-0.3, BETA_CAROTENE_FLUX)

    def test_the_typical_fold_error_is_the_log_error_exponentiated(self):
        """The unit a reader thinks in. 0.2002 in log space is 1.22x, not 20%.

        Reporting the log rmse alone invites it being read as a fraction, which understates
        the spread by about a tenth here and much more for the weaker genes.
        """
        assert BETA_CAROTENE_FLUX.typical_fold_error == math.exp(
            BETA_CAROTENE_FLUX.loso_rmse_log)

    def test_the_interval_brackets_the_point_prediction(self):
        predicted = predict_flux_from_expression(1.0, BETA_CAROTENE_FLUX)
        low, high = predicted.interval()

        assert low < predicted.flux_mmol_per_gdcw_h < high

    def test_the_interval_is_the_held_out_spread_and_not_a_fit_uncertainty(self):
        """What a user of this number faces is the spread of held-out predictions.

        A confidence interval on ``alpha`` would be several times narrower with six states,
        and would describe how well the scalar is pinned rather than how wrong the next
        strain is likely to be.
        """
        predicted = predict_flux_from_expression(1.0, BETA_CAROTENE_FLUX)
        low, high = predicted.interval()

        assert high / low == pytest.approx(BETA_CAROTENE_FLUX.typical_fold_error ** 2,
                                           rel=1e-12)

    def test_an_expression_outside_the_fitted_range_is_flagged(self):
        """Extrapolating a proportionality is not refused -- proportionality is the whole
        model, so there is no boundary to enforce -- but it must arrive with a note.

        The window is the calibration's OWN range, and that is the fix this test records.
        It used to be the literals ``0.2 <= e <= 5.0``, described in the note as "roughly
        0.5-3x", while the measured CrtE range is [0.245, 1.000] -- a 4.08-fold span whose
        maximum sits 5x below the literal ceiling. A guard whose bounds are unrelated to the
        data cannot fire where it matters, and an expression of 2.0 -- twice anything ever
        measured -- passed silently."""
        predicted = predict_flux_from_expression(2.0, BETA_CAROTENE_FLUX)

        assert any("outside the range the law was fitted over" in note
                   for note in predicted.notes)

    def test_the_window_is_the_measured_range_and_not_a_literal(self):
        low, high = BETA_CAROTENE_FLUX.expression_range

        assert (low, high) == pytest.approx((0.245160473, 1.0))

    def test_below_the_range_is_flagged_too(self):
        """Both directions. The flux law is fitted in log space, so an input a factor below
        the range is as much an extrapolation as one a factor above."""
        predicted = predict_flux_from_expression(0.1, BETA_CAROTENE_FLUX)

        assert any("outside the range" in note for note in predicted.notes)

    def test_the_note_says_how_far_outside(self):
        """A flag without a magnitude cannot be triaged."""
        predicted = predict_flux_from_expression(2.0, BETA_CAROTENE_FLUX)

        assert any("2.0x beyond it" in note for note in predicted.notes)

    def test_a_calibration_with_no_declared_range_says_that_instead(self):
        """Silence would read as "inside the range", which is the stronger claim."""
        import dataclasses

        unranged = dataclasses.replace(BETA_CAROTENE_FLUX, expression_range=None)

        predicted = predict_flux_from_expression(1.0, unranged)

        assert any("declares no expression range" in note for note in predicted.notes)

    def test_an_expression_inside_the_fitted_range_carries_no_note(self):
        predicted = predict_flux_from_expression(1.0, BETA_CAROTENE_FLUX)

        assert predicted.notes == ()


class TestTheLawIsAboutStrainsAndNotDilutionRates:
    def test_the_shipped_scalar_reproduces_every_measured_state_within_the_stated_fold(
            self, measurements):
        """The in-sample residuals have to sit inside the held-out spread it advertises.

        Not evidence on their own -- the held-out score above is that -- but a fit whose
        training residuals exceeded its own LOSO error would mean the two were computed on
        different quantities.
        """
        errors = [
            abs(math.log(predict_flux_from_expression(
                row.CrtE, BETA_CAROTENE_FLUX).flux_mmol_per_gdcw_h) - math.log(row.flux))
            for _, row in measurements.iterrows()
        ]

        assert max(errors) < 2 * BETA_CAROTENE_FLUX.loso_rmse_log

    def test_the_law_has_no_growth_rate_term(self, measurements):
        """Deliberate, and the reason is on both sides of the evidence.

        Fitting the growth exponent freely gives -0.16 with a 95% profile interval of
        [-0.52, +0.20], which contains zero and is too wide to act on with three paired
        contrasts. So the same expression at two dilution rates must give the same flux
        here -- if a mu term ever appears, this is what has to be argued with.
        """
        low = measurements[measurements.condition == "2D01"].iloc[0]
        high = measurements[measurements.condition == "2D025"].iloc[0]
        predicted = [predict_flux_from_expression(row.CrtE, BETA_CAROTENE_FLUX)
                     for row in (low, high)]

        assert (predicted[1].flux_mmol_per_gdcw_h / predicted[0].flux_mmol_per_gdcw_h
                == pytest.approx(high.CrtE / low.CrtE, rel=1e-12))

    def test_the_two_dilution_rates_really_are_far_apart(self, measurements):
        """Guards the test above: with one dilution rate it would assert nothing.

        The pair spans 2.5-fold, which is the whole range over which the absent growth term
        was tested.
        """
        low = measurements[measurements.condition == "2D01"].iloc[0]
        high = measurements[measurements.condition == "2D025"].iloc[0]

        assert high.mu_per_h > 2 * low.mu_per_h


class TestTrainingOnlyBaseline:
    def test_uninformative_expression_equals_the_training_mean_baseline(self):
        fit = fit_flux_law([1.0, 1.0, 1.0], [1.0, 10.0, 100.0],
                           ["a", "b", "c"], "entry", "synthetic")
        assert fit.loso_skill == pytest.approx(0.0, abs=1e-14)

    @pytest.mark.parametrize("invalid", [float("nan"), float("inf"), -float("inf")])
    def test_nonfinite_expression_is_not_a_prediction(self, invalid):
        with pytest.raises(ValueError, match="finite"):
            predict_flux_from_expression(invalid, BETA_CAROTENE_FLUX)

    @pytest.mark.parametrize("invalid", [float("nan"), float("inf"), -float("inf")])
    def test_nonfinite_training_target_is_refused(self, invalid):
        with pytest.raises(ValueError, match="finite"):
            fit_flux_law([1.0, 2.0], [1.0, invalid], ["a", "b"], "entry", "synthetic")
