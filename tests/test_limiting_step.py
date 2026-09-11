"""Which enzyme sets the flux, ranked from data rather than assumed.

The claim this module exists to make executable, on the shipped carotenoid panel: the entry
enzyme and the downstream cyclase are indistinguishable, and every native gene loses to the
flux's own mean. That says the fitted scalar reads CASSETTE EXPRESSION rather than the entry
enzyme's activity, which is why a capacity computed at the entry step would sit on the wrong
enzyme.

These tests pin the arithmetic, the refusals, and that the ranking can actually separate a
gene that predicts from one that does not -- a ranker that calls everything a winner would
report the same answer for a new product regardless of its data.
"""

from __future__ import annotations

import csv
import math

import pytest

from ystwin import paths
from ystwin.pathway.limiting_step import (
    MINIMUM_STATES,
    NotEnoughStates,
    rank_candidate_genes,
)

_MRNA = paths.data_dir() / "carotenoid" / "elizondo2025_relative_mrna.tsv"
_STATES = paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv"


def _panel():
    """The six shipped steady states, expression per gene against pathway flux."""
    with _MRNA.open() as handle:
        mrna = list(csv.DictReader(handle, delimiter="\t"))
    with _STATES.open() as handle:
        states = list(csv.DictReader(handle, delimiter="\t"))
    expr: dict[str, dict[str, float]] = {}
    for row in mrna:
        expr.setdefault(row["gene"], {})[row["condition"]] = float(row["rel_expression"])
    conditions = [s["condition"] for s in states]
    flux = [float(s["q_lycopene"]) + float(s["q_betacarotene"]) for s in states]
    return {g: [expr[g][c] for c in conditions] for g in expr}, flux


class TestTheArithmetic:
    def test_a_perfect_gene_fits_with_no_error(self):
        report = rank_candidate_genes({"perfect": [1.0, 2.0, 4.0]}, [2.0, 4.0, 8.0], "perfect")

        assert report.best.alpha == pytest.approx(2.0)
        assert report.best.rmse_log10 == pytest.approx(0.0, abs=1e-12)

    def test_the_fit_passes_through_the_origin(self):
        """No intercept, because the law's own shape is that a silent cassette makes no
        product. An intercept would let a gene score well while asserting otherwise."""
        report = rank_candidate_genes({"g": [1.0, 2.0, 3.0]}, [3.0, 6.0, 9.0], "g")

        assert report.best.alpha == pytest.approx(3.0)

    def test_genes_come_back_best_first(self):
        report = rank_candidate_genes(
            {"good": [1.0, 2.0, 4.0], "bad": [1.0, 1.0, 1.0]}, [2.0, 4.0, 8.0], "good")

        assert [f.gene for f in report.fits] == ["good", "bad"]
        assert report.best.gene == "good"


class TestItRefusesWhereItCannotAnswer:
    def test_two_states_are_refused_by_name(self):
        """With two states a one-parameter law has no residual and the ranking is noise."""
        with pytest.raises(NotEnoughStates, match=str(MINIMUM_STATES)):
            rank_candidate_genes({"g": [1.0, 2.0]}, [1.0, 2.0], "g")

    def test_the_refusal_names_the_measurement(self):
        with pytest.raises(NotEnoughStates, match="SAME steady states"):
            rank_candidate_genes({"g": [1.0, 2.0]}, [1.0, 2.0], "g")

    def test_a_length_mismatch_is_refused(self):
        with pytest.raises(ValueError, match="expression values against"):
            rank_candidate_genes({"g": [1.0, 2.0]}, [1.0, 2.0, 3.0], "g")

    @pytest.mark.parametrize("bad", [0.0, -1.0, math.inf, math.nan])
    def test_a_non_positive_or_non_finite_flux_is_refused(self, bad):
        with pytest.raises(ValueError, match="log space"):
            rank_candidate_genes({"g": [1.0, 2.0, 3.0]}, [1.0, 2.0, bad], "g")

    def test_no_candidates_is_refused(self):
        with pytest.raises(ValueError, match="no candidate genes"):
            rank_candidate_genes({}, [1.0, 2.0, 3.0], "g")


class TestTheShippedPanel:
    """The finding, as an executable claim about data in the repository."""

    def test_only_the_heterologous_genes_beat_the_flux_mean(self):
        expression, flux = _panel()

        report = rank_candidate_genes(expression, flux, entry_enzyme="CrtE")
        winners = {f.gene for f in report.fits if f.beats_the_mean}

        assert winners == {"CrtE", "CrtYB", "CrtI"}

    def test_the_entry_enzyme_is_indistinguishable_from_the_cyclase(self):
        """The whole point. If these two separated, the entry enzyme would be a mechanism;
        because they do not, it is a readout of cassette expression."""
        expression, flux = _panel()

        report = rank_candidate_genes(expression, flux, entry_enzyme="CrtE")
        by_gene = {f.gene: f for f in report.fits}

        assert report.entry_is_indistinguishable_from_best
        assert abs(by_gene["CrtE"].rmse_log10 - by_gene["CrtYB"].rmse_log10) < 0.01

    def test_every_native_gene_loses_to_the_mean(self):
        expression, flux = _panel()

        report = rank_candidate_genes(expression, flux, entry_enzyme="CrtE")
        native = [f for f in report.fits
                  if f.gene not in {"CrtE", "CrtYB", "CrtI"}]

        assert native, "the panel should carry native genes"
        assert all(f.r_squared_log < 0.0 for f in native)

    def test_the_fitted_alpha_is_recovered_for_the_shipped_entry_gene(self):
        """Cross-check against `calibrations.BETA_CAROTENE_FLUX`, fitted by a different
        script. Agreement to a few percent says both are fitting the same thing."""
        from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX
        expression, flux = _panel()

        report = rank_candidate_genes(expression, flux, entry_enzyme="CrtE")
        entry = report.entry

        assert entry.alpha == pytest.approx(BETA_CAROTENE_FLUX.alpha, rel=0.05)


class TestTheCheckCanActuallyFail:
    """A ranker that calls every gene a winner would say the same thing about any product."""

    def test_a_gene_uncorrelated_with_flux_is_reported_as_losing(self):
        report = rank_candidate_genes(
            {"flat": [1.0, 1.0, 1.0, 1.0]}, [1.0, 4.0, 2.0, 8.0], "flat")

        assert not report.best.beats_the_mean
        assert not report.any_gene_predicts
        assert "does not predict flux" in report.summary()

    def test_an_entry_enzyme_that_is_clearly_worse_is_said_to_be_worse(self):
        report = rank_candidate_genes(
            {"good": [1.0, 2.0, 4.0, 8.0], "entry": [1.0, 1.02, 0.98, 1.01]},
            [1.0, 2.0, 4.0, 8.0], entry_enzyme="entry")

        assert not report.entry_is_indistinguishable_from_best
        assert "wrong enzyme" in report.summary()

    def test_an_entry_enzyme_absent_from_the_panel_is_said_to_be_absent(self):
        report = rank_candidate_genes({"g": [1.0, 2.0, 3.0]}, [1.0, 2.0, 3.0], "missing")

        assert report.entry is None
        assert "not among the candidates" in report.summary()
