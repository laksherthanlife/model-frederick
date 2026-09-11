"""Which enzyme in a chain actually sets the flux, decided by measurement rather than assumed.

The flux law fits one scalar against ONE gene's expression, and which gene that is has to be
chosen before anything can be fitted. For beta-carotene the choice was `CrtE`, the entry
enzyme, and it works. This module exists because it works for a reason that is not the one
the name suggests, and the difference decides how a NEW product should be modelled.

WHAT THE BETA-CAROTENE DATA SAYS, and it is the reason this module was written. Regressing
the measured pathway flux on each gene's expression in turn, the entry enzyme `CrtE` and the
downstream cyclase `CrtYB` are indistinguishable -- their errors differ by less than the
spread of six states -- while every NATIVE gene in the same panel does worse than predicting
the mean. Two heterologous genes on one cassette both predict; the host's own genes do not.

So the fitted scalar is not reading the entry enzyme. It is reading CASSETTE EXPRESSION, for
which the entry enzyme happens to be one adequate readout among several. That matters
because the entry step is NOT the step that binds: its enzyme could carry orders of
magnitude more flux than is observed, so a capacity computed there is a ceiling with slack
under it, not a rate. The step that binds is downstream.

WHY THIS IS THE PIECE THAT TRANSFERS. For a product nobody has modelled, which enzyme binds
is exactly what cannot be guessed and exactly what decides the model's shape: a capacity
belongs on the limiting step and nowhere else. This module answers it from data, refuses
when the data cannot answer it, and says which of the two it is doing.

WHAT IT IS NOT. It is a diagnostic, not a fitter. It reports how well each candidate gene
would do and leaves the choosing to a caller who can see the numbers, because picking the
best of N genes by score and then reporting that score as an honest error is a selection
effect -- one this repository has already made once, and `scripts/fit_pathway_flux.py`
carries the note about it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "GeneFit",
    "LimitingStepReport",
    "NotEnoughStates",
    "rank_candidate_genes",
]

#: Below this many paired (expression, flux) states, ranking genes is not a measurement --
#: with two states any monotone gene fits perfectly and the ranking is noise. Three is the
#: minimum at which a residual exists at all for a one-parameter law.
MINIMUM_STATES = 3


class NotEnoughStates(ValueError):
    """Fewer paired states than a ranking can be computed from."""


@dataclass(frozen=True)
class GeneFit:
    """One candidate gene's one-parameter fit, ``flux = alpha * expression``."""

    gene: str
    alpha: float
    rmse_log10: float
    r_squared_log: float
    n_states: int

    @property
    def beats_the_mean(self) -> bool:
        """A negative ``R^2`` means the gene predicts worse than the flux's own average, and
        that is the honest bar -- a law that loses to a constant is not a law."""
        return self.r_squared_log > 0.0


@dataclass(frozen=True)
class LimitingStepReport:
    """Every candidate ranked, with the entry enzyme's standing called out."""

    fits: tuple[GeneFit, ...]
    entry_enzyme: str

    @property
    def best(self) -> GeneFit:
        return self.fits[0]

    @property
    def entry(self) -> GeneFit | None:
        for fit in self.fits:
            if fit.gene == self.entry_enzyme:
                return fit
        return None

    @property
    def entry_is_indistinguishable_from_best(self) -> bool:
        """True when the entry enzyme's error is within the fit's own scatter of the best.

        The comparison is against the BEST GENE'S OWN rmse rather than a fixed tolerance,
        because the question is whether the difference between two genes is larger than the
        error either of them carries. Where it is not, choosing between them on score is
        choosing on noise.
        """
        entry = self.entry
        if entry is None:
            return False
        return abs(entry.rmse_log10 - self.best.rmse_log10) < self.best.rmse_log10

    @property
    def any_gene_predicts(self) -> bool:
        return any(fit.beats_the_mean for fit in self.fits)

    def summary(self) -> str:
        if not self.any_gene_predicts:
            return (f"no candidate gene beats the flux's own mean; expression does not "
                    f"predict flux in these {self.best.n_states} states")
        best = self.best
        head = (f"best {best.gene} (rmse_log10 {best.rmse_log10:.4f}, "
                f"R2 {best.r_squared_log:.3f}) over {len(self.fits)} candidates on "
                f"{best.n_states} states")
        entry = self.entry
        if entry is None:
            return head + f"; entry enzyme {self.entry_enzyme!r} was not among the candidates"
        if best.gene == self.entry_enzyme:
            return head + "; the entry enzyme is the best candidate"
        if self.entry_is_indistinguishable_from_best:
            return (head + f"; the entry enzyme {entry.gene} is INDISTINGUISHABLE from it "
                    f"(rmse {entry.rmse_log10:.4f}), so expression is being read rather than "
                    f"that specific enzyme's activity")
        return (head + f"; the entry enzyme {entry.gene} is WORSE (rmse {entry.rmse_log10:.4f}), "
                f"so a capacity placed on the entry step would be on the wrong enzyme")


def rank_candidate_genes(expression: dict[str, list[float]],
                         flux: list[float],
                         entry_enzyme: str) -> LimitingStepReport:
    """Rank each gene by how well ``flux = alpha * expression`` fits, best first.

    Args:
        expression: Gene name to its relative expression, one value per state, in the same
            order as ``flux``.
        flux: Measured pathway flux per state, same order. Strictly positive: the fit and
            its error are computed in log space, because a flux spanning a decade is not
            well described by an additive residual.
        entry_enzyme: The gene the flux law currently uses, so its standing can be reported
            whether or not it wins.

    Raises:
        NotEnoughStates: with fewer than :data:`MINIMUM_STATES` paired states.
        ValueError: on a length mismatch, a non-positive flux, or no candidates.
    """
    if not expression:
        raise ValueError("no candidate genes supplied")
    n = len(flux)
    if n < MINIMUM_STATES:
        raise NotEnoughStates(
            f"ranking genes needs at least {MINIMUM_STATES} paired states and got {n}. "
            f"Below that a one-parameter law has no residual to be judged on and the "
            f"ranking reports noise. The measurement is expression and pathway flux from "
            f"the SAME steady states, which is the pairing a chemostat series gives")
    for gene, values in expression.items():
        if len(values) != n:
            raise ValueError(
                f"gene {gene!r} has {len(values)} expression values against {n} fluxes")
        if any(not math.isfinite(v) or v <= 0.0 for v in values):
            raise ValueError(f"gene {gene!r}: every expression must be finite and positive")
    if any(not math.isfinite(f) or f <= 0.0 for f in flux):
        raise ValueError("every flux must be finite and positive; the fit is in log space")

    log_flux = [math.log10(f) for f in flux]
    mean_log = sum(log_flux) / n
    total = sum((lf - mean_log) ** 2 for lf in log_flux)

    fits: list[GeneFit] = []
    for gene, values in expression.items():
        # Least squares through the origin, which is the law's own shape: no expression, no
        # flux. Fitting an intercept would let a gene score well while asserting product is
        # made when the cassette is silent.
        alpha = (sum(x * f for x, f in zip(values, flux))
                 / sum(x * x for x in values))
        residuals = [lf - math.log10(alpha * x) for lf, x in zip(log_flux, values)]
        rss = sum(r * r for r in residuals)
        fits.append(GeneFit(
            gene=gene, alpha=alpha,
            rmse_log10=math.sqrt(rss / n),
            r_squared_log=(1.0 - rss / total) if total > 0 else float("nan"),
            n_states=n))

    fits.sort(key=lambda f: f.rmse_log10)
    return LimitingStepReport(tuple(fits), entry_enzyme)
