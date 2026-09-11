"""Nulls for the claims this package makes, because none of them currently have one.

Three claims here are quantitative and none is falsifiable as stated.

**"A latent stress state transfers to an unseen stressor."** Reported as a held-out R2.
But a random subspace of the same dimension also predicts *something*, because the
channels are correlated and any low-rank projection captures shared variance. Without
knowing what a structureless model scores, a transfer R2 of 0.6 is a number and not
evidence.

**"These channels identify these modules."** Reported from the Fisher information of a
literature loading matrix. A random matrix of the same shape and scale has an
information matrix too, and the question is whether the real one does better.

**"Growth rate explains a median 0.82 of log(RFU/OD) variance."** This one is the most
exposed, and the reason is specific: `log(RFU/OD)` and `mu` are both smooth,
autocorrelated series *derived from the same optical density trace*. Two smooth
autocorrelated series correlate spuriously as a matter of course -- the effective number
of independent observations in a 25-point trace is far below 25 -- and here they are not
even independent by construction, since an OD error moves the numerator of one and the
derivative of the other. The right null preserves each series' own autocorrelation and
destroys only the relation between them, which is what phase randomisation does.

A null is not a formality. In the closest comparable published work a *random* smooth
auxiliary channel outperformed the real biosensor channels on one of six tasks, and that
was only visible because the random control was in the table.

The convention for reporting is the skill score: how much of the gap between a trivial
baseline and perfection a method actually closes. See :func:`skill_score`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "NullResult",
    "compare_to_null",
    "matched_marginals",
    "phase_randomised",
    "rotated_subspace",
    "shuffled_loadings",
    "skill_score",
    "skill_score_if_defined",
]

_MIN_DRAWS = 20


# ---------------------------------------------------------------------------
# the result object
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NullResult:
    """An observed statistic beside the distribution a structureless model gives it.

    Args:
        observed: The statistic on the real data.
        null: The statistic on each null draw.
        greater_is_better: Whether a larger statistic means a better model. R2 and
            power are; error and RMSE are not. The p-value's tail depends on it and
            getting it wrong silently inverts the conclusion.
        label: What was tested, for reporting.
    """

    observed: float
    null: np.ndarray = field(repr=False)
    greater_is_better: bool
    label: str = ""

    @property
    def n_draws(self) -> int:
        return int(np.sum(np.isfinite(self.null)))

    @property
    def null_median(self) -> float:
        finite = self.null[np.isfinite(self.null)]
        return float(np.median(finite)) if finite.size else float("nan")

    @property
    def p_value(self) -> float:
        """Fraction of null draws at least as extreme as the observed statistic.

        Uses the ``(hits + 1) / (draws + 1)`` form. A permutation p-value of exactly
        zero is not attainable -- the observed arrangement is itself one of the
        arrangements the null could have produced -- and reporting 0.0 from a finite
        number of draws claims a precision the procedure does not have. The floor is
        therefore ``1 / (draws + 1)``, which is the honest resolution limit.
        """
        finite = self.null[np.isfinite(self.null)]
        if finite.size == 0:
            return float("nan")
        if self.greater_is_better:
            hits = int(np.sum(finite >= self.observed))
        else:
            hits = int(np.sum(finite <= self.observed))
        return (hits + 1) / (finite.size + 1)

    @property
    def beats_null(self) -> bool:
        """Whether the observed statistic clears the null at the 5% level."""
        return bool(np.isfinite(self.p_value) and self.p_value < 0.05)

    @property
    def effect_over_null(self) -> float:
        """Observed minus the null median, signed so positive always means better."""
        gap = self.observed - self.null_median
        return float(gap if self.greater_is_better else -gap)

    def summary(self) -> str:
        verdict = "beats" if self.beats_null else "DOES NOT BEAT"
        name = f"{self.label}: " if self.label else ""
        return (
            f"{name}observed {self.observed:.4f}, null median {self.null_median:.4f} "
            f"({self.n_draws} draws) -- {verdict} the null, p = {self.p_value:.4f}, "
            f"gain {self.effect_over_null:+.4f}"
        )


def compare_to_null(
    statistic,
    data,
    surrogate,
    n_draws: int = 200,
    greater_is_better: bool = True,
    label: str = "",
    seed: int = 0,
) -> NullResult:
    """Score the real data, then score ``n_draws`` structureless surrogates of it.

    Args:
        statistic: Callable taking one dataset and returning a float. Must be
            deterministic given its input, or the null and the observation are not
            comparable -- seed anything stochastic inside it.
        data: The real dataset, in whatever form ``statistic`` and ``surrogate`` accept.
        surrogate: Callable ``(data, rng) -> data`` producing one null draw.
        n_draws: Number of null draws. The p-value cannot resolve below
            ``1 / (n_draws + 1)``, so 200 draws bottoms out at 0.005.
        greater_is_better: See :class:`NullResult`.
        label: What is being tested.
        seed: Random seed.

    Raises:
        ValueError: if ``n_draws`` is too small for the p-value to mean anything.
    """
    if n_draws < _MIN_DRAWS:
        raise ValueError(
            f"n_draws={n_draws} cannot resolve a p-value below "
            f"{1 / (n_draws + 1):.3f}; use at least {_MIN_DRAWS}"
        )
    rng = np.random.default_rng(seed)
    observed = float(statistic(data))
    draws = np.empty(n_draws, dtype=float)
    for i in range(n_draws):
        try:
            draws[i] = float(statistic(surrogate(data, rng)))
        except (ValueError, np.linalg.LinAlgError, ZeroDivisionError):
            # A surrogate can be degenerate -- a shuffle that leaves a channel
            # constant, say. Record it as missing rather than aborting the sweep or,
            # worse, substituting a value that would bias the null.
            draws[i] = np.nan
    return NullResult(observed, draws, greater_is_better, label)


# ---------------------------------------------------------------------------
# surrogates
# ---------------------------------------------------------------------------


def phase_randomised(series: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """A surrogate with the same power spectrum and destroyed phase relations.

    The right null for "does X explain Y" when both are smooth time series. Shuffling
    the samples of a smooth trace produces something visibly unlike it -- white noise
    with the same histogram -- and any statistic sensitive to smoothness then rejects
    the null for the wrong reason. Randomising the Fourier phases instead preserves the
    autocorrelation exactly, which is the property that generates spurious correlation,
    and removes only the alignment with the other series.

    This is the standard surrogate-data construction for testing against a linear
    Gaussian process with the observed spectrum (Theiler et al., 1992).

    Args:
        series: One real-valued series. Multi-column input is randomised per column
            with independent phases, which destroys cross-channel structure too.
        rng: Random generator.

    Returns:
        A surrogate of the same shape, real-valued, with the same mean and
        approximately the same autocorrelation.
    """
    x = np.asarray(series, dtype=float)
    if x.ndim == 1:
        return _phase_randomise_1d(x, rng)
    return np.column_stack([_phase_randomise_1d(x[:, j], rng) for j in range(x.shape[1])])


def _phase_randomise_1d(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n = x.size
    if n < 4:
        raise ValueError("need at least 4 points to randomise phases")
    spectrum = np.fft.rfft(x - x.mean())
    magnitude = np.abs(spectrum)
    phases = rng.uniform(0.0, 2.0 * np.pi, magnitude.size)
    # The DC term carries no phase, and for even n the Nyquist term must stay real or
    # the inverse transform is complex.
    phases[0] = 0.0
    if n % 2 == 0:
        phases[-1] = 0.0
    surrogate = np.fft.irfft(magnitude * np.exp(1j * phases), n=n)
    return surrogate + x.mean()


def shuffled_loadings(loadings: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Loadings with every entry permuted, preserving the multiset of values.

    Destroys which reporter reads which module while keeping the overall distribution
    of loading strengths, so a method that only exploits the *magnitude* of the entries
    scores the same on the null. What it cannot keep is the crosstalk *structure* --
    which is the thing the literature matrix is supposed to contribute.

    Args:
        loadings: ``(channels, modules)``.
        rng: Random generator.
    """
    matrix = np.asarray(loadings, dtype=float)
    flat = matrix.ravel().copy()
    rng.shuffle(flat)
    return flat.reshape(matrix.shape)


def rotated_subspace(readings: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Readings projected onto a random subspace of the same dimension and rescaled.

    The null for a transfer or dimension claim. A rank-*k* model of correlated channels
    captures shared variance whether or not the axes mean anything, so the honest
    question is whether the *fitted* axes beat *arbitrary* axes of the same count.
    A random orthonormal rotation of the data keeps the total variance and the number
    of directions and discards their identity.

    Args:
        readings: ``(observations, channels)``.
        rng: Random generator.
    """
    data = np.asarray(readings, dtype=float)
    if data.ndim != 2:
        raise ValueError("readings must be observations-by-channels")
    n_channels = data.shape[1]
    # A Haar-distributed orthogonal matrix from the QR of a Gaussian, with the sign
    # convention fixed so the distribution is genuinely uniform rather than biased by
    # LAPACK's choice of R's diagonal signs.
    gaussian = rng.normal(size=(n_channels, n_channels))
    q, r = np.linalg.qr(gaussian)
    q = q * np.sign(np.diag(r))
    return data @ q


def matched_marginals(readings: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Each channel independently resampled from its own observed values.

    Keeps every channel's marginal distribution exactly and removes all
    cross-channel correlation. The null for "these channels carry shared structure":
    if a method scores as well here, what it found was in the marginals.

    Unlike :func:`phase_randomised` this destroys autocorrelation as well, so it is
    the right null for observation-indexed data and the wrong one for a time course.

    Args:
        readings: ``(observations, channels)``.
        rng: Random generator.
    """
    data = np.asarray(readings, dtype=float)
    if data.ndim != 2:
        raise ValueError("readings must be observations-by-channels")
    out = np.empty_like(data)
    for j in range(data.shape[1]):
        column = data[:, j]
        finite = column[np.isfinite(column)]
        if finite.size == 0:
            raise ValueError(f"channel {j} has no finite observations to resample")
        out[:, j] = rng.choice(finite, size=column.size, replace=True)
        out[~np.isfinite(column), j] = np.nan  # keep the missingness pattern
    return out


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------


def skill_score_if_defined(observed: float, baseline: float,
                           perfect: float = 0.0) -> float | None:
    """:func:`skill_score`, or ``None`` where it is undefined.

    For a report that scores many partitions and should not lose the other rows because
    one baseline happened to be perfect. Code that must not carry on past an undefined
    comparison should call :func:`skill_score` and let it raise.
    """
    try:
        return skill_score(observed, baseline, perfect)
    except ValueError:
        return None


def skill_score(observed: float, baseline: float, perfect: float = 0.0) -> float:
    """Fraction of the achievable improvement over a baseline that was realised.

    ``(baseline - observed) / (baseline - perfect)``, the convention from forecast
    verification. One means the baseline's error was eliminated, zero means the method
    matched the baseline, and **negative means it did worse than the trivial
    prediction** -- which is the outcome a bare error figure hides and which a skill
    score makes impossible to miss.

    Defaults assume an error-like statistic where zero is perfect. For a statistic
    where larger is better, pass its perfect value and the sign works out.

    Args:
        observed: The method's score.
        baseline: The trivial predictor's score on the same data.
        perfect: The score a perfect prediction would get.

    Raises:
        ValueError: if the baseline already achieves the perfect score, which leaves
            no room to improve and no denominator.
    """
    denominator = baseline - perfect
    if abs(denominator) < 1e-15:
        raise ValueError(
            "the baseline already scores perfectly; a skill score is undefined and "
            "the comparison needs a harder baseline or a different statistic"
        )
    return float((baseline - observed) / denominator)
