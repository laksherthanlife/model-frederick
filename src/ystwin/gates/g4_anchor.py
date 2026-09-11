"""G4: does a biosensor track an independent anchor, or does it track growth rate?

PASS, REFUTED, or INCONCLUSIVE -- and keeping the third apart from the first is the point. The
usual failure is a flat anchor beside a moving reporter written up as agreement, so the anchor's
effect is compared against its own between-replicate spread before any correlation is computed.
A failing anchor QC forces INCONCLUSIVE: a reporter cannot be refuted by a broken measurement.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

__all__ = ["AnchorResult", "anchor_agreement"]

_MIN_REPLICATES = 2
_MIN_DOSES = 2


@dataclass(frozen=True)
class AnchorResult:
    """Outcome of the G4 comparison, with the numbers behind it."""

    verdict: str
    reason: str
    anchor_effect: float
    anchor_spread: float
    reporter_anchor_r: float
    growth_r2: float
    n_replicates: int
    n_doses: int
    replicates_needed: int

    def summary(self) -> str:
        return (
            f"{self.verdict}: {self.reason} "
            f"[anchor effect {self.anchor_effect:.2f} vs spread {self.anchor_spread:.2f}, "
            f"r(reporter, anchor) {self.reporter_anchor_r:+.2f}, "
            f"growth R^2 {self.growth_r2:.2f}, "
            f"{self.n_replicates} replicates x {self.n_doses} doses]"
        )


def _correlation(a, b, default: float = float("nan")) -> float:
    """Pearson r, returning ``default`` when either series has no variance."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.size < 2 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return default
    return float(np.corrcoef(a, b)[0, 1])


def _inconclusive(reason, **kw):
    base = dict(
        anchor_effect=float("nan"), anchor_spread=float("nan"),
        reporter_anchor_r=float("nan"), growth_r2=float("nan"),
        n_replicates=0, n_doses=0, replicates_needed=0,
    )
    base.update(kw)
    return AnchorResult(verdict="INCONCLUSIVE", reason=reason, **base)


def anchor_agreement(
    frame: pd.DataFrame,
    min_anchor_snr: float = 1.5,
    min_correlation: float = 0.7,
    growth_confound_r2: float = 0.9,
    target_snr: float = 2.0,
    anchor_qc_pass_rate: float | None = None,
    min_anchor_qc_pass_rate: float = 0.5,
) -> AnchorResult:
    """Compare a reporter's dose response with an independent anchor's.

    Args:
        frame: Columns ``replicate``, ``dose_mM``, ``anchor_fold_change``,
            ``reporter_response`` and ``growth_rate``.
        min_anchor_snr: Smallest ratio of anchor effect to between-replicate spread
            that counts as a detectable anchor response.
        min_correlation: Smallest reporter-anchor correlation counted as tracking.
        growth_confound_r2: Above this, the reporter's response is treated as
            explained by growth rate rather than by the anchor.
        target_snr: Signal-to-noise the replicate-count estimate aims for.
        anchor_qc_pass_rate: Fraction of the anchor's readings that cleared their
            no-reverse-transcriptase control. Below ``min_anchor_qc_pass_rate`` the
            verdict is forced to INCONCLUSIVE: a reporter cannot be refuted by an
            anchor that is measuring genomic DNA.
        min_anchor_qc_pass_rate: Threshold for the above.
    """
    required = {"replicate", "dose_mM", "anchor_fold_change", "reporter_response", "growth_rate"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"frame is missing columns {sorted(missing)}")

    n_reps = int(frame.replicate.nunique())
    n_doses = int(frame.dose_mM.nunique())
    if anchor_qc_pass_rate is not None and anchor_qc_pass_rate < min_anchor_qc_pass_rate:
        return _inconclusive(
            f"only {anchor_qc_pass_rate:.0%} of the anchor's readings cleared their "
            "no-reverse-transcriptase control, so the anchor is measuring largely "
            "genomic DNA; nothing about the reporter follows from it either way",
            n_replicates=n_reps, n_doses=n_doses,
        )
    if n_reps < _MIN_REPLICATES:
        return _inconclusive(
            f"only {n_reps} biological replicate; the anchor's reliability cannot "
            "be established without at least two",
            n_replicates=n_reps, n_doses=n_doses,
        )
    if n_doses < _MIN_DOSES:
        return _inconclusive(
            f"only {n_doses} dose level; a dose response needs at least two",
            n_replicates=n_reps, n_doses=n_doses,
        )

    per_dose = frame.groupby("dose_mM").agg(
        anchor_mean=("anchor_fold_change", "mean"),
        anchor_sd=("anchor_fold_change", "std"),
        reporter=("reporter_response", "mean"),
        growth=("growth_rate", "mean"),
    ).reset_index()

    # Signed per-replicate effect: an unsigned range scores reversals as signal.
    control_dose = float(per_dose.dose_mM.min())
    top_dose = float(per_dose.dose_mM.max())
    per_replicate = frame.pivot_table(
        index="replicate", columns="dose_mM", values="anchor_fold_change", aggfunc="mean"
    )
    effects = (per_replicate[top_dose] - per_replicate[control_dose]).to_numpy(dtype=float)
    anchor_effect = float(np.mean(effects))
    anchor_spread = float(np.std(effects, ddof=1)) if effects.size > 1 else 0.0
    snr = abs(anchor_effect) / anchor_spread if anchor_spread > 0 else np.inf

    reps_needed = n_reps
    if np.isfinite(snr) and snr > 0:
        # Spread of a mean falls as 1/sqrt(n), so the replicates needed to reach the
        # target SNR scale with the square of the shortfall.
        reps_needed = int(np.ceil(n_reps * (target_snr / snr) ** 2))

    if snr < min_anchor_snr:
        return _inconclusive(
            f"anchor moved {anchor_effect:+.2f} fold from {control_dose} to {top_dose} mM "
            f"but varies {anchor_spread:.2f} between replicates (SNR {snr:.2f} < "
            f"{min_anchor_snr}); about {reps_needed} replicates would be needed "
            "to settle it",
            anchor_effect=anchor_effect, anchor_spread=anchor_spread,
            reporter_anchor_r=float("nan"), growth_r2=float("nan"),
            n_replicates=n_reps, n_doses=n_doses, replicates_needed=reps_needed,
        )

    r = _correlation(per_dose.reporter, per_dose.anchor_mean)
    # A growth rate with no variance was not measured across doses. Reporting 0.00
    # would read as a cleared confound, so it stays NaN and the note says why.
    growth_r = _correlation(per_dose.reporter, per_dose.growth)
    growth_r2 = growth_r**2 if np.isfinite(growth_r) else float("nan")
    growth_note = (
        " (growth rate was not measured across doses, so the confound is untested)"
        if not np.isfinite(growth_r2) else ""
    )

    common = dict(
        anchor_effect=anchor_effect, anchor_spread=anchor_spread,
        reporter_anchor_r=r, growth_r2=growth_r2,
        n_replicates=n_reps, n_doses=n_doses, replicates_needed=reps_needed,
    )
    if np.isfinite(growth_r2) and growth_r2 >= growth_confound_r2:
        return AnchorResult(
            verdict="REFUTED",
            reason=(
                f"the reporter's dose response is explained by growth rate "
                f"(R^2 {growth_r2:.2f}); tracking the anchor as well does not "
                "distinguish stress from slower growth"
            ),
            **common,
        )
    if r < min_correlation:
        return AnchorResult(
            verdict="REFUTED",
            reason=(
                f"the anchor responded ({anchor_effect:+.2f} fold) but the reporter "
                f"did not track it (r {r:+.2f} < {min_correlation})" + growth_note
            ),
            **common,
        )
    return AnchorResult(
        verdict="PASS",
        reason=(
            f"the anchor responded {anchor_effect:+.2f} fold at SNR {snr:.1f} and the "
            f"reporter tracked it (r {r:+.2f})"
            + (growth_note if growth_note else f", not explained by growth (R^2 {growth_r2:.2f})")
        ),
        **common,
    )
