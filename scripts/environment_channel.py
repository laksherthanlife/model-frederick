"""The environment -> titre channel, measured: carbon source moving content at FIXED growth rate.

    python3 scripts/environment_channel.py

Writes ``outputs/environment_channel.csv`` (the raw contrasts) and
``outputs/environment_law_scores.csv`` (the model comparison). Reads only the vendored
``data/phb/kocharin2013_chemostat_states.tsv``; no GEM, no wet-lab data.

WHY. `scripts/product_environment_sweep.py` measured that content was a strict function of
mu: 62 environments at a held growth rate returned ONE content value, and 92 batch
environments re-predicted at their own mu disagreed by 6.4e-07 mg/gDCW. This script carries
the measurement that says yeast disagrees, and it exists so that every number in
`pathway/environment_flux.py`'s docstring has a cell to point at rather than living in prose
-- the same discipline that caught two wrong claims in `bridge/maintenance_calibration.py`.

WHAT IT ESTABLISHES, and the two halves must not be conflated:

  * THE CHANNEL, which needs no model. One genotype (SCKK006), one vessel, carbon matched at
    0.666 Cmol/L, growth rate identical by row: PHB content is 4.33 +/- 0.19 mg/gDW on
    glucose and 16.55 +/- 0.02 on ethanol at mu = 0.05 /h. That separation is tens of
    standard deviations and requires no fit.
  * THE LAW, which is fitted and unearned. Three carbon sources put the smallest attainable
    permutation p at 1/3! = 0.1667, and a relabelling test returns exactly that. The
    coefficients are reported with that floor beside them, never without.

AND THE EFFECT REVERSES WITH GROWTH RATE. Ethanol beats glucose below the crossover and
loses above it, so a single scalar coefficient is the wrong shape -- which is why the scored
model comparison below includes the interaction term and why it wins.
"""
from __future__ import annotations

import itertools
import math
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.pathway.environment_flux import carbon_descriptor

STATES = pathlib.Path(__file__).resolve().parents[1] / "data" / "phb" / \
    "kocharin2013_chemostat_states.tsv"


def _load() -> pd.DataFrame:
    frame = pd.read_csv(STATES, sep="\t")
    frame["z"] = [carbon_descriptor(g, e) for g, e
                  in zip(frame.feed_glucose_g_per_l, frame.feed_ethanol_g_per_l)]
    frame["log_content"] = np.log(frame.phb_mg_per_gdw)
    frame["log_mu"] = np.log(frame.mu_per_h)
    return frame


def _contrasts(frame: pd.DataFrame) -> pd.DataFrame:
    """Every non-glucose feed against glucose at the SAME growth rate. No model."""
    rows = []
    for mu in sorted(frame.mu_per_h.unique()):
        block = frame[frame.mu_per_h == mu]
        reference = block[block.carbon_source == "glucose"]
        if reference.empty:
            continue
        reference = reference.iloc[0]
        for _, state in block.iterrows():
            if state.carbon_source == "glucose":
                continue
            pooled = math.hypot(reference.phb_sd_mg_per_gdw, state.phb_sd_mg_per_gdw)
            rows.append({
                "growth_rate_per_h": mu,
                "carbon_source": state.carbon_source,
                "ethanol_carbon_fraction": round(state.z, 4),
                "glucose_mg_per_gdw": reference.phb_mg_per_gdw,
                "test_mg_per_gdw": state.phb_mg_per_gdw,
                "fold_vs_glucose": round(state.phb_mg_per_gdw / reference.phb_mg_per_gdw, 4),
                "separation_sd": (round((state.phb_mg_per_gdw - reference.phb_mg_per_gdw)
                                        / pooled, 2) if pooled > 0 else None),
                "ethanol_better": bool(state.phb_mg_per_gdw > reference.phb_mg_per_gdw),
            })
    return pd.DataFrame(rows)


def _held_out(frame: pd.DataFrame, columns: list[str],
              group: str | None = None) -> tuple[float, float]:
    """Median and worst held-out fold error, holding out one STATE or one CARBON SOURCE.

    One function for both, because they are the same computation over a different grouping
    and writing them twice is how the two got different code paths: the held-out-FEED numbers
    lived only inside `tests/test_environment_flux.py` and were quoted as eight literals in
    `pathway/environment_flux.py`'s docstring with no cell to point at.
    """
    design = np.column_stack([np.ones(len(frame))] + [frame[c].values for c in columns])
    target = frame.log_content.values
    groups = frame[group].to_numpy() if group else np.arange(len(frame))
    errors = []
    for value in np.unique(groups):
        mask = groups != value
        if mask.sum() < design.shape[1]:
            continue
        beta, *_ = np.linalg.lstsq(design[mask], target[mask], rcond=None)
        errors += list(abs(target[~mask] - design[~mask] @ beta))
    errors = np.asarray(errors)
    return math.exp(float(np.median(errors))), math.exp(float(errors.max()))


def _loo(frame: pd.DataFrame, columns: list[str]) -> tuple[float, float]:
    """Leave-one-STATE-out. Kept as a name because the permutation test reads it."""
    return _held_out(frame, columns, group=None)


def _permutation_p(frame: pd.DataFrame, columns: list[str], sources, zed) -> float:
    """Relabelling test for ONE law, against the mu-only baseline.

    PER LAW as of 2026-09-04. This used to be computed once, for the z-only law, and then
    BROADCAST onto all five rows of `outputs/environment_law_scores.csv` -- so the rows
    "content constant" and "mu only", which contain no z term at all, asserted a p-value from
    a different law's relabelling. Nothing quoted the wrong cells, because every marker
    happened to target row="law=z only"; but markers address this file by (column, row), so a
    marker for row="law=mu only" would have pinned prose to a number that law never produced.

    For a law with no z term the honest answer is 1.0 and it is returned without solving
    anything: permuting the z labels cannot change a prediction that does not read z, so every
    relabelling ties with the observed one.
    """
    if "z" not in columns and "log_mu_z" not in columns:
        return 1.0
    baseline = _loo(frame, ["log_mu"])[0]
    observed = _loo(frame, columns)[0]
    permutations = list(itertools.permutations(sources))
    wins = 0
    for permutation in permutations:
        shuffled = frame.copy()
        shuffled["z"] = frame.carbon_source.map(
            {a: zed[b] for a, b in zip(sources, permutation)})
        shuffled["log_mu_z"] = shuffled.log_mu * shuffled.z
        if (baseline - _loo(shuffled, columns)[0]) >= (baseline - observed):
            wins += 1
    return wins / len(permutations)


def main() -> int:
    frame = _load()
    frame["log_mu_z"] = frame.log_mu * frame.z
    out_dir = paths.outputs_dir()

    contrasts = _contrasts(frame)
    contrasts.to_csv(out_dir / "environment_channel.csv", index=False)
    print(contrasts.to_string(index=False))

    laws = [("content constant", []),
            ("mu only", ["log_mu"]),
            ("z only", ["z"]),
            ("mu + z", ["log_mu", "z"]),
            ("mu + z + mu*z", ["log_mu", "z", "log_mu_z"])]
    # The exchangeable unit is the CARBON SOURCE, so every relabelling test runs over the 3!
    # assignments of z to carbon source -- not over states, which would overstate the
    # evidence by treating two dilution rates of one feed as independent.
    sources = frame.carbon_source.unique()
    zed = {s: frame[frame.carbon_source == s].z.iloc[0] for s in sources}

    # EVERY COLUMN BELOW IS COMPUTED PER LAW. As of 2026-09-04: b_z, b_mu_z,
    # sign_flip_growth_rate_per_h and permutation_p were computed once and assigned to the
    # whole frame, so five rows carried one law's numbers. A law that has no z term now writes
    # NA rather than borrowing another law's coefficient, and a quantity a law cannot define
    # -- the sign flip needs an interaction -- is NA too.
    rows = []
    for name, columns in laws:
        median, worst = _held_out(frame, columns)
        feed_median, feed_worst = _held_out(frame, columns, group="carbon_source")
        design = np.column_stack(
            [np.ones(len(frame))] + [frame[c].values for c in columns])
        beta, *_ = np.linalg.lstsq(design, frame.log_content.values, rcond=None)
        coefficient = dict(zip(columns, beta[1:]))
        b_z = coefficient.get("z")
        b_mu_z = coefficient.get("log_mu_z")
        rows.append({
            "law": name,
            "free_parameters": len(columns) + 1,
            "loo_median_fold_error": round(median, 4),
            "loo_worst_fold_error": round(worst, 4),
            # Held out one whole CARBON SOURCE: the bar a new medium faces. These four pairs
            # were quoted as eight literals in `pathway/environment_flux.py` and computed
            # nowhere any of them could point at.
            "loco_median_fold_error": round(feed_median, 4),
            "loco_worst_fold_error": round(feed_worst, 4),
            "b_z": None if b_z is None else round(float(b_z), 4),
            "b_mu_z": None if b_mu_z is None else round(float(b_mu_z), 4),
            # Where the carbon source stops helping and starts hurting: b_z + b_muz*log(mu)
            # = 0. Only the interaction law can define it.
            "sign_flip_growth_rate_per_h": (
                None if b_z is None or b_mu_z is None
                else round(math.exp(-float(b_z) / float(b_mu_z)), 4)),
            "n_states": len(frame),
            "n_carbon_sources": len(sources),
            "permutation_p": round(_permutation_p(frame, columns, sources, zed), 4),
            "permutation_floor": round(1.0 / math.factorial(len(sources)), 4),
        })

    scores = pd.DataFrame(rows)
    scores.to_csv(out_dir / "environment_law_scores.csv", index=False)
    interaction = scores[scores.law == "mu + z + mu*z"].iloc[0]
    beta = [None, None, interaction.b_z, interaction.b_mu_z]
    flip = interaction.sign_flip_growth_rate_per_h
    wins = interaction.permutation_p * math.factorial(len(sources))
    permutations = list(itertools.permutations(sources))

    print(f"\n{scores[['law','free_parameters','loo_median_fold_error','loo_worst_fold_error']].to_string(index=False)}")
    print(f"\nfull fit: b_z = {beta[2]:+.4f}, b_mu_z = {beta[3]:+.4f}, "
          f"sign flip at mu = {flip:.4f} /h")
    print(f"carbon-source relabelling test: p = {wins / len(permutations):.4f} "
          f"(floor 1/{len(permutations)} = {1 / len(permutations):.4f}) -- the coefficients "
          f"cannot be established from {len(sources)} carbon sources")
    print(f"\nwrote {out_dir / 'environment_channel.csv'} and "
          f"{out_dir / 'environment_law_scores.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
