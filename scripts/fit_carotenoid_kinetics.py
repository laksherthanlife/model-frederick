"""Fit conditional branch diagnostics and separately validate nested forward model selection.

The forward entry score is nested by STRAIN: gene, entry law and branch family are
selected inside each outer training fold, with both beta-carotene and lycopene held out.
There are three independent strains and six condition predictions. Use --entry-only
for this validation without the older conditional kinetic diagnostics below; those
state-LOO tables condition on measured lycopene and retain sibling strain outcomes,
so their globally ranked winner is not an independent new-strain forward score.

This is where `kinetic/carotenoid.py` stops being parked. The question it answers is not
"what are the enzyme parameters" -- the dataset cannot answer that -- but the two questions
that come before it: which steps does this dataset identify at all, and does the surviving
form beat the trivial answer on a state it did not see.

Three things it is built to make hard to fudge.

**Counting.** Six steady states, and the free parameters are printed beside them for every
candidate. A four-parameter form on six points is reported as such, and it does not win.

**Every row varies ONE thing, and its name says which.** That rule was added on 2026-09-04,
because breaking it had already cost four documents a wrong sentence. The row now called
`vmax per strain, growth-independent` was called `vmax per strain` and varied two things:
per-strain capacities AND the loss of the growth term, which is the axis the sign test above
refutes before any fitting happens. Its LOO score of 0.5503 was then quoted in four places as
evidence that "a per-strain capacity scores -0.58 out of sample" -- attributing to gene
dosage a loss that growth-rate-independence was causing. The missing row, `capacity per
strain x mu`, costs the same four parameters on the same six points and scores 0.2054, ABOVE
the baseline. The pooled two-parameter fit still wins and no shipped constant moves; what
moves is the REASON, which is now a nested F test (F(2,2) = 0.5666, p = 0.6383) and a measured
capacity spread of 1.1494x across strains whose crtYB copy number is 1:2:3, rather than a
score that was measuring something else. Both are written to
``outputs/carotenoid_parsimony.csv``.

**A held-out score against a baseline that is allowed to win.** Leave-one-out over all six
states, against predicting the training mean and against carrying the mean *content*
forward. An RMSE with no baseline says nothing, and both baselines here beat several
candidates that look reasonable on the full-data fit.

**Refutation reported as loudly as fit.** The candidate this module was originally written
around -- a Michaelis-Menten cyclase whose vmax does not depend on growth rate -- is refuted
by the data, and the script prints the sign test that refutes it before it prints any fitted
number. Deleting that section would make the surviving fit look like a confirmation of the
model rather than a replacement for it.

The prediction target is deliberately not the one that would flatter the fit. Predicting
beta-carotene from the *total* pathway flux is partly circular, because that total contains
beta-carotene; here the predictor sees only the lycopene rate and the growth rate, which are
a different HPLC peak and a pump setting.

Usage: python scripts/fit_carotenoid_kinetics.py
"""
from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import f as f_dist
from scipy.stats import t as student_t

from ystwin import paths
from ystwin.kinetic.carotenoid import calibration_states
from ystwin.pathway.flux import (
    carotenoid_measurements,
    score_product_validation,
    summarize_product_validation,
)
from ystwin.pathway.spec import load_pathway

STRAINS = ("b-car2", "b-car3", "b-car4")

REFERENCE_GROWTH_RATE = 0.100987987
"""The lower of Elizondo's two dilution rates, 1/h, as the constant the growth term divides by.

Fixed here rather than taken from whatever frame is in hand, and that is not a style
preference. A leave-one-out *test* frame is a single row, so ``frame.mu_per_h.min()``
evaluates to that row's own growth rate, the growth term collapses to exactly 1, and every
held-out state at the high dilution rate is under-predicted by the full 2.5-fold ratio. The
first draft did that and scored the growth-power candidate at RMSE(log) 0.66 -- worse than
the training mean -- when its real out-of-sample error is 0.21. ``main`` asserts this
constant against the data rather than trusting the comment.
"""


# --------------------------------------------------------------------------- #
# candidates
# --------------------------------------------------------------------------- #
#
# Every candidate predicts q_betacarotene from (lycopene content, growth rate) as
#
#     q_bcar = vmax_lcy * L / (km_lcy + L)
#
# and differs only in what vmax_lcy is allowed to depend on. Parameters are carried as
# logs so the optimiser works on a scale where 1e-4 and 1e-3 are one step apart.

def _shared(t, d):
    return np.full(len(d), np.exp(t[0])), np.exp(t[1])


def _by_mrna(t, d):
    return np.exp(t[0]) * d.crtyb_mrna.to_numpy(), np.exp(t[1])


def _by_strain(t, d):
    """One capacity per strain, and NO growth term. Both, which is the point of the name.

    Named `vmax per strain, growth-independent` in `CANDIDATES` rather than `vmax per
    strain`, because every other candidate in that table varies exactly one thing against
    `vmax shared` and this one varies two. Reading it as a test of gene dosage attributes to
    per-strain capacity a loss that growth-rate-independence is causing -- the axis the
    script's own sign test refutes before any fitting happens. See `_by_strain_times_growth`,
    which is the same dosage hypothesis with the growth term the winner carries.
    """
    return np.exp(np.array([t[STRAINS.index(s)] for s in d.strain])), np.exp(t[3])


def _by_strain_times_growth(t, d):
    """One capacity per strain, times growth rate. The dosage hypothesis, stated properly.

    ADDED 2026-09-04. Its absence is why four documents recorded that "a per-strain capacity
    scores -0.58 out of sample": the only per-strain row in the table also dropped the growth
    term, so its score measured the missing mu and not the per-strain capacities. This row
    costs the same four parameters and the same six points, and it lands above the baseline.
    The pooled fit still wins -- so no shipped constant moves -- but it wins on a nested F
    test rather than on a number that was measuring something else.
    """
    capacity = np.exp(np.array([t[STRAINS.index(s)] for s in d.strain]))
    return capacity * d.mu_per_h.to_numpy(), np.exp(t[3])


def _growth_proportional(t, d):
    return np.exp(t[0]) * d.mu_per_h.to_numpy(), np.exp(t[1])


def _growth_power(t, d):
    scaled = (d.mu_per_h.to_numpy() / REFERENCE_GROWTH_RATE) ** t[2]
    return np.exp(t[0]) * scaled, np.exp(t[1])


def _mrna_times_growth_power(t, d):
    scaled = (d.mu_per_h.to_numpy() / REFERENCE_GROWTH_RATE) ** t[2]
    return np.exp(t[0]) * d.crtyb_mrna.to_numpy() * scaled, np.exp(t[1])


def _growth_proportional_unsaturated(t, d):
    """First-order coefficient in 1/h; None distinguishes absence of saturation."""
    return np.exp(t[0]) * d.mu_per_h.to_numpy(), None


CANDIDATES = {
    "vmax shared":                (_shared, [-8.2, -7.5], "the module's original form"),
    "vmax ~ CrtYB mRNA":          (_by_mrna, [-7.7, -7.5], "measured enzyme proxy"),
    "vmax per strain, growth-independent":
                                  (_by_strain, [-8.7, -8.3, -8.1, -7.5],
                                   "gene dosage, unconstrained, and no growth term"),
    "capacity per strain x mu":   (_by_strain_times_growth, [-6.0, -6.0, -6.0, -7.5],
                                   "gene dosage with the winner's growth term"),
    "vmax ~ mu":                  (_growth_proportional, [-6.0, -7.5], "growth-rate-proportional capacity"),
    "vmax ~ mu^beta":             (_growth_power, [-8.2, -7.5, 1.0], "the above with the exponent freed"),
    "vmax ~ mRNA x mu^beta":      (_mrna_times_growth_power, [-7.7, -7.5, 1.0], "both, together"),
    "vmax ~ mu, no saturation":   (_growth_proportional_unsaturated, [-6.0], "first order in lycopene"),
}


def predict(name: str, theta: np.ndarray, frame: pd.DataFrame) -> np.ndarray:
    """Specific beta-carotene rate, mmol/gDCW/h, under one candidate."""
    vmax, km = CANDIDATES[name][0](theta, frame)
    lyc = frame.lycopene_content.to_numpy()
    if km is None:
        return vmax * lyc
    return vmax * lyc / (km + lyc)


def fit(name: str, frame: pd.DataFrame, restarts: int = 80):
    """Least squares on log rate, multi-start because the surface has flat arms.

    Multi-start is not decoration. Several candidates have a ridge along which km falls to
    zero while vmax compensates, and a single start from one guess lands on it and reports a
    saturation constant of 1e-19 as though it meant something.
    """
    target = np.log(frame.q_betacarotene.to_numpy())
    start = np.asarray(CANDIDATES[name][1], dtype=float)
    rng = np.random.default_rng(0)

    def residual(theta):
        return np.log(np.maximum(predict(name, theta, frame), 1e-18)) - target

    best = None
    with np.errstate(all="ignore"):
        for attempt in range(restarts):
            guess = start if attempt == 0 else start + rng.normal(0.0, 2.0, start.size)
            try:
                found = least_squares(
                    residual, guess, method="lm",
                    xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000,
                )
            except (ValueError, np.linalg.LinAlgError):
                continue
            if not np.all(np.isfinite(found.fun)):
                continue
            if best is None or np.sum(found.fun ** 2) < np.sum(best.fun ** 2):
                best = found
    if best is None:
        raise RuntimeError(f"no start converged for candidate {name!r}")
    return best


# --------------------------------------------------------------------------- #
# the refutation, which comes first
# --------------------------------------------------------------------------- #

def sign_test_against_fixed_vmax(frame: pd.DataFrame) -> pd.DataFrame:
    """Within each strain, does more lycopene buy more beta-carotene?

    Michaelis-Menten with a fixed vmax says the rate is non-decreasing in the substrate
    concentration, whatever vmax and km are. Each strain was run at two dilution rates, so
    each strain is a paired two-point test of exactly that, with the enzyme complement held
    as fixed as an experiment can hold it. No fitting is involved and no parameter can
    rescue a violation, which is why this runs before anything is fitted.
    """
    rows = []
    for strain, group in frame.groupby("strain"):
        low, high = group.sort_values("lycopene_content").iloc[0], group.sort_values("lycopene_content").iloc[-1]
        rows.append({
            "strain": strain,
            "lower [lycopene] at mu": low.mu_per_h,
            "higher [lycopene] at mu": high.mu_per_h,
            "[lycopene] ratio": high.lycopene_content / low.lycopene_content,
            "q_bcar ratio": high.q_betacarotene / low.q_betacarotene,
            "MM violated": bool(high.q_betacarotene < low.q_betacarotene),
        })
    return pd.DataFrame(rows)


def upstream_is_unidentifiable(frame: pd.DataFrame) -> pd.DataFrame:
    """Show that two very different desaturases reproduce the data equally well.

    Not an argument, a demonstration. The desaturase flux is pinned by the measurements at
    ``q_lycopene + q_betacarotene``, but its substrate -- phytoene -- was never measured, so
    for any (vmax, km) with vmax above that flux there is a phytoene pool that delivers it
    exactly. This tabulates the pool each candidate would need. Nothing observable
    distinguishes them, which is why `calibrated_kinetics` refuses to pick one.
    """
    flux = frame.desaturase_flux.to_numpy()
    rows = []
    for vmax_multiple in (1.5, 10.0, 1000.0):
        vmax = vmax_multiple * flux.max()
        for km in (1e-5, 1e-3):
            # v = vmax * p / (km + p)  =>  p = km * v / (vmax - v)
            pool = km * flux / (vmax - flux)
            rows.append({
                "vmax_crti": vmax, "km_crti": km,
                "implied [phytoene] min": pool.min(), "implied [phytoene] max": pool.max(),
                "reproduces every state": bool(np.allclose(vmax * pool / (km + pool), flux)),
            })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #

def leave_one_out(name: str, frame: pd.DataFrame) -> np.ndarray:
    """Predicted q_betacarotene for each state, from a fit that never saw it."""
    predicted = []
    for i in range(len(frame)):
        train, test = frame.drop(frame.index[i]), frame.iloc[[i]]
        if name == "training mean rate":
            predicted.append(train.q_betacarotene.mean())
        elif name == "training mean content":
            predicted.append(train.betacarotene_content.mean() * test.mu_per_h.to_numpy()[0])
        elif name == "training mean ratio":
            ratio = (train.q_betacarotene / train.q_lycopene).mean()
            predicted.append(ratio * test.q_lycopene.to_numpy()[0])
        else:
            predicted.append(predict(name, fit(name, train).x, test)[0])
    return np.asarray(predicted, dtype=float)


BASELINES = ("training mean rate", "training mean content", "training mean ratio")


def score_all(frame: pd.DataFrame) -> pd.DataFrame:
    observed = frame.q_betacarotene.to_numpy()
    rows = []
    for name in (*BASELINES, *CANDIDATES):
        predicted = leave_one_out(name, frame)
        log_error = np.log(predicted) - np.log(observed)
        if name in BASELINES:
            free, rmse_train = 1, np.nan
        else:
            found = fit(name, frame)
            free = found.x.size
            rmse_train = np.sqrt(np.sum(found.fun ** 2) / len(frame))
        rows.append({
            "model": name,
            "free params": free,
            "RMSE(log) train": rmse_train,
            "RMSE(log) LOO": np.sqrt(np.mean(log_error ** 2)),
            "MAPE LOO %": 100 * np.mean(np.abs(predicted - observed) / observed),
            "worst rel. err %": 100 * np.max(np.abs(predicted - observed) / observed),
        })
    table = pd.DataFrame(rows)
    best_baseline = table.loc[table.model.isin(BASELINES), "RMSE(log) LOO"].min()
    table["skill vs best baseline"] = 1 - (table["RMSE(log) LOO"] / best_baseline) ** 2
    return table


def nested_f_test(restricted: str, full: str, frame: pd.DataFrame) -> dict:
    """Does the larger model buy enough fit to pay for its extra parameters?

    ADDED 2026-09-04, and it is what the four documents quoting "-0.58" should have been
    quoting. `capacity per strain x mu` nests `vmax ~ mu` exactly -- set the three log
    capacities equal and it IS the smaller model -- so the comparison is an F test and not a
    matter of preference, and it does not depend on which held-out split the LOO happened to
    draw. With six points and four parameters the denominator has two degrees of freedom,
    which is a weak test; that weakness is the honest reason to pool, and it is a different
    reason from "the per-strain model scores worse".
    """
    rows = {}
    for label, name in (("restricted", restricted), ("full", full)):
        found = fit(name, frame)
        rows[label] = {"model": name, "params": int(found.x.size),
                       "rss": float(np.sum(found.fun ** 2))}
    n = len(frame)
    df_num = rows["full"]["params"] - rows["restricted"]["params"]
    df_den = n - rows["full"]["params"]
    if df_num <= 0 or df_den <= 0:
        raise ValueError(f"{full!r} does not nest {restricted!r} with spare degrees of freedom")
    statistic = (((rows["restricted"]["rss"] - rows["full"]["rss"]) / df_num)
                 / (rows["full"]["rss"] / df_den))
    return {**rows, "df_num": df_num, "df_den": df_den, "f": float(statistic),
            "p": float(f_dist.sf(statistic, df_num, df_den))}


def per_strain_capacities(frame: pd.DataFrame) -> pd.DataFrame:
    """The fitted capacity per strain, beside the measured CrtYB transcript for that strain.

    The measurement that replaces "both cannot be true, flagged not modelled". If crtYB
    dosage set the branch capacity, these three would spread like the dosage does; what they
    do instead is a number, and it is small.
    """
    theta = fit("capacity per strain x mu", frame).x
    rows = []
    for index, strain in enumerate(STRAINS):
        group = frame[frame.strain == strain]
        rows.append({"strain": strain,
                     "capacity_mmol_per_gdcw": float(np.exp(theta[index])),
                     "mean_crtyb_mrna": float(group.crtyb_mrna.mean())})
    return pd.DataFrame(rows)


def uncertainty(name: str, frame: pd.DataFrame, draws: int = 2000) -> dict:
    """Two independent uncertainty estimates for the winning candidate, plus the LOO spread.

    They answer different questions and the wider one is the one to quote. The asymptotic
    interval propagates the *scatter of the fit* -- how badly the two-parameter law misses
    six points -- and the Monte Carlo propagates the *measurement error* Elizondo reports.
    The first is roughly twice the second here, which says the residual is model error
    rather than measurement error, and that is worth knowing before anyone tightens the
    confidence interval by measuring more replicates.
    """
    found = fit(name, frame)
    n, k = len(frame), found.x.size
    rss = float(np.sum(found.fun ** 2))
    covariance = (rss / (n - k)) * np.linalg.inv(found.jac.T @ found.jac)
    se_log = np.sqrt(np.diag(covariance))
    crit = student_t.ppf(0.975, n - k)
    asymptotic = [
        (float(np.exp(v)), float(np.exp(v - crit * s)), float(np.exp(v + crit * s)))
        for v, s in zip(found.x, se_log)
    ]

    sd_lyc = (frame.q_lycopene_hi95 - frame.q_lycopene_lo95) / (2 * 1.959964)
    sd_bcar = (frame.q_betacarotene_hi95 - frame.q_betacarotene_lo95) / (2 * 1.959964)
    rng = np.random.default_rng(7)
    sampled = []
    for _ in range(draws):
        perturbed = frame.copy()
        perturbed["q_lycopene"] = np.maximum(rng.normal(frame.q_lycopene, sd_lyc), 1e-9)
        perturbed["q_betacarotene"] = np.maximum(rng.normal(frame.q_betacarotene, sd_bcar), 1e-9)
        perturbed["lycopene_content"] = perturbed.q_lycopene / perturbed.mu_per_h
        perturbed["desaturase_flux"] = perturbed.q_lycopene + perturbed.q_betacarotene
        sampled.append(np.exp(fit(name, perturbed, restarts=8).x))
    sampled = np.asarray(sampled)

    dropped = []
    for i in range(len(frame)):
        dropped.append(np.exp(fit(name, frame.drop(frame.index[i])).x))
    dropped = np.asarray(dropped)

    return {
        "asymptotic": asymptotic,
        "monte_carlo": np.percentile(sampled, [2.5, 50, 97.5], axis=0),
        "leave_one_out": dropped,
        "residual_sd_log": float(np.sqrt(rss / (n - k))),
        "log_correlation": float(covariance[0, 1] / np.prod(se_log)),
    }


def score_entry_models(frame: pd.DataFrame | None = None, *, draws: int = 0,
                       seed: int = 0) -> pd.DataFrame:
    """Shared nested forward validation, not the conditional state-LOO branch diagnostic."""
    states = carotenoid_measurements() if frame is None else frame
    return score_product_validation(load_pathway("beta_carotene"), states, draws=draws, seed=seed)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--csv", action="store_true",
                        help="accepted and ignored; the tables are always written now")
    parser.add_argument("--entry-only", action="store_true",
                        help="run only nested forward gene/model validation")
    parser.add_argument("--entry-draws", type=int, default=0,
                        help="conditional two-stage draws for the nested entry score")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=pathlib.Path)
    args = parser.parse_args()
    directory = args.output_dir if args.output_dir is not None else paths.outputs_dir()
    directory.mkdir(parents=True, exist_ok=True)
    forward = score_entry_models(draws=args.entry_draws, seed=args.seed)
    forward_summary = summarize_product_validation(forward)
    print("Nested forward gene AND model selection: three strains, six condition predictions, both channels")
    print(forward[["state", "strain", "selected_model", "status"]].to_string(index=False))
    print(forward_summary[["model", "n_independent_strains", "n_condition_predictions",
                           "n_failed_predictions", "product_rmse_log", "lycopene_rmse_log",
                           "joint_rmse_log"]].to_string(index=False))
    # This forward validation is the same call fit_pathway_flux makes; its four tables were
    # byte-identical twins. They live at outputs/pathway_flux_{nested,scores,comparisons}.csv.
    if args.entry_only:
        print(f"wrote nested entry validation to {directory}")
        return 0

    frame = calibration_states()
    mrna = pd.read_csv(
        paths.data_dir() / "carotenoid" / "elizondo2025_relative_mrna.tsv", sep="\t"
    )
    crtyb = mrna[mrna.gene == "CrtYB"].set_index("condition").rel_expression
    frame["crtyb_mrna"] = frame.condition.map(crtyb)
    if frame.crtyb_mrna.isna().any():
        raise SystemExit("CrtYB rows missing from elizondo2025_relative_mrna.tsv")
    if not np.isclose(frame.mu_per_h.min(), REFERENCE_GROWTH_RATE, rtol=1e-6):
        raise SystemExit(
            f"REFERENCE_GROWTH_RATE is {REFERENCE_GROWTH_RATE} but the calibration set's "
            f"lowest dilution rate is {frame.mu_per_h.min()}. The growth term would be "
            "normalised against a rate nothing was measured at."
        )

    scale = 1e6  # mmol -> nmol, the units the paper prints
    pd.set_option("display.width", 160)

    print("The six steady states (rates nmol/gDCW/h, contents nmol/gDCW)")
    shown = frame[[
        "condition", "strain", "mu_per_h", "q_glucose", "q_ethanol",
        "q_lycopene", "q_betacarotene", "lycopene_content", "betacarotene_content",
        "desaturase_flux", "cyclase_fraction", "crtyb_mrna",
    ]].copy()
    for column in ("q_lycopene", "q_betacarotene", "lycopene_content",
                   "betacarotene_content", "desaturase_flux"):
        shown[column] *= scale
    print(shown.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    print("\nRefutation first: is q_betacarotene non-decreasing in lycopene content?")
    print(sign_test_against_fixed_vmax(frame).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("Michaelis-Menten with a growth-rate-independent vmax requires all three to read False.")

    print("\nWhy the upstream steps are not fitted: every row below reproduces all six states")
    print(upstream_is_unidentifiable(frame).to_string(index=False, float_format=lambda v: f"{v:.3g}"))

    print("\nExploratory conditional branch diagnostics: state-LOO, not new-strain validation")
    table = score_all(frame)
    print(table.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print("Every row varies ONE thing against 'vmax shared', and the row name says which. "
          "'vmax per strain, growth-independent' varies two and is named for both -- reading "
          "it as a test of gene dosage attributes its loss to the wrong factor.")

    print("\nIs the per-strain capacity worth its parameters? Nested F, full data")
    nested = nested_f_test("vmax ~ mu", "capacity per strain x mu", frame)
    print(f"  {nested['restricted']['model']!r}  {nested['restricted']['params']} params, "
          f"RSS {nested['restricted']['rss']:.6f}")
    print(f"  {nested['full']['model']!r}  {nested['full']['params']} params, "
          f"RSS {nested['full']['rss']:.6f}")
    print(f"  F({nested['df_num']}, {nested['df_den']}) = {nested['f']:.4f}, "
          f"p = {nested['p']:.4f} -- the extra capacities do not pay for themselves")

    strains = per_strain_capacities(frame)
    spread = strains.capacity_mmol_per_gdcw.max() / strains.capacity_mmol_per_gdcw.min()
    dosage_corr = float(np.corrcoef(np.log(strains.mean_crtyb_mrna),
                                    np.log(strains.capacity_mmol_per_gdcw))[0, 1])
    print("\nThe per-strain capacities themselves, against the measured CrtYB transcript")
    print(strains.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print(f"  spread {spread:.4f}x across strains whose crtYB copy number is 1:2:3, "
          f"and Pearson(log mRNA, log capacity) = {dosage_corr:+.4f}")

    winner = table.loc[table["RMSE(log) LOO"].idxmin(), "model"]
    if winner in BASELINES:
        print(f"\nNo candidate beats the baselines. Winner is {winner!r}; nothing is fitted.")
        return 1

    print(f"\nWinning candidate: {winner!r} -- {CANDIDATES[winner][2]}")
    stats = uncertainty(winner, frame)
    names = ["capacity (nmol/gDCW)", "km (nmol/gDCW)", "exponent"][: len(stats["asymptotic"])]
    print(f"  residual sd in log rate: {stats['residual_sd_log']:.4f}")
    print(f"  correlation of the two log parameters: {stats['log_correlation']:+.3f}")
    for i, label in enumerate(names):
        point, lo, hi = stats["asymptotic"][i]
        mc_lo, mc_mid, mc_hi = stats["monte_carlo"][:, i]
        loo_lo, loo_hi = stats["leave_one_out"][:, i].min(), stats["leave_one_out"][:, i].max()
        unit = scale if "nmol" in label else 1.0
        print(f"  {label}: {point * unit:.1f}")
        print(f"      fit scatter, 95%:        [{lo * unit:.1f}, {hi * unit:.1f}]")
        print(f"      measurement error, 95%:  [{mc_lo * unit:.1f}, {mc_hi * unit:.1f}]  (median {mc_mid * unit:.1f})")
        print(f"      across the six LOO refits: [{loo_lo * unit:.1f}, {loo_hi * unit:.1f}]")

    predicted = leave_one_out(winner, frame)
    held = pd.DataFrame({
        "condition": frame.condition,
        "held-out q_bcar": frame.q_betacarotene * scale,
        "predicted": predicted * scale,
        "rel. err %": 100 * (predicted - frame.q_betacarotene) / frame.q_betacarotene,
        "baseline": leave_one_out("training mean rate", frame) * scale,
    })
    print("\nConditional state-LOO diagnostics: sibling strains retained, winner selected on these data")
    print(held.to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    # ALWAYS WRITTEN, as of 2026-09-04. This table is quoted by docs/research/KINETIC_FIT.md,
    # docs/MEASUREMENTS_NEEDED.md, docs/WHAT_IS_LEFT.md and tests/test_kinetic_carotenoid.py,
    # and while it was a `--csv` opt-in none of those four had a cell to be checked against --
    # which is how "a per-strain capacity scores -0.58" survived in all four after the row it
    # described turned out to be measuring growth-rate-independence instead.
    destination = directory / "carotenoid_fit.csv"
    table.to_csv(destination, index=False)
    held.to_csv(directory / "carotenoid_heldout.csv", index=False)
    parsimony = pd.DataFrame([
        {"quantity": "nested F statistic", "value": round(nested["f"], 4),
         "note": f"{nested['full']['model']} against {nested['restricted']['model']}"},
        {"quantity": "nested F p", "value": round(nested["p"], 4),
         "note": f"F({nested['df_num']}, {nested['df_den']})"},
        {"quantity": "rss restricted", "value": round(nested["restricted"]["rss"], 6),
         "note": f"{nested['restricted']['model']}, {nested['restricted']['params']} params"},
        {"quantity": "rss full", "value": round(nested["full"]["rss"], 6),
         "note": f"{nested['full']['model']}, {nested['full']['params']} params"},
        {"quantity": "capacity spread across strains", "value": round(float(spread), 4),
         "note": "max/min of the three fitted per-strain capacities"},
        {"quantity": "pearson log mrna vs log capacity", "value": round(dosage_corr, 4),
         "note": "measured CrtYB transcript against fitted capacity, n = 3"},
    ])
    parsimony.to_csv(directory / "carotenoid_parsimony.csv", index=False)
    strains.to_csv(directory / "carotenoid_strain_capacity.csv", index=False)
    print(f"\nwrote {destination}, carotenoid_heldout.csv, carotenoid_parsimony.csv "
          f"and carotenoid_strain_capacity.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
