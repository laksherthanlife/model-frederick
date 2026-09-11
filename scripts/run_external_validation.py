"""Does the latent stress state transfer on real data? Gasch 2000, GEO GSE18.

Every transfer number in this repository was fitted to `generator/stress_panel.py` and
scored on the same generator. This script runs the same test on measurements: 81 arrays
across 8 stressors from the Gasch et al. 2000 environmental stress compendium, with the
reporter channels read off the arrays as the marker genes the panel names, and the module
activities read off as regulon means under a binding-only, expression-independent
regulatory map from SGD.

Three things are scored, in this order:

1. **Positive controls.** Does each channel rise under the stressor it is supposed to
   read, and does each module's regulon? A transfer number computed on a mapping that
   fails this is a number about nothing, so it is printed first and the failures are
   named rather than dropped.
2. **Transfer.** Leave one stressor out, fit on the rest, score the held-out stressor's
   own target modules with the fixed `module_transfer` metric -- against
   `rotated_subspace` and `matched_marginals` nulls, because a rank-k projection of
   correlated channels transfers something whether or not its axes mean anything.
3. **The redundancy law.** Simulation says recovery tracks how many *other* stressors
   drive the same module (0 peers -> 0%, 6+ peers -> 84%, Spearman +0.59). Real data
   either corroborates that or refutes it.

Plus one test the repository has never had a null for at all: the known-loadings readout
`x = pinv(L) y` against `shuffled_loadings(L)`.

Circularity: `stress_panel.py` cites Gasch 2000 for `diamide` (agent, EC50, ESR arm),
`sorbitol` (agent, EC50) and the `ESR` module. `--independent` drops both stressors, so
the panel used carries nothing traceable to the dataset except the ESR arm, which is
reported separately either way. See docs/EXTERNAL_VALIDATION.md.

Usage:
  python scripts/fetch_external_validation_data.py     # first, once
  python scripts/run_external_validation.py [--independent] [--draws N] [--width K]
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import dataclasses
import gzip
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ystwin import paths  # noqa: E402
from ystwin.analysis.latent import select_dimension  # noqa: E402
from ystwin.analysis.nulls import (  # noqa: E402
    _MIN_DRAWS,
    compare_to_null,
    matched_marginals,
    rotated_subspace,
    shuffled_loadings,
    skill_score_if_defined,
)
from ystwin.analysis.validation import validate_modules  # noqa: E402
# The module -> transcription factor table lives in bridge/regulation.py, which is the
# only definition. It was duplicated here, and a divergence between the two would be
# invisible: both would run and quietly disagree about which genes a module owns. It sat
# 45 lines further down beside the table it replaced, which was the one place in this file
# an import was not with the others.
from ystwin.bridge.regulation import MODULE_FACTORS  # noqa: E402
from ystwin.generator.panel_experiment import PanelDataset  # noqa: E402
from ystwin.generator.stress_panel import (  # noqa: E402
    MODULES,
    STRESSORS,
    reporter_loadings,
)

GASCH = ROOT / "data" / "external" / "gasch2000"
SGD = ROOT / "data" / "external" / "sgd" / "regulation"
OUT = paths.outputs_dir()

# GPL63 carries one array (a rescan of H2O2 40 min) and GEO publishes no annotation table
# for it, so its probes cannot be mapped to genes.
PLATFORMS = ("51", "52", "53", "54", "55", "56", "57")
PLATFORMS_EXTRA = ("64",)

# Sample titles are the only condition labels GSE18 carries, and they are free text. Each
# rule below is a regex over the title, matched in order, case-insensitive. Everything
# unmatched is dropped and listed with a reason in outputs/external_conditions.csv.
CONDITION_RULES = (
    ("menadione", r"menadione", "1 mM menadione time course"),
    ("diamide", r"diamide", "1.5 mM diamide time course"),
    ("DTT", r"^\s*2\.5mM DTT|^\s*dtt \d+ min", "two DTT time courses"),
    ("H2O2", r"constant 0\.32 mM H2O2", "0.32 mM peroxide time course"),
    ("sorbitol", r"^1M sorbitol - \d+ min$", "1 M sorbitol hyper-osmotic shock"),
    ("heat", r"^heat shock \d+ min|^Heat Shock \d+|^heat shock \d+ to 37",
     "25->37 C shock, plus the variable pre-shift ladder"),
    ("rapamycin", r"^Nitrogen Depletion", "nitrogen depletion; matched on module, not agent"),
    ("glucose_starvation", r"^YPD \d+ [hd] \(25C\)|^YPD_\d+_[hd]_30C",
     "growth to stationary phase, i.e. glucose exhaustion"),
)

DROP_REASONS = (
    (r"^Hypo-osmotic", "no panel stressor: hypo-osmotic shock has no module"),
    (r"^Amino acid \+ adenine", "no panel module: general amino-acid starvation is Gcn4"),
    (r"vs reference pool", "reference is a pool, not a paired unstressed control"),
    (r"^steady state|deg growth", "adapted steady state, not a shift; cold has no module"),
    (r"^25C --|^33C |^29C ", "mild or combined shift with no matching panel dose"),
    (r"^1M sorbitol vs\. YPD", "adapted steady state, not a shock"),
    (r"centrifugation problem", "flagged by the depositors"),
    (r"^DBY|^FY3|^yap", "genotype perturbation, not a stressor"),
)


# The marker gene each transcriptional reporter's promoter comes from, as written in
# REPORTERS. A slash in the panel means either gene would serve, so both are read and
# averaged. CUP1 is a tandem duplication and SGD names the copies CUP1-1 and CUP1-2.
CHANNEL_GENES = {
    "STRE-general": ("HSP12", "CTT1"), "UPRE-ER": ("KAR2",),
    "TRX2-oxidative": ("TRX2",), "HSE-heat": ("HSP104",), "STRE-osmotic": ("GPD1",),
    "PACE-proteasome": ("RPT1",), "FeRE-iron": ("FIT3",),
    "RLM1-cellwall": ("PIR3", "CRH1"), "Xbox-dna": ("RNR3",),
    "CDRE-calcium": ("PMC1", "CMK2"), "CSRE-carbon": ("ADH2", "ICL1"),
    "GATA-nitrogen": ("DAL5", "MEP2"), "AR1-hypoxia": ("ANB1", "DAN1"),
    "CuRE-copper": ("CUP1-1", "CUP1-2"), "ZRE-zinc": ("ZRT1",),
    "METbox-sulfur": ("MET17",), "PDRE-xenobiotic": ("PDR5",),
    "Rbox-retrograde": ("CIT2",), "RIM101-alkaline": ("ENA1",),
}

# RNR3 is absent from 55% of these arrays and RPT1 from 16%, both because whole print runs
# omit the probe. `train_stress_model` centres with a plain mean, so one missing entry
# poisons a whole channel; dropping the two channels keeps 81 of 96 arrays where keeping
# them would keep 35. Recorded here rather than in a filter, because losing the proteasome
# channel is a real cost and the reader should see it named.
SPARSE_CHANNELS = ("Xbox-dna", "PACE-proteasome")

# MacIsaac et al. 2006, PMID 16522208 -- conserved binding motifs over the Harbison and
# Lee ChIP data. Chosen over "every binding record SGD holds" because the wider set runs
# to 500+ promiscuous ChIP hits per factor and its regulon mean is then the array mean:
# under the wide set the oxidative and DNA-damage modules read 0.01 and -0.03 under
# peroxide. Chosen over anything with expression evidence because expression-derived
# labels would partly be Gasch's own clustering. Gasch 2000 (PMID 11102521) contributes
# zero records to any of these 32 factors, and the script asserts that.
REGULON_SOURCE_PMID = 16522208
GASCH_PMID = 11102521
MIN_REGULON = 5

# Excluded from the circularity-restricted panel. stress_panel.py cites Gasch 2000 for
# both: diamide's agent, its 1.5 mM EC50 and its 0.60 ESR arm, and sorbitol's agent and
# its 1.0 M EC50. Both EC50s are literally the dose Gasch dosed.
GASCH_DERIVED_STRESSORS = ("diamide", "sorbitol")


def rule(title: str) -> None:
    print(f"\n{title}\n{'-' * len(title)}")


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------


def _annotation(gpl: str) -> pd.DataFrame:
    rows, inside = [], False
    with gzip.open(GASCH / f"GPL{gpl}.annot.gz", "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                inside = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if inside:
                rows.append(line.rstrip("\n").split("\t"))
    frame = pd.DataFrame(rows[1:], columns=rows[0])
    return frame[["ID", "Gene symbol", "Platform_ORF"]]


def _series_matrix(gpl: str):
    titles = accessions = None
    rows, inside = [], False
    with gzip.open(GASCH / f"GSE18-GPL{gpl}_series_matrix.txt.gz", "rt",
                   errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accessions = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                inside = True
            elif line.startswith("!series_matrix_table_end"):
                break
            elif inside:
                rows.append(line.rstrip("\n").split("\t"))
    header = [h.strip('"') for h in rows[0]]
    frame = pd.DataFrame(rows[1:], columns=header)
    frame[header[0]] = frame[header[0]].str.strip('"')
    for column in header[1:]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return titles, accessions, frame.rename(columns={header[0]: "ID"})


def load_expression() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    """Gene-by-array log2 ratios, the array titles, and a gene-symbol-to-ORF map."""
    blocks, meta, symbol_to_orf = [], [], {}
    for gpl in PLATFORMS + PLATFORMS_EXTRA:
        annotation = _annotation(gpl)
        for _, row in annotation.iterrows():
            if row["Gene symbol"] and row["Platform_ORF"]:
                symbol_to_orf.setdefault(row["Gene symbol"], row["Platform_ORF"])
        titles, accessions, frame = _series_matrix(gpl)
        merged = frame.merge(annotation[["ID", "Platform_ORF"]], on="ID", how="left")
        merged = merged[merged["Platform_ORF"].notna() & (merged["Platform_ORF"] != "")]
        # Several spots per gene per array; the median of them is the gene's ratio.
        merged = merged.drop(columns=["ID"]).groupby("Platform_ORF").median()
        merged.columns = accessions
        blocks.append(merged)
        meta += [{"gsm": a, "title": t, "platform": f"GPL{gpl}"}
                 for t, a in zip(titles, accessions)]
    expression = pd.concat(blocks, axis=1)
    expression = expression.loc[:, ~expression.columns.duplicated()]
    frame = pd.DataFrame(meta).drop_duplicates("gsm").set_index("gsm")
    return expression, frame, symbol_to_orf


def assign_conditions(meta: pd.DataFrame) -> pd.DataFrame:
    """Attach a panel stressor to every array, or the reason it was dropped."""
    stressors, notes = [], []
    for title in meta["title"]:
        match = next((name for name, pattern, _ in CONDITION_RULES
                      if pd.Series([title]).str.contains(pattern, case=False,
                                                         regex=True).iloc[0]), None)
        stressors.append(match)
        if match is not None:
            notes.append(next(note for name, _, note in CONDITION_RULES if name == match))
        else:
            notes.append(next((reason for pattern, reason in DROP_REASONS
                               if pd.Series([title]).str.contains(pattern, case=False,
                                                                  regex=True).iloc[0]),
                              "unmatched"))
    out = meta.copy()
    out["stressor"] = stressors
    out["note"] = notes
    return out


def load_regulons(expression: pd.DataFrame, excluded: set[str]) -> dict[str, list[str]]:
    """Module gene sets from one binding-only source, with the channel genes removed.

    Removing every channel gene from every regulon -- not just its own module's -- is what
    stops a channel predicting a module because the module's truth contains that channel.

    Raises:
        RuntimeError: if the dataset's own paper turns out to contribute records, which
            would make the labels partly the measurements. A real raise rather than an
            assert, because the claim that it contributes none is load-bearing and
            `python -O` would silently drop the check.
    """
    regulons = {}
    for module, factors in MODULE_FACTORS.items():
        targets = set()
        for factor in factors:
            records = json.loads((SGD / f"{factor}.json").read_text())
            gasch = [r for r in records
                     if r["reference"].get("pubmed_id") == GASCH_PMID]
            if gasch:
                raise RuntimeError(
                    f"{factor}: {len(gasch)} regulation records come from PMID "
                    f"{GASCH_PMID}, the paper supplying the measurements. The module "
                    "labels would not be independent of them.")
            for record in records:
                if record["regulation_of"] != "transcription":
                    continue
                if record["reference"].get("pubmed_id") != REGULON_SOURCE_PMID:
                    continue
                # SGD's per-gene file also holds records where this gene is somebody
                # else's TARGET. Taking locus2 unconditionally then puts the factor into
                # its own regulon: CAT8.json has two conserved-motif records and CAT8 is
                # the target in both, which yielded a one-gene regulon of CAT8 itself.
                if record.get("locus1", {}).get("display_name") != factor:
                    continue
                targets.add(record["locus2"]["format_name"])
        regulons[module] = sorted((targets - excluded) & set(expression.index))
    return regulons


def build(independent: bool):
    """Everything downstream needs: the dataset, the mapping tables and the channel list."""
    expression, meta, symbol_to_orf = load_expression()
    meta = assign_conditions(meta)

    channels = {name: [symbol_to_orf[s] for s in genes
                       if s in symbol_to_orf and symbol_to_orf[s] in expression.index]
                for name, genes in CHANNEL_GENES.items()
                if name not in SPARSE_CHANNELS}
    missing = [name for name, orfs in channels.items() if not orfs]
    if missing:
        raise RuntimeError(f"no array probe for {missing}")
    excluded = {orf for genes in CHANNEL_GENES.values() for symbol in genes
                for orf in [symbol_to_orf.get(symbol)] if orf}
    regulons = load_regulons(expression, excluded)

    meta = meta.copy()
    meta["excluded"] = ""
    unmapped = meta["stressor"].isna()
    meta.loc[unmapped, "excluded"] = meta.loc[unmapped, "note"]
    if independent:
        restricted = meta["stressor"].isin(GASCH_DERIVED_STRESSORS)
        meta.loc[restricted, "excluded"] = "stress_panel.py takes this agent from Gasch 2000"

    candidate = meta[meta["excluded"] == ""]
    block = expression[list(candidate.index)]
    readings = pd.DataFrame({name: block.loc[orfs].mean(axis=0)
                             for name, orfs in channels.items()}, index=candidate.index)
    incomplete = readings.index[~readings.notna().all(axis=1)]
    meta.loc[incomplete, "excluded"] = "a channel gene has no probe on this array's print run"
    kept = meta[meta["excluded"] == ""]
    readings = readings.loc[kept.index]
    block = expression[list(kept.index)]

    measurable = [m for m in MODULES if len(regulons.get(m, [])) >= MIN_REGULON]
    truth = pd.DataFrame(
        {m: (block.loc[regulons[m]].mean(axis=0) if m in measurable
             else pd.Series(0.0, index=kept.index)) for m in MODULES},
        index=kept.index)

    dataset = PanelDataset(
        readings=readings.to_numpy(dtype=float),
        labels=kept["stressor"].to_numpy(),
        doses=np.zeros(len(kept)),
        modules=truth.to_numpy(dtype=float),
        reporters=list(readings.columns))
    return dataset, kept, meta, readings, truth, regulons, measurable


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------


def own_targets(stressor: str, measurable: list[str]) -> list[str]:
    """Modules the panel says this agent acts on directly, and that this data can see."""
    return [m for m in STRESSORS[stressor].targets if m in measurable]


def peer_count(module: str, stressor: str, panel: list[str]) -> int:
    """How many OTHER stressors in this panel drive the same module, per the panel."""
    return sum(1 for other in panel
               if other != stressor and module in STRESSORS[other].targets)


# A module counts as driven when its regulon mean moves this far in log2, i.e. about 19%.
# Used for the data-defined peer count, which does not depend on the panel's asserted
# target lists -- and those turn out to be incomplete, which is the point of having it.
DRIVEN_LOG2 = 0.25


def peer_count_measured(dataset, module: str, stressor: str, panel: list[str]) -> int:
    """The same count, read off the data instead of off the panel's target lists.

    Counts only the other stressors, all of which are in training when this one is held
    out, so nothing here looks at the held-out block.
    """
    index = list(MODULES).index(module)
    hits = 0
    for other in panel:
        if other == stressor:
            continue
        block = dataset.modules[dataset.mask(other), index]
        if block.size and abs(float(np.mean(block))) > DRIVEN_LOG2:
            hits += 1
    return hits


def _validate_workers(workers) -> None:
    if (not isinstance(workers, (int, np.integer)) or isinstance(workers, (bool, np.bool_))
            or workers < 1):
        raise ValueError("workers must be a positive integer")


def transfer_pairs(dataset, panel, measurable, width, seed=0, *, workers: int = 1) -> pd.DataFrame:
    _validate_workers(workers)

    def fold_rows(held):
        validation = validate_modules(dataset, held, width=width, seed=seed)
        scores = validation.scores
        block = dataset.modules[validation.fold.test]
        rows = []
        for module in own_targets(held, measurable):
            index = list(MODULES).index(module)
            rows.append({
                "held_out": held, "module": module, "score": scores[module],
                "n_states": validation.n_states,
                "n_train": len(validation.fold.train), "n_test": len(validation.fold.test),
                "selection": validation.selection,
                "baseline": "outer training module mean",
                "score_method": "1 - SSE / training-mean SSE; negative scores retained",
                "peers": peer_count(module, held, panel),
                "peers_measured": peer_count_measured(dataset, module, held, panel),
                "truth_sd": float(np.std(block[:, index])),
                "truth_span": float(np.ptp(block[:, index])),
            })
        return rows

    if workers == 1 or len(panel) < 2:
        return pd.DataFrame([row for held in panel for row in fold_rows(held)])
    with ThreadPoolExecutor(max_workers=min(workers, len(panel))) as executor:
        return pd.DataFrame([row for rows in executor.map(fold_rows, panel) for row in rows])


def median_pair_score(dataset, panel, measurable, width, seed=0, *, workers: int = 1) -> float:
    table = transfer_pairs(dataset, panel, measurable, width, seed, workers=workers)
    return float(np.nanmedian(table["score"])) if len(table) else float("nan")


def loadings_readout(readings: np.ndarray, loadings: np.ndarray) -> np.ndarray:
    """Path A: modules from the known loading matrix, no fitting.

    ``y = x L'`` with L (channels, modules), so the minimum-norm inverse is
    ``x = y (L+)'``. L has 17 rows against 24 columns, so 7 modules are unreachable by
    construction -- which is the Ledermann-adjacent point `IDENTIFIABILITY.md` makes, met
    here on real readings.
    """
    return readings @ np.linalg.pinv(loadings).T


def loadings_agreement(readings, truth, loadings, measurable) -> float:
    """Median Spearman correlation between the known-loadings estimate and the truth."""
    estimate = loadings_readout(readings, loadings)
    if len(estimate) < 2:
        return float("nan")
    order = list(MODULES)
    values = []
    for module in measurable:
        column = order.index(module)
        if np.ptp(estimate[:, column]) == 0 or np.ptp(truth[:, column]) == 0:
            continue
        rho = stats.spearmanr(estimate[:, column], truth[:, column]).statistic
        if np.isfinite(rho):
            values.append(rho)
    return float(np.median(values)) if values else float("nan")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main(independent: bool, draws: int, width: int | None, output_dir=None,
         *, workers: int = 1) -> None:
    if not isinstance(draws, (int, np.integer)) or draws < _MIN_DRAWS:
        raise ValueError(f"null comparisons require at least {_MIN_DRAWS} draws")
    _validate_workers(workers)
    out = pathlib.Path(output_dir) if output_dir is not None else OUT
    out.mkdir(parents=True, exist_ok=True)
    dataset, kept, meta, readings, truth, regulons, measurable = build(independent)
    panel = sorted(set(dataset.labels))

    rule("what the dataset supplies")
    print(f"  {len(meta)} arrays in GSE18; {len(kept)} used, over {len(panel)} stressors")
    print("  " + kept.groupby("stressor").size().to_string().replace("\n", "\n  "))
    print(f"  channels: {len(dataset.reporters)} of {len(CHANNEL_GENES)} "
          f"(dropped {list(SPARSE_CHANNELS)}: probe absent from whole print runs)")
    print(f"  modules measurable: {len(measurable)} of {len(MODULES)}")
    print("  unmeasurable: " + ", ".join(m for m in MODULES if m not in measurable))
    if independent:
        print(f"  CIRCULARITY-RESTRICTED: {list(GASCH_DERIVED_STRESSORS)} dropped entirely")

    suffix = "_independent" if independent else ""
    conditions = meta.copy()
    conditions["used"] = conditions["excluded"] == ""
    conditions.to_csv(out / f"external_conditions{suffix}.csv")

    module_map = pd.DataFrame([
        {"module": m, "factors": "/".join(MODULE_FACTORS.get(m, ())),
         "regulon_genes": len(regulons.get(m, [])),
         "measurable": m in measurable} for m in MODULES])
    module_map.to_csv(out / f"external_module_map{suffix}.csv", index=False)

    rule("positive controls: does each channel rise under the stressor it reads?")
    channel_table = readings.copy()
    channel_table["stressor"] = kept["stressor"].to_numpy()
    per_channel = channel_table.groupby("stressor").mean().T
    per_channel["argmax"] = per_channel.idxmax(axis=1)
    print(per_channel.round(2).to_string())

    rule("positive controls: does each module's regulon rise under its own stressor?")
    module_table = truth[measurable].copy()
    module_table["stressor"] = kept["stressor"].to_numpy()
    per_module = module_table.groupby("stressor").mean().T
    per_module["argmax"] = per_module.idxmax(axis=1)
    print(per_module.round(2).to_string())
    per_channel.to_csv(out / f"external_channel_response{suffix}.csv")
    per_module.to_csv(out / f"external_module_response{suffix}.csv")

    rule("latent width")
    centred = dataset.readings - dataset.readings.mean(axis=0)
    exploratory_width = select_dimension(centred, max_states=min(8, len(dataset.reporters)))
    chosen = width
    share = np.linalg.svd(centred, compute_uv=False) ** 2
    print(f"  variance per component: {np.round(share / share.sum(), 3)[:6]}")
    print(f"  full-data reconstruction width: {exploratory_width} (descriptive only)")
    sweep = pd.DataFrame([
        {"width": w, "median_pair_score": median_pair_score(dataset, panel, measurable, w,
                                                           workers=workers)}
        for w in range(1, min(6, len(dataset.reporters)) + 1)])
    print(sweep.round(4).to_string(index=False))
    print("  This sweep is descriptive, not a selector. Negative skill is not clipped.")
    print("  Width is selected within each outer training fold unless --width was prespecified.")

    rule(f"transfer: leave one stressor out, width {chosen if chosen is not None else 'nested'}, own-target modules only")
    table = transfer_pairs(dataset, panel, measurable, chosen, workers=workers)
    print(table.round(4).to_string(index=False))
    by_stressor = table.groupby("held_out").agg(
        pairs=("score", "size"), median=("score", "median"), best=("score", "max"),
        recovered=("score", lambda s: int((s > 0.25).sum())))
    print("\n" + by_stressor.round(4).to_string())
    observed = float(np.nanmedian(table["score"]))
    print(f"\n  median over all {len(table)} (held-out stressor, own target) pairs: "
          f"{observed:+.4f}")
    print(f"  pairs clearing the 0.25 recovery threshold: "
          f"{int((table['score'] > 0.25).sum())} of {len(table)}")

    rule("and is that better than no structure at all?")
    def statistic(data) -> float:
        return median_pair_score(data, panel, measurable, chosen, workers=workers)

    def rotate(data, rng):
        return dataclasses.replace(data, readings=rotated_subspace(data.readings, rng))

    def scramble(data, rng):
        return dataclasses.replace(data, readings=matched_marginals(data.readings, rng))

    null_rows = []
    for surrogate, label in ((rotate, "rotated_subspace"),
                             (scramble, "matched_marginals")):
        result = compare_to_null(statistic, dataset, surrogate, n_draws=draws,
                                 greater_is_better=True, label=label, seed=0)
        print("  " + result.summary())
        null_rows.append({"test": "transfer (median own-target score)", "null": label,
                          "observed": result.observed, "null_median": result.null_median,
                          "p_value": result.p_value, "draws": result.n_draws,
                          "beats_null": result.beats_null})
    rotated = null_rows[0]["null_median"]
    skill = skill_score_if_defined(observed, rotated, perfect=1.0)
    print(f"  skill over the rotated null: "
          f"{'undefined' if skill is None else f'{skill:+.3f}'} "
          "(1.0 closes the gap to perfect, negative is worse than no structure)")

    rule("the redundancy law: does peer count predict recovery on real data?")
    finite = table.dropna(subset=["score"])
    summaries = []
    for column, tag in (("peers", "panel target lists"),
                        ("peers_measured", "measured regulon movement")):
        bins = pd.cut(table[column], [-1, 0, 2, 5, 99], labels=["0", "1-2", "3-5", "6+"])
        summary = table.groupby(bins, observed=False).agg(
            pairs=("score", "size"), mean_score=("score", "mean"),
            recovered=("score", lambda s: float((s > 0.25).mean()) if len(s) else np.nan))
        rho = stats.spearmanr(finite[column], finite["score"])
        print(f"\n  peers from {tag}:")
        print("  " + summary.round(3).to_string().replace("\n", "\n  "))
        print(f"  Spearman = {rho.statistic:+.3f}, p = {rho.pvalue:.3f}, n = {len(finite)}")
        summary["peer_definition"] = tag
        summaries.append(summary)

        def permuted(data, rng, key=column):
            shuffled = data.copy()
            shuffled[key] = rng.permutation(shuffled[key].to_numpy())
            return shuffled

        law = compare_to_null(
            lambda d, key=column: float(stats.spearmanr(d[key], d["score"]).statistic),
            finite, permuted, n_draws=max(draws, 200), greater_is_better=True,
            label=f"peers shuffled ({tag})", seed=0)
        print("  " + law.summary())
        null_rows.append({"test": f"redundancy law, peers from {tag}",
                          "null": "peers shuffled", "observed": law.observed,
                          "null_median": law.null_median, "p_value": law.p_value,
                          "draws": law.n_draws, "beats_null": law.beats_null})

    # The sharper question. A permutation null says the relation is not chance. It does
    # not say the relation is evidence for the FITTED state, because the transfer scores
    # themselves do not beat arbitrary axes. So run the law on rotated readings: if the
    # relation survives there, it is a fact about which modules the panel covers twice,
    # not about the latent basis.
    def law_statistic(data) -> float:
        pairs = transfer_pairs(data, panel, measurable, chosen, workers=workers).dropna(subset=["score"])
        return float(stats.spearmanr(pairs["peers"], pairs["score"]).statistic)

    law_rotated = compare_to_null(law_statistic, dataset, rotate, n_draws=draws,
                                  greater_is_better=True,
                                  label="the law under rotated axes", seed=0)
    print("\n  " + law_rotated.summary())
    null_rows.append({"test": "redundancy law, peers from panel target lists",
                      "null": "rotated_subspace", "observed": law_rotated.observed,
                      "null_median": law_rotated.null_median,
                      "p_value": law_rotated.p_value, "draws": law_rotated.n_draws,
                      "beats_null": law_rotated.beats_null})
    pd.concat(summaries).to_csv(out / f"external_redundancy{suffix}.csv")

    rule("the other null the repository has never run: do these channels identify "
         "these modules?")
    loadings = reporter_loadings(dataset.reporters)
    print(f"  L is {loadings.shape[0]} channels x {loadings.shape[1]} modules, "
          f"rank {np.linalg.matrix_rank(loadings)}")
    agreement = compare_to_null(
        lambda matrix: loadings_agreement(dataset.readings, dataset.modules, matrix,
                                          measurable),
        loadings, shuffled_loadings, n_draws=max(draws, 200), greater_is_better=True,
        label="shuffled_loadings", seed=0)
    print("  " + agreement.summary())
    null_rows.append({"test": "known-loadings readout (median Spearman)",
                      "null": "shuffled_loadings", "observed": agreement.observed,
                      "null_median": agreement.null_median, "p_value": agreement.p_value,
                      "draws": agreement.n_draws, "beats_null": agreement.beats_null})

    nulls = pd.DataFrame(null_rows)
    nulls.to_csv(out / f"external_nulls{suffix}.csv", index=False)
    table.to_csv(out / f"external_transfer_pairs{suffix}.csv", index=False)
    print(f"\nwrote {out}/external_*{suffix}.csv: conditions, module_map, "
          "channel_response, module_response, transfer_pairs, redundancy, nulls")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--independent", action="store_true",
                        help="drop the stressors stress_panel.py takes from Gasch 2000")
    parser.add_argument("--draws", type=int, default=200)
    parser.add_argument("--width", type=int, default=None)
    parser.add_argument("--output-dir", type=pathlib.Path, default=None)
    parser.add_argument("--workers", type=int, default=1,
                        help="maximum threads for independent outer folds (default: 1)")
    args = parser.parse_args()
    main(args.independent, args.draws, args.width, args.output_dir, workers=args.workers)
