"""Every product in this repository, scored -- or the exact reason it cannot be.

    python3 scripts/score_all_products.py

"Does it work on more than one compound" is the question that decides whether this is a
model or a curve fit, and the answer has two halves that are usually reported separately.
This prints both in one place: what is scored, and for everything else the specific column
that is missing. A negative that can be re-run is worth more than a paragraph, which is the
same argument `tests/test_tal_pathway.py` makes for pinning the TAL rejection as tests.

**What a real score requires**, and why so few clear it. The chain is
``genotype -> flux -> X = v_in/mu -> product``, so scoring it needs three columns measured
on the same steady state: the entry enzyme's **expression**, the product's **q**, and
**mu**. Titre will not substitute for q -- it is time-integrated and confounded by biomass
and duration -- and mu is the denominator, not a covariate. `docs/SECOND_DATASET_HUNT.md`
searched eight product classes and found no second dataset with all three.

**The trap this script exists to avoid.** PHB looks like an eleven-state validation across
three carbon sources. It is not: its `q` column is `content x mu` to five decimal places, so
`content = q/mu` recovers the input and would report a perfect fit on arithmetic. That check
is run below and printed as the refusal it is, because a fit that cannot fail is the most
convincing wrong answer available here.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from ystwin import paths
from ystwin.pathway.spec import load_pathway

REPO = paths.REPO_ROOT


def scored_products() -> list[dict]:
    """The products `scripts/predict_product.py` scores leave-one-strain-out."""
    table = pd.read_csv(paths.outputs_dir() / "product_prediction.csv")
    rows = []
    for label, measured, predicted in (
            ("beta_carotene", "product_measured", "product_predicted"),
            ("lycopene", "lycopene_measured", "lycopene_predicted")):
        if measured not in table.columns:
            continue
        fold = table[predicted] / table[measured]
        rows.append({
            "product": label, "status": "SCORED",
            "n_states": len(table),
            "median_abs_error_pct": float(np.median(np.abs(fold - 1.0)) * 100),
            "worst_fold": float(max(fold.max(), 1.0 / fold.min())),
            "detail": "leave-one-strain-out; no product measurement of the held-out "
                      "strain enters anywhere",
        })
    return rows


def phb_is_circular() -> dict:
    """PHB has content, q and mu on eleven states -- and its q is derived from its content."""
    frame = pd.read_csv(REPO / "data" / "phb" / "kocharin2013_chemostat_states.tsv", sep="\t")
    frame = frame[frame.phb_mmol_per_gdcw.notna() & frame.q_phb_mmol_per_gdcw_h.notna()]
    recovered = frame.q_phb_mmol_per_gdcw_h / frame.mu_per_h
    deviation = float(np.max(np.abs(recovered / frame.phb_mmol_per_gdcw - 1.0)))
    return {
        "product": "phb", "status": "CANNOT SCORE", "n_states": len(frame),
        "median_abs_error_pct": float("nan"), "worst_fold": float("nan"),
        "detail": (f"q_phb is content x mu to {deviation:.1e} on all {len(frame)} states, so "
                   "content = q/mu returns its own input. No expression series exists either, "
                   "which is why spec.calibratable is False"),
    }


def spec_refusals() -> list[dict]:
    """Products whose spec loads but cannot yet produce a number, with the refusal's words.

    `glycogen` moved from REFUSED AT LOAD to REFUSED AT SOLVE on 2026-09-01, when
    `solve_pathway` gained a degradation outlet. The distinction is the point: it was
    reported as outside the arithmetic, and it is one measured rate constant short of it.
    """
    from ystwin.pathway.solve import solve_pathway
    rows = []
    for name in ("glycogen",):
        try:
            spec = load_pathway(name)
        except Exception as exc:
            rows.append({"product": name, "status": "REFUSED AT LOAD", "n_states": 0,
                         "median_abs_error_pct": float("nan"), "worst_fold": float("nan"),
                         "detail": " ".join(str(exc).split())[:200]})
            continue
        try:
            solve_pathway(spec, 0.5, 0.1, {})
            status, detail = "SCORED", "solves"
        except Exception as exc:
            status = "REFUSED AT SOLVE"
            detail = " ".join(str(exc).split())[:200]
        rows.append({"product": name, "status": status, "n_states": 0,
                     "median_abs_error_pct": float("nan"), "worst_fold": float("nan"),
                     "detail": detail})
    return rows


#: Products investigated in full and rejected, with the blocking column named. Each is
#: recorded in `docs/SECOND_DATASET_HUNT.md` and, for TAL, in `tests/test_tal_pathway.py`.
NOT_VENDORED = [
    ("triacetic acid lactone", "K. marxianus, PMID 32995271",
     "product axis is a bar chart (18 states, one published number); expression numeric at "
     "30 C only, and the authors report the rank correspondence is lost at 37 and 41 C. TAL "
     "is assayed in the supernatant, so its fate is SECRETED -- since 2026-09-01 that is a "
     "refusal at SOLVE naming the missing export capacity rather than at load, and it was "
     "never the binding blocker here: the product axis would still be a bar chart"),
    ("caffeic acid", "S. cerevisiae, PMID 35300487",
     "titre only, no growth rate; RT-qPCR and titre are both bitmap figures"),
    ("(2S)-naringenin", "E. coli, PMID 35346204", "titre only, n = 1 per titre, no replicates"),
    ("alpha-santalene", "S. cerevisiae, PMC3527295",
     "q is figure-only (Table 1 prints a volumetric Totsant); farnesol has no printed "
     "number anywhere; and there is no expression measurement. It is also SECRETED -- "
     "captured in a dodecane overlay and assayed in the organic phase -- which since "
     "2026-09-01 refuses at SOLVE for want of an export capacity rather than at load. "
     "Three blockers precede that one. Verified against the source 2026-09-01"),
    ("resveratrol", "S. cerevisiae, PMID 26369953",
     "15 chemostats with q per steady state -- one genotype, and the microarray carries no "
     "probes for the heterologous PAL/4CL/VST1"),
]


def main() -> int:
    rows = scored_products() + [phb_is_circular()] + spec_refusals()
    table = pd.DataFrame(rows)

    print("\nProducts this repository can score, and what it costs to add one\n")
    print(table[["product", "status", "n_states", "median_abs_error_pct",
                 "worst_fold", "detail"]]
          .to_string(index=False, max_colwidth=64,
                     float_format=lambda v: "" if pd.isna(v) else f"{v:.1f}"))

    print("\nInvestigated, not vendored -- the blocking column, per candidate:")
    for name, where, why in NOT_VENDORED:
        print(f"  {name} ({where})")
        print(f"      {why}")

    scored = table[table.status == "SCORED"]
    print(f"\n{len(scored)} product(s) scored, both in one chain, one organism, one paper.")
    print("What would add a third, unchanged since docs/SECOND_DATASET_HUNT.md:")
    print("  RT-qPCR of the pathway entry enzyme on >= 5 producing strains held at ONE")
    print("  dilution rate in glucose-limited chemostat, sampled from the same steady state")
    print("  as the product measurement. Every candidate above already has the product side.")

    out = paths.outputs_dir() / "product_scoreboard.csv"
    table.to_csv(out, index=False)
    print(f"\n-> {paths.display_path(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
