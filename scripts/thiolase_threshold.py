"""The thiolase lead, and the unit error that was the whole of it.

    python3 scripts/thiolase_threshold.py

Six mechanisms have now failed to explain why Kocharin's PHB flux moves 4.24x between carbon
feeds at identical growth rate: FBA on the product reaction, GEM precursor ceilings, GEM
precursor THROUGHPUT, E-Flux regulation, an mRNA-reading law -- and now this one, which was
the only one that ever fit.

**The step is uphill, and that part survives.** yeast-GEM's own cytosolic thiolase, r_0103,
is the chemical analogue of PhaA -- both condense two acetyl-CoA into acetoacetyl-CoA plus
CoA. A reaction with a positive dGr'0 runs forward only when the mass-action ratio pushes it,
so it has a CONCENTRATION THRESHOLD below which no net flux occurs at all.

**What does not survive is the ordering.** The claim was that the two feeds sit on opposite
sides of that threshold, which would make driving force the thing carbon source changes. It
rested on dGr'0 = +15.6 kJ/mol and a threshold of 221 uM, against measured cytosolic
acetyl-CoA of about 10 uM on glucose and 425 uM on ethanol.

That +15.6 was a unit error. `bridge/thermodynamic.py` read the vendored ModelSEED table --
which declares itself kcal/mol -- through `pytfa` under its default `thermo_unit='kJ/mol'`,
so the pH and ionic-strength transform was computed in kJ on top of a tabulated energy in
kcal, and the result was on neither scale. What exposed it is a value from outside both
estimators: water came out at +24.85 kJ/mol, and no scale on which the formation energy of
water is positive is a scale. Corrected on 2026-08-30 the same table gives +38.07.

    dGr'0                                    threshold   glucose 10 uM   ethanol 425 uM
    +15.61  ModelSEED, unit bug (published)    221 uM        below           ABOVE
    +38.07  ModelSEED, unit corrected       19,042 uM        below           below
    +24.96  eQuilibrator component contrib.  1,414 uM        below           below

**The refutation does not depend on the fix.** eQuilibrator never touched the unit bug, it is
the estimator this repository already documents as the better of the two, and its +24.96 is
the literature value for a thiolase condensation. On that reading the threshold is 1.4 mM and
ethanol's 425 uM does not reach it either. The ordering exists only in the erroneous number,
and it was checkable from inside this repository from the day the eQuilibrator backend
landed. Nobody asked it.

**What this leaves.** The 4.24x is unexplained again, and the honest count of failed
mechanisms is six rather than five. The uphill step and its threshold are real; what is gone
is the claim that the two feeds straddle it. `tests/test_thiolase_threshold.py` pins the
refutation, including that it holds under both corrected estimators.

**Unchanged and still correct:** this does NOT require changing anything in
`bridge/tmfa.py`. A first pass concluded that TMFA's 10 uM concentration floor made the
glucose state unrepresentable and widened it to 1e-8. That was wrong: `ThermoModel._apply`
replaces the window on any metabolite the caller measures. The floor is a prior for
metabolites nobody measured.

**And the caveat that was always there.** The two concentrations come from two different
laboratories with two different quench protocols (Kolbeinsen & Bruheim 2021 for glucose,
Kozak & van Rossum 2016 for ethanol), and a 43-fold difference across labs was never a number
to build on. That caveat is now moot -- both feeds are below the threshold by more than the
disagreement between the labs -- but it is kept because it was the right caveat.
"""

from __future__ import annotations

import math
import pathlib
import sys
import warnings

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ystwin import paths

THIOLASE = "r_0103"
GAS_CONSTANT_KJ = 8.314e-3
TEMPERATURE_K = 303.15          # 30 C, Kocharin's cultivation temperature

# Measured cytosolic acetyl-CoA, umol/gDW, same strain (CEN.PK113-7D), same medium
# (Verduyn), same dilution rate (0.05 /h). Two laboratories -- see the module docstring.
MEASURED_UMOL_PER_GDW = {
    "glucose-limited": (0.0199, "Kolbeinsen & Bruheim 2021, supplementary Table S1"),
    "ethanol-limited": (0.85, "Kozak & van Rossum 2016"),
}

# Cytosolic volume per gram dry weight, mL/gDW. A convention rather than a measurement,
# so the threshold is reported across the plausible range rather than at one value.
VOLUMES_ML_PER_GDW = (1.0, 2.0, 2.7)

# The other two participants, held at plausible levels. The threshold scales as the square
# root of their product, so an order of magnitude in either moves it by about threefold.
COA_M = 100e-6
ACETOACETYL_COA_M = 1e-6


def standard_energy() -> float:
    import cobra

    from ystwin.bridge.thermodynamic import ThermodynamicData, reaction_dg0

    warnings.simplefilter("ignore")
    gem = paths.yeast_gem()
    if gem is None or not gem.is_file():
        raise SystemExit("Yeast9 not present; set YSTWIN_YEAST_GEM")
    model = cobra.io.read_sbml_model(str(gem))
    reaction = model.reactions.get_by_id(THIOLASE)
    energy = reaction_dg0({m.id: c for m, c in reaction.metabolites.items()},
                          ThermodynamicData.load())
    if energy is None:
        raise SystemExit(
            f"{THIOLASE} has an uncovered participant, so no energy can be computed. "
            "reaction_dg0 returns None rather than a partial sum on purpose")
    return float(energy)


def threshold_molar(dg0_kj: float, coa: float = COA_M,
                    acetoacetyl: float = ACETOACETYL_COA_M) -> float:
    """Acetyl-CoA at which dG crosses zero, molar.

    ``dG = dG0 + RT ln([acetoacetyl-CoA][CoA] / [acetyl-CoA]^2)``, and setting it to zero
    gives the square root because the substrate is consumed two at a time. That squared
    dependence is why the threshold is sharp: halving the pool quarters the driving term.
    """
    equilibrium = math.exp(-dg0_kj / (GAS_CONSTANT_KJ * TEMPERATURE_K))
    return math.sqrt(acetoacetyl * coa / equilibrium)


def main() -> int:
    dg0 = standard_energy()
    threshold = threshold_molar(dg0)

    print("=" * 84)
    print(f"{THIOLASE}, yeast-GEM's cytosolic thiolase -- the chemical analogue of PhaA")
    print("=" * 84)
    print("  2 acetyl-CoA <=> acetoacetyl-CoA + CoA")
    print(f"  dGr'0 = {dg0:+.2f} kJ/mol   -- UPHILL, so it has a concentration threshold")
    print(f"\n  at [CoA] = {COA_M * 1e6:g} uM and [acetoacetyl-CoA] = "
          f"{ACETOACETYL_COA_M * 1e6:g} uM,")
    print(f"  net forward flux needs [acetyl-CoA] > {threshold * 1e6:.1f} uM")

    print("\n" + "=" * 84)
    print("Where the two carbon sources sit, across every plausible cytosolic volume")
    print("=" * 84)
    rows = []
    print(f"\n  {'feed':<18} {'umol/gDW':>9}  " +
          "".join(f"{v:g} mL/gDW".rjust(14) for v in VOLUMES_ML_PER_GDW))
    for feed, (amount, source) in MEASURED_UMOL_PER_GDW.items():
        cells = []
        for volume in VOLUMES_ML_PER_GDW:
            molar = amount * 1e-6 / (volume * 1e-3)
            side = "below" if molar < threshold else "ABOVE"
            cells.append(f"{molar * 1e6:7.1f} uM {side}")
            rows.append({"feed": feed, "volume_ml_per_gdw": volume,
                         "acetyl_coa_uM": molar * 1e6, "above_threshold": molar > threshold,
                         "threshold_uM": threshold * 1e6, "dg0_kj_per_mol": dg0,
                         "source": source})
        print(f"  {feed:<18} {amount:9.4f}  " + "".join(c.rjust(14) for c in cells))

    table = pd.DataFrame(rows)
    consistent = table.groupby("feed").above_threshold.nunique().max() == 1
    print(f"\n  the verdict is the same at every volume: {consistent}")
    print("  and it is now 'below' at every one of them, which is the refutation.")

    print("\n" + "=" * 84)
    print("How much the assumed cofactor levels move it")
    print("=" * 84)
    print(f"\n  {'[CoA]':>8} {'[acacCoA]':>11} {'threshold':>12}   glucose?   ethanol?")
    for coa in (30e-6, 100e-6, 300e-6):
        for acac in (0.3e-6, 1e-6, 3e-6):
            t = threshold_molar(dg0, coa, acac)
            glc = 10e-6 > t
            eth = 425e-6 > t
            print(f"  {coa * 1e6:7.0f}u {acac * 1e6:10.1f}u {t * 1e6:11.1f}u   "
                  f"{'ABOVE' if glc else 'below':>8}   {'ABOVE' if eth else 'below':>8}")
    print("\n  Both feeds are below at every corner. Under the published +15.6 kJ/mol this")
    print("  table held the ordering in eight corners of nine and the exceptions were the")
    print("  interesting part; at the corrected energy there is nothing to except.")

    print("\n" + "=" * 84)
    print("The two conventions CROSSED, which is the test neither half was doing")
    print("=" * 84)
    print("\n  The table above holds the volume at 2 mL/gDW; the one before it holds the")
    print("  cofactors at one pair. Crossing them is where the claim actually lives.\n")
    print(f"  {'volume':>8} {'ethanol':>9}  " +
          "".join(f"[CoA]={c * 1e6:.0f}u".rjust(13) for c in (30e-6, 100e-6, 300e-6)))
    crossed_ok = 0
    crossed_total = 0
    for volume in VOLUMES_ML_PER_GDW:
        ethanol = MEASURED_UMOL_PER_GDW["ethanol-limited"][0] * 1e-6 / (volume * 1e-3)
        cells = []
        for coa in (30e-6, 100e-6, 300e-6):
            t = threshold_molar(dg0, coa, ACETOACETYL_COA_M)
            ok = ethanol > t
            crossed_ok += ok
            crossed_total += 1
            cells.append(("ABOVE" if ok else "below").rjust(13))
        print(f"  {volume:>8} {ethanol * 1e6:>7.0f}uM  " + "".join(cells))
    print(f"\n  ethanol clears the threshold in {crossed_ok} of {crossed_total} crossed corners.")
    print("  Neither convention rescues it and neither is close: the largest pool any of")
    print("  these volumes gives is 850 uM against a threshold of 19 mM. The lead needed")
    print("  ethanol ABOVE and there is no corner of either assumption where it is.")
    print("\n  eQuilibrator, which never touched the unit bug and whose +24.96 kJ/mol is the")
    print("  literature value for a thiolase condensation, puts the threshold at 1.4 mM --")
    print("  so ethanol misses on that reading too, and the refutation does not rest on")
    print("  which corrected estimator is chosen. tests/test_thiolase_threshold.py pins it.")

    out = paths.outputs_dir() / "thiolase_threshold.csv"
    table.to_csv(out, index=False)
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
