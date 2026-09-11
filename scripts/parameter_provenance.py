"""Every constant the prediction rests on, where it came from, and what would upgrade it.

    python3 scripts/parameter_provenance.py
    python3 scripts/parameter_provenance.py --grade measured   # filter to one grade

The discipline already exists everywhere in this package -- `calibrated_kinetics` raises and
names the missing HPLC peak, `calibrate.py` reports `gain x gdcw_per_od` as one product and
lists the pair as unidentifiable, `solve_pathway` refuses a degraded node without its rate.
What did not exist is a single place that answers **"what do we actually know, and what
would one experiment fix"**, which is the question anybody deciding what to run next is
holding.

Five grades, and the distinction between the middle three is the useful part:

* ``MEASURED``  -- fitted or read from data this repository holds, with the n it rests on.
* ``BOUNDED``   -- not measured, but constrained from data already in hand rather than
  assumed. `minimum_consistent_kdeg` is the model: a free lower bound from an ordinary
  growth trace, and the right prior until a chase is run.
* ``BORROWED``  -- measured, by somebody else, on a culture that is not this one. Numerically
  real and possibly wrong here, which is worse than a guess because it looks solid.
* ``ASSERTED``  -- a value chosen because the arithmetic needs one. Tier 0.
* ``REFUSED``   -- no value, and the code raises rather than inventing one.

A parameter's grade is not its quality. `biomass_yield` is BORROWED and precise; the
carotenoid ceiling is MEASURED and refuted by 63x. Read the note, not the badge.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ystwin import paths
from ystwin.generator.plate import PlateConditions
from ystwin.kinetic.carotenoid import UNIDENTIFIABLE_STEPS
from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS

MEASURED, BOUNDED, BORROWED, ASSERTED, REFUSED = (
    "MEASURED", "BOUNDED", "BORROWED", "ASSERTED", "REFUSED")


def inventory() -> pd.DataFrame:
    """One row per constant, with its live value read from the code where it has one."""
    flux = BETA_CAROTENE_FLUX
    cyclase = BETA_CAROTENE_KINETICS["lycopene"]
    ceiling = cyclase.vmax_per_growth * 536.87
    rows = [
        dict(parameter="alpha (expression -> flux)", value=f"{flux.alpha:.6g}",
             unit="mmol gDCW-1 h-1 per unit", grade=MEASURED,
             source=f"fitted on {flux.n_states} chemostat states, {flux.entry_enzyme}, "
                    "Elizondo 2025 PMID 40891387",
             upgrade="6 producing strains at one dilution rate with entry-enzyme qPCR; "
                     "size computed in scripts/design_flux_experiment.py"),
        dict(parameter="leave-one-strain-out spread", value=f"{flux.loso_rmse_log:.4g}",
             unit="log", grade=MEASURED,
             source="the same fit, held out by strain",
             upgrade="nothing: it is a property of the six states, not a knob. It shrinks only if the product assay gets more precise -- see scripts/error_budget.py"),
        dict(parameter="cyclase vmax/mu", value=f"{cyclase.vmax_per_growth:.6g}",
             unit="mmol gDCW-1", grade=MEASURED,
             source="fitted on the same six states",
             upgrade="content on 3+ strains differing in crtYB dosage alone"),
        dict(parameter="cyclase Km", value=f"{cyclase.km:.6g}", unit="mmol gDCW-1",
             grade=MEASURED, source="fitted on the same six states",
             upgrade="the same crtYB dosage series as the vmax above -- the two are fitted "
                     "together and neither moves without it"),
        dict(parameter="content ceiling", value=f"{ceiling:.4f}", unit="mg gDCW-1",
             grade=MEASURED,
             source="falls out of vmax/mu; mu cancels, so it binds every strain at every rate",
             upgrade="REFUTED, not upgradeable: Arhar 2024 PMID 39215465 measured 79 mg/gDCW "
                     "by HPLC, about 63x this. The fit saw only entry-enzyme variation"),
        dict(parameter=f"upstream kinetics ({', '.join(UNIDENTIFIABLE_STEPS)})", value="none",
             unit="-", grade=REFUSED,
             source="calibrated_kinetics() raises: phytoene and GGPP were never measured, so "
                    "any (vmax, km) pair reproduces the data",
             upgrade="one phytoene peak on HPLC runs that already exist -- M8, run order 1"),
        dict(parameter="growth rate mu", value="per condition", unit="h-1", grade=MEASURED,
             source="the pump in a chemostat; d(ln OD)/dt in a flask, via "
                    "growth.specific_growth_rate",
             upgrade="nothing to buy: it is measurable on any plate you already run. The constraint is not knowing mu, it is landing inside 0.101-0.254 /h"),
        dict(parameter="reporter loss k_deg", value="lower bound", unit="h-1",
             grade=BOUNDED,
             source="minimum_consistent_kdeg: the smallest loss under which inferred "
                    "activity is never negative, from an ordinary growth trace",
             upgrade="a cycloheximide chase gives the value rather than the bound; "
                     "docs/PROTOCOLS.md carries it"),
        dict(parameter="gdcw_per_od", value=f"{PlateConditions(seed=0).gdcw_per_od}",
             unit="g L-1 per OD", grade=ASSERTED,
             source="working value; calibrate.py fits gain x gdcw_per_od as ONE product and "
                    "lists the pair as unidentifiable without a dry weight",
             upgrade="one gravimetric dry-weight measurement on a known OD"),
        dict(parameter="biomass yield Ysx", value="0.45-0.49", unit="g g-1",
             grade=BORROWED,
             source="van Hoek 1998 PMID 9797269, GLUCOSE-LIMITED chemostat",
             upgrade="the calibration strains were glucose-EXCESS at 4.1-6.4x that uptake; "
                     "measure dry weight on your own culture"),
        dict(parameter="glucose uptake q_glc", value="0.3-11.1", unit="mmol gDCW-1 h-1",
             grade=BORROWED, source="van Hoek 1998, same table, same caveat",
             upgrade="as above -- it sets the yield column and nothing else"),
        dict(parameter="molar mass", value="536.87", unit="g mol-1", grade=MEASURED,
             source="chemistry; beta-carotene C40H56",
             upgrade="nothing to buy: it is exact, and the only way it moves is if the molecule being reported is not the one named"),
        dict(parameter="inner-filter coefficient", value="none", unit="-", grade=REFUSED,
             source="ReporterOptics refuses the correction when it is uncalibrated rather "
                    "than treating it as zero",
             upgrade="spike purified carotenoid into a non-producing strain"),
        dict(parameter="autofluorescence a", value="~249", unit="RFU per OD",
             grade=BOUNDED,
             source="one recovered control plate; splices two runs eleven days apart and "
                    "the control has no medium-only fluorescing well, so it is an "
                    "over-estimate",
             upgrade="a per-construct value from the recovered channel -- a data-reduction "
                     "task, not a plate to run"),
        dict(parameter="ATP per unit activity", value="envelope", unit="mmol gDW-1",
             grade=ASSERTED,
             source="Tier 0; was found wrong by a factor of ten thousand and is now a range "
                    "rather than a number",
             upgrade="the latent branch it feeds is refused by G4 anyway"),
    ]
    return pd.DataFrame(rows)


#: Which constants each reported quantity actually leans on. Written down because the
#: alternative is a reader inferring it from the arithmetic, and the whole point of the
#: inventory is that nobody should have to.
DEPENDS_ON = {
    "rate": ["alpha (expression -> flux)", "leave-one-strain-out spread",
             "cyclase vmax/mu", "cyclase Km", "growth rate mu", "molar mass"],
    "content": ["alpha (expression -> flux)", "leave-one-strain-out spread",
                "cyclase vmax/mu", "cyclase Km", "growth rate mu", "molar mass",
                "content ceiling"],
    "yield": ["alpha (expression -> flux)", "leave-one-strain-out spread",
              "cyclase vmax/mu", "cyclase Km", "growth rate mu", "molar mass",
              "glucose uptake q_glc"],
    "titre": ["alpha (expression -> flux)", "leave-one-strain-out spread",
              "cyclase vmax/mu", "cyclase Km", "growth rate mu", "molar mass",
              "content ceiling", "biomass yield Ysx"],
}


def provenance_for(quantity: str) -> pd.DataFrame:
    """The graded constants one reported quantity rests on.

    This is the link the inventory was missing. `ProductPrediction` carries `layers`, which
    says which LAYERS ran; it says nothing about which of its CONSTANTS were measured and
    which were asserted, and a caller reading a titre had no way to find out except by
    reading the arithmetic. The van Hoek caveat reaches the notes only because somebody
    typed it at that call site -- so a new borrowed constant would arrive silently, which is
    the failure this whole table exists to prevent.

    Raises:
        KeyError: naming the quantities that are mapped, rather than returning an empty
            frame that would read as "this rests on nothing".
    """
    if quantity not in DEPENDS_ON:
        raise KeyError(
            f"no dependency list for {quantity!r}; mapped quantities are "
            f"{sorted(DEPENDS_ON)}. Add one rather than reading an empty frame as 'no "
            "dependencies'")
    table = inventory()
    return table[table.parameter.isin(DEPENDS_ON[quantity])]


def provenance_summary(quantity: str) -> str:
    """One line a caller can print beside a number: what it rests on, worst grade first."""
    rows = provenance_for(quantity)
    counts = rows.grade.value_counts()
    parts = [f"{int(counts[g])} {g.lower()}"
             for g in (REFUSED, ASSERTED, BORROWED, BOUNDED, MEASURED) if g in counts]
    weakest = next((g for g in (REFUSED, ASSERTED, BORROWED, BOUNDED, MEASURED)
                    if g in counts), MEASURED)
    return f"{quantity}: {', '.join(parts)} — weakest link {weakest}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grade", default=None,
                        help="show one grade only: measured, bounded, borrowed, asserted, refused")
    args = parser.parse_args()

    table = inventory()
    if args.grade:
        table = table[table.grade == args.grade.upper()]

    order = [MEASURED, BOUNDED, BORROWED, ASSERTED, REFUSED]
    counts = table.grade.value_counts()
    print("\nWhat the prediction rests on\n")
    for grade in order:
        rows = table[table.grade == grade]
        if rows.empty:
            continue
        print(f"  {grade}  ({len(rows)})")
        for r in rows.itertuples():
            print(f"    {r.parameter:34s} {r.value:>12s} {r.unit}")
            print(f"      from:  {r.source}")
            print(f"      fix:   {r.upgrade}")
        print()

    total = len(table)
    solid = int(counts.get(MEASURED, 0))
    print(f"  {solid} of {total} constants are measured on data this repository holds.")
    print("  The rest are bounded, borrowed from another culture, asserted, or refused --")
    print("  and every one names the measurement that would move it.")

    print("\n  what each reported quantity rests on:")
    for quantity in DEPENDS_ON:
        print(f"    {provenance_summary(quantity)}")

    out = paths.outputs_dir() / "parameter_provenance.csv"
    table.to_csv(out, index=False)
    print(f"\n-> {paths.display_path(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
