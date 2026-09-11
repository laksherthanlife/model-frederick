"""Enzyme abundance from the measured proteome, so a capacity needs no per-product fit.

`pathway/enzyme_capacity.py` computes ``vmax = kcat * [E]`` and refuses without ``[E]``,
because the enzyme's share of the proteome moves the answer by the width of its own window
and a default would BE the result. This module supplies that term for NATIVE enzymes from a
vendored measurement -- PaxDb's integrated *S. cerevisiae* whole-organism dataset, 6,351
proteins, abundance in molar ppm.

    [E] mmol/gDCW  =  ppm/1e6  x  (total protein g/gDCW / average protein MW)  x  1000

THIS MODULE OWNS THE UNIT MAP, both halves of it, and each half is named for its own
convention. :func:`enzyme_content_mmol_per_gdcw` takes a MOLAR ppm, which is what PaxDb
reports. :func:`enzyme_content_from_mass_fraction` takes a MASS fraction, which is what a
targeted PRM assay reports, and divides by the enzyme's own molar mass instead. They differ
by ``AVERAGE_PROTEIN_MW_G_PER_MOL / MW_enzyme``. Until 2026-09-04 the mass map lived in
`pathway/enzyme_capacity.py` behind an argument called `proteome_fraction` -- a name that
says nothing about which convention it carries -- and this module's own
:func:`proteome_fraction` advertised the molar number as being in that form. One conversion
per convention, in one module, named: see that module's docstring for the size of the error
that arrangement produced.

WHAT THIS BUYS AND WHAT IT COSTS. It removes the last per-product measurement from the
capacity: a kcat is a property of a protein and an abundance is a property of the host, so
neither is fitted to the product being predicted. What it costs is accuracy, and the cost
was measured rather than guessed -- see the table in `data/proteome/SOURCE.md`. Against
native enzymes whose flux, kcat and abundance are each known independently, the derived
capacity lands within 6% for one enzyme and falls short by up to twenty-fold for another.
**The error is about an order of magnitude and it is ONE-SIDED, biased low.**

THAT BIAS IS A FINDING, NOT NOISE. A capacity below the flux the network carries is
infeasible -- an enzyme cannot pass more than it can catalyse -- so where the ratio comes out
under one, the in vitro kcat understates the in vivo rate. The enzyme it fails worst on is
squalene epoxidase, the classical rate-limiting step of sterol synthesis, which is where an
in vitro number is most likely to be an underestimate.

WHAT IT CANNOT DO, and this is the boundary a caller must not cross by accident. PaxDb
measures a WILD-TYPE proteome, so **no heterologous protein is in it**. The Crt enzymes, EEVS
and every other cassette gene are absent by construction, and :func:`abundance_ppm` returns
``None`` for them rather than a number. The promoter anchors in
:data:`PROMOTER_ANCHOR_PPM` exist for that case and are CEILINGS: a cassette expressed from
a promoter cannot exceed that promoter's own native product by much, and will usually sit
well below it, because a heterologous protein folds and is degraded differently.
"""

from __future__ import annotations

import csv
import functools

from .. import paths

__all__ = [
    "AVERAGE_PROTEIN_MW_G_PER_MOL",
    "PROMOTER_ANCHOR_PPM",
    "TOTAL_PROTEIN_G_PER_GDCW",
    "abundance_ppm",
    "enzyme_content_from_mass_fraction",
    "enzyme_content_mmol_per_gdcw",
    "proteome_fraction",
]

#: Total protein per gram dry cell weight, g/gDCW. A host physiological constant.
TOTAL_PROTEIN_G_PER_GDCW = 0.49974

#: Average protein molar mass used to convert a MOLAR ppm into a molar content, g/mol.
#: **This is an assumption and it scales every derived content linearly.** It is declared
#: here rather than buried in the arithmetic so that a reader can see what a factor-of-two
#: error in it would do: a factor of two in every capacity this module feeds.
AVERAGE_PROTEIN_MW_G_PER_MOL = 50_000.0

#: Native abundance of the proteins whose promoters heterologous cassettes usually borrow,
#: molar ppm, read from the vendored PaxDb table. **Upper bounds, not estimates.** A cassette
#: on the TDH3 promoter does not reach TDH3's own abundance -- transcription is only one of
#: the steps between a promoter and a folded protein -- so using one of these as ``[E]``
#: yields a capacity CEILING and must be labelled as such by the caller.
PROMOTER_ANCHOR_PPM = {"TDH3": 23256.0, "PGK1": 13767.0, "TEF1": 8658.0}

_TABLE = paths.data_dir() / "proteome" / "paxdb_scerevisiae_integrated.tsv"


@functools.lru_cache(maxsize=1)
def _abundances() -> dict[str, float]:
    """Gene name to molar ppm, read once.

    Raises:
        FileNotFoundError: if the vendored table is absent, rather than returning an empty
            map -- a silent empty map would make every lookup return ``None`` and read as
            "this protein was not measured" instead of "the measurement is missing".
    """
    if not _TABLE.is_file():
        raise FileNotFoundError(
            f"{_TABLE} is not vendored. It is the PaxDb S. cerevisiae integrated "
            f"whole-organism table; data/proteome/SOURCE.md records where it came from and "
            f"the one check that says the column still means molar ppm")
    with _TABLE.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {r["gene_name"]: float(r["abundance_ppm"]) for r in rows if r["abundance_ppm"]}


def abundance_ppm(gene: str) -> float | None:
    """Molar ppm for a gene, or ``None`` where it is not in the measured proteome.

    ``None`` is the answer for every heterologous enzyme, and it is the right answer rather
    than an inconvenience: the cassette was not present in the cells PaxDb measured.
    """
    return _abundances().get(gene)


def proteome_fraction(gene: str) -> float | None:
    """The gene's share of protein MOLECULES, 0-1. A MOLAR fraction, not a mass fraction.

    PaxDb's ppm counts molecules, so this is ``ppm/1e6`` and nothing else. Its consumer is
    :func:`enzyme_content_mmol_per_gdcw` -- pass the ppm, or multiply this by 1e6.

    Until 2026-09-04 this docstring said the number was "in the form `enzyme_capacity`
    takes", which was false and pointed callers down that module's MASS-fraction path. The
    error is exactly ``AVERAGE_PROTEIN_MW_G_PER_MOL / MW_enzyme`` -- 0.6690x for CrtYB. That
    path no longer exists: :func:`ystwin.pathway.enzyme_capacity.derived_capacity` now takes
    a content in mmol/gDCW and this module owns both conversions into it.
    """
    ppm = abundance_ppm(gene)
    return None if ppm is None else ppm / 1e6


def enzyme_content_mmol_per_gdcw(
        ppm: float,
        total_protein_g_per_gdcw: float = TOTAL_PROTEIN_G_PER_GDCW,
        average_protein_mw_g_per_mol: float = AVERAGE_PROTEIN_MW_G_PER_MOL) -> float:
    """Convert a molar ppm to an enzyme content.

    The conversion is molar throughout: PaxDb's ppm counts MOLECULES, so it is divided into
    the host's total protein MOLES rather than its grams. Mixing the two is a silent error
    of about the ratio of the enzyme's mass to the average, which for a large enzyme is a
    factor of two and looks like an ordinary answer.
    """
    if ppm < 0:
        raise ValueError(f"ppm must be non-negative, got {ppm}")
    total_moles = total_protein_g_per_gdcw / average_protein_mw_g_per_mol
    return ppm / 1e6 * total_moles * 1000.0


def enzyme_content_from_mass_fraction(
        mass_fraction: float,
        molar_mass_g_per_mol: float,
        total_protein_g_per_gdcw: float = TOTAL_PROTEIN_G_PER_GDCW) -> float:
    """Convert a MASS fraction of total protein to an enzyme content, mmol/gDCW.

    The other conversion in this module. This one is for the measurement the enzyme-capacity
    promotion criterion actually asks for -- targeted PRM/SRM with a heavy peptide, which
    reports the analyte as a share of protein MASS -- so the enzyme's own molar mass is the
    divisor and the host's average protein mass never enters.

    The two functions differ by ``AVERAGE_PROTEIN_MW_G_PER_MOL / molar_mass_g_per_mol``,
    which for a 75 kDa enzyme is 0.67x and for a 25 kDa one is 2x. That is why the
    convention is in the function NAME: from 2026-09-02 to 2026-09-04 it lived in a
    docstring instead, and a molar fraction was routed down the mass path.

    Raises:
        ValueError: on a negative fraction, a non-positive mass, or a fraction above 1.0 --
            the last because ``0.5`` meaning half the proteome is almost always ``0.5%``
            mistyped, and the error is a hundredfold in the dominant term.
    """
    if mass_fraction < 0:
        raise ValueError(f"mass_fraction must be non-negative, got {mass_fraction}")
    if mass_fraction > 1.0:
        raise ValueError(
            f"mass_fraction is a FRACTION of total protein, not a percentage, and must be "
            f"<= 1.0; got {mass_fraction}. A cassette at one percent of protein is 0.01 here")
    if molar_mass_g_per_mol <= 0 or total_protein_g_per_gdcw <= 0:
        raise ValueError(
            f"molar_mass_g_per_mol and total_protein_g_per_gdcw must be positive, got "
            f"{molar_mass_g_per_mol} and {total_protein_g_per_gdcw}")
    return total_protein_g_per_gdcw * mass_fraction / molar_mass_g_per_mol * 1000.0
