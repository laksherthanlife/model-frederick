"""Enzyme abundance from the measured proteome, and the boundary it must not cross.

Two things are pinned here and the second matters more. The conversion from PaxDb's molar
ppm to mmol/gDCW has to stay molar end to end, because mixing molar and mass units is an
error of order the enzyme's mass ratio and comes out looking like an ordinary number. And
every HETEROLOGOUS enzyme must come back `None`: PaxDb measured a wild-type proteome, so a
cassette gene is absent by construction, and a lookup that quietly returned something for
`CrtE` would put a wild-type number where a producing strain's should be.
"""

from __future__ import annotations

import csv

import pytest

from ystwin import paths
from ystwin.pathway.proteome import (
    AVERAGE_PROTEIN_MW_G_PER_MOL,
    PROMOTER_ANCHOR_PPM,
    TOTAL_PROTEIN_G_PER_GDCW,
    abundance_ppm,
    enzyme_content_from_mass_fraction,
    enzyme_content_mmol_per_gdcw,
    proteome_fraction,
)

CRTYB_MW = 74736.0

_TABLE = paths.data_dir() / "proteome" / "paxdb_scerevisiae_integrated.tsv"

pytestmark = pytest.mark.skipif(not _TABLE.is_file(), reason="PaxDb table not vendored")


class TestTheVendoredTableIsWhatItClaims:
    def test_the_abundances_sum_to_a_million(self):
        """The check that says the column is a MOLAR FRACTION and not an intensity in
        arbitrary units. If a refresh breaks this, the column has changed meaning and the
        conversion downstream is invalid -- see data/proteome/SOURCE.md."""
        with _TABLE.open() as handle:
            total = sum(float(r["abundance_ppm"]) for r in csv.DictReader(handle, delimiter="\t"))

        assert total == pytest.approx(1e6, rel=1e-4)

    def test_it_covers_most_of_the_proteome(self):
        with _TABLE.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

        assert len(rows) > 5000

    def test_every_abundance_is_positive(self):
        with _TABLE.open() as handle:
            assert all(float(r["abundance_ppm"]) > 0
                       for r in csv.DictReader(handle, delimiter="\t"))


class TestNativeEnzymesResolve:
    @pytest.mark.parametrize("gene", ["ERG9", "ERG1", "ERG20", "BTS1", "TDH3"])
    def test_a_native_enzyme_has_an_abundance(self, gene):
        assert abundance_ppm(gene) > 0.0

    def test_the_fraction_is_the_ppm_over_a_million(self):
        assert proteome_fraction("ERG9") == pytest.approx(abundance_ppm("ERG9") / 1e6)

    def test_the_glycolytic_workhorse_is_the_most_abundant_thing_here(self):
        """A sanity check on the whole table rather than on one lookup: if TDH3 is not at
        the top, the file being read is not a yeast proteome."""
        with _TABLE.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        top = max(rows, key=lambda r: float(r["abundance_ppm"]))

        assert top["gene_name"] == "TDH3"


class TestHeterologousEnzymesDoNotResolve:
    """The boundary. PaxDb measured wild-type cells; a cassette gene was not in them."""

    @pytest.mark.parametrize("gene", ["CrtE", "CrtYB", "CrtI", "EEVS", "MT-Ox"])
    def test_a_cassette_gene_is_absent_rather_than_defaulted(self, gene):
        assert abundance_ppm(gene) is None
        assert proteome_fraction(gene) is None


class TestTheConversion:
    def test_it_is_molar_end_to_end(self):
        """ppm counts MOLECULES, so it divides into total protein MOLES. Dividing into grams
        instead is wrong by the enzyme's mass ratio and looks like an ordinary answer."""
        ppm = 1000.0
        expected = ppm / 1e6 * (TOTAL_PROTEIN_G_PER_GDCW / AVERAGE_PROTEIN_MW_G_PER_MOL) * 1000.0

        assert enzyme_content_mmol_per_gdcw(ppm) == pytest.approx(expected, rel=1e-12)

    def test_it_is_linear_in_the_abundance(self):
        assert (enzyme_content_mmol_per_gdcw(200.0)
                == pytest.approx(2 * enzyme_content_mmol_per_gdcw(100.0)))

    def test_the_assumed_average_mass_scales_the_answer_linearly(self):
        """Stated as a test because the assumption is invisible otherwise: halving the
        assumed average protein mass doubles every capacity this module feeds."""
        half = enzyme_content_mmol_per_gdcw(100.0, average_protein_mw_g_per_mol=25_000.0)

        assert half == pytest.approx(2 * enzyme_content_mmol_per_gdcw(100.0))

    def test_a_negative_ppm_is_refused(self):
        with pytest.raises(ValueError, match="non-negative"):
            enzyme_content_mmol_per_gdcw(-1.0)


class TestThePromoterAnchorsAreCeilings:
    def test_each_anchor_matches_the_vendored_table(self):
        """They are read from the measurement, not typed in. If the table is refreshed and
        an anchor stops matching, this fails rather than the constant silently aging."""
        for gene, ppm in PROMOTER_ANCHOR_PPM.items():
            assert abundance_ppm(gene) == pytest.approx(ppm)

    def test_an_anchor_is_far_above_a_typical_metabolic_enzyme(self):
        """Which is why using one as an abundance gives a CEILING: a cassette does not reach
        its promoter's own native product."""
        assert PROMOTER_ANCHOR_PPM["TDH3"] > 100 * abundance_ppm("ERG9")


class TestTheMassFractionConversionIsADifferentFunction:
    """The other half of the unit map, named for its own convention.

    A targeted PRM assay -- the measurement `pathway/enzyme_capacity.py` names as its
    promotion criterion -- reports a share of protein MASS, so the enzyme's own molar mass
    is the divisor and the host's assumed average never enters. Until 2026-09-04 this
    conversion lived inside `derived_capacity` behind an argument called `proteome_fraction`,
    which said nothing about which convention it carried.
    """

    def test_it_divides_by_the_enzymes_own_mass_and_not_the_average(self):
        assert enzyme_content_from_mass_fraction(0.005, 100_000.0) == pytest.approx(
            TOTAL_PROTEIN_G_PER_GDCW * 0.005 / 100_000.0 * 1000.0, rel=1e-12)

    def test_the_two_conversions_differ_by_exactly_the_mass_ratio(self):
        """The size of the error the old seam produced, stated as an identity rather than
        as a remembered number: feeding a MOLAR fraction down the MASS path is wrong by
        `AVERAGE_PROTEIN_MW / MW_enzyme`."""
        ppm = 5000.0
        molar = enzyme_content_mmol_per_gdcw(ppm)
        mass = enzyme_content_from_mass_fraction(ppm / 1e6, CRTYB_MW)

        assert mass / molar == pytest.approx(
            AVERAGE_PROTEIN_MW_G_PER_MOL / CRTYB_MW, rel=1e-12)

    def test_a_percentage_mistaken_for_a_fraction_is_refused(self):
        with pytest.raises(ValueError, match="FRACTION"):
            enzyme_content_from_mass_fraction(5.0, CRTYB_MW)

    def test_a_negative_fraction_is_refused(self):
        with pytest.raises(ValueError, match="non-negative"):
            enzyme_content_from_mass_fraction(-0.001, CRTYB_MW)


class TestItFeedsTheCapacityLayerWithNothingFitted:
    """The seam. `derived_capacity` takes a CONTENT, and this module is what produces one.

    This class replaces `test_a_native_capacity_needs_no_per_product_number`, which asked
    only that the answer was positive and that a fraction survived a round trip -- and it
    asked it about ERG9, at 51000 g/mol against an assumed average of 50000. That is a
    0.9804x ratio, so the test passed under EITHER unit convention and could not see that
    the two modules disagreed. The cases below are chosen for the opposite property: a molar
    mass at least 1.4x away from the average, where the two conventions cannot be mistaken
    for one another.
    """

    #: (label, molar ppm, molar mass g/mol). CrtYB is heterologous, so its abundance is the
    #: TDH3 promoter anchor and the number is a CEILING, which is all this test needs. FAS1's
    #: mass is UniProt P07149 (2051 aa).
    CASES = [("crtYB at the TDH3 anchor", PROMOTER_ANCHOR_PPM["TDH3"], 74736.0),
             ("FAS1", 1284.0, 228_689.0)]

    @pytest.mark.parametrize("label, ppm, molar_mass", CASES)
    def test_the_mass_ratio_is_far_enough_from_one_to_tell_the_conventions_apart(
            self, label, ppm, molar_mass):
        """The property that makes the rest of this class a real cross-check. If a future
        edit swaps these for enzymes near the average mass, this fails and says why."""
        ratio = AVERAGE_PROTEIN_MW_G_PER_MOL / molar_mass

        assert ratio > 1.4 or ratio < 1 / 1.4, (
            f"{label} is only {ratio:.4g}x from the average protein mass, so it cannot "
            f"distinguish the molar conversion from the mass one")

    @pytest.mark.parametrize("label, ppm, molar_mass", CASES)
    def test_the_content_reaching_the_capacity_layer_is_the_molar_conversion(
            self, label, ppm, molar_mass):
        from ystwin.pathway.enzyme_capacity import derived_capacity

        content = enzyme_content_mmol_per_gdcw(ppm)
        derived = derived_capacity(label, kcat_per_s=1.65, enzyme_mmol_per_gdcw=content)

        assert derived.enzyme_mmol_per_gdcw == content
        assert derived.vmax_mmol_per_gdcw_h == pytest.approx(1.65 * 3600.0 * content)

    @pytest.mark.parametrize("label, ppm, molar_mass", CASES)
    def test_and_the_mass_path_would_have_given_a_different_number(
            self, label, ppm, molar_mass):
        """What the old seam did. This is the assertion the ERG9 test could not make."""
        molar = enzyme_content_mmol_per_gdcw(ppm)
        mass = enzyme_content_from_mass_fraction(ppm / 1e6, molar_mass)

        assert mass != pytest.approx(molar, rel=0.05)

    def test_a_native_abundance_still_needs_no_per_product_number(self):
        from ystwin.pathway.enzyme_capacity import derived_capacity

        derived = derived_capacity(
            "ERG9", kcat_per_s=1.65,
            enzyme_mmol_per_gdcw=enzyme_content_mmol_per_gdcw(abundance_ppm("ERG9")))

        assert derived.vmax_mmol_per_gdcw_h > 0.0

    def test_there_is_exactly_one_total_protein_constant(self):
        """`enzyme_capacity` carried its own copy of 0.49974 until 2026-09-04. Two literals
        that must agree and nothing that makes them is the shape of the bug this class
        exists to prevent, so the check is that the second copy is GONE rather than that the
        two match."""
        from ystwin.pathway import enzyme_capacity

        assert not hasattr(enzyme_capacity, "TOTAL_PROTEIN_G_PER_GDCW")
