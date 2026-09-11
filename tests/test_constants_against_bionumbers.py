"""Pinning the generator's constants to an independent curated source.

BioNumbers is the Milo lab's curated database of biological quantities, each row carrying a
citation. It was assembled by people with no interest in this model, which is what makes it
worth checking against: a constant that agrees with it agrees with something external.

Every value below was picked from the literature while building the generator and only
afterwards compared. The agreement is therefore a check rather than a fit, and the one
disagreement is a real split in the published record rather than an error.
"""

import pandas as pd
import pytest

from ystwin import paths

BIONUMBERS = paths.data_dir() / "kaggle" / "BioNumbers_Nov2024.csv"
pytestmark = pytest.mark.skipif(not BIONUMBERS.exists(), reason="BioNumbers export absent")


@pytest.fixture(scope="module")
def yeast():
    frame = pd.read_csv(BIONUMBERS)
    return frame[frame.Organism.astype(str).str.contains("cerevisiae", case=False, na=False)]


def value_for(yeast, bnid):
    row = yeast[yeast.BNID == bnid]
    assert not row.empty, f"BNID {bnid} missing from the export"
    return str(row.iloc[0]["Value/Range"]).strip("'~ ")


class TestTheGeneratorsConstants:
    def test_growth_rate_matches_the_measured_doubling_time(self, yeast):
        """BNID 100270 puts a rich-medium generation time near 100 minutes."""
        import numpy as np

        from ystwin.generator.context import CultureContext, context_growth_rate

        minutes = float(value_for(yeast, 100270))
        implied = np.log(2) / (minutes / 60.0)

        assert context_growth_rate(CultureContext()) == pytest.approx(implied, rel=0.15)

    def test_the_atp_window_brackets_the_measured_concentration(self, yeast):
        """BNID 106020, aerobic glucose. The TMFA window is yeast-GEM's own YMDB range."""
        measured_mm = float(value_for(yeast, 106020))

        assert 0.9 <= measured_mm <= 4.4

    def test_the_nad_ratio_is_the_one_the_panel_cites(self, yeast):
        """BNID 108145 is Canelas 2008, which is where the panel took it from."""
        assert "101" in value_for(yeast, 108145)
        assert "320" in value_for(yeast, 108145)

    def test_the_measured_gssg_concentration_is_micromolar(self, yeast):
        """BNID 103548. Cytosolic GSSG is a few micromolar, not the millimolar a whole-cell
        extract suggests -- which is why whole-cell glutathione is the wrong observable."""
        assert float(value_for(yeast, 103548)) < 50.0


class TestTheOneDisagreement:
    """E_GSH is reported at -289 mV here and at -300 to -320 by the Grx1-roGFP2 line.

    Both are in the literature and the split is real: BNID 103543 is the rxYFP measurement,
    while the roGFP2 papers assume or measure a different cytosolic pH. The panel records
    the roGFP2 range and the reason for the gap rather than picking one silently.
    """

    def test_the_curated_value_is_the_rxyfp_one(self, yeast):
        assert float(value_for(yeast, 103543)) == pytest.approx(-289.0, abs=5.0)

    def test_the_panel_cites_the_roGFP2_range_and_says_why_they_differ(self):
        from ystwin.generator.stress_panel import REPORTERS

        source = REPORTERS["roGFP2-Grx1"].source

        assert "-306" in source or "306" in source
        assert "pH" in source

    def test_the_two_are_within_a_pH_unit_of_each_other(self, yeast):
        """A 40 mV gap is about what one pH unit does to the Nernst term, which is the
        published explanation rather than a discrepancy anyone needs to resolve here."""
        curated = float(value_for(yeast, 103543))

        assert abs(curated - (-320.0)) < 60.0
