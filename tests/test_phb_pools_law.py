"""`X = v_in/mu` on PHB: holds on ethanol, fails on glucose.

The terminal DILUTED line says content falls as growth rises, because washout is the only
outlet. Kocharin 2013's eleven chemostat states are already vendored, so the law is
falsifiable by one multiplication -- and it splits by carbon source:

    glucose   mu 0.05 -> 0.20   content 4.33 -> 5.59   predicted ratio 4.00, measured 0.775
    ethanol   mu 0.05 -> 0.15   content 16.55 -> 5.30  predicted 3.00, measured 3.12
    mixed     mu 0.05 -> 0.20   content 18.34 -> 4.56  predicted 4.00, measured 4.02

On glucose the pool goes UP with growth. The arithmetic is not wrong: `X = v_in/mu` is an
identity given the balance. What fails is the assumption underneath every prediction this
package makes -- that `v_in` does not depend on `mu`. On glucose `v_in = X*mu` rises 5x
across the range; on ethanol and mixed feed it is flat.
"""

from __future__ import annotations

import csv

import pytest

from ystwin import paths

_TABLE = paths.data_dir() / "phb" / "kocharin2013_chemostat_states.tsv"


def _states():
    """(carbon source, mu, content mg/gDW), from the vendored table."""
    with _TABLE.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    out: dict[str, list[tuple[float, float]]] = {}
    for row in rows:
        out.setdefault(row["carbon_source"], []).append(
            (float(row["mu_per_h"]), float(row["phb_mg_per_gdw"])))
    return {k: sorted(v) for k, v in out.items()}


@pytest.mark.skipif(not _TABLE.is_file(), reason="phb chemostat table not vendored")
class TestThePoolsLawSplitsByCarbonSource:
    @pytest.mark.parametrize("source", ["ethanol", "glucose_ethanol_1to2"])
    def test_content_falls_with_growth_as_the_law_says(self, source):
        """Where `v_in` is flat, `X = v_in/mu` predicts the ratio to within 4%."""
        states = _states()[source]
        (low_mu, low_x), (high_mu, high_x) = states[0], states[-1]

        assert low_x / high_x == pytest.approx(high_mu / low_mu, rel=0.05)

    def test_on_glucose_the_content_rises_instead(self):
        """The refutation. The law predicts 4.0x; the measurement is 0.775x, the wrong way."""
        states = _states()["glucose"]
        (low_mu, low_x), (high_mu, high_x) = states[0], states[-1]

        assert low_x < high_x
        assert low_x / high_x == pytest.approx(0.775, abs=0.01)
        assert high_mu / low_mu == pytest.approx(4.0)

    def test_the_broken_assumption_is_that_the_entry_flux_is_constant(self):
        """`v_in = X*mu` is flat on ethanol and rises 5x on glucose. That is the whole story:
        the balance is an identity, and what varies is what feeds it."""
        flux = {source: [x * mu for mu, x in states]
                for source, states in _states().items()}

        assert max(flux["glucose"]) / min(flux["glucose"]) > 5.0
        assert max(flux["ethanol"]) / min(flux["ethanol"]) < 2.0

    def test_the_table_still_holds_the_states_this_rests_on(self):
        """An exemption-free check that the numbers above came from somewhere."""
        states = _states()

        assert len(states["glucose"]) == 4
        assert states["glucose"][0] == (0.05, 4.33)
        assert states["ethanol"][0] == (0.05, 16.55)
