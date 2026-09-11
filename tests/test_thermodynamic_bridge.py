"""Coupling a measured cofactor ratio to the metabolic model, the way the literature allows.

The obvious coupling -- a sensor reads a pool, the pool sets a flux -- is refuted by the two
best-characterised cases. Larsson 1997 (PMID 9393686) found glycolytic flux correlates
NEGATIVELY with intracellular ATP in yeast chemostats and not at all with the ATP/ADP ratio.
Vemuri 2007 (PMID 17287356) drained cytosolic NADH enough to abolish 80% of glycerol
production and the critical dilution rate did not move; Agrimi 2011 (PMID 21335394) found
perturbing the mitochondrial NAD pool shifts D_crit while leaving the critical glucose flux
invariant. Cytosolic redox is a consequence of the flux, not its control variable.

What has precedent is narrower and real. Kummel 2006 (PMID 16788595) put the adenylate
energy charge and the NADH/NAD ratio on iND750, and Martinez 2014 (PMID 25028891) did the
same on Yeast 5, as bounds on Gibbs energies -- constraining reaction DIRECTIONALITY, not
flux magnitude. Both set those bounds deliberately loose so they would not bind.

So the coupling built here is a directionality constraint, and the test that matters is
whether it narrows anything at all. A measurement that leaves the feasible flux space
untouched has told the model nothing, and this suite is willing to report that.
"""

import numpy as np
import pytest

from ystwin.bridge.thermodynamic import (
    GAS_CONSTANT_KJ,
    ThermodynamicData,
    dg_prime,
    reaction_dg0,
)

from ystwin import paths

pytest.importorskip("pytfa")

pytestmark = pytest.mark.skipif(
    paths.yeast_gem() is None or paths.thermo_dir() is None,
    reason="needs yeast-GEM and the thermodynamic tables; see docs/REPRODUCING.md",
)


@pytest.fixture(scope="module")
def thermo(thermo_data):
    return thermo_data


class TestTheThermodynamicData:
    def test_it_maps_yeast_metabolites_onto_formation_energies(self, thermo):
        assert thermo.dgf("s_0434") is not None

    def test_the_cofactors_the_sensors_read_are_all_covered(self, thermo):
        for met in ("s_0434", "s_0394", "s_1198", "s_1203", "s_0750", "s_0754"):
            assert thermo.dgf(met) is not None, met

    def test_an_unmapped_metabolite_returns_nothing_rather_than_guessing(self, thermo):
        assert thermo.dgf("s_not_a_real_metabolite") is None

    def test_it_reports_how_much_of_the_model_it_covers(self, thermo):
        assert 0.0 < thermo.coverage < 1.0

    def test_atp_hydrolysis_lands_where_the_textbook_puts_it(self, thermo):
        """The sanity check on the whole table, and the reason the pH transform is applied.

        Untransformed the tabulated energies give -6.8 kJ/mol for ATP hydrolysis against a
        textbook -30. Transformed they give about -40, which is what the textbook figure
        becomes without the magnesium it assumes -- so this range is right and the raw sum
        was not.
        """
        dg0 = reaction_dg0({"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1}, thermo)

        assert dg0 is not None
        assert -50.0 < dg0 < -25.0

    def test_the_ph_it_was_transformed_to_is_recorded(self, thermo):
        assert thermo.ph == pytest.approx(7.5)

    def test_a_different_ph_gives_different_energies(self, thermo_data):
        """Not a detail: the 40 mV disagreement between the two E_GSH literatures is
        largely one assuming pH 7.0 where the other measured 7.5."""
        seven = ThermodynamicData.load(ph=7.0)
        seven_five = ThermodynamicData.load(ph=7.5)

        assert seven.dgf("s_0434") != pytest.approx(seven_five.dgf("s_0434"))


class TestGibbsEnergyUnderMeasuredConcentrations:
    def test_it_reduces_to_the_standard_energy_at_unit_activity(self):
        assert dg_prime(-30.0, {}, {}) == pytest.approx(-30.0)

    def test_a_product_accumulating_makes_a_reaction_less_favourable(self):
        loaded = dg_prime(-30.0, {"P": 1}, {"P": 10.0})
        empty = dg_prime(-30.0, {"P": 1}, {"P": 0.1})

        assert loaded > empty

    def test_a_substrate_accumulating_makes_it_more_favourable(self):
        rich = dg_prime(-30.0, {"S": -1}, {"S": 10.0})
        poor = dg_prime(-30.0, {"S": -1}, {"S": 0.1})

        assert rich < poor

    def test_the_shift_is_rt_ln_q(self):
        shift = dg_prime(0.0, {"P": 1, "S": -1}, {"P": 10.0, "S": 1.0})

        assert shift == pytest.approx(GAS_CONSTANT_KJ * 303.15 * np.log(10.0), rel=1e-6)

    def test_temperature_enters(self):
        cold = dg_prime(0.0, {"P": 1}, {"P": 10.0}, temperature_k=273.15)
        warm = dg_prime(0.0, {"P": 1}, {"P": 10.0}, temperature_k=310.15)

        assert warm > cold

    def test_a_missing_concentration_is_refused_rather_than_assumed(self):
        with pytest.raises(KeyError, match="no concentration"):
            dg_prime(-30.0, {"S": -1}, {})


class TestWhatTheMeasuredRatioDoes:
    def test_a_measured_atp_ratio_shifts_the_energy_of_a_coupled_reaction(self, thermo):
        """The whole mechanism: the ratio enters Q, and Q moves dG'."""
        stoich = {"s_0434": -1, "s_0394": 1}
        high = dg_prime(0.0, stoich, {"s_0434": 4.4e-3, "s_0394": 0.3e-3})
        low = dg_prime(0.0, stoich, {"s_0434": 0.9e-3, "s_0394": 2.2e-3})

        assert low > high

    def test_the_measured_yeast_range_is_wide_enough_to_matter(self, thermo):
        """yeast-GEM ships measured ATP 0.9-4.4 mM and ADP 0.3-2.2 mM, then discards them."""
        stoich = {"s_0434": -1, "s_0394": 1}
        span = (dg_prime(0.0, stoich, {"s_0434": 0.9e-3, "s_0394": 2.2e-3})
                - dg_prime(0.0, stoich, {"s_0434": 4.4e-3, "s_0394": 0.3e-3}))

        assert span > 5.0
