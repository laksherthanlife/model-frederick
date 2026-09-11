"""The first mechanism that explains why carbon source moves PHB flux, and why nothing else could.

Five mechanisms failed before this one, and every failure had the same shape: FBA on the
product reaction, GEM precursor ceilings, GEM precursor THROUGHPUT, E-Flux regulation, and an
mRNA-reading law are all arguments about how much the network CAN carry. A rate cap is blind
to a driving force, and a driving force is what changes here.

yeast-GEM's own cytosolic thiolase `r_0103` is the chemical analogue of PhaA -- both condense
two acetyl-CoA into acetoacetyl-CoA plus CoA. It is uphill, so it has a concentration
threshold below which no net forward flux occurs. That much survives.

**THE LEAD IS REFUTED, AND THIS FILE IS NOW THE RECORD OF HOW.** The claim was that the two
carbon sources sit on OPPOSITE sides of that threshold -- glucose-limited about 10 uM,
ethanol-limited about 425 uM, against a threshold of 221 uM -- which would make the driving
force the thing carbon source changes. It rested on dGr'0 = +15.6 kJ/mol, and that number
was an artefact of a unit error: `bridge/thermodynamic.py` read a kcal/mol table through
`pytfa` under its default `thermo_unit='kJ/mol'`, so the pH transform landed in kJ on top of
a tabulated energy in kcal. Corrected on 2026-08-30, the same table gives +38.07.

Both feeds are below the threshold on every defensible reading of the energy:

    dGr'0                                    threshold   glucose 10 uM   ethanol 425 uM
    +15.61  ModelSEED, unit bug (published)    221 uM        below           ABOVE
    +38.07  ModelSEED, unit corrected       19,042 uM        below           below
    +24.96  eQuilibrator component contrib.  1,414 uM        below           below

The third row is the one that settles it, because eQuilibrator never touched the unit bug
and it is the estimator this repository already documents as the more trustworthy of the two
(`tests/test_equilibrator_backend.py`). Its +24.96 is the literature value for a thiolase
condensation (~+26). Under the only reading that matches literature, ethanol does not clear
the threshold either -- so the ordering the lead was built on exists only in the erroneous
number, and this was checkable from inside the repository the whole time.

**What these tests pin.** The arithmetic, the refutation, and that the refutation does not
depend on which corrected estimator is chosen.
"""

from __future__ import annotations

import math

import pytest

from ystwin import paths

pytestmark = pytest.mark.integration

THIOLASE = "r_0103"
GAS_CONSTANT_KJ = 8.314e-3
TEMPERATURE_K = 303.15

GLUCOSE_UM = 10.0        # 0.0199 umol/gDW at 2 mL/gDW cytosolic volume
ETHANOL_UM = 425.0       # 0.85   umol/gDW at 2 mL/gDW


def _threshold_um(dg0_kj, coa_m=100e-6, acetoacetyl_m=1e-6):
    equilibrium = math.exp(-dg0_kj / (GAS_CONSTANT_KJ * TEMPERATURE_K))
    return math.sqrt(acetoacetyl_m * coa_m / equilibrium) * 1e6


@pytest.fixture(scope="module")
def energy():
    """dGr'0 for the thiolase, from this repository's own vendored tables."""
    import warnings

    import cobra

    from ystwin.bridge.thermodynamic import ThermodynamicData, reaction_dg0

    gem = paths.yeast_gem()
    if gem is None or not gem.is_file():
        pytest.skip("Yeast9 not present; set YSTWIN_YEAST_GEM")
    warnings.simplefilter("ignore")
    model = cobra.io.read_sbml_model(str(gem))
    reaction = model.reactions.get_by_id(THIOLASE)
    value = reaction_dg0({m.id: c for m, c in reaction.metabolites.items()},
                         ThermodynamicData.load())
    if value is None:
        pytest.skip("thermodynamic tables do not cover every participant")
    return float(value)


class TestTheStepIsUphill:
    def test_the_standard_energy_is_positive(self, energy):
        """Which is what gives it a threshold at all. A downhill step runs at any
        concentration and could not produce a carbon-source switch."""
        assert energy > 0

    def test_and_it_is_about_thirty_eight_kilojoules(self, energy):
        """+15.6 until 2026-08-30, which was the unit error, not the chemistry."""
        assert energy == pytest.approx(38.07, abs=0.5)

    def test_the_reaction_is_the_one_claimed(self):
        """Two acetyl-CoA in, acetoacetyl-CoA and CoA out -- PhaA's chemistry. If yeast-GEM
        ever renumbers this, the energy above would be some other reaction's."""
        import warnings

        import cobra

        gem = paths.yeast_gem()
        if gem is None or not gem.is_file():
            pytest.skip("Yeast9 not present")
        warnings.simplefilter("ignore")
        reaction = cobra.io.read_sbml_model(str(gem)).reactions.get_by_id(THIOLASE)
        stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}

        assert stoichiometry == {"s_0373": -2.0, "s_0367": 1.0, "s_0529": 1.0}


class TestBothFeedsSitBelowTheThreshold:
    """The refutation. The lead needed ethanol ABOVE, and it is not."""

    def test_glucose_sits_below_it(self, energy):
        assert GLUCOSE_UM < _threshold_um(energy)

    def test_ethanol_also_sits_below_it(self, energy):
        """This test asserted the opposite until 2026-08-30. The inversion is the record."""
        assert ETHANOL_UM < _threshold_um(energy)

    def test_the_threshold_is_where_it_is_claimed(self, energy):
        assert _threshold_um(energy) == pytest.approx(19042.0, rel=0.05)

    @pytest.mark.parametrize("volume_ml_per_gdw", [1.0, 2.0, 2.7])
    def test_no_plausible_cytosolic_volume_puts_ethanol_above_it(self, energy,
                                                                 volume_ml_per_gdw):
        """The umol/gDW to molar conversion is a convention, not a measurement. It was not
        what carried the old ordering and it cannot rescue it: the largest pool any of these
        volumes gives is 850 uM against a threshold of 19 mM."""
        threshold = _threshold_um(energy)
        glucose = 0.0199e-6 / (volume_ml_per_gdw * 1e-3) * 1e6
        ethanol = 0.85e-6 / (volume_ml_per_gdw * 1e-3) * 1e6

        assert glucose < ethanol < threshold

    def test_no_corner_of_the_cofactor_grid_restores_the_ordering(self, energy):
        """The threshold scales as the square root of [CoA]*[acetoacetyl-CoA], so the
        assumed levels move it -- and they used to move it enough to matter, holding the
        ordering in eight of nine corners. At the corrected energy none of the nine puts
        ethanol above, so the refutation does not rest on the cofactor assumptions either."""
        held = 0
        for coa in (30e-6, 100e-6, 300e-6):
            for acetoacetyl in (0.3e-6, 1e-6, 3e-6):
                threshold = _threshold_um(energy, coa, acetoacetyl)
                if GLUCOSE_UM < threshold < ETHANOL_UM:
                    held += 1

        assert held == 0

    def test_the_independent_estimator_refutes_it_too(self):
        """The load-bearing test. eQuilibrator never touched the unit bug, and its energy
        is the literature value for a thiolase condensation. If the refutation held only
        under the corrected ModelSEED table it would be one estimator's opinion."""
        pytest.importorskip("equilibrator_api")
        import cobra

        from ystwin.bridge.equilibrator import EquilibratorData

        gem = paths.yeast_gem()
        if gem is None or not gem.is_file():
            pytest.skip("Yeast9 not present; set YSTWIN_YEAST_GEM")
        model = cobra.io.read_sbml_model(str(gem))
        reaction = model.reactions.get_by_id(THIOLASE)
        value = EquilibratorData.load(model=model).reaction_dg0(
            {m.id: c for m, c in reaction.metabolites.items()})

        assert value is not None
        assert value == pytest.approx(24.96, abs=1.0)
        assert ETHANOL_UM < _threshold_um(float(value))


class TestWhyEveryEarlierMechanismMissedIt:
    def test_the_threshold_is_a_concentration_and_fva_returns_rates(self):
        """The one-line reason, kept as a test so it is not lost as prose. Everything tried
        before this bounded a RATE. Nothing about a rate bound can express that a reaction
        stops running when a pool drops, because that is a statement about mass action."""
        from ystwin.fba.fva import FluxRange

        assert "minimum" in FluxRange.__dataclass_fields__
        assert not any("concentration" in f for f in FluxRange.__dataclass_fields__)

    def test_and_the_thermodynamics_layer_that_can_express_it_is_now_called(self):
        """This test used to assert the opposite, and the inversion is the record.

        It read ``assert "thermodynamic" not in chain`` and said of itself: "780 lines,
        tested, and reached from nothing in the prediction chain. This test will start
        failing the day that changes, which is the intended signal." It fired. Wiring
        `pathway/thermo_gate.py` into `predict.py` is what made it fail, so the assertion
        is turned around rather than deleted -- a tripwire that is removed once it goes off
        leaves nothing behind saying the thing it guarded ever happened.

        What it guards now is the other direction: the gate reaching the chain is the only
        way the threshold above can be expressed on a real solve, and a refactor that
        quietly unhooks it would put this repository back where it was.
        """
        import pathlib

        chain = (pathlib.Path(__file__).resolve().parents[1] / "src" / "ystwin"
                 / "predict.py").read_text()

        assert "from .pathway.thermo_gate import" in chain
        assert "gate_across_volumes(" in chain

    def test_and_the_step_it_gates_is_the_one_this_module_is_about(self):
        """The wiring is only worth anything if it gates r_0103 specifically. `phb.toml`
        declares that reaction on its entry node, so the chain and this module are talking
        about the same chemistry rather than two things that happen to share a name."""
        from ystwin.pathway.spec import load_pathway

        entry = load_pathway("phb").node("acetoacetyl_coa")

        assert entry.reaction == THIOLASE
        assert entry.metabolite == "s_0367"
