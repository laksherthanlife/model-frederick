"""The two terms the Gibbs energies were missing, and why each one matters.

Both were documented as limits rather than fixed, and both change answers.

**Magnesium.** ATP, ADP and phosphate all bind Mg2+, and at the roughly 1 mM free Mg2+ of a
yeast cytosol most ATP is present as MgATP. Binding stabilises the more highly charged
species preferentially, so it shifts hydrolysis by several kJ/mol. The textbook -30.5 for
ATP hydrolysis is quoted at 1 mM Mg; without it the same reaction is nearer -37, which is
what an unmodelled table returns and what looked like an error.

**Membrane potential.** Moving a charged species across a membrane costs zF*dPsi whatever the
chemistry does. The yeast mitochondrial inner membrane sits near -160 mV, so a proton
entering the matrix releases about 15 kJ/mol of electrical work that a chemistry-only
energy cannot see. Without it proton pumping looks thermodynamically impossible, which is
exactly what oxidative phosphorylation did.

Compartments also differ in pH, and the formation energies are transformed at a pH. Using
the cytosol's 7.5 everywhere puts the matrix and the vacuole at the wrong protonation.
"""

import warnings

import pytest

from ystwin.bridge.thermodynamic import (
    COMPARTMENTS,
    FARADAY_KJ_PER_V,
    FREE_MAGNESIUM_M,
    ThermodynamicData,
    electrical_work,
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


class TestCompartments:
    def test_the_ones_the_model_uses_are_described(self):
        for compartment in ("c", "m", "e", "v"):
            assert compartment in COMPARTMENTS

    def test_the_cytosol_carries_the_measured_ph(self):
        assert COMPARTMENTS["c"].ph == pytest.approx(7.5)

    def test_the_matrix_is_more_alkaline_than_the_cytosol(self):
        assert COMPARTMENTS["m"].ph > COMPARTMENTS["c"].ph

    def test_the_vacuole_is_acidic(self):
        assert COMPARTMENTS["v"].ph < COMPARTMENTS["c"].ph

    def test_the_matrix_is_electrically_negative_to_the_cytosol(self):
        assert COMPARTMENTS["m"].potential_v < COMPARTMENTS["c"].potential_v

    def test_the_cytosol_is_the_reference(self):
        assert COMPARTMENTS["c"].potential_v == 0.0

    def test_every_compartment_states_a_source(self):
        assert all(c.source for c in COMPARTMENTS.values())


class TestElectricalWork:
    def test_moving_nothing_charged_costs_nothing(self):
        assert electrical_work({}) == pytest.approx(0.0)

    def test_a_reaction_inside_one_compartment_costs_nothing(self):
        """**This asserted nothing until 2026-08-31.** It read
        ``{("c", 1): -1, ("c", 1): 1}`` -- a duplicate dict key, so the literal was a single
        term, and it passed only because the cytosol's potential is the 0 V reference. If
        the cancellation in `electrical_work` broke it would still have passed. Written
        against a compartment with a real potential, it now tests what it says."""
        assert electrical_work({("m", 1): -1, ("m", 2): 0.5}) == pytest.approx(0.0)
        assert electrical_work({("m", 1): -1}) != pytest.approx(0.0)

    def test_a_proton_entering_the_matrix_releases_work(self):
        """The matrix is negative, so a positive charge falling into it is downhill."""
        assert electrical_work({("c", 1): -1, ("m", 1): 1}) < 0

    def test_pumping_one_out_costs_the_same_amount(self):
        into = electrical_work({("c", 1): -1, ("m", 1): 1})
        out = electrical_work({("m", 1): -1, ("c", 1): 1})

        assert out == pytest.approx(-into)

    def test_it_scales_with_charge(self):
        one = electrical_work({("c", 1): -1, ("m", 1): 1})
        two = electrical_work({("c", 2): -1, ("m", 2): 1})

        assert two == pytest.approx(2 * one)

    def test_the_size_is_what_the_potential_implies(self):
        """A proton crossing 160 mV is about 15 kJ/mol, which is not a rounding error."""
        work = abs(electrical_work({("c", 1): -1, ("m", 1): 1}))

        assert 10.0 < work < 20.0
        assert work == pytest.approx(
            FARADAY_KJ_PER_V * abs(COMPARTMENTS["m"].potential_v), rel=1e-6)


class TestMagnesium:
    def test_free_magnesium_is_a_stated_concentration(self):
        assert 1e-4 < FREE_MAGNESIUM_M < 1e-2

    def test_binding_makes_atp_hydrolysis_less_negative(self, thermo):
        """Mg stabilises the more charged species, so it costs the reaction some drive."""
        hydrolysis = {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1}
        without = ThermodynamicData.load(free_magnesium_m=0.0)

        assert reaction_dg0(hydrolysis, thermo) > reaction_dg0(hydrolysis, without)

    def test_it_lands_near_the_textbook_value(self, thermo):
        """-30.5 kJ/mol is quoted at 1 mM Mg, pH 7, and is what this should approach."""
        hydrolysis = {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1}

        assert -40.0 < reaction_dg0(hydrolysis, thermo) < -28.0

    def test_more_magnesium_shifts_it_further(self, thermo_data):
        # thermo_data is taken only for its skip: this needs two non-default loads.
        hydrolysis = {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1}
        low = ThermodynamicData.load(free_magnesium_m=1e-4)
        high = ThermodynamicData.load(free_magnesium_m=1e-2)

        assert reaction_dg0(hydrolysis, high) > reaction_dg0(hydrolysis, low)

    def test_a_metabolite_that_does_not_bind_is_untouched(self, thermo):
        without = ThermodynamicData.load(free_magnesium_m=0.0)

        assert thermo.dgf("s_0750") == pytest.approx(without.dgf("s_0750"))


class TestTheElectricalTermReachesReactionEnergies:
    """Available is not the same as applied. It has to be in the number TMFA uses."""

    def test_a_transport_reaction_gets_an_electrical_term(self, thermo):
        from ystwin.bridge.thermodynamic import _default_model, reaction_energy

        model = _default_model()
        proton_pump = next(r for r in model.reactions
                           if r.id == "r_1110")                    # ADP/ATP transporter
        chemistry, electrical = reaction_energy(proton_pump, thermo, split=True)

        assert electrical != pytest.approx(0.0)

    def test_a_reaction_inside_one_compartment_gets_none(self, thermo):
        from ystwin.bridge.thermodynamic import _default_model, reaction_energy

        model = _default_model()
        enolase = model.reactions.get_by_id("r_0366")
        _, electrical = reaction_energy(enolase, thermo, split=True)

        assert electrical == pytest.approx(0.0)

    def test_the_total_is_the_two_parts(self, thermo):
        from ystwin.bridge.thermodynamic import _default_model, reaction_energy

        reaction = _default_model().reactions.get_by_id("r_1110")
        total = reaction_energy(reaction, thermo)
        chemistry, electrical = reaction_energy(reaction, thermo, split=True)

        assert total == pytest.approx(chemistry + electrical)

    def test_an_uncovered_reaction_returns_nothing(self, thermo):
        from ystwin.bridge.thermodynamic import _default_model, reaction_energy

        model = _default_model()
        uncovered = next(r for r in model.reactions
                         if any(thermo.dgf(m.id) is None for m in r.metabolites))

        assert reaction_energy(uncovered, thermo) is None


class TestChargedTransportWithoutAMechanismCannotBeEstimated:
    """The electrical term is only right if the stoichiometry shows every charge that moves.

    yeast-GEM writes many carriers as a bare uniport of a charged species. Charging those the
    full zF*dPsi gives +77 kJ/mol for PRPP transport and +124 for propionyl-CoA, which no
    carrier does -- the real mechanisms are symports and antiports whose counter-ion the
    stoichiometry omits. Applied anyway it forbids directions the model needs and growth
    stops.

    So a charged transport step with no proton in its own equation has an energy that cannot
    be estimated, and is left unconstrained rather than given a fabricated one. This is what
    Henry 2007 does with any reaction whose energy is unreliable, and it is the difference
    between declining to answer and answering wrongly.
    """

    def test_a_bare_uniport_of_a_highly_charged_species_returns_nothing(self, thermo):
        from ystwin.bridge.thermodynamic import _default_model, reaction_energy

        model = _default_model()
        prpp_transport = model.reactions.get_by_id("r_2023")

        assert reaction_energy(prpp_transport, thermo) is None

    def test_a_proton_coupled_carrier_is_still_estimated(self, thermo):
        from ystwin.bridge.thermodynamic import _default_model, reaction_energy

        model = _default_model()
        leak = model.reactions.get_by_id("r_2129")

        assert reaction_energy(leak, thermo) is not None

    def test_an_electrogenic_antiport_is_still_estimated(self, thermo):
        """The ADP/ATP exchanger swaps ATP(4-) out for ADP(3-) in: net one charge, which is
        a real mechanism the stoichiometry states correctly."""
        from ystwin.bridge.thermodynamic import _default_model, reaction_energy

        model = _default_model()
        carrier = model.reactions.get_by_id("r_1110")

        assert reaction_energy(carrier, thermo) is not None

    def test_an_uncharged_transport_step_is_fine(self, thermo):
        from ystwin.bridge.thermodynamic import _default_model, reaction_energy

        model = _default_model()
        water = model.reactions.get_by_id("r_2096")

        assert reaction_energy(water, thermo) is not None

    def test_a_reaction_inside_one_compartment_is_never_refused_for_this(self, thermo):
        from ystwin.bridge.thermodynamic import _default_model, reaction_energy

        model = _default_model()
        enolase = model.reactions.get_by_id("r_0366")

        assert reaction_energy(enolase, thermo) is not None


def _thermodb():
    """The vendored database. The bundled pickle unpacks through a numpy call deprecated
    in 2.4, which the suite turns into an error -- `_load` suppresses it the same way."""
    from pytfa.io import load_thermoDB

    from ystwin.bridge.thermodynamic import _table

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_thermoDB(str(_table("thermo_data.thermodb")))


class TestTheTableIsReadInItsOwnUnit:
    """The database is kcal/mol and says so; ``pytfa`` defaults to kJ/mol.

    Until 2026-08-30 the unit was not passed, so the pH and ionic-strength transform was
    computed in kJ on top of a tabulated energy in kcal and every reported energy was on a
    scale that was neither. What settles it is a value from outside both estimators, which
    is what the first test here is.
    """

    #: Transformed standard formation energy of water at pH 7.0, I = 0.25 M, kJ/mol.
    #: Alberty, *Thermodynamics of Biochemical Reactions* (2003), table 4.2.
    ALBERTY_WATER_KJ = -155.66

    def test_water_lands_on_the_published_value(self):
        """Reading the table in its own unit puts water within a kJ of the reference.
        Not passing the unit puts it at +24.85, and no scale on which the formation
        energy of water is positive is a scale."""
        from ystwin.bridge.thermodynamic import KJ_PER_KCAL

        entry = _thermodb()["metabolites"]
        from pytfa.thermo.metabolite import MetaboliteThermo

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            water = MetaboliteThermo(entry["cpd00001"], pH=7.0, ionicStr=0.25,
                                     temperature=298.15,
                                     thermo_unit="kcal/mol").deltaGf_tr

        assert water * KJ_PER_KCAL == pytest.approx(self.ALBERTY_WATER_KJ, abs=1.0)

    def test_the_database_declares_the_unit_this_reads_it_in(self):
        """If a future vendored copy ships in kJ/mol this test fails rather than the
        numbers silently moving by 4.184."""
        assert _thermodb()["units"] == "kcal/mol"

    @pytest.mark.integration
    def test_an_independent_estimator_prefers_the_corrected_reading(self):
        """The reason the unit was corrected rather than merely noticed.

        eQuilibrator's component contribution is a different method from a different
        group. Scored against it over every single-compartment reaction all three cover,
        the corrected reading is closer on every summary: this is what makes the change a
        fix rather than a preference.
        """
        pytest.importorskip("equilibrator_api")
        import dataclasses

        import numpy as np

        from ystwin.bridge.equilibrator import EquilibratorData
        from ystwin.bridge.thermodynamic import (
            GAS_CONSTANT_KJ,
            STANDARD_TEMPERATURE_K,
            _LOG_K_MAGNESIUM,
            _default_model,
        )

        model = _default_model()
        if model is None:
            pytest.skip("needs yeast-GEM; see docs/REPRODUCING.md")
        fixed = ThermodynamicData.load(model=model)
        reference = EquilibratorData.load(model=model)
        database = _thermodb()["metabolites"]

        # Rebuild the pre-2026-08-30 table: the transform in kJ over a table in kcal.
        from pytfa.thermo.metabolite import MetaboliteThermo

        mixed = {}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for metabolite, seed in fixed.seed_by_metabolite.items():
                value = MetaboliteThermo(database[seed], pH=fixed.ph,
                                         ionicStr=fixed.ionic_strength,
                                         temperature=STANDARD_TEMPERATURE_K).deltaGf_tr
                log_k = _LOG_K_MAGNESIUM.get(seed)
                if log_k is not None and fixed.free_magnesium_m > 0:
                    bound = 10.0**log_k * fixed.free_magnesium_m
                    value -= GAS_CONSTANT_KJ * STANDARD_TEMPERATURE_K * np.log1p(bound)
                mixed[metabolite] = float(value)
        before = dataclasses.replace(fixed, by_metabolite=mixed)

        rows = []
        for reaction in model.reactions:
            if len({m.compartment for m in reaction.metabolites}) != 1:
                continue
            stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
            new = reaction_dg0(stoichiometry, fixed)
            old = reaction_dg0(stoichiometry, before)
            if new is None or old is None:
                continue
            try:
                theirs = reference.reaction_dg0(stoichiometry)
            except Exception:
                continue
            if theirs is None:
                continue
            rows.append((old, new, theirs))

        assert len(rows) > 500, "too few shared reactions to compare on"
        old, new, theirs = (np.array([r[i] for r in rows]) for i in range(3))

        assert np.corrcoef(new, theirs)[0, 1] > np.corrcoef(old, theirs)[0, 1]
        assert np.sqrt(np.mean((new - theirs) ** 2)) < np.sqrt(np.mean((old - theirs) ** 2))
        assert np.median(np.abs(new - theirs)) < np.median(np.abs(old - theirs))
