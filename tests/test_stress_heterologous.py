"""PHB is installed and balanced, export costs ATP, and the precursor budget is not vacuous.

Every number below was measured on the vendored yeast-GEM v9.0.2 and is stated with the
regime it was measured in, because the growth cost and the precursor floor both move an
order of magnitude with the glucose bound.
"""

import cobra
import pytest

from ystwin.fba import stress_heterologous as sh
from ystwin.fba.audit import precursor_floor
from ystwin.fba.carotenoid import add_beta_carotene_pathway
from ystwin.fba.fva import product_flux_range
from ystwin.fba.physiology import cap_uptake
from ystwin.fba.solver import FVA_PROCESSES

pytestmark = pytest.mark.integration

GLUCOSE, BIOMASS = "r_1714", "r_2111"


@pytest.fixture(scope="module")
def _phb_model_template(yeast_gem_factory):
    # Opt in, so the tests below still pin what phaA does. The DEFAULT omits it.
    return sh.add_phb_pathway(yeast_gem_factory(), heterologous_thiolase=True)


@pytest.fixture
def phb_model(_phb_model_template, model_copy):
    return model_copy(_phb_model_template)


@pytest.fixture(scope="module")
def _phb_default_template(yeast_gem_factory):
    return sh.add_phb_pathway(yeast_gem_factory())


@pytest.fixture
def phb_default(_phb_default_template, model_copy):
    return model_copy(_phb_default_template)


def test_the_default_install_omits_the_duplicate_thiolase(phb_default):
    """phaA is off by default: it duplicates reversible r_0103 and only adds a cycle."""
    assert "PHAA" not in {r.id for r in phb_default.reactions}
    for rid in ("PHAB", "PHAC", sh.PHB_DEMAND_ID):
        assert rid in {r.id for r in phb_default.reactions}


def test_omitting_phaa_leaves_the_native_thiolase_bounded(yeast_gem, phb_default):
    """r_0103 keeps its native -14.358 floor; with phaA installed it runs to -1000."""
    def window(model):
        with model as m:
            m.reactions.get_by_id("r_1714").lower_bound = -10.0
            m.reactions.get_by_id("r_1992").lower_bound = -1000.0
            fva = cobra.flux_analysis.flux_variability_analysis(
                m, reaction_list=["r_0103"], fraction_of_optimum=0.0, processes=FVA_PROCESSES)
        return float(fva.loc["r_0103", "minimum"]), float(fva.loc["r_0103", "maximum"])

    shipped_low, _ = window(yeast_gem)
    low, high = window(phb_default)
    assert low == pytest.approx(shipped_low, rel=1e-9)
    assert high == pytest.approx(10.960766, rel=1e-6)


def test_omitting_phaa_costs_nothing_at_the_phb_ceiling(phb_model, phb_default):
    """The ceiling is identical with phaA and without it, so phaA is a strict negative."""
    def ceiling(model):
        with model as m:
            m.reactions.get_by_id("r_1714").lower_bound = -10.0
            m.reactions.get_by_id("r_1992").lower_bound = -1000.0
            m.objective = sh.PHB_DEMAND_ID
            return float(m.slim_optimize())

    assert ceiling(phb_default) == pytest.approx(ceiling(phb_model), rel=1e-12)


@pytest.fixture(scope="module")
def _carotenoid_template(yeast_gem_factory):
    return add_beta_carotene_pathway(yeast_gem_factory())


@pytest.fixture
def carotenoid(_carotenoid_template, model_copy):
    return model_copy(_carotenoid_template)


@pytest.fixture(scope="module")
def _insulin_model_template(yeast_gem_factory):
    from ystwin.fba.insulin import add_insulin_precursor_pathway

    return add_insulin_precursor_pathway(yeast_gem_factory())


@pytest.fixture
def insulin_model(_insulin_model_template, model_copy):
    return model_copy(_insulin_model_template)


def _max_growth(model, uptake=10.0):
    with model as m:
        cap_uptake(m, GLUCOSE, uptake)
        m.objective = BIOMASS
        return m.slim_optimize()


def _exported(model, cargo, **kwargs):
    """A copy of ``model`` with ``cargo`` given a priced export, and that reaction."""
    out = model.copy()
    kwargs.setdefault("provenance", "yeast-GEM v9.0.2 ABC-system convention")
    return out, sh.add_atp_coupled_export(
        out, cargo, atp_per_molecule=sh.ATP_PER_ABC_EXPORT, **kwargs)


# --------------------------------------------------------------------------- #
# PHB
# --------------------------------------------------------------------------- #


def test_every_added_phb_reaction_is_mass_and_charge_balanced(phb_model):
    """cobra 0.32.1 puts charge in the same dict, so empty covers both -- pinned below."""
    for rid in sh.PHB_REACTION_IDS:
        rxn = phb_model.reactions.get_by_id(rid)
        if rxn.boundary:
            continue
        balance = rxn.check_mass_balance()
        assert balance == {}, f"{rid}: {balance}"
        assert sum(c * (m.charge or 0) for m, c in rxn.metabolites.items()) == 0


def test_check_mass_balance_really_does_report_charge(yeast_gem):
    """Otherwise every "balanced" assertion above is silent about half of what it claims."""
    with yeast_gem as m:
        rxn = m.reactions.get_by_id("r_0103")
        rxn.add_metabolites({m.metabolites.get_by_id("s_0794"): 1.0})

        assert rxn.check_mass_balance() == {"charge": 1.0, "H": 1.0}


def test_the_thioester_formula_is_derived_and_agrees_with_yeast9s_own_species(yeast_gem,
                                                                              phb_model):
    """Solved from phaB's balance against the host's NADPH/NADP+, never typed.

    On yeast-GEM v9.0.2 that lands exactly on the model's own peroxisomal s_2902, which is
    the cross-check; the derivation, not the copy, is what put it there.
    """
    hbcoa = phb_model.metabolites.get_by_id(sh.HBCOA_ID)
    reference = yeast_gem.metabolites.get_by_id("s_2902")

    assert (hbcoa.formula, hbcoa.charge) == ("C25H38N7O18P3S", -4)
    assert (hbcoa.formula, hbcoa.charge) == (reference.formula, reference.charge)


def test_copying_that_species_would_have_been_wrong_on_the_gecko_export(ec_yeast_gem):
    """ecYeastGEM 8.3.4 writes s_2902 as the neutral acid and every cytosolic thioester as
    the anion skeleton -- four hydrogens apart, in one model. Deriving is what survives it."""
    aaccoa = ec_yeast_gem.metabolites.get_by_id("s_0367[c]")
    reference = ec_yeast_gem.metabolites.get_by_id("s_2902[p]")

    assert (aaccoa.formula, reference.formula) == ("C25H36N7O18P3S", "C25H42N7O18P3S")


def test_the_pathway_installs_balanced_on_the_gecko_export_too(ec_yeast_gem):
    # Opt in here: on an enzyme-constrained host phaA carries its own cost, so it is real.
    installed = sh.add_phb_pathway(ec_yeast_gem, heterologous_thiolase=True)

    for rid in ("PHAA", "PHAB", "PHAC"):
        assert installed.reactions.get_by_id(rid).check_mass_balance() == {}
    assert installed.metabolites.get_by_id(sh.HBCOA_ID).formula == "C25H38N7O18P3S"


def test_the_repeat_unit_is_the_monomer_and_not_the_polymer(phb_model):
    """(C4H6O2)n, KEGG C06143; the polymer has no molar mass, so PHB is per repeat unit."""
    phb = phb_model.metabolites.get_by_id(sh.PHB_REPEAT_UNIT_ID)

    assert (phb.formula, phb.charge) == ("C4H6O2", 0)


def test_the_pathway_draws_carbon_from_the_native_cytosolic_acetyl_coa_pool(phb_model):
    assert phb_model.metabolites.get_by_id("s_0373") in \
        phb_model.reactions.get_by_id("PHAA").reactants


def test_the_reduction_is_nadph_dependent_as_ec_1_1_1_36(phb_model):
    by_name = {m.name: c for m, c in phb_model.reactions.get_by_id("PHAB").metabolites.items()}

    assert by_name["NADPH"] == -1.0 and by_name["NADP(+)"] == 1.0
    assert "NADH" not in by_name


def test_reactions_carry_the_real_construct_gene_names(phb_model):
    rules = {rid: phb_model.reactions.get_by_id(rid).gene_reaction_rule
             for rid in ("PHAA", "PHAB", "PHAC")}

    assert rules == {"PHAA": "phaA", "PHAB": "phaB", "PHAC": "phaC"}


def test_adding_a_drain_pathway_does_not_raise_maximum_growth(yeast_gem, phb_model):
    """0.887685358559 /h at glucose <= 10 with oxygen free, before and after."""
    assert _max_growth(phb_model) == pytest.approx(_max_growth(yeast_gem), rel=1e-9)


def test_every_added_phb_reaction_can_carry_flux(phb_model):
    """FVA at fraction_of_optimum=0, glucose <= 10: all four reach 10.960766 loopless.

    PHAA's unbounded maximum is 1000 and that is the thermodynamic loop with r_0103, so it
    is measured here with the native thiolase made irreversible -- see the test below.
    """
    with phb_model as m:
        cap_uptake(m, GLUCOSE, 10.0)
        m.reactions.get_by_id("r_0103").bounds = (0.0, 1000.0)
        fva = cobra.flux_analysis.flux_variability_analysis(
            m, reaction_list=list(sh.PHB_REACTION_IDS), fraction_of_optimum=0.0,
            processes=FVA_PROCESSES)

    for rid in sh.PHB_REACTION_IDS:
        assert fva.loc[rid, "maximum"] == pytest.approx(10.960766, rel=1e-6)


def test_phb_can_carry_flux_and_the_ceiling_is_359x_the_derived_productivity(phb_model):
    """1.1664 mmol/gDCW/h at glucose <= 10, mu held at 90% of max.

    Against 3.2524 umol/gDW/h, the SCKK006 rate that section 4b of
    data/phb/SOURCE_kocharin2012_genotype.md derives from Table 3 and cross-checks three
    ways. Section 4d of that file REFUSES the text's own 99.3 umol/gDW/h -- it is 30.5x
    the paper's tables -- so the ceiling is 359x the usable number, not 11.7x the refused
    one.
    """
    envelope = product_flux_range(phb_model, sh.PHB_DEMAND_ID, glucose_uptake=10.0)

    assert envelope.maximum == pytest.approx(1.1664023423, rel=1e-6)
    assert envelope.maximum / 0.0032524 == pytest.approx(358.6, abs=0.5)


def test_the_heterologous_thiolase_adds_no_capacity_over_the_native_one(phb_model):
    """phaA duplicates the reversible native r_0103, as crtE duplicates BTS1."""
    with phb_model as m:
        cap_uptake(m, GLUCOSE, 10.0)
        m.objective = sh.PHB_DEMAND_ID
        both = m.slim_optimize()
        m.reactions.get_by_id("PHAA").bounds = (0.0, 0.0)
        native_only = m.slim_optimize()

    assert native_only == pytest.approx(both, rel=1e-9)


def test_phaas_unbounded_fva_maximum_is_a_loop_with_r_0103_not_a_capacity(phb_model):
    """PHAA reaches 1000 only by cycling against the reversible native thiolase."""
    with phb_model as m:
        cap_uptake(m, GLUCOSE, 10.0)
        looped = cobra.flux_analysis.flux_variability_analysis(
            m, reaction_list=["PHAA"], fraction_of_optimum=0.0, processes=FVA_PROCESSES)
        m.reactions.get_by_id("r_0103").bounds = (0.0, 1000.0)
        loopless = cobra.flux_analysis.flux_variability_analysis(
            m, reaction_list=["PHAA"], fraction_of_optimum=0.0, processes=FVA_PROCESSES)

    assert looped.loc["PHAA", "maximum"] == pytest.approx(1000.0)
    assert loopless.loc["PHAA", "maximum"] == pytest.approx(10.960766, rel=1e-6)


def test_that_loop_also_makes_the_NATIVE_thiolase_unbounded(yeast_gem, phb_model):
    """r_0103 goes from [-14.358, 0] to [-1000, 10.961] the moment PHAA is installed.

    Nothing is created -- with every boundary shut to uptake the cycle carries no ATP, no
    NADPH and no biomass -- but any FVA read on ERG10 after this install is the cycle.
    """
    def envelope(model):
        with model as m:
            cap_uptake(m, GLUCOSE, 10.0)
            fva = cobra.flux_analysis.flux_variability_analysis(
                m, reaction_list=["r_0103"], fraction_of_optimum=0.0, processes=FVA_PROCESSES)
        return float(fva.loc["r_0103", "minimum"]), float(fva.loc["r_0103", "maximum"])

    assert envelope(yeast_gem) == pytest.approx((-14.35816, 0.0), abs=1e-4)
    assert envelope(phb_model) == pytest.approx((-1000.0, 10.960766), rel=1e-6)


def test_the_loop_makes_nothing_from_nothing(phb_model):
    """With uptake shut off, the futile cycle spins but yields no ATP, NADPH or biomass."""
    with phb_model as m:
        for r in m.reactions:
            if r.boundary:
                r.bounds = (0.0, 1000.0)
        for rid in ("r_4046", "r_4041"):
            m.reactions.get_by_id(rid).lower_bound = 0.0
        free_atp = cobra.Reaction("FREE_ATP", lower_bound=0.0, upper_bound=1000.0)
        m.add_reactions([free_atp])
        free_atp.add_metabolites({m.metabolites.get_by_id(i): c for i, c in (
            ("s_0434", -1), ("s_0803", -1), ("s_0394", 1), ("s_1322", 1), ("s_0794", 1))})
        results = {}
        for rid in ("r_2111", "FREE_ATP", "DM_phb_c", "PHAA"):
            m.objective = rid
            results[rid] = m.slim_optimize()

    assert results["r_2111"] == pytest.approx(0.0, abs=1e-9)
    assert results["FREE_ATP"] == pytest.approx(0.0, abs=1e-9)
    assert results["DM_phb_c"] == pytest.approx(0.0, abs=1e-9)
    assert results["PHAA"] == pytest.approx(1000.0)


def test_installing_phb_is_what_connects_the_dead_end_acetoacetyl_coa_pool(yeast_gem, phb_model):
    """s_0367 has exactly one reaction able to consume it until phaB arrives."""
    before = sh.native_drain(yeast_gem, "s_0367", glucose_uptake=10.0)
    after = sh.native_drain(phb_model, "s_0367", glucose_uptake=10.0)

    assert before.consumers == ("r_0103",) and before.vacuous
    assert after.consumers == ("PHAB", "r_0103") and not after.vacuous


def test_a_second_install_is_refused(phb_model):
    with pytest.raises(ValueError, match="already installed"):
        sh.add_phb_pathway(phb_model)


def test_a_host_missing_a_precursor_is_named_rather_than_silently_skipped(yeast_gem):
    with pytest.raises(KeyError, match="s_9999"):
        sh.add_phb_pathway(yeast_gem, aaccoa_id="s_9999")


def test_a_host_whose_thioester_bookkeeping_does_not_yield_the_repeat_unit_is_refused(yeast_gem):
    """Point coa_id at water and the derivation lands on something that is not C4H6O2."""
    with pytest.raises(ValueError, match="C06143"):
        sh.add_phb_pathway(yeast_gem, coa_id="s_0803")


# --------------------------------------------------------------------------- #
# the export cost
# --------------------------------------------------------------------------- #


def test_the_abc_rate_is_read_off_the_host_model_not_remembered(yeast_gem):
    """All twelve reactions yeast-GEM names "via ABC system" charge exactly one ATP."""
    atp = yeast_gem.metabolites.get_by_id("s_0434")
    coefficients = {r.id: r.metabolites.get(atp)
                    for r in yeast_gem.reactions if "abc system" in (r.name or "").lower()}

    assert len(coefficients) == 12
    assert set(coefficients.values()) == {-sh.ATP_PER_ABC_EXPORT}


def test_export_without_a_cost_is_refused_rather_than_made_free(carotenoid):
    with pytest.raises(sh.SecretionCostUnknown, match="will not invent"):
        sh.add_atp_coupled_export(carotenoid.copy(), "betacarotene_c",
                                  atp_per_molecule=None, provenance="anything")


def test_a_cost_without_provenance_is_refused(carotenoid):
    with pytest.raises(ValueError, match="provenance is required"):
        sh.add_atp_coupled_export(carotenoid.copy(), "betacarotene_c",
                                  atp_per_molecule=1.0, provenance="   ")


def test_the_export_reaction_is_mass_and_charge_balanced(carotenoid):
    model, rxn = _exported(carotenoid, "betacarotene_c")

    assert rxn.check_mass_balance() == {}
    assert sum(c * (m.charge or 0) for m, c in rxn.metabolites.items()) == 0
    assert {m.id for m in rxn.reactants} == {"betacarotene_c", "s_0434", "s_0803"}
    assert {m.id for m in rxn.products} == {"betacarotene_e", "s_0394", "s_1322", "s_0794"}


def test_the_export_is_a_transport_reaction_and_the_exchange_is_the_separate_one(carotenoid):
    """`EX_` means a boundary exchange in every model here, so the transport is `SEC_`."""
    model, rxn = _exported(carotenoid, "betacarotene_c")

    assert (rxn.id, len(rxn.metabolites)) == ("SEC_betacarotene", 7)
    assert model.reactions.get_by_id("EX_betacarotene_e").boundary


def test_the_new_exchange_cannot_import_the_product_from_a_medium_that_never_had_it(carotenoid):
    """cobra defaults an exchange to (-1000, 1000); a heterologous product needs (0, 1000)."""
    model, _rxn = _exported(carotenoid, "betacarotene_c")

    assert model.reactions.get_by_id("EX_betacarotene_e").bounds == (0.0, 1000.0)


def test_the_exported_id_strips_the_cargos_compartment_rather_than_appending(insulin_model):
    """`insulin_precursor_er` exports to `insulin_precursor_e`, not `insulin_precursor_er_e`."""
    from ystwin.fba.insulin import FOLDED_METABOLITE_ID

    model, rxn = _exported(insulin_model, FOLDED_METABOLITE_ID)

    assert rxn.id == "SEC_insulin_precursor"
    assert model.metabolites.has_id("insulin_precursor_e")


def test_export_can_carry_flux_and_does_not_change_maximum_growth(yeast_gem, carotenoid):
    """FVA at fraction_of_optimum=0, glucose <= 10: [0, 0.87562]. Growth is untouched."""
    model, rxn = _exported(carotenoid, "betacarotene_c")

    with model as m:
        cap_uptake(m, GLUCOSE, 10.0)
        fva = cobra.flux_analysis.flux_variability_analysis(
            m, reaction_list=[rxn.id], fraction_of_optimum=0.0, processes=FVA_PROCESSES)

    assert fva.loc[rxn.id, "minimum"] == pytest.approx(0.0)
    assert fva.loc[rxn.id, "maximum"] == pytest.approx(0.87562, rel=1e-4)
    assert _max_growth(model) == pytest.approx(_max_growth(yeast_gem), rel=1e-9)


def test_the_provenance_and_the_itemised_cost_are_stamped_onto_the_reaction(carotenoid):
    model, rxn = _exported(carotenoid, "betacarotene_c", leader_residues=85,
                           provenance="ASSERTED: no measurement supports this leader")

    assert rxn.notes["secretion_cost_provenance"].startswith("ASSERTED")
    assert rxn.notes["secretion_cost_atp"] == {"transport": 1.0, "leader_translation": 340.0,
                                               "leader_residues": 85}


def test_charging_the_transporter_lowers_the_ceiling_by_half_a_percent(carotenoid):
    """0.0979986 -> 0.0974989 mmol/gDCW/h at glucose <= 10, mu held at 90% of max.

    One ATP per C40 molecule is nearly free. That is the answer, not a failure to find an
    effect: a small molecule's transport cost is small, and saying so is what makes the
    protein number below mean something.
    """
    model, rxn = _exported(carotenoid, "betacarotene_c")

    plain = product_flux_range(carotenoid, "DM_betacarotene_c", glucose_uptake=10.0)
    charged = product_flux_range(model, rxn.id, glucose_uptake=10.0)

    assert plain.maximum == pytest.approx(0.0979985514, rel=1e-6)
    assert charged.maximum == pytest.approx(0.0974988783, rel=1e-6)
    assert (plain.maximum - charged.maximum) / plain.maximum == pytest.approx(0.0051, abs=1e-4)


def test_an_85_residue_leader_costs_64_percent_of_the_ceiling(carotenoid):
    """341 ATP-equivalents per molecule: 1 for the pump, 340 for a leader that is discarded.

    The leader length is the CALLER's, deliberately, so this is the arithmetic of the
    surcharge and not a claim about any construct.
    """
    model, rxn = _exported(carotenoid, "betacarotene_c", leader_residues=85,
                           provenance="ASSERTED length")

    plain = product_flux_range(carotenoid, "DM_betacarotene_c", glucose_uptake=10.0)
    charged = product_flux_range(model, rxn.id, glucose_uptake=10.0)

    assert rxn.metabolites[model.metabolites.get_by_id("s_0434")] == -341.0
    assert charged.maximum == pytest.approx(0.0356670331, rel=1e-6)
    assert (plain.maximum - charged.maximum) / plain.maximum == pytest.approx(0.636, abs=1e-3)


def test_the_leader_rate_is_the_full_four_atp_not_the_elongation_half():
    """4 per residue (Martin 2025, PMID 40562331), because a leader charged as pure ATP
    hydrolysis never draws on the host's aminoacyl-tRNA pool and so never pays the
    synthetase half through the model's own reactions."""
    from ystwin.fba.insulin import GTP_PER_RESIDUE

    assert sh.ATP_EQUIVALENTS_PER_RESIDUE == 4.0
    assert sh.leader_peptide_atp_equivalents(85) == 340.0
    assert sh.ATP_EQUIVALENTS_PER_RESIDUE == 2 * GTP_PER_RESIDUE


def test_the_mfalpha1_leader_length_is_the_one_number_here_with_a_sequence_behind_it():
    """UniProt P01149: SIGNAL 1..19, PROPEP 20..89, Kex2 site Lys84-Arg85, EAEA spacer 86-89."""
    assert sh.MFALPHA1_PREPRO_RESIDUES == 85
    assert sh.MFALPHA1_PREPRO_WITH_SPACER_RESIDUES == 89
    assert sh.leader_peptide_atp_equivalents(sh.MFALPHA1_PREPRO_RESIDUES) == 340.0


def test_pricing_the_export_moves_the_published_insulin_ceiling_by_15_percent(insulin_model):
    """0.0098978 -> 0.0098925 (pump) -> 0.0083701 (pump + 85-residue leader) mmol/gDCW/h.

    Shipped default, growth off -- the regime `fba/insulin.py` reports 0.009898 under. That
    module says leaving the product on a demand "is the largest omission here" and that it
    "biases the ceiling UPWARD". Measured, it does, by 15.4%, and 99.7% of that is the
    leader rather than the pump.
    """
    from ystwin.fba.insulin import FOLDED_METABOLITE_ID, PRODUCT_DEMAND_ID

    pumped, pump_rxn = _exported(insulin_model, FOLDED_METABOLITE_ID)
    secreted, sec_rxn = _exported(insulin_model, FOLDED_METABOLITE_ID,
                                  leader_residues=sh.MFALPHA1_PREPRO_RESIDUES,
                                  provenance="ABC convention; leader from UniProt P01149")

    def ceiling(model, reaction):
        with model as m:
            m.objective = reaction
            return m.slim_optimize()

    demand = ceiling(insulin_model, PRODUCT_DEMAND_ID)

    assert demand == pytest.approx(0.0098978104, rel=1e-6)
    assert ceiling(pumped, pump_rxn.id) == pytest.approx(0.0098925155, rel=1e-6)
    assert ceiling(secreted, sec_rxn.id) == pytest.approx(0.0083701023, rel=1e-6)
    assert (demand - ceiling(secreted, sec_rxn.id)) / demand == pytest.approx(0.1543, abs=1e-4)


def test_the_oxidative_folding_term_is_not_duplicated_here(carotenoid):
    """`fba/secretion.py` already charges one O2 in and one H2O2 out per disulfide bond
    (PMID 16407158). This module prices the export and touches neither species."""
    from ystwin.fba.secretion import DISULFIDE_O2_PER_BOND

    _model, rxn = _exported(carotenoid, "betacarotene_c", leader_residues=85,
                            provenance="ASSERTED length")

    assert DISULFIDE_O2_PER_BOND == 1.0
    assert not {"s_1275", "s_1276", "s_0837"} & {m.id for m in rxn.metabolites}


def test_what_is_not_charged_is_enumerated_with_a_reason_each():
    assert len(sh.SECRETION_NOT_CHARGED) >= 5
    for what, why in sh.SECRETION_NOT_CHARGED:
        assert what and why


# --------------------------------------------------------------------------- #
# the precursor budget
# --------------------------------------------------------------------------- #


def test_the_declared_beta_carotene_node_cannot_have_a_floor_at_all(yeast_gem):
    """GGPP has one producer and one consumer, so the zero is topology, not biology.

    `fba/audit.py` reports 0.0 here and its docstring illustrates the check with FPP --
    ergosterol and dolichol -- which is a different metabolite, one step upstream.
    """
    drain = sh.native_drain(yeast_gem, "s_0189", glucose_uptake=10.0)

    assert drain.vacuous
    assert drain.consumers == ("r_0461",)
    assert precursor_floor(yeast_gem, "s_0189") == 0.0


def test_the_whole_ggpp_branch_is_dispensable_so_the_zero_is_real(yeast_gem):
    base = _max_growth(yeast_gem)

    with yeast_gem as m:
        cap_uptake(m, GLUCOSE, 10.0)
        m.reactions.get_by_id("r_0461").bounds = (0.0, 0.0)
        m.objective = BIOMASS
        blocked = m.slim_optimize()

    assert blocked == pytest.approx(base, rel=1e-9)


def test_the_branch_point_one_step_up_has_a_floor_20_7x_the_shipped_number(yeast_gem):
    """FPP: 0.00284914 shipped -> 0.0589248 at glucose <= 10. 10.34x of that is the medium.

    The remaining 2.00x is squalene synthase's -2 coefficient, which `precursor_floor`
    never applies -- NOT the two reversible consumers it drops, which add 8.0e-7 between
    them. Both causes are isolated below so the factor cannot be misattributed again.
    """
    shipped = precursor_floor(yeast_gem, "s_0190")
    drain = sh.native_drain(yeast_gem, "s_0190", glucose_uptake=10.0)
    with yeast_gem as m:
        cap_uptake(m, GLUCOSE, 10.0)
        regime_only = precursor_floor(m, "s_0190")

    assert shipped == pytest.approx(0.00284914, rel=1e-4)
    assert regime_only == pytest.approx(0.0294620012, rel=1e-6)
    assert drain.floor == pytest.approx(0.0589248, rel=1e-4)
    assert not drain.vacuous
    assert set(drain.consumers) - {"r_0373", "r_1012", "r_4768"} == {"r_1766", "r_4604"}
    fpp = yeast_gem.metabolites.get_by_id("s_0190")
    assert yeast_gem.reactions.get_by_id("r_1012").metabolites[fpp] == -2.0
    assert drain.floor - 2 * regime_only == pytest.approx(7.99e-07, rel=1e-2)


def test_the_phb_precursor_floor_is_10_3x_the_shipped_number_and_the_gap_is_the_medium(yeast_gem):
    """acetyl-CoA: 0.0983969 shipped -> 1.01749 at glucose <= 10. All of it is the regime.

    The four reversible consumers `precursor_floor` drops -- r_0103 among them, the very
    thiolase phaA duplicates -- add exactly nothing here, and reporting that is the point.
    """
    shipped = precursor_floor(yeast_gem, "s_0373")
    drain = sh.native_drain(yeast_gem, "s_0373", glucose_uptake=10.0)

    with yeast_gem as m:
        cap_uptake(m, GLUCOSE, 10.0)
        regime_only = precursor_floor(m, "s_0373")

    assert shipped == pytest.approx(0.0983969, rel=1e-4)
    assert regime_only == pytest.approx(1.01749, rel=1e-4)
    assert drain.floor == pytest.approx(regime_only, rel=1e-6)
    assert "r_0103" in drain.consumers


def test_the_floor_refuses_to_be_read_without_a_medium(yeast_gem):
    with pytest.raises(ValueError, match="required rather than defaulted"):
        sh.native_drain(yeast_gem, "s_0373", glucose_uptake=0.0)


def test_an_absent_precursor_is_named_rather_than_returning_zero(yeast_gem):
    with pytest.raises(KeyError, match="s_9999"):
        sh.native_drain(yeast_gem, "s_9999", glucose_uptake=10.0)
