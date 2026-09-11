"""Protein synthesis as a priced sink, and what that price does and does not reproduce.

The numbers pinned here are the point of the module. Charging amino acids and translation ATP
puts the burden of 15% of total protein at 6.8% of growth, against a measured limit whose
proportional bound is 15% and whose growth-law bound is 17.7% -- so the chemistry is roughly
half the cost. Taking the same protein out of the GECKO enzyme pool instead gives 15.1%. That
gap is a finding rather than a bug, so it is a test rather than a comment.

The GSMM fixtures are `tests/conftest.py`'s, which resolve through `ystwin.paths` and skip
when the models are absent.
"""

import math
from copy import deepcopy

import cobra
import numpy
import pytest
from cobra.flux_analysis import flux_variability_analysis

from ystwin.fba.carotenoid import add_beta_carotene_pathway
from ystwin.fba.stress_proteostasis import (
    ATP_EQUIVALENTS_PER_PEPTIDE_BOND,
    CHARGING_ATP_EQUIVALENTS,
    CRT_ENZYME_MASSES_KDA,
    CRT_ENZYME_OF_REACTION,
    CRT_KCAT_PROVENANCE,
    EGUCHI_BURDEN_LIMIT,
    ELONGATION_ATP_EQUIVALENTS,
    HAC1_CHAPERONE_ORFS,
    HSF1_CHAPERONE_ORFS,
    METZL_RAZ_RIBOSOME_RESERVE,
    METZL_RAZ_RIBOSOME_SLOPE,
    METZL_RAZ_RIBOSOME_SLOPE_PER_GENERATION,
    add_chaperone_sink,
    add_heterologous_protein_sink,
    add_protein_sink,
    attach_crt_enzyme_demands,
    attach_enzyme_demand,
    burden_against_literature,
    burden_curve,
    charge_balance,
    couple_sink_to_growth,
    formula_weight,
    grams_per_gdcw_at_protein_fraction,
    mass_balance,
    protein_pool_burden_curve,
    read_protein_composition,
    ribosome_allocation_penalty,
    scale_protein_pool,
)

pytestmark = pytest.mark.integration

GLUCOSE, OXYGEN, BIOMASS = "r_1714", "r_1992", "r_2111"
CRT_PLACEHOLDER = "PLACEHOLDER: no kcat is published for any X. dendrorhous crt enzyme"


def carbon_limited(model):
    """Glucose 10 mmol/gDW/h, oxygen free: the audit's reference regime, base mu 0.887685 /h."""
    model.reactions.get_by_id(GLUCOSE).lower_bound = -10.0
    model.reactions.get_by_id(OXYGEN).lower_bound = -1000.0
    return model


def oxygen_limited(model):
    """The repo's REFERENCE_AEROBIC_BATCH: glucose 11.1, oxygen 3.7, base mu 0.345984 /h."""
    model.reactions.get_by_id(GLUCOSE).lower_bound = -11.1
    model.reactions.get_by_id(OXYGEN).lower_bound = -3.7
    return model


@pytest.fixture(scope="module")
def _sink_model_template(yeast_gem_factory):
    return add_heterologous_protein_sink(yeast_gem_factory())


@pytest.fixture
def sink_model(_sink_model_template, model_copy):
    return model_copy(_sink_model_template)


@pytest.fixture(scope="module")
def _chaperone_model_template(yeast_gem_factory):
    return add_chaperone_sink(yeast_gem_factory())


@pytest.fixture
def chaperone_model(_chaperone_model_template, model_copy):
    return model_copy(_chaperone_model_template)


@pytest.fixture(scope="module")
def _composition_template(yeast_gem_factory):
    return read_protein_composition(yeast_gem_factory())


@pytest.fixture
def composition(_composition_template):
    return deepcopy(_composition_template)


def test_the_composition_is_the_models_own_not_one_we_typed_in(yeast_gem, composition):
    """Every coefficient must be traceable to r_4047."""
    source = yeast_gem.reactions.get_by_id("r_4047")

    assert composition.source_reaction == "r_4047"
    assert composition.aa_coefficients["s_0404"] == pytest.approx(0.527012401964609, rel=1e-12)
    assert composition.aa_coefficients["s_1561"] == pytest.approx(0.303939602869103, rel=1e-12)
    assert len(composition.aa_coefficients) == 20
    for mid, coef in composition.aa_coefficients.items():
        assert -source.metabolites[yeast_gem.metabolites.get_by_id(mid)] == pytest.approx(coef)


def test_the_polymer_mass_and_residue_count_are_physical(composition):
    """0.46 g protein/gDCW and a ~109 g/mol mean residue are what any yeast proteome gives."""
    assert composition.grams_per_unit == pytest.approx(0.459757, abs=1e-5)
    assert composition.residues_per_unit == pytest.approx(4.233559, abs=1e-5)
    assert composition.residues_per_gram == pytest.approx(9.20825, abs=1e-4)
    assert 105.0 < composition.mean_residue_mass < 115.0
    # Yeast proteins are net acidic; the sign is read from the model's own tRNA charges.
    assert composition.polymer_charge_per_gram < 0.0


def test_the_two_model_versions_agree_on_the_mean_residue(composition, ec_yeast_gem):
    ec_composition = read_protein_composition(ec_yeast_gem)

    assert ec_composition.residues_per_gram == pytest.approx(composition.residues_per_gram, rel=1e-5)


def test_the_atp_charge_is_the_two_elongation_gtps_only(sink_model, composition):
    """Charging is already paid by the aminoacyl-tRNA synthetases; only elongation is added."""
    assert CHARGING_ATP_EQUIVALENTS + ELONGATION_ATP_EQUIVALENTS == ATP_EQUIVALENTS_PER_PEPTIDE_BOND
    rxn = sink_model.reactions.get_by_id("PROTSYN_hetprot")
    atp = sink_model.metabolites.get_by_id("s_0434")

    assert rxn.metabolites[atp] == pytest.approx(
        -ELONGATION_ATP_EQUIVALENTS * composition.residues_per_gram
    )


def test_the_model_already_pays_the_other_two_equivalents(yeast_gem, composition):
    """Read out of the model, because that is what makes charging only two defensible.

    Every one of the 20 synthetases runs ATP -> AMP + PPi (two high-energy bonds), and r_0568
    hydrolyses the PPi so it cannot be recovered.
    """
    atp = yeast_gem.metabolites.get_by_id("s_0434")
    amp = yeast_gem.metabolites.get_by_id("s_0423")
    ppi = yeast_gem.metabolites.get_by_id("s_0633")
    charging = []
    for mid in composition.aa_coefficients:
        met = yeast_gem.metabolites.get_by_id(mid)
        charging.append(
            any(
                rxn.metabolites[met] > 0
                and rxn.metabolites.get(atp, 0) < 0
                and rxn.metabolites.get(amp, 0) > 0
                and rxn.metabolites.get(ppi, 0) > 0
                for rxn in met.reactions
                if rxn.id != "r_4047"
            )
        )

    assert all(charging) and len(charging) == 20
    assert yeast_gem.reactions.get_by_id("r_0568").metabolites[ppi] == -1.0


def test_the_maintenance_atp_convention_is_read_not_written(yeast_gem, ec_yeast_gem, sink_model):
    """r_4041 carries 55.3 ATP in Yeast9 and 56.6883 in the ec model; the sink uses neither.

    It uses the *shape* -- one ATP, one water, one ADP, one phosphate, one proton -- so the two
    models' different formula conventions both come out balanced.
    """
    atp = yeast_gem.reactions.get_by_id("r_4041").metabolites[
        yeast_gem.metabolites.get_by_id("s_0434")]
    ec_atp = ec_yeast_gem.reactions.get_by_id("r_4041").metabolites[
        ec_yeast_gem.metabolites.get_by_id("s_0434[c]")]

    assert atp == pytest.approx(-55.3)
    assert ec_atp == pytest.approx(-56.6883)
    assert mass_balance(sink_model.reactions.get_by_id("PROTSYN_hetprot")) == {}


def _forced_hydrolysis_cost(model, stoichiometry, rate):
    """Growth when a bare hydrolysis is forced at ``rate`` mmol/gDCW/h."""
    with model as scoped:
        carbon_limited(scoped)
        rxn = cobra.Reaction("FORCE_HYDROLYSIS", lower_bound=rate, upper_bound=rate)
        scoped.add_reactions([rxn])
        rxn.add_metabolites(
            {scoped.metabolites.get_by_id(k): v for k, v in stoichiometry.items()}
        )
        balance = mass_balance(rxn)
        return balance, scoped.slim_optimize()


def test_gtp_is_exactly_atp_equivalent_here_so_charging_atp_is_not_an_approximation(yeast_gem):
    """r_0800 runs ATP + GDP -> ADP + GTP, so the two currencies cost the same."""
    rate = 18.4164973365255
    atp_balance, atp_growth = _forced_hydrolysis_cost(
        yeast_gem,
        {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_0794": 1, "s_1322": 1},
        rate,
    )
    gtp_balance, gtp_growth = _forced_hydrolysis_cost(
        yeast_gem,
        {"s_0785": -1, "s_0803": -1, "s_0739": 1, "s_0794": 1, "s_1322": 1},
        rate,
    )

    assert atp_balance == {} and gtp_balance == {}
    assert atp_growth == pytest.approx(0.802192, abs=1e-5)
    assert gtp_growth == pytest.approx(atp_growth, rel=1e-9)


def test_the_ec_model_cannot_balance_a_gtp_hydrolysis_but_can_balance_an_atp_one(ec_yeast_gem):
    """Why the cost is charged in ATP: yeast 8.3.4 protonates GTP and not GDP."""
    atp_balance, _ = _forced_hydrolysis_cost(
        ec_yeast_gem,
        {"s_0434[c]": -1, "s_0803[c]": -1, "s_0394[c]": 1, "s_0794[c]": 1, "s_1322[c]": 1},
        0.0,
    )
    gtp_balance, _ = _forced_hydrolysis_cost(
        ec_yeast_gem,
        {"s_0785[c]": -1, "s_0803[c]": -1, "s_0739[c]": 1, "s_0794[c]": 1, "s_1322[c]": 1},
        0.0,
    )

    assert atp_balance == {}
    assert gtp_balance == {"H": -1.0}


def test_the_sink_is_mass_and_charge_balanced(sink_model, chaperone_model):
    """r_4047 itself is not; this is stricter than the model it was read from."""
    for model, rid in ((sink_model, "PROTSYN_hetprot"), (chaperone_model, "PROTSYN_chaperone")):
        rxn = model.reactions.get_by_id(rid)

        assert mass_balance(rxn) == {}
        assert charge_balance(rxn) == pytest.approx(0.0, abs=1e-9)


def test_one_unit_of_the_protein_metabolite_weighs_one_gram(sink_model):
    protein = sink_model.metabolites.get_by_id("hetprot_c")

    assert formula_weight(protein) == pytest.approx(1000.0, rel=1e-6)


def test_both_sinks_carry_flux(sink_model, chaperone_model):
    """FVA at fraction_of_optimum=0: the floor is zero, the ceiling is ~0.98 g/gDCW/h."""
    for model, sink, demand in (
        (sink_model, "PROTSYN_hetprot", "DM_hetprot_c"),
        (chaperone_model, "PROTSYN_chaperone", "DM_chaperone_c"),
    ):
        with model as scoped:
            carbon_limited(scoped)
            fva = flux_variability_analysis(
                scoped, [sink, demand], fraction_of_optimum=0.0, processes=1
            )

        assert fva.loc[sink, "minimum"] == pytest.approx(0.0, abs=1e-9)
        assert fva.loc[sink, "maximum"] == pytest.approx(0.97845, abs=1e-4)
        assert fva.loc[demand, "maximum"] == pytest.approx(fva.loc[sink, "maximum"], rel=1e-9)


def test_installing_the_sink_does_not_move_growth(yeast_gem, sink_model):
    """A drain can only cost carbon, never create it, and an unused drain costs nothing."""
    for constrain, expected in ((carbon_limited, 0.887685), (oxygen_limited, 0.345984)):
        with yeast_gem as model:
            before = constrain(model).slim_optimize()
        with sink_model as model:
            after = constrain(model).slim_optimize()

        assert before == pytest.approx(expected, abs=1e-5)
        assert after == pytest.approx(before, rel=1e-9)


def test_growth_falls_monotonically_with_the_protein_carried(sink_model):
    frame = burden_curve(
        sink_model, [0.0, 0.01, 0.05, 0.1, 0.2, 0.4598], constrain=carbon_limited
    )

    assert (frame["status"] == "optimal").all()
    assert frame["growth"].is_monotonic_decreasing
    assert frame.iloc[-1]["relative_growth"] == pytest.approx(0.7063, abs=1e-3)


def test_the_slope_is_the_same_under_two_carbon_limited_regimes(sink_model):
    """Relative burden is regime-invariant while carbon is the binding constraint."""
    default = burden_curve(sink_model, [0.0, 0.05])
    limited = burden_curve(sink_model, [0.0, 0.05], constrain=carbon_limited)

    assert default.iloc[1]["relative_growth"] == pytest.approx(
        limited.iloc[1]["relative_growth"], rel=1e-5
    )
    assert (limited.iloc[1]["relative_growth"] - 1.0) / 0.05 == pytest.approx(-0.865, abs=0.01)


def test_the_slope_is_not_regime_invariant_once_oxygen_binds(sink_model):
    """Under REFERENCE_AEROBIC_BATCH the same gram costs a quarter less. State the regime."""
    limited = burden_curve(sink_model, [0.0, 0.05], constrain=carbon_limited)
    oxygen = burden_curve(sink_model, [0.0, 0.05], constrain=oxygen_limited)

    assert (oxygen.iloc[1]["relative_growth"] - 1.0) / 0.05 == pytest.approx(-0.650, abs=0.01)
    assert oxygen.iloc[1]["relative_growth"] > limited.iloc[1]["relative_growth"]


def test_the_gem_price_is_far_below_the_measured_burden_limit(sink_model, composition):
    """Eguchi 2018 puts the limit at 15% of total protein. The GEM charges 6.8%, not 15%."""
    load = grams_per_gdcw_at_protein_fraction(composition, EGUCHI_BURDEN_LIMIT)
    frame = burden_curve(sink_model, [0.0, load], constrain=carbon_limited)
    gem_loss = 1.0 - frame.iloc[1]["relative_growth"]

    assert load == pytest.approx(0.081134, abs=1e-5)
    assert gem_loss == pytest.approx(0.06832, abs=1e-4)
    assert gem_loss < EGUCHI_BURDEN_LIMIT
    assert gem_loss / EGUCHI_BURDEN_LIMIT == pytest.approx(0.455, abs=0.01)


def test_the_literature_table_reports_both_numbers_rather_than_reconciling_them(sink_model):
    """The comparison the module exists to make, and it is a disagreement."""
    frame = burden_against_literature(sink_model, constrain=carbon_limited).set_index("source")
    gem = frame.loc["this GEM: amino acids + translation ATP"]
    proportional = frame.loc[
        "proportional allocation, Eguchi 2018 eLife 7:e34595 limit"]
    growth_law = frame.loc["growth law, Metzl-Raz 2017 (PMID 28857745)"]
    kafri = frame.loc["Kafri 2016 (PMID 26725116)"]

    assert gem["relative_growth_loss"] == pytest.approx(0.06832, abs=1e-4)
    assert proportional["relative_growth_loss"] == pytest.approx(0.15)
    assert growth_law["relative_growth_loss"] == pytest.approx(0.17677, abs=1e-4)
    assert proportional["gem_over_reference"] == pytest.approx(0.4555, abs=1e-3)
    assert growth_law["gem_over_reference"] == pytest.approx(0.3865, abs=1e-3)
    # Kafri measured that no single coefficient exists, so the row carries no number.
    assert math.isnan(kafri["relative_growth_loss"])
    assert "low phosphate" in kafri["note"]


def test_the_growth_law_bound_collapses_to_the_naive_one_without_the_reserve():
    """Metzl-Raz's intercept is the whole difference between the two literature bounds."""
    naive = ribosome_allocation_penalty(0.15, 0.887685, reserve=0.0)
    measured = ribosome_allocation_penalty(0.15, 0.887685)
    slower = ribosome_allocation_penalty(0.15, 0.345984)

    assert naive == pytest.approx(0.15)
    assert measured == pytest.approx(0.17677, abs=1e-4)
    assert slower == pytest.approx(0.21869, abs=1e-4)
    assert METZL_RAZ_RIBOSOME_RESERVE == 0.08


def test_the_metzl_raz_slope_is_converted_off_the_papers_generations_per_hour_axis():
    """Fig. 2A's x-axis is generations/h, not mu, and the paper says so in the legend.

    Reading 0.35 as a per-mu slope overstates the growth-law bound by 1.2 points (18.9% vs
    17.7% at mu = 0.887685). The paper states the converted slope itself, as
    "Delta_r/Delta_mu = 21/ln(2) [min]", which is what this must reproduce.
    """
    assert METZL_RAZ_RIBOSOME_SLOPE_PER_GENERATION == 0.35
    assert METZL_RAZ_RIBOSOME_SLOPE == pytest.approx(21.0 / math.log(2.0) / 60.0, rel=1e-12)
    assert METZL_RAZ_RIBOSOME_SLOPE == pytest.approx(0.504943, abs=1e-6)
    # The uncorrected slope is what the 18.9% came from; it must no longer be reachable.
    assert ribosome_allocation_penalty(0.15, 0.887685, slope=0.35) == pytest.approx(0.18862, abs=1e-4)
    assert ribosome_allocation_penalty(0.15, 0.887685) == pytest.approx(0.17677, abs=1e-4)


def test_the_growth_law_bound_refuses_impossible_inputs():
    with pytest.raises(ValueError, match="fraction"):
        ribosome_allocation_penalty(1.0, 0.5)
    with pytest.raises(ValueError, match="growth rate"):
        ribosome_allocation_penalty(0.15, 0.0)


def test_the_burden_is_amino_acid_carbon_not_translation_atp(yeast_gem, composition):
    """Zeroing the translation ATP moves the penalty by 0.6 points out of 6.8; the full
    Stouthamer four moves it by 0.6 the other way."""
    load = grams_per_gdcw_at_protein_fraction(composition, EGUCHI_BURDEN_LIMIT)
    free_atp = add_protein_sink(yeast_gem, "hetprot", atp_per_residue=0.0)
    full_atp = add_protein_sink(
        yeast_gem, "hetprot4", atp_per_residue=ATP_EQUIVALENTS_PER_PEPTIDE_BOND)
    carbon_only = burden_curve(free_atp, [0.0, load], constrain=carbon_limited)
    everything = burden_curve(
        full_atp, [0.0, load], sink_reaction="PROTSYN_hetprot4", constrain=carbon_limited)

    assert 1.0 - carbon_only.iloc[1]["relative_growth"] == pytest.approx(0.0623, abs=1e-3)
    assert 1.0 - everything.iloc[1]["relative_growth"] == pytest.approx(0.0743, abs=1e-3)


def test_the_chaperone_sink_carries_the_regulon_no_model_has(yeast_gem, chaperone_model):
    """Zero chaperone genes exist in any yeast model; the GPR is the only handle Hsf1 gets."""
    rxn = chaperone_model.reactions.get_by_id("PROTSYN_chaperone")
    orfs = HSF1_CHAPERONE_ORFS + HAC1_CHAPERONE_ORFS

    assert {gene.id for gene in rxn.genes} == set(orfs)
    assert not any(yeast_gem.genes.has_id(orf) for orf in orfs)


def test_a_second_install_is_refused_by_name(sink_model):
    with pytest.raises(ValueError, match="already installed"):
        add_heterologous_protein_sink(sink_model)


def test_a_negative_translation_cost_is_refused(yeast_gem):
    with pytest.raises(ValueError, match="atp_per_residue"):
        add_protein_sink(yeast_gem, "bad", atp_per_residue=-1.0)


def test_a_protein_fraction_of_one_is_refused(composition):
    assert grams_per_gdcw_at_protein_fraction(composition, 0.0) == 0.0
    assert grams_per_gdcw_at_protein_fraction(composition, 0.5) == pytest.approx(
        composition.grams_per_unit
    )
    with pytest.raises(ValueError, match="fraction"):
        grams_per_gdcw_at_protein_fraction(composition, 1.0)


def test_a_burden_curve_without_an_unloaded_baseline_is_refused(sink_model, ec_yeast_gem):
    """Silently normalising to the lightest load reported zero burden AT a nonzero load."""
    with pytest.raises(ValueError, match="must include 0.0"):
        burden_curve(sink_model, [0.05, 0.10], constrain=carbon_limited)
    with pytest.raises(ValueError, match="must include 0.0"):
        protein_pool_burden_curve(ec_yeast_gem, [0.05, 0.15])


def test_a_model_without_a_protein_pseudoreaction_is_refused(yeast_gem):
    with pytest.raises(ValueError, match="charged tRNAs"):
        read_protein_composition(yeast_gem, "r_4048")


def test_couple_sink_to_growth_is_undone_by_the_model_context(sink_model):
    with sink_model as model:
        carbon_limited(model)
        model.objective = BIOMASS
        couple_sink_to_growth(model, "PROTSYN_hetprot", 0.4598)
        loaded = model.slim_optimize()
    with sink_model as model:
        released = carbon_limited(model).slim_optimize()

    assert loaded < released
    assert released == pytest.approx(0.887685, abs=1e-5)


# --- GECKO: the proteome constraint prices what the chemistry does not -------------------------


def test_taking_the_protein_out_of_the_enzyme_pool_costs_proportionally_by_construction(
    ec_yeast_gem,
):
    """15% of the proteome costs 15.1% of growth here and 6.8% in Yeast9 -- but see below.

    The 15.1% is NOT an independent reproduction of Eguchi's limit. Growth is exactly affine
    in the pool cap, so the loss is a fixed multiple of the fraction at every fraction; the
    next test pins that. The real content is the 2.2x gap against the Yeast9 chemistry.
    """
    frame = protein_pool_burden_curve(ec_yeast_gem, [0.0, 0.05, 0.15, 0.30])

    assert (frame["status"] == "optimal").all()
    assert frame.iloc[0]["growth"] == pytest.approx(0.376827, abs=1e-5)
    assert frame.iloc[0]["pool_cap"] == pytest.approx(0.10372, abs=1e-5)
    assert frame.iloc[2]["relative_growth_loss"] == pytest.approx(0.15077, abs=1e-4)
    assert frame.iloc[2]["relative_growth_loss"] / EGUCHI_BURDEN_LIMIT == pytest.approx(
        1.005, abs=0.01
    )
    assert frame.iloc[3]["relative_growth_loss"] == pytest.approx(0.30153, abs=1e-4)


def test_the_pool_price_is_proportional_at_every_fraction_so_it_predicts_nothing(ec_yeast_gem):
    """Why the agreement with a 15% proportional bound is arithmetic and not a result.

    loss/fraction is one constant, 1.005106, across two orders of magnitude of load, and
    growth is affine in the pool cap to machine precision. A model that returns 1.0051*f for
    every f has not measured f; it has restated the displacement assumption.
    """
    fractions = [0.01, 0.05, 0.15, 0.30, 0.50, 0.90]
    frame = protein_pool_burden_curve(ec_yeast_gem, [0.0] + fractions)
    ratios = (
        frame.loc[frame["fraction_of_total_protein"] > 0, "relative_growth_loss"]
        / frame.loc[frame["fraction_of_total_protein"] > 0, "fraction_of_total_protein"]
    )

    assert ratios.min() == pytest.approx(ratios.max(), rel=1e-9)
    assert ratios.iloc[0] == pytest.approx(1.005106, abs=1e-5)
    slope, intercept = numpy.polyfit(frame["pool_cap"], frame["growth"], 1)
    residual = numpy.abs(frame["growth"] - (slope * frame["pool_cap"] + intercept)).max()
    assert residual < 1e-12


def test_scaling_the_pool_leaves_the_original_model_alone(ec_yeast_gem):
    scaled = scale_protein_pool(ec_yeast_gem, 0.15)

    assert scaled.reactions.get_by_id("prot_pool_exchange").upper_bound == pytest.approx(
        0.0881620, abs=1e-6
    )
    assert ec_yeast_gem.reactions.get_by_id("prot_pool_exchange").upper_bound == pytest.approx(
        0.1037200, abs=1e-6
    )
    assert scaled.slim_optimize() == pytest.approx(0.320016, abs=1e-5)


def test_a_model_with_no_pool_cannot_be_charged_a_proteome_fraction(yeast_gem):
    """Yeast9 has no proteome budget, and saying so beats silently doing nothing."""
    with pytest.raises(ValueError, match="prot_pool_exchange"):
        scale_protein_pool(yeast_gem, 0.15)
    with pytest.raises(ValueError, match="fraction"):
        scale_protein_pool(yeast_gem, 1.0)


# --- GECKO: the installed carotenoid pathway carries no protein cost at all --------------------


@pytest.fixture(scope="module")
def _ec_carotenoid_template(ec_yeast_gem_factory):
    return add_beta_carotene_pathway(ec_yeast_gem_factory())


@pytest.fixture
def ec_carotenoid(_ec_carotenoid_template, model_copy):
    return model_copy(_ec_carotenoid_template)


@pytest.fixture(scope="module")
def _ec_charged_template(_ec_carotenoid_template, model_copy):
    """All four crt reactions charged at a PLACEHOLDER kcat of 1/s and UniProt masses."""
    return attach_crt_enzyme_demands(
        model_copy(_ec_carotenoid_template), {rid: 1.0 for rid in CRT_ENZYME_OF_REACTION},
        CRT_PLACEHOLDER
    )


@pytest.fixture
def ec_charged(_ec_charged_template, model_copy):
    return model_copy(_ec_charged_template)


def test_the_installed_carotenoid_reactions_are_free_in_gecko(ec_carotenoid):
    """The audit finding this module fixes."""
    for rid in CRT_KCAT_PROVENANCE:
        rxn = ec_carotenoid.reactions.get_by_id(rid)

        assert not any(met.id.startswith("prot_") for met in rxn.metabolites)


def test_no_kcat_is_shipped_for_any_crt_enzyme():
    """None is published, so the module refuses to carry one. The masses are not refused:
    they are UniProt sequence masses for the genes Verwaal 2007 expressed."""
    assert set(CRT_KCAT_PROVENANCE.values()) == {None}
    assert CRT_ENZYME_MASSES_KDA == {"crtE": 42.153, "crtYB": 74.736, "crtI": 65.091}


def test_attaching_an_enzyme_demand_couples_the_reaction_to_the_pool(ec_charged):
    rxn = ec_charged.reactions.get_by_id("CRTI")
    draw = ec_charged.reactions.get_by_id("draw_prot_crtI")
    pool = ec_charged.metabolites.get_by_id("prot_pool[c]")
    enzyme = ec_charged.metabolites.get_by_id("prot_crtI[c]")

    assert rxn.metabolites[enzyme] == pytest.approx(-1.0 / 3600.0)
    assert draw.metabolites[pool] == pytest.approx(-65.091)
    assert draw.notes["kcat_provenance"] == CRT_PLACEHOLDER
    # prot_* species have no formula, so the balance check sees the metabolic half only.
    assert mass_balance(rxn) == {}
    assert charge_balance(rxn) == pytest.approx(0.0, abs=1e-9)


def test_the_bifunctional_enzyme_is_drawn_once_and_not_twice(ec_charged):
    """crtYB runs both the synthase and the cyclase step, so it is one protein cost."""
    crtyb = ec_charged.metabolites.get_by_id("prot_crtYB[c]")
    users = {rxn.id for rxn in crtyb.reactions}

    assert users == {"draw_prot_crtYB", "CRTYB_PSY", "CRTYB_LCY"}
    assert len([r for r in ec_charged.reactions if r.id.startswith("draw_prot_crt")]) == 3


def test_the_enzyme_demand_changes_what_a_forced_product_flux_costs(ec_carotenoid, ec_charged):
    """Without it the pathway is free; with it, product flux buys proteome mass."""
    assert ec_charged.slim_optimize() == pytest.approx(ec_carotenoid.slim_optimize(), abs=1e-5)
    with ec_carotenoid as model:
        model.reactions.get_by_id("DM_betacarotene_c").bounds = (0.1, 0.1)
        free = model.slim_optimize()
    with ec_charged as model:
        model.reactions.get_by_id("DM_betacarotene_c").bounds = (0.1, 0.1)
        priced = model.slim_optimize()

    assert free == pytest.approx(0.197806, abs=1e-5)
    assert priced == pytest.approx(0.173666, abs=1e-5)
    assert 1.0 - priced / free == pytest.approx(0.1220, abs=1e-3)


def test_the_charged_pathway_carries_flux_and_the_draws_carry_it_with_it(ec_charged):
    with ec_charged as model:
        model.reactions.get_by_id("DM_betacarotene_c").bounds = (0.1, 0.1)
        fva = flux_variability_analysis(
            model,
            ["draw_prot_crtI", "draw_prot_crtYB", "draw_prot_crtE"],
            fraction_of_optimum=0.0,
            processes=1,
        )

    assert fva.loc["draw_prot_crtI", "minimum"] == pytest.approx(2.7778e-5, rel=1e-3)
    # crtYB runs two steps at 0.1 each, so it draws twice what crtI does.
    assert fva.loc["draw_prot_crtYB", "minimum"] == pytest.approx(5.5556e-5, rel=1e-3)
    # crtE is optional: native BTS1 makes the same GGPP, so its floor stays at zero.
    assert fva.loc["draw_prot_crtE", "minimum"] == pytest.approx(0.0, abs=1e-9)
    assert fva.loc["draw_prot_crtE", "maximum"] > 0.0


def test_installing_the_enzyme_demand_does_not_move_growth_on_its_own(ec_yeast_gem, ec_charged):
    """An unused enzyme draws nothing; the cost only appears when product flux is forced."""
    assert ec_charged.slim_optimize() == pytest.approx(ec_yeast_gem.slim_optimize(), abs=1e-5)


def test_an_enzyme_demand_without_provenance_is_refused(ec_carotenoid):
    with pytest.raises(ValueError, match="provenance"):
        attach_enzyme_demand(
            ec_carotenoid, "CRTI", "crtI", kcat_per_s=1.0, molecular_weight_kda=65.091,
            provenance="",
        )


def test_a_citation_shaped_provenance_is_refused_for_an_enzyme_with_no_published_kcat(
    ec_carotenoid,
):
    """The guard that matters: no kcat exists for these four, so a citation would be false."""
    with pytest.raises(ValueError, match="no kcat is published"):
        attach_enzyme_demand(
            ec_carotenoid, "CRTI", "crtI", kcat_per_s=1.0, molecular_weight_kda=65.091,
            provenance="Verwaal 2007",
        )


def test_a_non_carotenoid_reaction_is_refused_by_the_crt_helper(ec_carotenoid):
    with pytest.raises(ValueError, match="not carotenoid reactions"):
        attach_crt_enzyme_demands(ec_carotenoid, {"r_0005": 1.0}, CRT_PLACEHOLDER)


def test_a_model_with_no_protein_pool_is_refused(yeast_gem):
    """Yeast9 has no proteome to charge, and saying so beats silently doing nothing."""
    carotenoid = add_beta_carotene_pathway(yeast_gem)
    with pytest.raises(ValueError, match="prot_pool"):
        attach_enzyme_demand(
            carotenoid, "CRTI", "crtI", kcat_per_s=1.0, molecular_weight_kda=65.091,
            provenance="PLACEHOLDER",
        )


def test_charging_the_same_reaction_twice_is_refused(ec_charged):
    with pytest.raises(ValueError, match="already carries"):
        attach_enzyme_demand(
            ec_charged, "CRTE", "crtE2", kcat_per_s=1.0, molecular_weight_kda=42.153,
            provenance="PLACEHOLDER",
        )


def test_reusing_an_enzyme_id_with_a_different_mass_is_refused(ec_charged):
    """Reuse is what makes crtYB one protein; a disagreeing mass would make it two."""
    with pytest.raises(ValueError, match="already draws"):
        attach_enzyme_demand(
            ec_charged, "DM_betacarotene_c", "crtYB", kcat_per_s=1.0,
            molecular_weight_kda=40.0, provenance="PLACEHOLDER",
        )


def test_the_sink_installs_on_the_enzyme_constrained_model_too(ec_yeast_gem):
    model = add_heterologous_protein_sink(ec_yeast_gem)
    rxn = model.reactions.get_by_id("PROTSYN_hetprot")
    fva = flux_variability_analysis(
        model, ["PROTSYN_hetprot"], fraction_of_optimum=0.0, processes=1)

    assert mass_balance(rxn) == {}
    assert charge_balance(rxn) == pytest.approx(0.0, abs=1e-9)
    assert fva.loc["PROTSYN_hetprot", "maximum"] == pytest.approx(0.49253, abs=1e-4)
    assert model.slim_optimize() == pytest.approx(ec_yeast_gem.slim_optimize(), abs=1e-5)
