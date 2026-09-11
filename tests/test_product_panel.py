from __future__ import annotations

import ast
import csv
import inspect
import io
import math
from dataclasses import FrozenInstanceError

import cobra
import pandas as pd
import pytest

from ystwin import paths
from ystwin.fba import product_panel
from ystwin.fba.physiology import cap_uptake
from ystwin.fba.solver import load_model
from ystwin.pathway.proteome import enzyme_content_mmol_per_gdcw


MODEL_ID = "M_ecYeastGEM_batch_v8__46__3__46__4"
PRODUCTS = (
    ("squalene", "s_1447[c]", "C30H50", "MNXM292", "c", 30),
    ("glycogen", "s_0773[c]", "C6H12O6", "MNXM55375", "c", 6),
    ("trehalose", "s_1520[c]", "C12H22O11", "MNXM198", "c", 12),
    ("glutathione", "s_0750[c]", "C10H16N3O6S", "MNXM57", "c", 10),
    ("glycerol", "s_0766[e]", "C3H8O3", "MNXM89612", "e", 3),
)
ENZYMES = (
    ("ERG9", "YHR190W", "P29704", 165.0),
    ("GSY1", "YFR015C", "P23337", 30.0),
    ("GSY2", "YLR258W", "P27472", 86.5),
    ("GLG1", "YKR058W", "P36143", 2.27),
    ("GLG2", "YJL137C", "P47011", 5.13),
    ("TPS1", "YBR126C", "Q00764", 280.0),
    ("TPS2", "YDR074W", "P31688", 191.0),
    ("TPS3", "YMR261C", "P38426", 105.0),
    ("TSL1", "YML100W", "P38427", 130.0),
    ("GSH1", "YJL101C", "P32477", 31.1),
    ("GSH2", "YOL049W", "Q08220", 78.0),
    ("GPD1", "YDL022W", "Q00055", 436.0),
    ("GPD2", "YOL059W", "P41911", 230.0),
    ("GPP1", "YIL053W", "P41277", 2919.0),
    ("GPP2", "YER062C", "P40106", 461.0),
)


class MemoryTable:
    def __init__(self, rows):
        self.rows = rows

    def open(self, *args, **kwargs):
        stream = io.StringIO()
        writer = csv.DictWriter(
            stream, fieldnames=("gene_name", "string_id", "abundance_ppm"), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(self.rows)
        stream.seek(0)
        return stream


@pytest.fixture
def proteomics(monkeypatch):
    rows = [
        {"gene_name": gene, "string_id": f"4932.{orf}", "abundance_ppm": ppm}
        for gene, orf, _, ppm in ENZYMES
    ]
    monkeypatch.setattr(product_panel, "_PAXDB_PATH", MemoryTable(rows))
    return rows


@pytest.fixture
def panel_model():
    model = cobra.Model(MODEL_ID)
    model.solver = "glpk"
    pool = cobra.Metabolite("prot_pool[c]", compartment="c")
    model.add_metabolites([pool])
    supply = cobra.Reaction("prot_pool_exchange", upper_bound=1.0)
    supply.add_metabolites({pool: 1.0})
    model.add_reactions([supply])
    for _, orf, protein_id, _ in ENZYMES:
        protein = cobra.Metabolite(f"prot_{protein_id}[c]", compartment="c")
        draw = cobra.Reaction(f"draw_prot_{protein_id}", upper_bound=float("inf"))
        draw.add_metabolites({pool: -50.0, protein: 1.0})
        draw.gene_reaction_rule = orf
        consumption = cobra.Reaction(f"turnover_{orf}", upper_bound=float("inf"))
        consumption.add_metabolites({protein: -0.001})
        consumption.gene_reaction_rule = orf
        model.add_reactions([draw, consumption])
    shared = cobra.Reaction("second_TPS2_task", upper_bound=float("inf"))
    shared.add_metabolites({model.metabolites.get_by_id("prot_P31688[c]"): -0.001})
    shared.gene_reaction_rule = "YDR074W"
    model.add_reactions([shared])
    for _, mid, formula, identity, compartment, _ in PRODUCTS:
        metabolite = cobra.Metabolite(mid, formula=formula, charge=0, compartment=compartment)
        metabolite.annotation["metanetx.chemical"] = identity
        model.add_metabolites([metabolite])
    glycerol = cobra.Metabolite("s_0765[c]", formula="C3H8O3", charge=0, compartment="c")
    glycerol.annotation["metanetx.chemical"] = "MNXM89612"
    model.add_metabolites([glycerol])
    export = cobra.Reaction("r_1808", upper_bound=float("inf"))
    export.add_metabolites({model.metabolites.get_by_id("s_0766[e]"): -1.0})
    uptake = cobra.Reaction("r_1808_REV", upper_bound=0.0)
    uptake.add_metabolites({model.metabolites.get_by_id("s_0766[e]"): 1.0})
    native_turnover = cobra.Reaction("native_storage_turnover", upper_bound=3.0)
    native_turnover.add_metabolites({model.metabolites.get_by_id("s_0773[c]"): -1.0})
    model.add_reactions([export, uptake, native_turnover])
    model.objective = "turnover_YHR190W"
    return model


def snapshot(model):
    return (
        model.id,
        str(model.objective.expression),
        tuple(
            (r.id, r.bounds, r.gene_reaction_rule, tuple(sorted((m.id, c) for m, c in r.metabolites.items())))
            for r in model.reactions
        ),
        tuple((m.id, m.formula, m.compartment) for m in model.metabolites),
    )


def test_relative_abundance_changes_a_binding_enzyme_budget(panel_model, proteomics):
    prepared, inputs = product_panel.prepare_panel_model(panel_model)
    prepared.objective = "turnover_YDL022W"
    before = snapshot(prepared)
    baseline = prepared.slim_optimize(error_value=None)
    changed, applied = product_panel.apply_relative_abundances(
        prepared, inputs, {"YDL022W": 2.0})
    assert snapshot(prepared) == before
    assert changed.slim_optimize(error_value=None) == pytest.approx(2.0 * baseline)
    assert applied.iloc[0]["relative_abundance"] == 2.0
    assert changed.reactions.prot_pool_exchange.bounds == prepared.reactions.prot_pool_exchange.bounds
    assert changed.reactions.draw_prot_P41911.bounds == prepared.reactions.draw_prot_P41911.bounds


def test_relative_abundance_preserves_independent_hard_caps(panel_model, proteomics):
    content = enzyme_content_mmol_per_gdcw(436.0)
    panel_model.reactions.draw_prot_Q00055.upper_bound = 0.5 * content
    prepared, inputs = product_panel.prepare_panel_model(panel_model)
    changed, _ = product_panel.apply_relative_abundances(prepared, inputs, {"YDL022W": 0.75})
    assert changed.reactions.draw_prot_Q00055.upper_bound == pytest.approx(0.5 * content)
    changed, _ = product_panel.apply_relative_abundances(prepared, inputs, {"YDL022W": 0.25})
    assert changed.reactions.draw_prot_Q00055.upper_bound == pytest.approx(0.25 * content)


def test_relative_abundance_rejects_wrong_identity_and_ambiguous_reference(panel_model, proteomics):
    prepared, inputs = product_panel.prepare_panel_model(panel_model)
    bad = inputs.copy()
    bad.loc[bad.orf == "YDL022W", "draw_reaction_id"] = "draw_prot_P41911"
    with pytest.raises(ValueError, match="identity"):
        product_panel.apply_relative_abundances(prepared, bad, {"YDL022W": 2})
    prepared.reactions.draw_prot_Q00055.upper_bound *= 0.5
    with pytest.raises(ValueError, match="reference"):
        product_panel.apply_relative_abundances(prepared, inputs, {"YDL022W": 2})


def test_intracellular_glycerol_is_not_reported_as_secreted_product(panel_model):
    prepared, task = product_panel.install_product(panel_model, "glycerol", output_compartment="c")
    assert task.metabolite_id == "s_0765[c]"
    assert task.reaction_id != "r_1808"
    assert "intracellular" in task.output_basis
    assert prepared.reactions.get_by_id(task.reaction_id).metabolites == {
        prepared.metabolites.get_by_id("s_0765[c]"): -1.0}
    assert prepared.reactions.r_1808.bounds == panel_model.reactions.r_1808.bounds
    assert prepared.reactions.r_1808_REV.bounds == (0.0, 0.0)
    with pytest.raises(ValueError, match="compartment"):
        product_panel.install_product(panel_model, "squalene", output_compartment="e")


def test_product_names_are_frozen_in_the_requested_order():
    assert product_panel.product_names() == tuple(row[0] for row in PRODUCTS)
    assert isinstance(product_panel.product_names(), tuple)


@pytest.mark.parametrize("name,mid,formula,identity,compartment,carbons", PRODUCTS)
def test_product_task_uses_validated_species_and_net_product_sink(
    panel_model, name, mid, formula, identity, compartment, carbons
):
    before = snapshot(panel_model)
    prepared, task = product_panel.install_product(panel_model, name)
    assert snapshot(panel_model) == before
    assert prepared is not panel_model
    assert isinstance(task, product_panel.ProductTask)
    assert task.name == name
    assert task.metabolite_id == mid
    assert task.carbon_atoms == carbons
    assert task.molar_mass_g_per_mol == pytest.approx(prepared.metabolites.get_by_id(mid).formula_weight)
    assert math.isfinite(task.molar_mass_g_per_mol) and task.molar_mass_g_per_mol > 0
    reaction = prepared.reactions.get_by_id(task.reaction_id)
    assert reaction.boundary
    assert {m.id: c for m, c in reaction.metabolites.items()} == {mid: -1.0}
    assert reaction.lower_bound == 0.0
    assert isinstance(task.assumptions, tuple) and task.assumptions
    assert any(formula in assumption for assumption in task.assumptions)
    if name == "glycerol":
        assert task.reaction_id == "r_1808"
        assert "extracellular" in task.output_basis
        assert len(prepared.reactions) == len(panel_model.reactions)
    else:
        assert "extra intracellular net accumulation" in task.output_basis
        assert "biomass" in task.output_basis
        assert len(prepared.reactions) == len(panel_model.reactions) + 1
    with pytest.raises(FrozenInstanceError):
        task.name = "other"


def test_glycogen_mass_is_the_frozen_model_glucose_equivalent_repeat_unit(panel_model):
    _, task = product_panel.install_product(panel_model, "glycogen")
    text = " ".join((task.output_basis, *task.assumptions)).lower()
    assert "repeat unit" in text
    assert "glucose equivalent" in text
    assert task.molar_mass_g_per_mol == pytest.approx(180.15588)
    assert "starting" in text and "pool" in text


def test_installation_retains_competition_turnover_objective_and_existing_bounds(panel_model):
    prepared, task = product_panel.install_product(panel_model, "glycogen")
    assert prepared.reactions.get_by_id("native_storage_turnover").bounds == (0.0, 3.0)
    assert str(prepared.objective.expression) == str(panel_model.objective.expression)
    prepared.reactions.get_by_id("native_storage_turnover").upper_bound = 1.0
    prepared.metabolites.get_by_id(task.metabolite_id).name = "changed copy"
    assert panel_model.reactions.get_by_id("native_storage_turnover").upper_bound == 3.0
    assert panel_model.metabolites.get_by_id(task.metabolite_id).name != "changed copy"


def test_installation_is_idempotent_but_rejects_a_colliding_sink(panel_model):
    first, task = product_panel.install_product(panel_model, "squalene")
    second, repeated = product_panel.install_product(first, "squalene")
    assert second is not first
    assert repeated == task
    assert len(second.reactions) == len(first.reactions)
    sink = first.reactions.get_by_id(task.reaction_id)
    sink.add_metabolites({first.metabolites.get_by_id("s_1520[c]"): -1.0})
    with pytest.raises(ValueError, match="demand|sink|stoichiometry"):
        product_panel.install_product(first, "squalene")


def test_installation_rejects_unsupported_names(panel_model):
    with pytest.raises(ValueError, match="unsupported product.*gadusol"):
        product_panel.install_product(panel_model, "gadusol")


@pytest.mark.parametrize("change", ("missing", "formula", "identity", "compartment"))
def test_installation_refuses_missing_or_wrong_species(panel_model, change):
    metabolite = panel_model.metabolites.get_by_id("s_1447[c]")
    if change == "missing":
        panel_model.remove_metabolites([metabolite])
    elif change == "formula":
        metabolite.formula = "C40H56"
    elif change == "identity":
        metabolite.annotation["metanetx.chemical"] = "MNXM198"
    else:
        metabolite.compartment = "m"
    with pytest.raises((ValueError, KeyError), match="s_1447|squalene"):
        product_panel.install_product(panel_model, "squalene")


def test_species_identity_does_not_depend_on_names(panel_model):
    panel_model.metabolites.get_by_id("s_1447[c]").name = "trehalose misleading label"
    _, task = product_panel.install_product(panel_model, "squalene")
    assert task.metabolite_id == "s_1447[c]"


def test_glycerol_refuses_import_reexport_as_net_secretion(panel_model):
    panel_model.reactions.get_by_id("r_1808_REV").upper_bound = 1.0
    with pytest.raises(ValueError, match="r_1808_REV|reuptake|net"):
        product_panel.install_product(panel_model, "glycerol")


def test_glycerol_refuses_an_exchange_of_the_wrong_species(panel_model):
    exchange = panel_model.reactions.get_by_id("r_1808")
    exchange.add_metabolites({panel_model.metabolites.get_by_id("s_0766[e]"): 1.0})
    exchange.add_metabolites({panel_model.metabolites.get_by_id("s_1447[c]"): -1.0})
    with pytest.raises(ValueError, match="r_1808|glycerol"):
        product_panel.install_product(panel_model, "glycerol")


def test_preparation_caps_the_same_union_without_adding_enzyme_copies(panel_model, proteomics):
    before = snapshot(panel_model)
    prepared, provenance = product_panel.prepare_panel_model(panel_model)
    assert snapshot(panel_model) == before
    assert prepared is not panel_model
    assert isinstance(provenance, pd.DataFrame)
    assert set(provenance.orf) == {orf for _, orf, _, _ in ENZYMES}
    assert len(provenance) == 15
    assert provenance.draw_reaction_id.is_unique
    assert provenance.protein_id.is_unique
    assert len(prepared.reactions) == len(panel_model.reactions)
    assert len(prepared.metabolites) == len(panel_model.metabolites)
    for gene, orf, protein_id, ppm in ENZYMES:
        row = provenance.set_index("orf").loc[orf]
        assert row.gene == gene
        assert row.string_id == f"4932.{orf}"
        assert row.protein_id == f"prot_{protein_id}[c]"
        assert row.draw_reaction_id == f"draw_prot_{protein_id}"
        assert row.abundance_ppm == ppm
        assert row.enzyme_mmol_per_gdcw == pytest.approx(enzyme_content_mmol_per_gdcw(ppm))
        assert math.isinf(row.upper_bound_before)
        assert row.upper_bound_after == pytest.approx(row.enzyme_mmol_per_gdcw)
        assert prepared.reactions.get_by_id(row.draw_reaction_id).upper_bound == row.upper_bound_after
        assert "PaxDb" in row.source
        assert "model-derived" in row.turnover_source
        assert "not necessarily measured" in row.turnover_source
        assert isinstance(row.assumptions, tuple)
    for name in product_panel.product_names():
        target, _ = product_panel.install_product(prepared, name)
        assert {
            rid: target.reactions.get_by_id(rid).upper_bound for rid in provenance.draw_reaction_id
        } == {
            rid: prepared.reactions.get_by_id(rid).upper_bound for rid in provenance.draw_reaction_id
        }


def test_shared_enzyme_consumers_retain_a_single_total_budget(panel_model, proteomics):
    prepared, provenance = product_panel.prepare_panel_model(panel_model)
    protein = prepared.metabolites.get_by_id("prot_P31688[c]")
    assert sum(r.metabolites[protein] > 0 for r in protein.reactions) == 1
    prepared.objective = {
        prepared.reactions.get_by_id("turnover_YDR074W"): 1.0,
        prepared.reactions.get_by_id("second_TPS2_task"): 1.0,
    }
    solution = prepared.optimize()
    assert solution.status == "optimal"
    cap = provenance.set_index("orf").loc["YDR074W", "upper_bound_after"]
    assert solution.objective_value == pytest.approx(cap / 0.001)
    assert prepared.reactions.get_by_id("turnover_YDR074W").metabolites[protein] == -0.001
    assert prepared.reactions.get_by_id("second_TPS2_task").metabolites[protein] == -0.001


def test_mapping_uses_string_ids_and_draw_gprs_not_gene_or_protein_names(panel_model, proteomics):
    proteomics[0]["gene_name"] = "source_label_only"
    panel_model.genes.get_by_id("YHR190W").name = "misleading_gene_name"
    panel_model.metabolites.get_by_id("prot_P29704[c]").id = "structured_enzyme_identity"
    panel_model.reactions.get_by_id("draw_prot_P29704").id = "structured_draw_identity"
    _, provenance = product_panel.prepare_panel_model(panel_model)
    row = provenance.set_index("orf").loc["YHR190W"]
    assert row.gene == "source_label_only"
    assert row.protein_id == "structured_enzyme_identity"
    assert row.draw_reaction_id == "structured_draw_identity"


@pytest.mark.parametrize("multiplier", (0.0, 0.5, 2.0))
def test_multiplier_is_applied_once_and_recorded(panel_model, proteomics, multiplier):
    prepared, provenance = product_panel.prepare_panel_model(
        panel_model, abundance_multiplier=multiplier
    )
    assert (provenance.abundance_multiplier == multiplier).all()
    for row in provenance.itertuples():
        assert row.scaled_enzyme_mmol_per_gdcw == pytest.approx(row.enzyme_mmol_per_gdcw * multiplier)
        assert row.upper_bound_after == pytest.approx(row.scaled_enzyme_mmol_per_gdcw)
        assert prepared.reactions.get_by_id(row.draw_reaction_id).upper_bound == row.upper_bound_after
        assert any("multiplier" in text for text in row.assumptions)


@pytest.mark.parametrize("multiplier", (-1.0, float("nan"), float("inf"), -float("inf")))
def test_invalid_multiplier_is_rejected_without_mutation(panel_model, proteomics, multiplier):
    before = snapshot(panel_model)
    with pytest.raises(ValueError, match="abundance_multiplier"):
        product_panel.prepare_panel_model(panel_model, abundance_multiplier=multiplier)
    assert snapshot(panel_model) == before


def test_preexisting_tighter_caps_are_not_relaxed(panel_model, proteomics):
    panel_model.reactions.get_by_id("draw_prot_P29704").upper_bound = 1e-8
    prepared, provenance = product_panel.prepare_panel_model(panel_model, abundance_multiplier=2.0)
    row = provenance.set_index("orf").loc["YHR190W"]
    assert row.upper_bound_before == 1e-8
    assert row.upper_bound_after == 1e-8
    assert prepared.reactions.get_by_id(row.draw_reaction_id).upper_bound == 1e-8


def test_disabled_proteomics_preserves_bounds_and_records_the_ablation(panel_model, proteomics):
    prepared, provenance = product_panel.prepare_panel_model(
        panel_model, use_proteomics=False, abundance_multiplier=2.0
    )
    assert snapshot(prepared) == snapshot(panel_model)
    assert prepared is not panel_model
    assert not provenance.use_proteomics.any()
    assert (provenance.upper_bound_before == provenance.upper_bound_after).all()
    assert all(any("not applied" in text for text in assumptions) for assumptions in provenance.assumptions)


def test_prepared_model_and_provenance_do_not_mutate_the_source(panel_model, proteomics):
    prepared, provenance = product_panel.prepare_panel_model(panel_model)
    prepared.reactions.get_by_id("draw_prot_P29704").upper_bound = 0.0
    provenance.loc[:, "abundance_ppm"] = 0.0
    assert math.isinf(panel_model.reactions.get_by_id("draw_prot_P29704").upper_bound)
    _, fresh = product_panel.prepare_panel_model(panel_model)
    assert fresh.set_index("orf").loc["YHR190W", "abundance_ppm"] == 165.0


@pytest.mark.parametrize("failure", ("missing", "ambiguous", "complex", "unmapped_protein"))
def test_missing_or_ambiguous_draw_mapping_is_refused(panel_model, proteomics, failure):
    draw = panel_model.reactions.get_by_id("draw_prot_P29704")
    if failure == "missing":
        draw.gene_reaction_rule = ""
    elif failure == "ambiguous":
        duplicate = draw.copy()
        duplicate.id = "another_ERG9_draw"
        panel_model.add_reactions([duplicate])
    elif failure == "complex":
        draw.gene_reaction_rule = "YHR190W or YOL049W"
    else:
        draw.add_metabolites({panel_model.metabolites.get_by_id("prot_P29704[c]"): -1.0})
    before = snapshot(panel_model)
    with pytest.raises(ValueError, match="YHR190W|draw|protein"):
        product_panel.prepare_panel_model(panel_model)
    assert snapshot(panel_model) == before


@pytest.mark.parametrize("failure", ("missing", "wrong_taxon", "duplicate", "nan", "negative"))
def test_invalid_selected_proteomics_is_refused(panel_model, proteomics, failure):
    if failure == "missing":
        proteomics.pop(0)
    elif failure == "wrong_taxon":
        proteomics[0]["string_id"] = "9606.YHR190W"
    elif failure == "duplicate":
        proteomics.append(dict(proteomics[0]))
    elif failure == "nan":
        proteomics[0]["abundance_ppm"] = "nan"
    else:
        proteomics[0]["abundance_ppm"] = -1.0
    with pytest.raises(ValueError, match="PaxDb|YHR190W|abundance"):
        product_panel.prepare_panel_model(panel_model)


def test_panel_refuses_a_different_model_version(panel_model, proteomics):
    panel_model.id = "yeastGEM_v9__46__0__46__2"
    with pytest.raises(ValueError, match="8.3.4|frozen|model"):
        product_panel.prepare_panel_model(panel_model)
    with pytest.raises(ValueError, match="8.3.4|frozen|model"):
        product_panel.install_product(panel_model, "squalene")
    with pytest.raises(ValueError, match="8.3.4|frozen|model"):
        product_panel.biomass_pool_quota(panel_model, "glycogen")


def test_product_panel_imports_only_model_and_proteomics_inputs():
    tree = ast.parse(inspect.getsource(product_panel))
    allowed = {"__future__", "ast", "csv", "math", "dataclasses", "cobra", "pandas", "pathway.proteome", ""}
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
    assert set(imports) <= allowed
    assert not any("calibration" in name or "stress_pools" in name or "stress_biomass" in name for name in imports)


@pytest.fixture
def biomass_model(panel_model):
    panel_model.add_metabolites([
        cobra.Metabolite(mid, compartment="c") for mid in ("s_3718[c]", "s_0450[c]")
    ])
    for rid, stoichiometry in (
        ("r_4048", {"s_0773[c]": -0.6, "s_1520[c]": -0.4, "s_3718[c]": 2.0}),
        ("r_4041", {"s_3718[c]": -3.0, "s_0450[c]": 4.0}),
        ("r_2111", {"s_0450[c]": -1.0}),
        ("r_2111_REV", {"s_0450[c]": 1.0}),
    ):
        reaction = cobra.Reaction(rid, upper_bound=0.0 if rid == "r_2111_REV" else 1000.0)
        reaction.add_metabolites({
            panel_model.metabolites.get_by_id(mid): coefficient
            for mid, coefficient in stoichiometry.items()
        })
        panel_model.add_reactions([reaction])
    return panel_model


@pytest.mark.parametrize("product,expected", (("glycogen", 0.225), ("trehalose", 0.15)))
def test_biomass_pool_quota_uses_composition_and_normalizes_both_links(
    biomass_model, monkeypatch, product, expected
):
    monkeypatch.setattr(biomass_model, "optimize", None)
    monkeypatch.setattr(biomass_model, "slim_optimize", None)
    before = snapshot(biomass_model)
    assert "biomass_pool_quota" in product_panel.__all__
    assert product_panel.biomass_pool_quota(biomass_model, product) == pytest.approx(expected)
    assert snapshot(biomass_model) == before


@pytest.mark.parametrize("growth", (0.0, 0.05, 0.1))
@pytest.mark.parametrize("extra", (0.0, 0.2))
def test_biomass_pool_quota_ignores_growth_objective_and_extra_demand(biomass_model, growth, extra):
    model, task = product_panel.install_product(biomass_model, "glycogen")
    model.reactions.get_by_id("r_2111").bounds = (growth, growth)
    model.reactions.get_by_id(task.reaction_id).bounds = (extra, extra)
    model.objective = task.reaction_id
    before = snapshot(model)
    assert product_panel.biomass_pool_quota(model, "glycogen") == pytest.approx(0.225)
    assert snapshot(model) == before
    assert task.reaction_id not in biomass_model.reactions


@pytest.mark.parametrize("rid", ("r_4048", "r_4041", "r_2111", "r_2111_REV"))
def test_biomass_pool_quota_refuses_missing_assembly(biomass_model, rid):
    biomass_model.remove_reactions([rid])
    before = snapshot(biomass_model)
    with pytest.raises(ValueError, match="quota unavailable"):
        product_panel.biomass_pool_quota(biomass_model, "glycogen")
    assert snapshot(biomass_model) == before


@pytest.mark.parametrize("mid", ("s_3718[c]", "s_0450[c]"))
@pytest.mark.parametrize("coefficient", (-1.0, 1.0))
def test_biomass_pool_quota_refuses_ambiguous_pseudometabolite_links(biomass_model, mid, coefficient):
    alternative = cobra.Reaction("alternative_biomass_link")
    alternative.add_metabolites({biomass_model.metabolites.get_by_id(mid): coefficient})
    biomass_model.add_reactions([alternative])
    before = snapshot(biomass_model)
    with pytest.raises(ValueError, match="quota unavailable.*ambiguous"):
        product_panel.biomass_pool_quota(biomass_model, "glycogen")
    assert snapshot(biomass_model) == before


@pytest.mark.parametrize("failure", ("missing_pool", "produced_pool", "growth_import", "direct_pool"))
def test_biomass_pool_quota_does_not_assume_zero_or_ignore_altered_assembly(biomass_model, failure):
    glycogen = biomass_model.metabolites.get_by_id("s_0773[c]")
    if failure in ("missing_pool", "produced_pool"):
        biomass_model.reactions.get_by_id("r_4048").add_metabolites({
            glycogen: 0.6 if failure == "missing_pool" else 1.2
        })
    elif failure == "growth_import":
        biomass_model.reactions.get_by_id("r_2111_REV").upper_bound = 1.0
    else:
        biomass_model.reactions.get_by_id("r_4041").add_metabolites({glycogen: -0.1})
    before = snapshot(biomass_model)
    with pytest.raises(ValueError, match="quota unavailable"):
        product_panel.biomass_pool_quota(biomass_model, "glycogen")
    assert snapshot(biomass_model) == before


@pytest.mark.parametrize("product", ("gadusol", "squalene", "glutathione", "glycerol"))
def test_biomass_pool_quota_refuses_unsupported_products_instead_of_assuming_zero(biomass_model, product):
    with pytest.raises(ValueError, match="quota unavailable"):
        product_panel.biomass_pool_quota(biomass_model, product)


@pytest.fixture(scope="module")
def vendored_panel():
    model_path = paths.data_dir() / "gem" / "ecYeastGEM_batch.xml.gz"
    proteomics_path = paths.data_dir() / "proteome" / "paxdb_scerevisiae_integrated.tsv"
    if not model_path.is_file() or not proteomics_path.is_file():
        pytest.skip("frozen EC batch model or PaxDb table not vendored")
    model, _ = load_model(model_path)
    prepared, provenance = product_panel.prepare_panel_model(model)
    return model, prepared, provenance


@pytest.mark.integration
@pytest.mark.parametrize("name,mid,formula,identity,compartment,carbons", PRODUCTS)
def test_five_actual_model_products_have_supported_positive_net_flux(
    vendored_panel, name, mid, formula, identity, compartment, carbons
):
    original, prepared, provenance = vendored_panel
    model, task = product_panel.install_product(prepared, name)
    metabolite = model.metabolites.get_by_id(mid)
    assert metabolite.formula == formula
    assert metabolite.annotation["metanetx.chemical"] == identity
    assert task.molar_mass_g_per_mol == pytest.approx(metabolite.formula_weight)
    assert task.carbon_atoms == carbons
    assert len(provenance) == 15
    cap_uptake(model, "r_1714", 20.0 / 6.0)
    cap_uptake(model, "r_1992", 1000.0)
    model.reactions.get_by_id("r_2111").bounds = (0.05, 0.05)
    model.objective = task.reaction_id
    solution = model.optimize()
    assert solution.status == "optimal"
    assert solution.fluxes[task.reaction_id] > 1e-6
    assert original.reactions.get_by_id("r_1714_REV").upper_bound == 1000.0
    assert math.isinf(original.reactions.get_by_id("draw_prot_P29704").upper_bound)
    if name != "glycerol":
        assert task.reaction_id not in original.reactions
        assert task.reaction_id not in prepared.reactions


@pytest.mark.integration
def test_actual_tps_complexes_reuse_the_original_protein_budgets(vendored_panel):
    original, prepared, provenance = vendored_panel
    for rid in ("r_0195No1", "r_0195No2", "r_1051No1", "r_1051No2"):
        expected = {m.id: c for m, c in original.reactions.get_by_id(rid).metabolites.items()}
        assert {m.id: c for m, c in prepared.reactions.get_by_id(rid).metabolites.items()} == expected
    tps2 = provenance.set_index("orf").loc["YDR074W"]
    assert tps2.protein_id == "prot_P31688[c]"
    assert tps2.draw_reaction_id == "draw_prot_P31688"
    protein = prepared.metabolites.get_by_id(tps2.protein_id)
    assert {r.id for r in protein.reactions if r.metabolites[protein] > 0} == {tps2.draw_reaction_id}


@pytest.mark.integration
@pytest.mark.parametrize("product,expected,mg_per_gdcw", (
    ("glycogen", 0.330522, 59.54548176936),
    ("trehalose", 0.126456, 43.28544367488),
))
def test_actual_biomass_pool_quota_is_native_content_not_extra_accumulation(
    vendored_panel, product, expected, mg_per_gdcw
):
    original, prepared, _ = vendored_panel
    original_before, prepared_before = snapshot(original), snapshot(prepared)
    assert product_panel.biomass_pool_quota(original, product) == pytest.approx(expected)
    assert product_panel.biomass_pool_quota(prepared, product) == pytest.approx(expected)
    model, task = product_panel.install_product(prepared, product)
    assert expected * task.molar_mass_g_per_mol == pytest.approx(mg_per_gdcw)
    model.objective = task.reaction_id
    for glucose, oxygen, growth in ((1.5, 5.0, 0.025), (20.0 / 6.0, 1000.0, 0.05)):
        cap_uptake(model, "r_1714", glucose)
        cap_uptake(model, "r_1992", oxygen)
        model.reactions.get_by_id("r_2111").bounds = (growth, growth)
        for extra in (0.0, 0.005):
            model.reactions.get_by_id(task.reaction_id).bounds = (extra, extra)
            before = snapshot(model)
            assert product_panel.biomass_pool_quota(model, product) == pytest.approx(expected)
            assert snapshot(model) == before
            solution = model.optimize()
            assert solution.status == "optimal"
            assert solution.fluxes["r_2111"] == pytest.approx(growth)
            assert solution.fluxes[task.reaction_id] == pytest.approx(extra)
            incorporation = -model.reactions.get_by_id("r_4048").get_coefficient(task.metabolite_id)
            assert incorporation * solution.fluxes["r_4048"] / growth == pytest.approx(expected)
            assert product_panel.biomass_pool_quota(model, product) == pytest.approx(expected)
    assert snapshot(original) == original_before
    assert snapshot(prepared) == prepared_before
