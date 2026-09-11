from __future__ import annotations

import json
import math
import pickle
from dataclasses import FrozenInstanceError, asdict
from types import MappingProxyType

import cobra
import pytest
from cobra.util.solver import linear_reaction_coefficients

from ystwin import paths
from ystwin.fba.product_panel import biomass_pool_quota, install_product, prepare_panel_model
from ystwin.fba.physiology import cap_uptake
from ystwin.fba.solver import load_model
from ystwin.fba.storage import StorageBasis, separate_storage_biomass


MODEL_ID = "M_ecYeastGEM_batch_v8__46__3__46__4"
QUOTAS = {"glycogen": 0.330522, "trehalose": 0.126456}
MASSES = {"glycogen": 180.15588, "trehalose": 342.29648}
SPECIES = {"glycogen": "s_0773[c]", "trehalose": "s_1520[c]"}
UNSCALED_RATES = {"r_4048", "r_4041", "r_2111", "r_2111_REV"}
FRACTION = 1.0 - sum(QUOTAS[p] * MASSES[p] / 1000.0 for p in QUOTAS)


def stoichiometry(reaction):
    return {metabolite.id: coefficient for metabolite, coefficient in reaction.metabolites.items()}


def snapshot(model):
    return (
        model.id,
        str(model.objective.expression),
        model.objective.direction,
        repr(model.notes),
        repr(model.annotation),
        tuple(
            (r.id, r.bounds, r.gene_reaction_rule, tuple(sorted(stoichiometry(r).items())))
            for r in model.reactions
        ),
        tuple((m.id, m.formula, m.compartment, repr(m.annotation)) for m in model.metabolites),
    )


def install_reserve_sinks(model):
    for product in QUOTAS:
        model, _ = install_product(model, product)
    return model


@pytest.fixture
def reference_model():
    model = cobra.Model(MODEL_ID)
    model.solver = "glpk"
    species = {
        "s_0773[c]": ("C6H12O6", "MNXM55375"),
        "s_1520[c]": ("C12H22O11", "MNXM198"),
    }
    ids = (
        "precursor[c]", "wall[c]", "other_structure[c]", "probe[c]",
        "s_3718[c]", "s_0450[c]", "s_0434[c]", "s_0803[c]",
        "s_0394[c]", "s_0794[c]", "s_1322[c]", "prot_pool[c]", "enzyme[c]",
        *species,
    )
    for mid in ids:
        metabolite = cobra.Metabolite(mid, compartment="c")
        if mid in species:
            metabolite.formula, identity = species[mid]
            metabolite.annotation["metanetx.chemical"] = identity
        model.add_metabolites([metabolite])
    hydrolysis = {
        "s_0434[c]": -1.0, "s_0803[c]": -1.0,
        "s_0394[c]": 1.0, "s_0794[c]": 1.0, "s_1322[c]": 1.0,
    }
    reactions = (
        ("feed", {"precursor[c]": 1.0}, (0.0, 10.0)),
        ("prot_pool_exchange", {"prot_pool[c]": 1.0}, (0.0, 0.5)),
        ("enzyme_draw", {"prot_pool[c]": -50.0, "enzyme[c]": 1.0}, (0.0, 0.005)),
        ("wall_synthesis", {"precursor[c]": -1.0, "enzyme[c]": -0.001, "wall[c]": 1.0}, (0.0, math.inf)),
        ("other_synthesis", {"precursor[c]": -1.0, "other_structure[c]": 1.0}, (0.0, math.inf)),
        ("glycogen_synthesis", {"precursor[c]": -1.0, "enzyme[c]": -0.001, "s_0773[c]": 1.0}, (0.0, math.inf)),
        ("trehalose_synthesis", {"precursor[c]": -2.0, "enzyme[c]": -0.002, "s_1520[c]": 1.0}, (0.0, math.inf)),
        ("native_glycogen_turnover", {"s_0773[c]": -1.0, "precursor[c]": 1.0, "enzyme[c]": -0.001}, (0.0, 3.0)),
        ("native_trehalose_turnover", {"s_1520[c]": -1.0, "precursor[c]": 2.0, "enzyme[c]": -0.001}, (0.0, 2.0)),
        ("regenerate_energy", {**{mid: -2.0 * c for mid, c in hydrolysis.items()}, "precursor[c]": -1.0}, (0.0, math.inf)),
        ("r_4046", hydrolysis, (0.05, 0.05)),
        ("r_4048", {"wall[c]": -0.5, "s_0773[c]": -QUOTAS["glycogen"], "s_1520[c]": -QUOTAS["trehalose"], "s_3718[c]": 1.0}, (0.0, 2.0)),
        ("r_4041", {**{mid: 3.0 * c for mid, c in hydrolysis.items()}, "s_3718[c]": -1.0, "other_structure[c]": -0.4, "s_0450[c]": 1.0}, (0.0, 2.0)),
        ("r_2111", {"s_0450[c]": -1.0}, (0.0, 1.0)),
        ("r_2111_REV", {"s_0450[c]": 1.0}, (0.0, 0.0)),
        ("reversible_probe", {"precursor[c]": -1.0, "probe[c]": 1.0}, (-4.0, 5.0)),
    )
    for rid, coefficients, bounds in reactions:
        reaction = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
        reaction.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in coefficients.items()})
        model.add_reactions([reaction])
    model.objective = "r_2111"
    model.notes = {"source": {"labels": ["reference"]}}
    model.annotation = {"source": ["reference"]}
    return model


@pytest.fixture
def basis(reference_model):
    return separate_storage_biomass(reference_model)[1]


def test_basis_is_derived_from_quotas_and_explicit_model_species(reference_model):
    transformed, basis = separate_storage_biomass(reference_model)
    assert isinstance(basis, StorageBasis)
    assert basis.original_model_id == MODEL_ID
    assert basis.structural_mass_fraction == pytest.approx(FRACTION, abs=1e-14)
    assert basis.removed_quotas_mmol_per_reference_g == pytest.approx(QUOTAS)
    assert basis.molecular_masses == pytest.approx(MASSES)
    assert transformed.id != reference_model.id
    assert basis.metadata["transformed_model_id"] == transformed.id
    assert basis.metadata["species"]["glycogen"]["metabolite_id"] == "s_0773[c]"
    assert basis.metadata["species"]["glycogen"]["formula"] == "C6H12O6"
    text = json.dumps(basis.metadata).lower()
    assert "glucose equivalent" in text and "not" in text and "polymer" in text
    assert basis.metadata["inventory_remobilization_supported"] is False


def test_only_nonstorage_biomass_requirements_are_renormalized(reference_model):
    transformed, basis = separate_storage_biomass(reference_model)
    f = basis.structural_mass_fraction
    for original in reference_model.reactions:
        actual = stoichiometry(transformed.reactions.get_by_id(original.id))
        expected = stoichiometry(original)
        if original.id == "r_4048":
            expected = {
                mid: coefficient if mid == "s_3718[c]" else coefficient / f
                for mid, coefficient in expected.items() if mid not in SPECIES.values()
            }
        elif original.id == "r_4041":
            expected = {
                mid: coefficient if mid in {"s_3718[c]", "s_0450[c]"} else coefficient / f
                for mid, coefficient in expected.items()
            }
        assert actual == pytest.approx(expected)
    assert transformed.reactions.get_by_id("r_4048").get_coefficient("s_3718[c]") == 1.0
    assert transformed.reactions.get_by_id("r_4041").get_coefficient("s_3718[c]") == -1.0
    assert transformed.reactions.get_by_id("r_4041").get_coefficient("s_0450[c]") == 1.0


def test_quota_derivation_does_not_assume_unit_pseudometabolite_coefficients(reference_model):
    reference_model.reactions.get_by_id("r_4048").add_metabolites({
        reference_model.metabolites.get_by_id("s_0773[c]"): -0.6,
        reference_model.metabolites.get_by_id("s_1520[c]"): -0.4,
        reference_model.metabolites.get_by_id("s_3718[c]"): 2.0,
    }, combine=False)
    reference_model.reactions.get_by_id("r_4041").add_metabolites({
        reference_model.metabolites.get_by_id("s_3718[c]"): -3.0,
        reference_model.metabolites.get_by_id("s_0450[c]"): 4.0,
    }, combine=False)
    transformed, basis = separate_storage_biomass(reference_model)
    quotas = {"glycogen": 0.225, "trehalose": 0.15}
    assert basis.removed_quotas_mmol_per_reference_g == pytest.approx(quotas)
    assert basis.structural_mass_fraction == pytest.approx(
        1.0 - sum(quotas[p] * MASSES[p] / 1000.0 for p in quotas)
    )
    assert transformed.reactions.get_by_id("r_4048").get_coefficient("s_3718[c]") == 2.0
    assert transformed.reactions.get_by_id("r_4041").get_coefficient("s_3718[c]") == -3.0
    assert transformed.reactions.get_by_id("r_4041").get_coefficient("s_0450[c]") == 4.0


def test_only_selected_storage_is_excluded_from_structural_mass(reference_model):
    transformed, basis = separate_storage_biomass(reference_model, products=("glycogen",))
    f = 1.0 - QUOTAS["glycogen"] * MASSES["glycogen"] / 1000.0
    assert basis.structural_mass_fraction == pytest.approx(f)
    assert basis.removed_quotas_mmol_per_reference_g == {"glycogen": QUOTAS["glycogen"]}
    assert "s_0773[c]" not in stoichiometry(transformed.reactions.get_by_id("r_4048"))
    assert transformed.reactions.get_by_id("r_4048").get_coefficient("s_1520[c]") == pytest.approx(-QUOTAS["trehalose"] / f)
    with pytest.raises(ValueError, match="unknown|unsupported|separated"):
        basis.total_dry_mass(1.0, {"trehalose": 0.0})


def test_zero_storage_has_no_implicit_basal_floor(basis):
    assert basis.total_dry_mass(2.0, {}) == 2.0
    for product in QUOTAS:
        assert basis.content_mg_per_g_total(product, 2.0, {}) == 0.0
        assert basis.content_mg_per_g_total(product, 2.0, dict.fromkeys(QUOTAS, 0.0)) == 0.0


def test_positive_inventories_contribute_mass_and_share_the_total_denominator(basis):
    pools = MappingProxyType({"glycogen": 3.0, "trehalose": 0.7})
    total = 2.0 + sum(pools[p] * MASSES[p] / 1000.0 for p in pools)
    assert basis.total_dry_mass(2.0, pools) == pytest.approx(total)
    for product in pools:
        assert basis.content_mg_per_g_total(product, 2.0, pools) == pytest.approx(pools[product] * MASSES[product] / total)


def test_reconstructed_reference_composition_is_not_double_counted(basis):
    reference_biomass = 2.3
    structural = FRACTION * reference_biomass
    pools = {p: quota * reference_biomass for p, quota in QUOTAS.items()}
    assert basis.total_dry_mass(structural, pools) == pytest.approx(reference_biomass)
    for product in pools:
        assert basis.content_mg_per_g_total(product, structural, pools) == pytest.approx(QUOTAS[product] * MASSES[product])


def test_zero_total_mass_is_valid_but_its_content_is_undefined(basis):
    assert basis.total_dry_mass(0.0, {}) == 0.0
    with pytest.raises(ValueError, match="total.*mass|zero|positive"):
        basis.content_mg_per_g_total("glycogen", 0.0, {})
    assert basis.content_mg_per_g_total("glycogen", 0.0, {"glycogen": 1.0}) == pytest.approx(1000.0)


@pytest.mark.parametrize("value", (-11.1, -0.01, 0.0, 0.49974, 1.0e-8, 56.6883))
def test_specific_unit_conversion_roundtrip_preserves_volumetric_rates(basis, value):
    structural_specific = basis.to_structural_specific(value)
    assert structural_specific * FRACTION == pytest.approx(value)
    reference_biomass = 2.3
    assert structural_specific * (FRACTION * reference_biomass) == pytest.approx(value * reference_biomass)


def test_all_specific_bounds_including_enzyme_amounts_and_maintenance_are_converted(reference_model):
    reference_model = install_reserve_sinks(reference_model)
    reference_model.reactions.get_by_id("r_2111").bounds = (0.03, 0.12)
    transformed, basis = separate_storage_biomass(reference_model)
    for original in reference_model.reactions:
        expected = original.bounds if original.id in UNSCALED_RATES else tuple(b / FRACTION for b in original.bounds)
        assert transformed.reactions.get_by_id(original.id).bounds == pytest.approx(expected)
    assert set(basis.metadata["unchanged_rate_reaction_ids"]) == UNSCALED_RATES
    assert set(basis.metadata["converted_bound_reaction_ids"]) == set(reference_model.reactions.list_attr("id")) - UNSCALED_RATES
    assert basis.metadata["specific_amount_conversion_factor"] == pytest.approx(1.0 / FRACTION)
    assert basis.metadata["growth_rate_units"] == "1/h"
    assert "enzyme_draw" in basis.metadata["enzyme_amount_reaction_ids"]
    assert basis.metadata["protein_pool_supply_reaction_ids"] == ("prot_pool_exchange",)
    assert basis.metadata["accumulation_reaction_ids"] == {
        p: (f"DM_panel_{p}",) for p in QUOTAS
    }
    with pytest.raises(ValueError, match="frozen|already|model"):
        separate_storage_biomass(transformed)


def test_existing_linear_objective_is_reexpressed_in_compatible_units(reference_model):
    reference_model.objective = {
        reference_model.reactions.get_by_id("r_2111"): 2.0,
        reference_model.reactions.get_by_id("feed"): -0.1,
    }
    reference_model.objective.direction = "min"
    transformed, _ = separate_storage_biomass(reference_model)
    coefficients = {r.id: c for r, c in linear_reaction_coefficients(transformed).items()}
    assert coefficients == pytest.approx({"r_2111": 2.0, "feed": -0.1 * FRACTION})
    assert transformed.objective.direction == "min"


def test_transformation_is_opt_in_does_not_solve_and_leaves_original_untouched(reference_model, monkeypatch):
    before = snapshot(reference_model)
    monkeypatch.setattr(reference_model, "optimize", None)
    monkeypatch.setattr(reference_model, "slim_optimize", None)
    transformed, _ = separate_storage_biomass(reference_model)
    assert transformed is not reference_model
    assert snapshot(reference_model) == before
    transformed.notes["source"]["labels"].append("copy")
    transformed.annotation["source"].append("copy")
    transformed.reactions.get_by_id("feed").upper_bound = 0.0
    transformed.metabolites.get_by_id("s_0773[c]").annotation["metanetx.chemical"] = "changed"
    assert snapshot(reference_model) == before


def test_no_reaction_is_added_or_removed_and_no_signed_inventory_supply_is_opened(reference_model):
    reference_model = install_reserve_sinks(reference_model)
    transformed, _ = separate_storage_biomass(reference_model)
    assert transformed.reactions.list_attr("id") == reference_model.reactions.list_attr("id")
    assert transformed.metabolites.list_attr("id") == reference_model.metabolites.list_attr("id")
    for product in QUOTAS:
        sink = transformed.reactions.get_by_id(f"DM_panel_{product}")
        assert stoichiometry(sink) == {SPECIES[product]: -1.0}
        assert sink.lower_bound == 0.0
        turnover = f"native_{product}_turnover"
        assert stoichiometry(transformed.reactions.get_by_id(turnover)) == stoichiometry(reference_model.reactions.get_by_id(turnover))
    assert transformed.reactions.get_by_id("r_2111_REV").bounds == (0.0, 0.0)


def test_growth_can_be_feasible_with_zero_storage_synthesis(reference_model):
    model, _ = separate_storage_biomass(install_reserve_sinks(reference_model))
    for product in QUOTAS:
        model.reactions.get_by_id(f"{product}_synthesis").upper_bound = 0.0
        model.reactions.get_by_id(f"DM_panel_{product}").upper_bound = 0.0
    model.reactions.get_by_id("r_2111").bounds = (0.1, 0.1)
    solution = model.optimize()
    assert solution.status == "optimal"
    assert solution.fluxes["r_2111"] == pytest.approx(0.1)


def test_basis_and_its_accounting_maps_are_frozen_and_json_serializable(basis):
    with pytest.raises(FrozenInstanceError):
        basis.structural_mass_fraction = 1.0
    with pytest.raises(TypeError):
        basis.molecular_masses["glycogen"] = 162.14
    with pytest.raises(TypeError):
        basis.removed_quotas_mmol_per_reference_g["glycogen"] = 0.0
    with pytest.raises(TypeError):
        basis.metadata["species"]["glycogen"]["formula"] = "C6H10O5"
    payload = json.loads(json.dumps(asdict(basis)))
    assert payload["molecular_masses"] == pytest.approx(MASSES)
    assert isinstance(basis.molecular_masses, dict)


def test_basis_can_cross_a_worker_boundary_without_losing_immutability(basis):
    restored = pickle.loads(pickle.dumps(basis))
    assert restored == basis
    with pytest.raises(TypeError):
        restored.molecular_masses["glycogen"] = 0.0


def test_direct_basis_construction_detaches_the_input_maps():
    quotas, masses, metadata = dict(QUOTAS), dict(MASSES), {"source": ["model"]}
    basis = StorageBasis(MODEL_ID, FRACTION, quotas, masses, metadata)
    quotas["glycogen"] = 0.0
    masses["glycogen"] = 0.0
    metadata["source"].append("changed")
    assert basis.removed_quotas_mmol_per_reference_g == QUOTAS
    assert basis.molecular_masses == MASSES
    assert basis.metadata["source"] == ("model",)


@pytest.mark.parametrize("value", (-1.0, math.nan, math.inf, -math.inf, None, "unknown", True))
def test_invalid_structural_mass_is_refused(basis, value):
    with pytest.raises(ValueError, match="structural|mass"):
        basis.total_dry_mass(value, {})
    with pytest.raises(ValueError, match="structural|mass"):
        basis.content_mg_per_g_total("glycogen", value, {})


@pytest.mark.parametrize("pools", (
    {"glycogen": -1.0}, {"trehalose": math.nan}, {"glycogen": math.inf},
    {"trehalose": -math.inf}, {"glycogen": None}, {"glycogen": "unknown"},
    {"glycogen": True}, {"squalene": 0.0}, {"glycogen": 1.0, "typo": 0.0}, None, [],
))
def test_invalid_or_unknown_inventory_values_are_refused(basis, pools):
    with pytest.raises(ValueError, match="pool|inventory|unknown|supported|glycogen|trehalose"):
        basis.total_dry_mass(1.0, pools)
    with pytest.raises(ValueError):
        basis.content_mg_per_g_total("glycogen", 1.0, pools)


@pytest.mark.parametrize("product", ("squalene", "glycogen_typo", None))
def test_unknown_content_product_is_refused(basis, product):
    with pytest.raises(ValueError, match="product|unknown|separated"):
        basis.content_mg_per_g_total(product, 1.0, {})


@pytest.mark.parametrize("value", (math.nan, math.inf, -math.inf, None, "unknown", True))
def test_specific_conversion_refuses_nonfinite_or_invalid_values(basis, value):
    with pytest.raises(ValueError, match="finite|specific"):
        basis.to_structural_specific(value)


def test_finite_inputs_cannot_silently_overflow_accounting(basis):
    with pytest.raises(ValueError, match="finite|overflow"):
        basis.total_dry_mass(1.7e308, {"glycogen": 1e308, "trehalose": 1e308})
    assert math.isfinite(basis.content_mg_per_g_total("glycogen", 1.0, {"glycogen": 1e308}))
    with pytest.raises(ValueError, match="finite|overflow"):
        basis.to_structural_specific(1.7e308)


@pytest.mark.parametrize("products", ((), ("glycogen", "glycogen"), ("squalene",), "glycogen", None))
def test_invalid_product_selection_is_refused_without_mutation(reference_model, products):
    before = snapshot(reference_model)
    with pytest.raises(ValueError, match="product|storage|supported|duplicate"):
        separate_storage_biomass(reference_model, products=products)
    assert snapshot(reference_model) == before


@pytest.mark.parametrize("change", ("identity", "formula", "model", "growth_import", "mass_fraction"))
def test_invalid_reconstructions_are_refused_without_mutation(reference_model, change):
    glycogen = reference_model.metabolites.get_by_id("s_0773[c]")
    if change == "identity":
        glycogen.annotation["metanetx.chemical"] = "MNXM198"
    elif change == "formula":
        glycogen.formula = "C6H10O5"
    elif change == "model":
        reference_model.id = "yeast9"
    elif change == "growth_import":
        reference_model.reactions.get_by_id("r_2111_REV").upper_bound = 1.0
    else:
        reference_model.reactions.get_by_id("r_4048").add_metabolites({glycogen: -10.0}, combine=False)
    before = snapshot(reference_model)
    with pytest.raises(ValueError):
        separate_storage_biomass(reference_model)
    assert snapshot(reference_model) == before


@pytest.mark.parametrize("field,value", (
    ("structural_mass_fraction", 1.0), ("structural_mass_fraction", math.nan),
    ("structural_mass_fraction", 0.0),
    ("removed_quotas_mmol_per_reference_g", {"glycogen": -1.0, "trehalose": QUOTAS["trehalose"]}),
    ("molecular_masses", {"glycogen": math.inf, "trehalose": MASSES["trehalose"]}),
    ("molecular_masses", {"glycogen": 180.0}),
))
def test_inconsistent_direct_basis_construction_is_refused(field, value):
    arguments = dict(
        original_model_id=MODEL_ID,
        structural_mass_fraction=FRACTION,
        removed_quotas_mmol_per_reference_g=dict(QUOTAS),
        molecular_masses=dict(MASSES),
        metadata={},
    )
    arguments[field] = value
    with pytest.raises(ValueError):
        StorageBasis(**arguments)


@pytest.mark.parametrize("kind", ("constraint", "variable", "nonzero_balance", "matrix_edit", "variable_bound"))
def test_unsupported_solver_extensions_are_refused_instead_of_misconverted(reference_model, kind):
    if kind == "constraint":
        reference_model.add_cons_vars(reference_model.problem.Constraint(
            reference_model.reactions.get_by_id("feed").flux_expression,
            lb=0.0, ub=1.0, name="external_budget",
        ))
    elif kind == "variable":
        reference_model.add_cons_vars(reference_model.problem.Variable("external_amount", lb=0.0))
    elif kind == "variable_bound":
        reference_model.reactions.get_by_id("feed").forward_variable.ub = 1.0
    elif kind == "matrix_edit":
        reference_model.metabolites.get_by_id("wall[c]").constraint.set_linear_coefficients({
            reference_model.reactions.get_by_id("feed").forward_variable: 0.5,
        })
    else:
        constraint = reference_model.metabolites.get_by_id("wall[c]").constraint
        constraint.ub = 0.1
        constraint.lb = 0.1
    before = snapshot(reference_model)
    with pytest.raises(ValueError, match="constraint|variable|balance|steady"):
        separate_storage_biomass(reference_model)
    assert snapshot(reference_model) == before


@pytest.mark.parametrize("kind", ("signed", "supply", "nonunit_sink"))
def test_storage_boundaries_cannot_masquerade_as_inventory_remobilization(reference_model, kind):
    reference = install_reserve_sinks(reference_model)
    demand = reference.reactions.get_by_id("DM_panel_glycogen")
    if kind == "signed":
        demand.lower_bound = -1.0
    else:
        demand.add_metabolites({
            reference.metabolites.get_by_id("s_0773[c]"): 1.0 if kind == "supply" else -2.0,
        }, combine=False)
    before = snapshot(reference)
    with pytest.raises(ValueError, match="nonnegative|unit demand|inventory supply"):
        separate_storage_biomass(reference)
    assert snapshot(reference) == before


def test_an_objective_not_written_in_net_fluxes_is_refused(reference_model):
    reference_model.objective = reference_model.problem.Objective(
        reference_model.reactions.get_by_id("reversible_probe").forward_variable,
        direction="max",
    )
    with pytest.raises(ValueError, match="objective|net reaction"):
        separate_storage_biomass(reference_model)


@pytest.fixture(scope="module")
def frozen_model():
    model_path = paths.data_dir() / "gem" / "ecYeastGEM_batch.xml.gz"
    if not model_path.is_file():
        pytest.skip("vendored frozen EC batch model is absent")
    return load_model(model_path, tolerance=1e-9)[0]


@pytest.fixture(scope="module")
def frozen_prepared(frozen_model):
    if not (paths.data_dir() / "proteome" / "paxdb_scerevisiae_integrated.tsv").is_file():
        pytest.skip("vendored PaxDb enzyme budgets are absent")
    prepared, _ = prepare_panel_model(frozen_model)
    return install_reserve_sinks(prepared)


def maximum_balance_residual(model, fluxes):
    return max(abs(math.fsum(r.metabolites[m] * fluxes[r.id] for r in m.reactions)) for m in model.metabolites)


@pytest.mark.integration
def test_real_model_basis_energy_and_native_turnover_are_preserved(frozen_prepared):
    before = snapshot(frozen_prepared)
    transformed, basis = separate_storage_biomass(frozen_prepared)
    assert basis.removed_quotas_mmol_per_reference_g == pytest.approx(QUOTAS)
    assert basis.molecular_masses == pytest.approx(MASSES)
    assert basis.structural_mass_fraction == pytest.approx(0.89716907455576)
    assert transformed.reactions.get_by_id("r_4041").get_coefficient("s_0434[c]") == pytest.approx(-56.6883 / FRACTION)
    assert transformed.reactions.get_by_id("r_4046").bounds == pytest.approx((0.7 / FRACTION, 0.7 / FRACTION))
    assert transformed.reactions.get_by_id("prot_pool_exchange").upper_bound == pytest.approx(0.103720024310123 / FRACTION)
    for original in frozen_prepared.reactions:
        if original.id not in {"r_4048", "r_4041"}:
            assert stoichiometry(transformed.reactions.get_by_id(original.id)) == stoichiometry(original)
    for product in QUOTAS:
        assert biomass_pool_quota(frozen_prepared, product) == pytest.approx(QUOTAS[product])
    assert snapshot(frozen_prepared) == before


@pytest.mark.integration
@pytest.mark.parametrize("glucose,oxygen,proteomics", ((1.5, 5.0, True), (20.0, 1000.0, True), (20.0, 1000.0, False)))
def test_real_model_equivalent_balanced_growth_reference(frozen_model, frozen_prepared, glucose, oxygen, proteomics):
    reference = frozen_prepared.copy() if proteomics else install_reserve_sinks(frozen_model)
    glucose_id = cap_uptake(reference, "r_1714", glucose)
    reference.reactions.get_by_id(glucose_id).lower_bound = glucose
    cap_uptake(reference, "r_1992", oxygen)
    reference.objective = "r_2111"
    with reference:
        for product in QUOTAS:
            reference.reactions.get_by_id(f"DM_panel_{product}").bounds = (0.0, 0.0)
        old = reference.optimize()
        assert old.status == "optimal"
        assert old.fluxes["r_2111"] > 0.0
        assert maximum_balance_residual(reference, old.fluxes) < 1e-6
    transformed, basis = separate_storage_biomass(reference)
    f = basis.structural_mass_fraction
    growth = transformed.reactions.get_by_id("r_2111")
    for product, quota in basis.removed_quotas_mmol_per_reference_g.items():
        demand = transformed.reactions.get_by_id(f"DM_panel_{product}")
        transformed.add_cons_vars(transformed.problem.Constraint(
            demand.flux_expression - (quota / f) * growth.flux_expression,
            lb=0.0, ub=0.0, name=f"reference_inventory_{product}",
        ))
    new = transformed.optimize()
    assert new.status == "optimal"
    assert new.fluxes["r_2111"] == pytest.approx(old.fluxes["r_2111"], rel=1e-6, abs=1e-8)
    assert maximum_balance_residual(transformed, new.fluxes) < 1e-6
    forward = {
        r.id: old.fluxes[r.id] if r.id in UNSCALED_RATES else old.fluxes[r.id] / f
        for r in reference.reactions
    }
    reverse = {
        r.id: new.fluxes[r.id] if r.id in UNSCALED_RATES else new.fluxes[r.id] * f
        for r in reference.reactions
    }
    for product, quota in QUOTAS.items():
        rid = f"DM_panel_{product}"
        forward[rid] += quota * old.fluxes["r_2111"] / f
        reverse[rid] -= quota * new.fluxes["r_2111"]
        assert new.fluxes[rid] == pytest.approx(quota * new.fluxes["r_2111"] / f)
    assert maximum_balance_residual(transformed, forward) < 1e-6
    assert maximum_balance_residual(reference, reverse) < 1e-6
    for r in transformed.reactions:
        assert r.lower_bound - 1e-7 <= forward[r.id] <= r.upper_bound + 1e-7
    for r in reference.reactions:
        assert r.lower_bound - 1e-7 <= reverse[r.id] <= r.upper_bound + 1e-7
    new_total_mass_rate = new.fluxes["r_2111"] + sum(
        new.fluxes[f"DM_panel_{p}"] * MASSES[p] / 1000.0 for p in QUOTAS
    )
    assert new_total_mass_rate == pytest.approx(old.fluxes["r_2111"] / f, rel=1e-6)
    assert new_total_mass_rate / new.fluxes[glucose_id] == pytest.approx(
        old.fluxes["r_2111"] / old.fluxes[glucose_id], rel=1e-6,
    )
    pools = {p: quota * 1.7 for p, quota in QUOTAS.items()}
    assert basis.total_dry_mass(f * 1.7, pools) == pytest.approx(1.7)
    for product, quota in QUOTAS.items():
        assert new.fluxes[f"DM_panel_{product}"] * (f * 1.7) == pytest.approx(quota * old.fluxes["r_2111"] * 1.7, rel=1e-6)


@pytest.mark.integration
def test_real_model_supports_growth_without_a_hidden_storage_synthesis_floor(frozen_prepared):
    transformed, basis = separate_storage_biomass(frozen_prepared)
    glucose_id = cap_uptake(transformed, "r_1714", basis.to_structural_specific(1.5))
    cap_uptake(transformed, "r_1992", basis.to_structural_specific(5.0))
    transformed.objective = glucose_id
    transformed.objective.direction = "min"
    for product, mid in SPECIES.items():
        metabolite = transformed.metabolites.get_by_id(mid)
        for reaction in metabolite.reactions:
            if reaction.metabolites[metabolite] > 0.0:
                assert reaction.lower_bound == 0.0
                reaction.upper_bound = 0.0
        transformed.reactions.get_by_id(f"DM_panel_{product}").upper_bound = 0.0
    transformed.reactions.get_by_id("r_2111").bounds = (0.05, 0.05)
    solution = transformed.optimize()
    assert solution.status == "optimal"
    assert solution.fluxes["r_2111"] == pytest.approx(0.05)
    assert maximum_balance_residual(transformed, solution.fluxes) < 1e-6
