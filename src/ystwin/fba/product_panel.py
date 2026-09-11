from __future__ import annotations

import ast
import csv
import math
from dataclasses import dataclass

import cobra
import pandas as pd

from .. import paths
from ..pathway.proteome import (
    AVERAGE_PROTEIN_MW_G_PER_MOL,
    TOTAL_PROTEIN_G_PER_GDCW,
    enzyme_content_mmol_per_gdcw,
)


__all__ = ["ProductTask", "product_names", "prepare_panel_model", "apply_relative_abundances", "install_product", "biomass_pool_quota"]

_FROZEN_MODEL_ID = "M_ecYeastGEM_batch_v8__46__3__46__4"
_PROTEIN_POOL_ID = "prot_pool[c]"
_PAXDB_PATH = paths.data_dir() / "proteome" / "paxdb_scerevisiae_integrated.tsv"
_PAXDB_SOURCE = (
    "PaxDb S. cerevisiae whole-organism integrated molar ppm; "
    "data/proteome/paxdb_scerevisiae_integrated.tsv; gene mapping by string_id"
)
_TURNOVER_SOURCE = (
    "GECKO turnover parameters are model-derived, not necessarily measured; "
    "data/gem/ecYeastGEM_batch.xml.gz on yeast-GEM 8.3.4"
)
_SELECTED_ORFS = (
    "YHR190W",
    "YFR015C", "YLR258W", "YKR058W", "YJL137C",
    "YBR126C", "YDR074W", "YMR261C", "YML100W",
    "YJL101C", "YOL049W",
    "YDL022W", "YOL059W", "YIL053W", "YER062C",
)
_PRODUCTS = {
    "squalene": ("s_1447[c]", "C30H50", "MNXM292", "c"),
    "glycogen": ("s_0773[c]", "C6H12O6", "MNXM55375", "c"),
    "trehalose": ("s_1520[c]", "C12H22O11", "MNXM198", "c"),
    "glutathione": ("s_0750[c]", "C10H16N3O6S", "MNXM57", "c"),
    "glycerol": ("s_0766[e]", "C3H8O3", "MNXM89612", "e"),
}
_INTRACELLULAR_BASIS = (
    "extra intracellular net accumulation beyond the frozen biomass-incorporated pools"
)


@dataclass(frozen=True)
class ProductTask:
    name: str
    reaction_id: str
    metabolite_id: str
    molar_mass_g_per_mol: float
    carbon_atoms: float
    output_basis: str
    assumptions: tuple[str, ...]


def product_names() -> tuple[str, ...]:
    return tuple(_PRODUCTS)


def _require_frozen_model(model: cobra.Model) -> None:
    if getattr(model, "id", None) != _FROZEN_MODEL_ID:
        raise ValueError(
            f"the product panel requires the frozen EC batch model on yeast-GEM 8.3.4 "
            f"({_FROZEN_MODEL_ID!r}), got {getattr(model, 'id', None)!r}"
        )
    if _PROTEIN_POOL_ID not in model.metabolites:
        raise ValueError(f"the frozen EC model is missing protein pool {_PROTEIN_POOL_ID!r}")


def _read_proteomics() -> dict[str, tuple[str, str, float]]:
    selected = {f"4932.{orf}": orf for orf in _SELECTED_ORFS}
    records = {}
    with _PAXDB_PATH.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_name", "string_id", "abundance_ppm"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError(f"PaxDb table requires columns {sorted(required)}")
        for row in reader:
            string_id = (row.get("string_id") or "").strip()
            if string_id not in selected:
                continue
            orf = selected[string_id]
            if orf in records:
                raise ValueError(f"ambiguous PaxDb mapping: duplicate string_id {string_id!r}")
            gene = (row.get("gene_name") or "").strip()
            if not gene:
                raise ValueError(f"PaxDb gene_name is missing for {string_id!r}")
            try:
                ppm = float(row["abundance_ppm"])
            except (ValueError, TypeError) as exc:
                raise ValueError(f"PaxDb abundance is missing or invalid for {string_id!r}") from exc
            if not math.isfinite(ppm) or ppm < 0.0:
                raise ValueError(f"PaxDb abundance must be finite and non-negative for {string_id!r}")
            records[orf] = (gene, string_id, ppm)
    missing = sorted(set(_SELECTED_ORFS) - records.keys())
    if missing:
        raise ValueError(f"PaxDb string_id mapping is missing for selected yeast ORFs: {missing}")
    return records


def _protein_draws(model: cobra.Model) -> dict[str, tuple[cobra.Reaction, cobra.Metabolite]]:
    pool = model.metabolites.get_by_id(_PROTEIN_POOL_ID)
    selected = set(_SELECTED_ORFS)
    draws = {}
    proteins = set()
    for reaction in sorted(pool.reactions, key=lambda r: r.id):
        if reaction.metabolites[pool] >= 0.0:
            continue
        genes = {gene.id for gene in reaction.genes}
        if not genes.intersection(selected):
            continue
        rule = reaction.gpr.body
        if not isinstance(rule, ast.Name) or len(genes) != 1 or rule.id not in selected:
            raise ValueError(f"protein draw {reaction.id!r} has an ambiguous GPR: {reaction.gene_reaction_rule!r}")
        orf = rule.id
        others = [(met, coefficient) for met, coefficient in reaction.metabolites.items() if met != pool]
        if len(others) != 1 or others[0][1] != 1.0:
            raise ValueError(f"protein draw {reaction.id!r} for {orf} must produce exactly one enzyme molecule")
        protein = others[0][0]
        if not math.isfinite(reaction.metabolites[pool]) or reaction.lower_bound < 0.0:
            raise ValueError(f"protein draw {reaction.id!r} for {orf} has invalid mass or direction")
        if orf in draws or protein.id in proteins:
            raise ValueError(f"ambiguous protein draw mapping for {orf}: {reaction.id!r}, {protein.id!r}")
        suppliers = {
            r.id for r in protein.reactions
            if r.metabolites[protein] > 0.0
            or (r.metabolites[protein] < 0.0 and r.lower_bound < 0.0)
        }
        if suppliers != {reaction.id}:
            raise ValueError(f"protein {protein.id!r} for {orf} has multiple supply reactions: {sorted(suppliers)}")
        if not any(
            r.metabolites[protein] < 0.0 and orf in {gene.id for gene in r.genes}
            for r in protein.reactions
        ):
            raise ValueError(f"protein {protein.id!r} has no enzyme-consuming reaction mapped to {orf}")
        draws[orf] = (reaction, protein)
        proteins.add(protein.id)
    missing = sorted(selected - draws.keys())
    if missing:
        raise ValueError(f"missing native protein draw GPR mapping for selected ORFs: {missing}")
    return draws


def prepare_panel_model(
    model: cobra.Model,
    *,
    use_proteomics: bool = True,
    abundance_multiplier: float = 1.0,
) -> tuple[cobra.Model, pd.DataFrame]:
    try:
        multiplier = float(abundance_multiplier)
    except (TypeError, ValueError) as exc:
        raise ValueError("abundance_multiplier must be finite and non-negative") from exc
    if not math.isfinite(multiplier) or multiplier < 0.0:
        raise ValueError("abundance_multiplier must be finite and non-negative")
    _require_frozen_model(model)
    abundances = _read_proteomics()
    out = model.copy()
    draws = _protein_draws(out)
    assumptions = (
        "Integrated wild-type molar abundances are capacity proxies, not condition-matched active enzyme measurements.",
        f"The ppm conversion assumes total protein {TOTAL_PROTEIN_G_PER_GDCW:g} g/gDW and average protein mass {AVERAGE_PROTEIN_MW_G_PER_MOL:g} g/mol.",
        "The same union of 15 native enzyme budgets is used for every product; shared enzyme costs remain in the GEM.",
        "Only the selected protein draws are capped; alternative synthesis routes and other model enzyme budgets remain intact.",
        f"The global abundance multiplier is {multiplier:g}; it is an assumption, not fitted to product amounts.",
        "Proteomics caps are applied without relaxing tighter existing bounds."
        if use_proteomics else "Proteomics caps are not applied in this ablation; source values remain in the provenance.",
    )
    rows = []
    for orf in _SELECTED_ORFS:
        gene, string_id, ppm = abundances[orf]
        reaction, protein = draws[orf]
        content = enzyme_content_mmol_per_gdcw(ppm)
        scaled_content = content * multiplier
        if not math.isfinite(scaled_content):
            raise ValueError(f"abundance_multiplier produces a non-finite enzyme content for {orf}")
        before = float(reaction.upper_bound)
        after = min(before, scaled_content) if use_proteomics else before
        if after < reaction.lower_bound:
            raise ValueError(f"proteomics cap for {orf} on {reaction.id!r} is below its existing lower bound")
        if use_proteomics:
            reaction.upper_bound = after
        rows.append({
            "gene": gene,
            "orf": orf,
            "string_id": string_id,
            "draw_reaction_id": reaction.id,
            "protein_id": protein.id,
            "abundance_ppm": ppm,
            "enzyme_mmol_per_gdcw": content,
            "scaled_enzyme_mmol_per_gdcw": scaled_content,
            "abundance_multiplier": multiplier,
            "use_proteomics": bool(use_proteomics),
            "upper_bound_before": before,
            "upper_bound_after": after,
            "source": _PAXDB_SOURCE,
            "turnover_source": _TURNOVER_SOURCE,
            "model_id": model.id,
            "assumptions": assumptions,
        })
    return out, pd.DataFrame(rows)


def apply_relative_abundances(model, reference_inputs, abundance_ratios):
    _require_frozen_model(model)
    required = {"orf", "draw_reaction_id", "protein_id", "model_id", "use_proteomics",
                "scaled_enzyme_mmol_per_gdcw", "upper_bound_before", "upper_bound_after", "source"}
    if required - set(reference_inputs.columns):
        raise ValueError("complete reference enzyme-capacity provenance is required")
    out = model.copy()
    draws = _protein_draws(out)
    rows = []
    for orf, value in abundance_ratios.items():
        ratio = float(value)
        if isinstance(value, bool) or not math.isfinite(ratio) or ratio < 0:
            raise ValueError(f"relative abundance for {orf} must be finite and nonnegative")
        selected = reference_inputs.loc[reference_inputs.orf == orf]
        if len(selected) != 1 or orf not in draws:
            raise ValueError(f"a unique reference protein capacity is required for {orf}")
        reference = selected.iloc[0]
        reaction, protein = draws[orf]
        if (reference["model_id"] != model.id or reference["draw_reaction_id"] != reaction.id
                or reference["protein_id"] != protein.id):
            raise ValueError(f"enzyme identity does not match the reference capacity for {orf}")
        if reference["use_proteomics"] not in (True,):
            raise ValueError("relative abundance requires an applied reference abundance, not an uncapped ablation")
        if not math.isclose(reaction.upper_bound, float(reference["upper_bound_after"]), rel_tol=1e-9, abs_tol=1e-14):
            raise ValueError(f"current capacity for {orf} differs from the supplied reference; do not compound controls or mix dry-mass bases")
        content = float(reference["scaled_enzyme_mmol_per_gdcw"])
        hard_cap = float(reference["upper_bound_before"])
        requested = content * ratio
        if not math.isfinite(content) or content < 0 or not math.isfinite(requested) or math.isnan(hard_cap) or hard_cap < 0:
            raise ValueError(f"invalid reference abundance or hard capacity for {orf}")
        cap = min(hard_cap, requested)
        if cap < reaction.lower_bound:
            raise ValueError(f"relative abundance for {orf} conflicts with a required enzyme lower bound")
        reaction.upper_bound = cap
        rows.append({
            "orf": orf, "draw_reaction_id": reaction.id, "relative_abundance": ratio,
            "reference_abundance_mmol_per_gdcw": content,
            "reference_effective_cap": float(reference["upper_bound_after"]),
            "existing_hard_cap": hard_cap, "applied_cap_mmol_per_gdcw": cap,
            "reference_source": reference["source"],
            "assumption": "Caller-provided relative protein abundance scales available catalytic capacity; shared protein-pool costs and independent hard bounds remain unchanged.",
        })
    return out, pd.DataFrame(rows)


def _validated_species(
    model: cobra.Model, mid: str, formula: str, identity: str, compartment: str
) -> cobra.Metabolite:
    if mid not in model.metabolites:
        raise ValueError(f"required product species {mid!r} is missing from the frozen model")
    metabolite = model.metabolites.get_by_id(mid)
    aliases = metabolite.annotation.get("metanetx.chemical", ())
    if isinstance(aliases, str):
        aliases = (aliases,)
    if not isinstance(aliases, (tuple, list, set)) or identity not in aliases:
        raise ValueError(f"product species {mid!r} must have metanetx.chemical identity {identity!r}")
    if metabolite.compartment != compartment:
        raise ValueError(f"product species {mid!r} must be in compartment {compartment!r}")
    if metabolite.formula != formula:
        raise ValueError(f"product species {mid!r} has formula {metabolite.formula!r}, expected frozen-model formula {formula!r}")
    if not math.isfinite(metabolite.formula_weight) or metabolite.formula_weight <= 0.0:
        raise ValueError(f"product species {mid!r} has no positive finite formula mass")
    return metabolite


def biomass_pool_quota(model: cobra.Model, product: str) -> float:
    unavailable = f"biomass pool quota unavailable for {product!r}"
    if product not in ("glycogen", "trehalose"):
        raise ValueError(f"{unavailable}: only frozen carbohydrate storage pools are supported")
    _require_frozen_model(model)
    metabolite = _validated_species(model, *_PRODUCTS[product])
    try:
        carbohydrate = model.metabolites.get_by_id("s_3718[c]")
        biomass = model.metabolites.get_by_id("s_0450[c]")
        carbohydrate_assembly = model.reactions.get_by_id("r_4048")
        biomass_assembly = model.reactions.get_by_id("r_4041")
        growth = model.reactions.get_by_id("r_2111")
        reverse = model.reactions.get_by_id("r_2111_REV")
        coefficients = (
            -carbohydrate_assembly.metabolites[metabolite],
            carbohydrate_assembly.metabolites[carbohydrate],
            -biomass_assembly.metabolites[carbohydrate],
            biomass_assembly.metabolites[biomass],
        )
    except KeyError as exc:
        raise ValueError(f"{unavailable}: missing frozen biomass composition") from exc
    if (
        carbohydrate.reactions != {carbohydrate_assembly, biomass_assembly}
        or biomass.reactions != {biomass_assembly, growth, reverse}
        or set(carbohydrate_assembly.products) != {carbohydrate}
        or metabolite in biomass_assembly.metabolites
        or growth.metabolites != {biomass: -1.0}
        or reverse.metabolites != {biomass: 1.0}
        or reverse.bounds != (0.0, 0.0)
        or any(r.lower_bound < 0.0 for r in (carbohydrate_assembly, biomass_assembly, growth))
    ):
        raise ValueError(f"{unavailable}: ambiguous frozen biomass assembly")
    if any(not math.isfinite(c) or c <= 0.0 for c in coefficients):
        raise ValueError(f"{unavailable}: incorporation and biomass yields must be positive and finite")
    incorporation, carbohydrate_yield, carbohydrate_need, biomass_yield = coefficients
    quota = incorporation / carbohydrate_yield * carbohydrate_need / biomass_yield
    if not math.isfinite(quota) or quota <= 0.0:
        raise ValueError(f"{unavailable}: biomass quota must be positive and finite")
    return float(quota)


def _validate_sink(reaction: cobra.Reaction, metabolite: cobra.Metabolite) -> None:
    if not reaction.boundary or reaction.metabolites != {metabolite: -1.0}:
        raise ValueError(f"product demand or exchange {reaction.id!r} has incorrect sink stoichiometry for {metabolite.id!r}")
    if reaction.lower_bound < 0.0:
        raise ValueError(f"product sink {reaction.id!r} permits import, so its forward flux is not net production")


def install_product(model: cobra.Model, name: str, *, output_compartment: str | None = None) -> tuple[cobra.Model, ProductTask]:
    if name not in _PRODUCTS:
        raise ValueError(f"unsupported product {name!r}; supported products are {product_names()}")
    species = _PRODUCTS[name]
    if output_compartment is not None and output_compartment != species[3]:
        if name != "glycerol" or output_compartment != "c":
            raise ValueError(f"unsupported output compartment {output_compartment!r} for {name}")
        species = ("s_0765[c]", "C3H8O3", "MNXM89612", "c")
    _require_frozen_model(model)
    out = model.copy()
    metabolite = _validated_species(out, *species)
    assumptions = (
        f"Mass and carbon count use the frozen EC species {metabolite.id} formula {metabolite.formula}, without neutral-species or polymer-mass corrections.",
        "Native competing pathways, biomass incorporation, and storage turnover are retained unchanged.",
        "No starting product pool is supplied or inferred from product amount measurements.",
        "Installation supplies a net product task, not a realized production rate or an allocation rule.",
    )
    if name == "glycerol":
        _validated_species(out, "s_0765[c]", "C3H8O3", "MNXM89612", "c")
        extracellular = _validated_species(out, *_PRODUCTS["glycerol"])
        if "r_1808" not in out.reactions or "r_1808_REV" not in out.reactions:
            raise ValueError("native glycerol exchange r_1808 and its closed uptake arm r_1808_REV are required")
        reverse = out.reactions.get_by_id("r_1808_REV")
        if reverse.bounds != (0.0, 0.0) or reverse.metabolites != {extracellular: 1.0}:
            raise ValueError("r_1808_REV must remain closed to distinguish net glycerol production from import and reexport")
        _validate_sink(out.reactions.get_by_id("r_1808"), extracellular)
        assumptions += (
            "Native glycerol transport is reused; its uptake arm is closed and no new export capacity or energetic cost is invented.",
            "Alternative glycerol synthesis routes remain available; GPP capacity alone is not a total-product ceiling.",
        )
    if name == "glycerol" and metabolite.compartment == "e":
        reaction_id = "r_1808"
        output_basis = "extracellular net glycerol secretion through the native r_1808 exchange"
    else:
        reaction_id = f"DM_panel_{name}"
        output_basis = _INTRACELLULAR_BASIS
        if reaction_id not in out.reactions:
            out.add_boundary(metabolite, type="demand", reaction_id=reaction_id, lb=0.0, ub=float("inf"))
        assumptions += (
            "The demand measures extra net accumulation, not gross biosynthesis or the basal amount already incorporated into biomass.",
            "Added intracellular product mass must not be double counted as biomass; reservoir remobilization and time-dependent loss kinetics are not inferred here.",
        )
        if name == "glycogen":
            output_basis += "; glycogen repeat units reported as model glucose equivalents"
            assumptions += (
                "The frozen EC 8.3.4 model assigns C6H12O6 to a glycogen repeat unit; mass is therefore a glucose equivalent, not an anhydroglucose C6H10O5 polymer residue.",
            )
    _validate_sink(out.reactions.get_by_id(reaction_id), metabolite)
    task = ProductTask(
        name=name,
        reaction_id=reaction_id,
        metabolite_id=metabolite.id,
        molar_mass_g_per_mol=float(metabolite.formula_weight),
        carbon_atoms=float(metabolite.elements["C"]),
        output_basis=output_basis,
        assumptions=assumptions,
    )
    return out, task
