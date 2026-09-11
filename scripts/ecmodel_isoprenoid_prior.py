from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import pathlib
import sys
import tempfile
import xml.etree.ElementTree as ET
import zlib
from dataclasses import asdict, dataclass

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.pathway import enzyme_capacity, proteome

_SBML = "{http://www.sbml.org/sbml/level3/version1/core}"
_FBC = "{http://www.sbml.org/sbml/level3/version1/fbc/version2}"
_FLUX_UNIT_ID = "mmol_per_gDW_per_hr"
BOUNDARY_BASIS = (
    "GECKO2: prot_* reactant coefficients interpreted as hours (1/kcat_h); "
    "enzyme mmol/gDCW; reaction flux mmol/gDCW/h; declared convention, "
    "not inferred from unannotated protein species or flux bounds; no medium or flux bounds applied"
)
BRANCH_GROUPS = {
    "MVA": ("ERG10", "ERG13", "HMG1", "HMG2", "ERG12", "ERG8", "MVD1", "IDI1"),
    "prenyl_transferase": ("ERG20", "BTS1"),
    "sterol": ("ERG9", "ERG1", "ERG7", "ERG11", "ERG24", "ERG25", "ERG26", "ERG27",
               "ERG6", "ERG2", "ERG3", "ERG5", "ERG4"),
}
SWEEP_MASS_FRACTIONS = (0.0005, 0.001, 0.005, 0.01, 0.05)
DEFAULT_ENZYME = "crtYB"
DEFAULT_ENZYME_MOLAR_MASS_G_PER_MOL = 74736.0
DEFAULT_GROWTH_RATE_PER_H = 0.101
INVENTORY_FIELDS = (
    "reaction_id", "declared_gpr_genes", "gpr_genes", "branch", "protein_species_id",
    "protein_coefficient_h", "kcat_per_s", "is_reverse_arm", "included_in_forward_span",
    "association_basis", "interpretation", "source_model_sha256", "source_sbml_sha256",
    "boundary_basis",
)
SWEEP_FIELDS = (
    "enzyme", "kcat_per_s", "mass_fraction_of_protein", "enzyme_molar_mass_g_per_mol",
    "total_protein_g_per_gdcw", "growth_rate_per_h", "enzyme_mmol_per_gdcw",
    "vmax_mmol_per_gdcw_h", "capacity_mmol_per_gdcw", "reference_content_source",
    "reference_content_mg_per_gdcw", "reference_product_molar_mass_g_per_mol",
    "reference_content_mmol_per_gdcw", "fold_over_reference_content",
    "reference_capacity_source", "reference_capacity_mmol_per_gdcw",
    "fold_over_reference_capacity", "reference_basis", "interpretation", "kcat_basis",
    "source_model_sha256", "source_sbml_sha256", "boundary_basis",
)


def _positive(value, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f"{name} must be finite and positive, got {value!r}") from None
    if isinstance(value, bool) or not math.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be finite and positive, got {value!r}")
    return number


def _label(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty label")
    return value


def _axis(values, name: str) -> tuple[float, ...]:
    result = tuple(_positive(value, name) for value in values)
    if not result:
        raise ValueError(f"{name} must contain at least one value")
    return result


@dataclass(frozen=True)
class TurnoverInventory:
    rows: list[dict]
    reactions_with_gpr: int
    reactions_with_kcat: int
    source_path: str
    source_model_sha256: str
    source_sbml_sha256: str
    model_id: str


@dataclass(frozen=True)
class ContentReference:
    source: str
    content_mg_per_gdcw: float
    product_molar_mass_g_per_mol: float

    def __post_init__(self):
        _label(self.source, "reference_content_source")
        _positive(self.content_mg_per_gdcw, "reference_content_mg_per_gdcw")
        _positive(self.product_molar_mass_g_per_mol, "reference_product_molar_mass_g_per_mol")
        _positive(self.mmol_per_gdcw, "reference_content_mmol_per_gdcw")

    @property
    def mmol_per_gdcw(self) -> float:
        return self.content_mg_per_gdcw / self.product_molar_mass_g_per_mol


@dataclass(frozen=True)
class CapacityReference:
    source: str
    capacity_mmol_per_gdcw: float

    def __post_init__(self):
        _label(self.source, "reference_capacity_source")
        _positive(self.capacity_mmol_per_gdcw, "reference_capacity_mmol_per_gdcw")


def _index(elements, attribute: str, name: str) -> dict[str, ET.Element]:
    result = {}
    for element in elements:
        key = element.get(attribute)
        if not key or key in result:
            raise ValueError(f"missing or duplicate {name}: {key!r}")
        result[key] = element
    return result


def _verify_unit_basis(model: ET.Element) -> None:
    definition = model.find(
        f"{_SBML}listOfUnitDefinitions/{_SBML}unitDefinition[@id='{_FLUX_UNIT_ID}']")
    if definition is None:
        raise ValueError(f"GECKO2 basis requires the {_FLUX_UNIT_ID} unit definition")
    units = [
        (unit.get("kind"), float(unit.get("exponent", "nan")),
         float(unit.get("scale", "nan")), float(unit.get("multiplier", "nan")))
        for unit in definition.iter(_SBML + "unit")
    ]
    expected = [("mole", 1, -3, 1), ("gram", -1, 0, 1), ("second", -1, 0, 3600)]
    if len(units) != len(expected) or set(units) != set(expected):
        raise ValueError(f"unsupported {_FLUX_UNIT_ID} definition; expected mmol/gDW/hour")
    if any(model.get(field) is not None for field in ("timeUnits", "substanceUnits", "extentUnits")):
        raise ValueError("explicit model unit overrides are unsupported by this GECKO2 convention")


def _protein_demands(reaction: ET.Element, species: dict[str, ET.Element]):
    reaction_id = reaction.get("id")
    coefficients = {}
    for ref in reaction.findall(f"{_SBML}listOfReactants/{_SBML}speciesReference"):
        species_id = ref.get("species", "")
        if not species_id.startswith("prot_") or species_id.startswith("prot_pool"):
            continue
        name = f"{reaction_id} {species_id}"
        if species_id not in species or species_id in coefficients:
            raise ValueError(f"{name}: unknown or duplicate protein species")
        if species[species_id].get("substanceUnits") is not None:
            raise ValueError(f"{name}: explicit protein substanceUnits require a different unit map")
        if ref.get("constant") not in (None, "true", "1") or len(ref):
            raise ValueError(f"{name}: protein demand must be a constant literal stoichiometry")
        coefficient = _positive(ref.get("stoichiometry"), f"{name} coefficient")
        kcat = _positive(1.0 / (coefficient * 3600.0), f"{name} kcat_per_s")
        coefficients[species_id] = coefficient, kcat
    if coefficients and reaction.get("reversible") not in ("false", "0"):
        raise ValueError(f"{reaction_id}: GECKO2 protein-demand arms must be explicitly irreversible")
    return coefficients


def parse_kcats(path: pathlib.Path) -> TurnoverInventory:
    path = pathlib.Path(path)
    source = path.read_bytes()
    sbml = gzip.decompress(source) if source.startswith(b"\x1f\x8b") else source
    root = ET.fromstring(sbml)
    model = root.find(_SBML + "model")
    if model is None:
        raise ValueError("expected an SBML level 3 version 1 model")
    model_id = _label(model.get("id"), "model id")
    _verify_unit_basis(model)
    genes = _index(model.findall(f"{_FBC}listOfGeneProducts/{_FBC}geneProduct"),
                   _FBC + "id", "FBC gene product")
    labels = {key: _label(gene.get(_FBC + "label"), f"{key} gene label") for key, gene in genes.items()}
    species = _index(model.findall(f"{_SBML}listOfSpecies/{_SBML}species"), "id", "species")
    reactions = _index(model.findall(f"{_SBML}listOfReactions/{_SBML}reaction"), "id", "reaction")
    declared = set(enzyme_capacity.ISOPRENOID_BRANCH_GENES)
    if {gene for group in BRANCH_GROUPS.values() for gene in group} != declared:
        raise ValueError("branch groups do not match enzyme_capacity.ISOPRENOID_BRANCH_GENES")
    source_sha = hashlib.sha256(source).hexdigest()
    sbml_sha = hashlib.sha256(sbml).hexdigest()
    rows = []
    with_gpr = with_kcat = 0
    for reaction_id, reaction in reactions.items():
        if reaction_id.startswith(("arm_", "draw_prot_")):
            continue
        association = reaction.find(_FBC + "geneProductAssociation")
        if association is None:
            continue
        refs = {ref.get(_FBC + "geneProduct") for ref in association.iter(_FBC + "geneProductRef")}
        if not refs or not refs <= labels.keys():
            raise ValueError(f"{reaction_id}: empty GPR or unknown FBC gene product")
        with_gpr += 1
        coefficients = _protein_demands(reaction, species)
        with_kcat += bool(coefficients)
        gpr_genes = {labels[ref] for ref in refs}
        selected = gpr_genes & declared
        if not selected:
            continue
        reverse = "_REV" in reaction_id
        for species_id, (coefficient, kcat) in coefficients.items():
            rows.append({
                "reaction_id": reaction_id,
                "declared_gpr_genes": ";".join(sorted(selected)),
                "gpr_genes": ";".join(sorted(gpr_genes)),
                "branch": ";".join(name for name, group in BRANCH_GROUPS.items() if selected.intersection(group)),
                "protein_species_id": species_id,
                "protein_coefficient_h": coefficient,
                "kcat_per_s": kcat,
                "is_reverse_arm": reverse,
                "included_in_forward_span": not reverse,
                "association_basis": "reaction_GPR_membership_not_gene_to_protein_assignment",
                "interpretation": "model_encoded_prior_not_measurement",
                "source_model_sha256": source_sha,
                "source_sbml_sha256": sbml_sha,
                "boundary_basis": BOUNDARY_BASIS,
            })
    if not with_gpr:
        raise ValueError("no non-arm, non-draw GPR reactions found")
    if not rows:
        raise ValueError("no encoded protein turnover for the declared isoprenoid genes")
    rows.sort(key=lambda row: (row["kcat_per_s"], row["declared_gpr_genes"],
                               row["reaction_id"], row["protein_species_id"]))
    return TurnoverInventory(rows, with_gpr, with_kcat, str(path.resolve()),
                             source_sha, sbml_sha, model_id)


def forward_span(rows) -> tuple[dict, dict]:
    forward = [row for row in rows if row["included_in_forward_span"]]
    if not forward:
        raise ValueError("no forward declared-branch protein turnovers found")
    return min(forward, key=lambda row: row["kcat_per_s"]), max(forward, key=lambda row: row["kcat_per_s"])


def sweep_capacities(
        kcats_per_s, mass_fractions=SWEEP_MASS_FRACTIONS, *, enzyme: str = DEFAULT_ENZYME,
        enzyme_molar_mass_g_per_mol: float | None = None,
        total_protein_g_per_gdcw: float = proteome.TOTAL_PROTEIN_G_PER_GDCW,
        growth_rate_per_h: float = DEFAULT_GROWTH_RATE_PER_H,
        content_reference: ContentReference | None = None,
        capacity_reference: CapacityReference | None = None) -> list[dict]:
    kcats = _axis(kcats_per_s, "kcats_per_s")
    fractions = _axis(mass_fractions, "mass_fractions")
    if any(value > 1 for value in fractions):
        raise ValueError("mass_fractions must be fractions of total protein mass in (0, 1], not percentages")
    enzyme = _label(enzyme, "enzyme")
    if enzyme_molar_mass_g_per_mol is None:
        if enzyme != DEFAULT_ENZYME:
            raise ValueError("enzyme_molar_mass_g_per_mol is required for a nondefault enzyme")
        enzyme_molar_mass_g_per_mol = DEFAULT_ENZYME_MOLAR_MASS_G_PER_MOL
    mass = _positive(enzyme_molar_mass_g_per_mol, "enzyme_molar_mass_g_per_mol")
    protein = _positive(total_protein_g_per_gdcw, "total_protein_g_per_gdcw")
    if protein > 1:
        raise ValueError("total_protein_g_per_gdcw cannot exceed one gram per gram of dry biomass")
    growth = _positive(growth_rate_per_h, "growth_rate_per_h")
    reference_content = content_reference.mmol_per_gdcw if content_reference else None
    reference_capacity = capacity_reference.capacity_mmol_per_gdcw if capacity_reference else None
    rows = []
    for kcat in kcats:
        for fraction in fractions:
            content = proteome.enzyme_content_from_mass_fraction(fraction, mass, protein)
            _positive(content, "derived enzyme_mmol_per_gdcw")
            derived = enzyme_capacity.derived_capacity(enzyme, kcat, content)
            _positive(derived.vmax_mmol_per_gdcw_h, "derived vmax_mmol_per_gdcw_h")
            capacity = _positive(derived.capacity_mmol_per_gdcw(growth), "derived capacity_mmol_per_gdcw")
            content_fold = (_positive(capacity / reference_content, "fold_over_reference_content")
                            if reference_content is not None else None)
            capacity_fold = (_positive(capacity / reference_capacity, "fold_over_reference_capacity")
                             if reference_capacity is not None else None)
            rows.append({
                "enzyme": enzyme, "kcat_per_s": kcat, "mass_fraction_of_protein": fraction,
                "enzyme_molar_mass_g_per_mol": mass, "total_protein_g_per_gdcw": protein,
                "growth_rate_per_h": growth, "enzyme_mmol_per_gdcw": content,
                "vmax_mmol_per_gdcw_h": derived.vmax_mmol_per_gdcw_h,
                "capacity_mmol_per_gdcw": capacity,
                "reference_content_source": content_reference.source if content_reference else None,
                "reference_content_mg_per_gdcw": content_reference.content_mg_per_gdcw if content_reference else None,
                "reference_product_molar_mass_g_per_mol": content_reference.product_molar_mass_g_per_mol if content_reference else None,
                "reference_content_mmol_per_gdcw": reference_content,
                "fold_over_reference_content": content_fold,
                "reference_capacity_source": capacity_reference.source if capacity_reference else None,
                "reference_capacity_mmol_per_gdcw": reference_capacity,
                "fold_over_reference_capacity": capacity_fold,
                "reference_basis": "user_declared_comparison_only_not_calibration",
                "interpretation": "conditional_sensitivity_not_product_prediction",
            })
    return rows


def _csv_bytes(rows: list[dict], fields: tuple[str, ...]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def _source_record(path) -> dict:
    path = pathlib.Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _write_outputs(output_dir: pathlib.Path, payloads: dict[str, bytes], sources: set[pathlib.Path]):
    for name in payloads:
        target = output_dir / name
        if target.resolve() in sources:
            raise ValueError(f"refusing to overwrite an input source: {target}")
        if target.exists() and not target.is_file():
            raise ValueError(f"output target is not a regular file: {target}")
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".ecmodel-prior-", dir=output_dir) as temporary:
        staging = pathlib.Path(temporary)
        for name, payload in payloads.items():
            (staging / name).write_bytes(payload)
        for name in payloads:
            (staging / name).replace(output_dir / name)


def _optional_reference(factory, values, name):
    if all(value is None for value in values):
        return None
    if any(value is None for value in values):
        raise ValueError(f"{name}: supply all reference values, units and source label together")
    return factory(*values)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=("Export GECKO2 isoprenoid turnover priors and a conditional capacity sweep. "
                     "Model priors are not enzyme measurements; no calibration or solver is run."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--model", type=pathlib.Path, help="GECKO2 SBML (.xml or gzip); otherwise use YSTWIN_EC_YEAST_GEM / paths.ec_yeast_gem")
    parser.add_argument("--output-dir", type=pathlib.Path, help="destination; otherwise use YSTWIN_OUTPUTS / paths.outputs_dir")
    parser.add_argument("--kcats-per-s", type=float, nargs="+", help="declared sensitivity axis in 1/s; otherwise every unique forward encoded turnover")
    parser.add_argument("--mass-fractions", type=float, nargs="+", default=SWEEP_MASS_FRACTIONS, help="sensitivity axis as MASS fractions of total protein, not molar fractions or percentages")
    parser.add_argument("--enzyme", default=DEFAULT_ENZYME, help="scenario label, not a lookup of measured enzyme properties")
    parser.add_argument("--enzyme-molar-mass-g-per-mol", type=float, help=f"declared enzyme g/mol; default {DEFAULT_ENZYME_MOLAR_MASS_G_PER_MOL:g} only for {DEFAULT_ENZYME}; required for another label")
    parser.add_argument("--total-protein-g-per-gdcw", type=float, default=proteome.TOTAL_PROTEIN_G_PER_GDCW, help="declared total protein g/gDCW")
    parser.add_argument("--growth-rate-per-h", type=float, default=DEFAULT_GROWTH_RATE_PER_H, help="declared mu in 1/h; capacity is vmax/mu, not predicted product content")
    parser.add_argument("--reference-content-mg-per-gdcw", type=float, help="optional comparison content; never fitted or selected against")
    parser.add_argument("--reference-product-molar-mass-g-per-mol", type=float, help="product g/mol for converting the optional reference mg/gDCW to mmol/gDCW")
    parser.add_argument("--reference-content-source", help="source label for the user-declared content comparator; not independently verified")
    parser.add_argument("--reference-capacity-mmol-per-gdcw", type=float, help="optional capacity comparator in mmol/gDCW; never a calibration input")
    parser.add_argument("--reference-capacity-source", help="source label for the user-declared capacity comparator; not independently verified")
    args = parser.parse_args(argv)
    try:
        content_reference = _optional_reference(ContentReference, (
            args.reference_content_source, args.reference_content_mg_per_gdcw,
            args.reference_product_molar_mass_g_per_mol), "content reference")
        capacity_reference = _optional_reference(CapacityReference, (
            args.reference_capacity_source, args.reference_capacity_mmol_per_gdcw), "capacity reference")
        model_path = args.model if args.model is not None else paths.require(
            paths.ec_yeast_gem(), "GECKO2 ecYeastGEM", "YSTWIN_EC_YEAST_GEM")
        inventory = parse_kcats(model_path)
        low, high = forward_span(inventory.rows)
        kcats = args.kcats_per_s
        kcat_basis = "user_declared_sensitivity_axis_not_measurement"
        if kcats is None:
            kcats = sorted({row["kcat_per_s"] for row in inventory.rows if row["included_in_forward_span"]})
            kcat_basis = "forward_model_encoded_turnovers_as_sensitivity_priors"
        sweep = sweep_capacities(
            kcats, args.mass_fractions, enzyme=args.enzyme,
            enzyme_molar_mass_g_per_mol=args.enzyme_molar_mass_g_per_mol,
            total_protein_g_per_gdcw=args.total_protein_g_per_gdcw,
            growth_rate_per_h=args.growth_rate_per_h,
            content_reference=content_reference, capacity_reference=capacity_reference,
        )
        for row in sweep:
            row.update(kcat_basis=kcat_basis, source_model_sha256=inventory.source_model_sha256,
                       source_sbml_sha256=inventory.source_sbml_sha256, boundary_basis=BOUNDARY_BASIS)
        forward_genes = {gene for row in inventory.rows if row["included_in_forward_span"]
                         for gene in row["declared_gpr_genes"].split(";")}
        source_code = [_source_record(path) for path in (__file__, proteome.__file__, enzyme_capacity.__file__)]
        report = {
            "schema_version": 1, "calibration_performed": False, "solver_used": False,
            "boundary_basis": BOUNDARY_BASIS,
            "source_model": {
                "path": inventory.source_path, "sha256": inventory.source_model_sha256,
                "sbml_sha256": inventory.source_sbml_sha256, "model_id": inventory.model_id,
                "sha256_basis": "exact_input_file_bytes", "sbml_sha256_basis": "uncompressed_SBML_bytes",
                "declared_flux_unit_definition": _FLUX_UNIT_ID,
            },
            "source_code": source_code,
            "inventory": {
                "declared_genes": list(enzyme_capacity.ISOPRENOID_BRANCH_GENES),
                "branch_groups": BRANCH_GROUPS,
                "direction_rule": "exclude arm_* and draw_prot_*; report _REV arms but exclude from forward span and default sweep",
                "protein_rule": "all prot_* reactant coefficients except prot_pool*; no first-protein shortcut",
                "association_basis": "reaction_GPR_membership_not_gene_to_protein_assignment",
                "interpretation": "model_encoded_prior_not_measurement",
                "reactions_with_gpr": inventory.reactions_with_gpr,
                "reactions_with_kcat": inventory.reactions_with_kcat,
                "coverage_fraction": inventory.reactions_with_kcat / inventory.reactions_with_gpr,
                "row_count": len(inventory.rows),
                "forward_row_count": sum(row["included_in_forward_span"] for row in inventory.rows),
                "forward_span": {"minimum": low, "maximum": high},
                "genes_without_forward_turnovers": sorted(set(enzyme_capacity.ISOPRENOID_BRANCH_GENES) - forward_genes),
            },
            "sweep": {
                "enzyme": args.enzyme, "kcats_per_s": kcats, "kcat_basis": kcat_basis,
                "mass_fractions_of_protein": args.mass_fractions,
                "enzyme_molar_mass_g_per_mol": sweep[0]["enzyme_molar_mass_g_per_mol"],
                "total_protein_g_per_gdcw": args.total_protein_g_per_gdcw,
                "growth_rate_per_h": args.growth_rate_per_h,
                "input_basis": "declared_sensitivity_scenario_not_strain_measurements",
                "content_reference": asdict(content_reference) if content_reference else None,
                "capacity_reference": asdict(capacity_reference) if capacity_reference else None,
                "reference_basis": "user_declared_comparison_only_not_calibration",
                "interpretation": "conditional_sensitivity_not_product_prediction",
                "row_count": len(sweep),
            },
        }
        payloads = {
            "ecmodel_isoprenoid_kcats.csv": _csv_bytes(inventory.rows, INVENTORY_FIELDS),
            "enzyme_capacity_sweep.csv": _csv_bytes(sweep, SWEEP_FIELDS),
        }
        report["outputs"] = {name: {"sha256": hashlib.sha256(payload).hexdigest()}
                             for name, payload in payloads.items()}
        payloads["ecmodel_isoprenoid_prior.json"] = (json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
        output_dir = args.output_dir if args.output_dir is not None else paths.outputs_dir()
        sources = {pathlib.Path(record["path"]).resolve() for record in source_code}
        sources.add(pathlib.Path(inventory.source_path))
        _write_outputs(output_dir, payloads, sources)
    except (OSError, EOFError, ET.ParseError, ValueError, OverflowError, zlib.error) as failure:
        print(f"ecmodel_isoprenoid_prior: refused: {failure}", file=sys.stderr)
        return 2
    print("Encoded model priors, not measurements; capacity sensitivity, not product prediction.")
    print(f"GPR coverage: {inventory.reactions_with_kcat}/{inventory.reactions_with_gpr} "
          f"= {100 * report['inventory']['coverage_fraction']:.4f}% (non-arm, non-draw reactions)")
    print(f"Forward coefficient span: {low['declared_gpr_genes']} {low['reaction_id']} "
          f"{low['kcat_per_s']:.6g} /s to {high['declared_gpr_genes']} {high['reaction_id']} "
          f"{high['kcat_per_s']:.6g} /s")
    print(f"Inventory: {len(inventory.rows)} protein coefficients; "
          f"{len(inventory.rows) - report['inventory']['forward_row_count']} reverse coefficients excluded from span")
    if report["inventory"]["genes_without_forward_turnovers"]:
        print("No forward turnover for declared genes: " + ", ".join(report["inventory"]["genes_without_forward_turnovers"]))
    print(f"Sweep: {len(sweep)} cells; enzyme {args.enzyme}, "
          f"MW {sweep[0]['enzyme_molar_mass_g_per_mol']:g} g/mol, "
          f"protein {args.total_protein_g_per_gdcw:g} g/gDCW, mu {args.growth_rate_per_h:g} /h")
    for column in ("fold_over_reference_content", "fold_over_reference_capacity"):
        values = [row[column] for row in sweep if row[column] is not None]
        if values:
            print(f"Comparison only, {column}: {min(values):.6g} to {max(values):.6g}; no calibration performed")
    print(f"Wrote two CSVs and ecmodel_isoprenoid_prior.json to {paths.display_path(output_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
