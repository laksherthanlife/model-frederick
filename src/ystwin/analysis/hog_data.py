from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from io import BytesIO
from math import isfinite
from numbers import Real
from pathlib import Path

import pandas as pd

from .. import paths

__all__ = ["HogDataset", "load_hog_data", "load_native_hog_data"]

_DOI = "10.1371/journal.pcbi.1003084"
_SPLITS = ("fit", "heldout_dose", "heldout_assay", "unused")
_ASSETS = {
    "observations": (
        "observations.xls", "s002", 198656,
        "87854ff96bf22943f1347abdf5636e261ddfcdec041ea78fe498d62d924a967e",
    ),
    "western_blots": (
        "western_blots.xls", "s004", 91136,
        "9fdb4e3a9b35fce549f5635adc17c5e511ea057dc429de7b808601d0d32d5ef6",
    ),
    "methods": (
        "methods.pdf", "s030", 257003,
        "59c57d8816f5101dc5a6f2552df81eedecffe69e8d0e9848ffb9af5710c45344",
    ),
}
_S2_SHEETS = (
    "Readme", "Numbers", "Numbers 2", "glycerol_i", "glycerol_e", "glucose_i", "glucose_e",
    "trehalose_i", "trehalose_e", "ethanol_i", "ethanol_e", "acetate_i", "acetate_e",
    "Hog1PP_ALL", "Gpd1_ALL",
)
_S4_LOW = "wild type pfk26Δ27Δ 0.4M"
_S4_HIGH = "wild type pfk26Δ27Δ 0.8M"
_S4_SHEETS = ("Readme", "wild type gpd2Δ gpd1Δ", _S4_LOW, _S4_HIGH)
_S2_PROCESSING = (
    "experimental data, processed. All concentrations in mol/l, stress (0.4M NaCl) "
    "added at t=3600. Strains as described in main text"
)
_S4_NOT_FITTED = (
    "Data and has not been used for fitting but refers to the Preliminary Dataset "
    "mentioned in the main text."
)
_GENOTYPES = {
    202: "wild_type", 1586: "pfk2627_del", 1064: "gpd1_del", 115: "fps1_open",
    2293: "stl1_del", 444: "hog1_del",
}
_STRAIN_LABELS = (
    (10, 115, "fps1-Δ1"), (11, 202, "wild type (W303)"), (12, 444, "hog1Δ"),
    (13, 1064, "gpd1Δ"), (14, 1586, "pfk26/27Δ"), (15, 2293, "stl1Δ"),
)
_HPLC_SERIES = (
    ("A", "B", "fps1-Δ1", 115, "fps1_open", 0.4, None),
    ("C", "D", "gpd1Δ", 1064, "gpd1_del", 0.4, None),
    ("E", "F", "hog1-att", None, "hog1_att", 0.4, None),
    ("G", "H", "hog1Δ", 444, "hog1_del", 0.4, None),
    ("I", "J", "No Stress2", 202, "wild_type", 0.0, "2"),
    ("K", "L", "pfk26/27Δ", 1586, "pfk2627_del", 0.4, None),
    ("M", "N", "Wild Type4", 202, "wild_type", 0.4, "4"),
)
_PROTEIN_BLOCKS = (
    ("Hog1PP_ALL", "Hog1PP_measured", 21, 26, (202, 1586, 1064, 115, 2293)),
    ("Gpd1_ALL", "Gpd1_measured", 17, 19, (202, 1586, 115, 444, 2293)),
)
_S4_BLOCKS = (
    (_S4_LOW, 0.4, "Measurement1", "C", 1, 1, 2, 14, "D", "E", "wild type (1)", "(1)", 1),
    (_S4_LOW, 0.4, "Measurement1", "C", 1, 1, 2, 14, "D", "F", "wild type (2)", "(2)", 1),
    (_S4_LOW, 0.4, "Measurement2", "A", 20, 20, 21, 33, "C", "D", "wild type", None, 1),
    (_S4_LOW, 0.4, "Measurement3", "A", 40, 41, 42, 54, "C", "D", "wild type", None, 100),
    (_S4_HIGH, 0.8, "Measurement1", "A", 21, 23, 24, 36, "D", "E", "wild type", None, 1),
    (_S4_HIGH, 0.8, "Measurement2", "A", 42, 42, 43, 55, "D", "E", 202, None, 1),
    (_S4_HIGH, 0.8, "Measurement4", "A", 82, 82, 83, 95, "D", "E", "wt", None, 100),
)
_NORMALIZATIONS = {
    "Hog1PP_measured": "literature_molar_scale_assumed_90_percent_wt_peak",
    "Gpd1_measured": "literature_molar_scale_prestress_abundance",
    "glycerol_measured": "hplc_molar_constant_reference_volume",
    "glycerol_e": "hplc_molar_conversion",
}


@dataclass(frozen=True)
class HogDataset:
    observations: pd.DataFrame = field(repr=False)
    metadata: dict = field(repr=False)

    def split(self) -> dict[str, pd.DataFrame]:
        return {
            name: self.observations.loc[self.observations["split"] == name].copy()
            for name in _SPLITS
        }


def _verified_inputs(directory: Path, asset_keys=None) -> tuple[dict[str, bytes], dict]:
    keys = tuple(_ASSETS) if asset_keys is None else tuple(asset_keys)
    if not keys or set(keys) - set(_ASSETS):
        raise ValueError("Source assets must be an explicit nonempty registered selection")
    manifest_bytes = (directory / "sources.json").read_bytes()
    try:
        manifest = json.loads(manifest_bytes)
    except (ValueError, UnicodeError) as exc:
        raise ValueError("Invalid HOG source manifest: sources.json") from exc
    if not isinstance(manifest, dict) or not isinstance(manifest.get("publication"), dict):
        raise ValueError("Invalid HOG source manifest publication")
    if manifest["publication"].get("doi") != _DOI:
        raise ValueError("HOG source manifest publication DOI mismatch")
    assets = manifest.get("assets")
    if not isinstance(assets, dict):
        raise ValueError("Invalid HOG source manifest assets")
    payloads = {}
    sources = {}
    for key in keys:
        filename, supplement, size, digest = _ASSETS[key]
        expected = {
            "filename": filename,
            "sha256": digest,
            "bytes": size,
            "url": (
                "https://journals.plos.org/ploscompbiol/article/file"
                f"?id={_DOI}.{supplement}&type=supplementary"
            ),
        }
        record = assets.get(key)
        if not isinstance(record, dict) or any(record.get(k) != v for k, v in expected.items()):
            raise ValueError(f"HOG source manifest mismatch for {key}")
        payload = (directory / filename).read_bytes()
        if hashlib.sha256(payload).hexdigest() != digest:
            raise ValueError(f"SHA-256 mismatch for {filename}")
        payloads[key] = payload
        sources[key] = expected
    metadata = {
        "schema_version": 1,
        "source_doi": _DOI,
        "sources": sources,
        "source_manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "publication": manifest["publication"],
        "license": manifest.get("license"),
    }
    return payloads, metadata


def _workbook(payload: bytes, expected_sheets: tuple[str, ...], selected: tuple[str, ...]) -> dict:
    with pd.ExcelFile(BytesIO(payload), engine="xlrd") as workbook:
        if tuple(workbook.sheet_names) != expected_sheets:
            raise ValueError(f"Unexpected HOG workbook sheet names: {workbook.sheet_names!r}")
        return {sheet: workbook.parse(sheet, header=None) for sheet in selected}


def _column_index(column: str) -> int:
    index = 0
    for character in column:
        index = index * 26 + ord(character) - ord("A") + 1
    return index - 1


def _cell(frame: pd.DataFrame, sheet: str, column: str, row: int):
    try:
        return frame.iat[row - 1, _column_index(column)]
    except IndexError as exc:
        raise ValueError(f"Missing source cell {sheet}!{column}{row}") from exc


def _expect(frame: pd.DataFrame, sheet: str, column: str, row: int, expected) -> None:
    actual = _cell(frame, sheet, column, row)
    if pd.isna(actual) or actual != expected:
        raise ValueError(f"{sheet}!{column}{row}: expected {expected!r}, got {actual!r}")


def _number(frame: pd.DataFrame, sheet: str, column: str, row: int) -> float:
    value = _cell(frame, sheet, column, row)
    if isinstance(value, bool) or not isinstance(value, Real) or not isfinite(float(value)):
        raise ValueError(f"{sheet}!{column}{row}: expected a finite numeric value, got {value!r}")
    return float(value)


def _context(
    supplement: str, observable: str, experiment: str, block: str,
    strain: int | None, genotype: str, dose: float, replicate: str | None = None,
) -> dict:
    preliminary = supplement == "s004"
    protein = observable in ("Hog1PP_measured", "Gpd1_measured")
    return {
        "source_doi": _DOI,
        "supplement_id": supplement,
        "experiment_id": experiment,
        "measurement_block": block,
        "replicate_label": replicate,
        "strain_id": strain,
        "genotype": genotype,
        "background": "W303" if strain == 202 else None,
        "medium": "YPD",
        "nacl_molar": dose,
        "time_unit": "minutes" if preliminary else "seconds",
        "time_unit_evidence": (
            "caller_supplied_inference_unlabelled_workbook" if preliminary
            else "explicit_workbook_header"
        ),
        "observable_id": observable,
        "unit": "relative_intensity" if preliminary else "mol/L",
        "raw_unit": "relative_intensity" if preliminary else "mol/l",
        "measurement_provenance": (
            "preliminary_western_blot" if preliminary else
            "literature_scaled_western_blot" if protein else "processed_hplc_measurement"
        ),
        "normalization": "reported_peak_scale_1" if preliminary else _NORMALIZATIONS[observable],
        "original_fit_status": (
            "not_used_in_published_parameter_fit" if preliminary else "published_fitting_dataset"
        ),
        "split": "unused",
        "statistic": "mean" if protein and not preliminary else "reported_value",
    }


def _series(
    frame: pd.DataFrame, sheet: str, rows: range, time_column: str, value_column: str,
    context: dict, *, sd_column: str | None = None, divisor: float = 1,
) -> list[dict]:
    records = []
    for row in rows:
        if pd.isna(_cell(frame, sheet, value_column, row)):
            continue
        raw_value = _number(frame, sheet, value_column, row)
        native_time = _number(frame, sheet, time_column, row)
        relative_time = native_time * 60 if context["time_unit"] == "minutes" else native_time - 3600
        sd = float("nan")
        if sd_column is not None and pd.notna(_cell(frame, sheet, sd_column, row)):
            sd = _number(frame, sheet, sd_column, row) / divisor
            if sd < 0:
                raise ValueError(f"{sheet}!{sd_column}{row}: negative standard deviation")
        records.append({
            **context,
            "source_sheet": sheet,
            "source_cell": f"{value_column}{row}",
            "source_time_cell": f"{time_column}{row}",
            "sd_source_cell": f"{sd_column}{row}" if sd_column is not None else None,
            "time_native": native_time,
            "time_relative_s": relative_time,
            "time_model_s": relative_time + 3600,
            "value": raw_value / divisor,
            "raw_value": raw_value,
            "sd": sd,
        })
    return records


def _s2_rows(sheets: dict) -> list[dict]:
    readme = sheets["Readme"]
    for row in (2, 3):
        _expect(readme, "Readme", "H", row, _S2_PROCESSING)
    for row, strain, label in _STRAIN_LABELS:
        _expect(readme, "Readme", "A", row, strain)
        _expect(readme, "Readme", "B", row, label)
    records = []
    for sheet, observable, unit_row, first, strains in _PROTEIN_BLOCKS:
        frame = sheets[sheet]
        if sheet == "Hog1PP_ALL":
            _expect(frame, sheet, "A", 20, "Hog1PP")
            _expect(frame, sheet, "A", 25, "time (min)")
            _expect(frame, sheet, "B", 25, "time (seconds, -60=0)")
        else:
            _expect(frame, sheet, "B", 17, "s")
            _expect(frame, sheet, "B", 18, "time in sec, 0=-60")
        for row in range(first, first + 13):
            seconds = _number(frame, sheet, "B", row)
            _expect(frame, sheet, "A", row, (seconds - 3600) / 60)
        for column, sd_column, strain in zip("CEGIK", "DFHJL", strains):
            _expect(frame, sheet, column, unit_row, "mol/l")
            _expect(frame, sheet, sd_column, unit_row, "mol/l")
            if sheet == "Hog1PP_ALL":
                _expect(frame, sheet, column, 23, strain)
                _expect(frame, sheet, column, 24, "avg")
                _expect(frame, sheet, sd_column, 24, "stdev")
            else:
                _expect(frame, sheet, column, 18, f"avg {strain}")
                _expect(frame, sheet, sd_column, 18, "stdev")
            context = _context(
                "s002", observable, f"s002:{sheet}:{strain}", "molar",
                strain, _GENOTYPES[strain], 0.4,
            )
            if observable == "Hog1PP_measured" and strain == 202:
                context["split"] = "fit"
            records.extend(_series(
                frame, sheet, range(first, first + 13), "B", column, context,
                sd_column=sd_column,
            ))
    for sheet, observable in (("glycerol_i", "glycerol_measured"), ("glycerol_e", "glycerol_e")):
        frame = sheets[sheet]
        for time_column, column, label, strain, genotype, dose, replicate in _HPLC_SERIES:
            _expect(frame, sheet, column, 1, label)
            _expect(frame, sheet, time_column, 2, "time (s)")
            header = "glycerol extra" if sheet == "glycerol_e" else (
                "glycerol_intra" if strain == 1064 else "Glycerol_int"
            )
            _expect(frame, sheet, column, 2, header)
            context = _context(
                "s002", observable, f"s002:HPLC:{label}", "HPLC",
                strain, genotype, dose, replicate,
            )
            records.extend(_series(frame, sheet, range(3, 18), time_column, column, context))
    return records


def _s4_rows(sheets: dict) -> list[dict]:
    _expect(sheets["Readme"], "Readme", "A", 4, _S4_NOT_FITTED)
    _expect(sheets["Readme"], "Readme", "A", 8, 202)
    _expect(sheets["Readme"], "Readme", "B", 8, "wild type (W303)")
    for sheet, dose in ((_S4_LOW, "0.4M NaCl"), (_S4_HIGH, "0.8M NaCl")):
        frame = sheets[sheet]
        _expect(frame, sheet, "A", 1, "Western blot")
        _expect(frame, sheet, "A", 3, "Medium: ")
        _expect(frame, sheet, "A", 4, "YPD")
        _expect(frame, sheet, "A", 6, "Stress: ")
        _expect(frame, sheet, "A", 7, dose)
        _expect(frame, sheet, "A", 12, "wt 202")
    _expect(sheets[_S4_HIGH], _S4_HIGH, "A", 19, "very noisy data")
    records = []
    for (
        sheet, dose, block, label_column, label_row, header_row, first, last,
        time_column, column, header, replicate, divisor,
    ) in _S4_BLOCKS:
        frame = sheets[sheet]
        _expect(frame, sheet, label_column, label_row, f"{block}:")
        _expect(frame, sheet, column, header_row, header)
        context = _context(
            "s004", "Hog1PP_measured", f"s004:{sheet}:{block}", block,
            202, "wild_type", dose, replicate,
        )
        context["split"] = "heldout_dose" if dose == 0.8 else "heldout_assay"
        if divisor == 100:
            context["raw_unit"] = "percent_of_peak"
            context["normalization"] = "reported_peak_scale_100_divided_by_100"
        records.extend(_series(
            frame, sheet, range(first, last + 1), time_column, column, context, divisor=divisor,
        ))
    return records


def _provenance(preliminary_time_unit: str) -> dict:
    return {
        "preliminary_time_unit": {
            "unit": preliminary_time_unit,
            "explicit_workbook_label": False,
            "basis": (
                "S4 time columns are unlabelled. The caller explicitly accepts minutes inferred "
                "from the article's relative-minute convention and sampling grids, not an XLS label."
            ),
            "conversion": "time_relative_s = 60 * time_native; time_model_s = time_relative_s + 3600",
        },
        "s2_time": {
            "unit": "seconds", "stress_onset_s": 3600,
            "evidence": "S2 Readme!H2:H3 and selected time-column headers",
            "conversion": "time_relative_s = time_native - 3600; time_model_s = time_native",
        },
        "protein_units": {
            "reported_unit": "mol/L",
            "direct_absolute_concentration": False,
            "hog1_peak_fraction_assumed": 0.9,
            "model_unit_conversion_applied": False,
            "basis": (
                "Western intensities normalized against WT and literature-scaled using "
                "Ghaemmaghami et al. 2003. Hog1PP assumes a WT peak of 90% of constant total "
                "Hog1; Gpd1 assumes literature prestress abundance. Reported molar values and "
                "SDs are preserved, not equated with model-native protein units."
            ),
            "evidence": "Text S1 p.4; S2 Hog1PP_ALL!O5:O14 and Gpd1_ALL!O15:O24",
        },
        "intracellular_volume_basis": {
            "comparison": "Use model measurement observables, not instantaneous intracellular concentrations.",
            "equation": "C_measured(t) = C_intracellular(t) * Vos(t) / Vos(0)",
            "evidence": "Text S1 p.13, equation 10",
            "source_discrepancy": (
                "Methods describe 50 fL while S2 describes roughly half as accessible and "
                "protein conversion cells use 0.5 * 50 fL. No concentrations are recomputed."
            ),
        },
        "s4_normalization": {
            "unit": "relative_intensity",
            "basis": "Numerical peak scales of 1 or 100; detailed preliminary correction protocol is not supplied.",
            "percent_blocks": [f"{_S4_LOW}:Measurement3", f"{_S4_HIGH}:Measurement4"],
            "baseline_subtraction_applied": False,
            "observation_gain_fitted": False,
            "noise_warning": {"source_cell": f"{_S4_HIGH}!A19", "text": "very noisy data"},
        },
        "published_fit": {
            "s002": {
                "status": "published_fitting_dataset",
                "includes_wild_type_and_mutants": True,
                "evidence": "S2 Readme!C3; Text S1 pp.4,13",
                "scope": "A fitting-data source, not a claim that every imported series entered the final objective.",
            },
            "s004": {
                "used_for_parameter_fitting": False,
                "part_of_model_development": True,
                "evidence": "S4 Readme!A4",
                "source_text": _S4_NOT_FITTED,
                "scope": "Not used in published parameter fitting, but preliminary data were part of model development.",
            },
        },
        "replication": {
            "biological_independence_established": False,
            "basis": (
                "Preserve Measurement groups and source replicate labels without assigning biological "
                "or technical independence. S2 supplies means and SDs without replicate counts. "
                "Assay-specific S2 series IDs do not establish independence between assays."
            ),
        },
        "condition_evidence": {
            "s2": "Primary article Experimental methods; S2 Readme strain IDs and experiment labels; No Stress2 is the unstressed WT control.",
            "s4": "Readme!A8:B8 plus each imported sheet's A4, A7 and A12 and exact WT column headers.",
            "mutant_backgrounds": "Not populated when unspecified; no W303 background imputed for mutant rows.",
        },
        "split_policy": {
            "name": "wt_hog1_0.4M_s2_fit_s4_whole_measurement_holdouts_v1",
            "group_column": "experiment_id",
            "fit": "Only S2 WT202/W303/YPD Hog1PP_measured, 13 rows.",
            "heldout_dose": "All imported S4 WT 0.8 M Measurement1/2/4 groups, 38 rows.",
            "heldout_assay": "All imported S4 WT 0.4 M Measurement1/2/3 groups, 52 rows.",
            "unused": "All other imported S2 observations, diagnostic only.",
            "heldout_fitting_allowed": False,
            "restriction": "No scale, offset or parameter fitting on heldouts; no random timepoint split.",
        },
        "missing_observables": {
            "gpd1mRNA_measured": "Absent from S2/S4; Text S1 p.5 says the fitting series came from Klipp et al. 2005 and was scaled.",
            "relVM": "Absent from S2/S4; optical density is not relative single-cell volume.",
        },
        "not_imported": [
            {
                "supplement_id": "s002", "source_sheet": "Numbers",
                "scope": "Entire sheet", "reason": "Wide-table observations are not duplicated; additional WT/control experiments and non-target metabolites are outside this selection.",
            },
            {
                "supplement_id": "s002", "source_sheet": "Numbers 2",
                "scope": "Entire sheet", "reason": "Convenient glycerol sheets already supply the selected HPLC observations.",
            },
            {
                "supplement_id": "s002", "source_sheet": "Hog1PP_ALL",
                "scope": "Rows 6:18 and 42:54", "reason": "Normalized duplicates and inferred unphosphorylated Hog1 are not independent measurements.",
            },
            {
                "supplement_id": "s002", "source_sheet": "Gpd1_ALL",
                "scope": "Rows 3:15", "reason": "Normalized-intensity duplicate of the imported molar block.",
            },
            {
                "supplement_id": "s002", "source_sheet": "glycerol_i",
                "scope": "D3:D17", "reason": "gpd1 deletion intracellular measurements are absent; inferred values elsewhere are not substituted.",
            },
            *[
                {
                    "supplement_id": "s002", "source_sheet": sheet,
                    "scope": "Entire sheet", "reason": "Non-target metabolite outside the requested observation mapping.",
                }
                for sheet in _S2_SHEETS[5:13]
            ],
            {
                "supplement_id": "s004", "source_sheet": "wild type gpd2Δ gpd1Δ",
                "scope": "Entire sheet", "reason": "YNB and mixed or unspecified backgrounds, including MATalpha WT689, are not WT202/W303/YPD; scaled duplicates are also present.",
            },
            {
                "supplement_id": "s004", "source_sheet": _S4_LOW,
                "scope": "Measurement1 G:H; Measurement2 E:G; Measurement3 E",
                "reason": "Mutant curves are outside the WT202/W303/YPD selection.",
            },
            {
                "supplement_id": "s004", "source_sheet": _S4_HIGH,
                "scope": "Measurement1 F:H; Measurement2 F:H; Measurement3 E:H; Measurement4 F",
                "reason": "Mutant curves are outside the WT selection; Measurement3 E/F are original and renormalized versions of one pfk26 deletion trace, not independent replicates.",
            },
        ],
    }


def load_hog_data(directory=None, *, preliminary_time_unit: str) -> HogDataset:
    if preliminary_time_unit != "minutes":
        raise ValueError("preliminary_time_unit must be explicitly supplied as 'minutes'")
    directory = paths.data_dir() / "hog2013" if directory is None else Path(directory)
    payloads, metadata = _verified_inputs(directory)
    s2 = _workbook(
        payloads["observations"], _S2_SHEETS,
        ("Readme", "Hog1PP_ALL", "Gpd1_ALL", "glycerol_i", "glycerol_e"),
    )
    s4 = _workbook(payloads["western_blots"], _S4_SHEETS, ("Readme", _S4_LOW, _S4_HIGH))
    metadata.update(_provenance(preliminary_time_unit))
    return _dataset(_s2_rows(s2) + _s4_rows(s4), metadata)


def _dataset(rows, metadata):
    observations = pd.DataFrame.from_records(rows)
    observations["strain_id"] = pd.array(observations["strain_id"], dtype="Int64")
    if observations.duplicated(["supplement_id", "source_sheet", "source_cell"]).any():
        raise ValueError("Duplicate HOG source observations")
    if not observations.groupby("experiment_id")["split"].nunique().eq(1).all():
        raise ValueError("HOG experiment groups cross split boundaries")
    metadata["row_counts"] = {
        name: int(observations["split"].eq(name).sum()) for name in _SPLITS
    }
    return HogDataset(observations=observations, metadata=metadata)


def load_native_hog_data(directory=None) -> HogDataset:
    directory = paths.data_dir() / "hog2013" if directory is None else Path(directory)
    payloads, metadata = _verified_inputs(directory, ("observations", "methods"))
    workbook = _workbook(
        payloads["observations"], _S2_SHEETS,
        ("Readme", "Hog1PP_ALL", "Gpd1_ALL", "glycerol_i", "glycerol_e"),
    )
    metadata.update(
        native_only=True,
        excluded_sources={"s004": "not opened by this native-only reader"},
        normalization_policy=dict(_NORMALIZATIONS),
        original_fit_status="Published source priors were fitted using these native data; they are not independent validation of those priors",
        role_policy="Original row roles are preserved; new training selections must be supplied by an explicit protocol",
        time_unit_evidence="S2 model time is explicitly seconds with the salt shock at 3600; relative times subtract that offset",
    )
    return _dataset(_s2_rows(workbook), metadata)
