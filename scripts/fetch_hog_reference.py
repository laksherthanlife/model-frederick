from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import urllib.request
from dataclasses import dataclass

import libsbml

from ystwin import paths


@dataclass(frozen=True)
class SourceAsset:
    filename: str
    supplement: str
    kind: str
    purpose: str
    model_id: str | None = None

    @property
    def url(self):
        return f"https://journals.plos.org/ploscompbiol/article/file?id={self.supplement}&type=supplementary"


ASSETS = {
    "model_wt": SourceAsset(
        "model_wt.xml", "10.1371/journal.pcbi.1003084.s023", "sbml",
        "Published fitted wild-type kinetic model; imported structure/parameters are not independent validation.",
        "PetelenzKuehn_osmoadaptation_WT"),
    "model_pfk2627_del": SourceAsset(
        "model_pfk2627_del.xml", "10.1371/journal.pcbi.1003084.s024", "sbml",
        "Published pfk26/27 deletion model; parameters participate in the authors' joint fit.",
        "PetelenzKuehn_osmoadaptation_pfk2627D"),
    "model_gpd1_del": SourceAsset(
        "model_gpd1_del.xml", "10.1371/journal.pcbi.1003084.s027", "sbml",
        "Published gpd1 deletion model; inferred intracellular glycerol is not an independent measurement.",
        "PetelenzKuehn_osmoadaptation_gpd1D"),
    "model_hog1_del": SourceAsset(
        "model_hog1_del.xml", "10.1371/journal.pcbi.1003663.s001", "sbml",
        "Corrected hog1 deletion identity; includes an empirical time-dependent Fps1 reopening term.",
        "PetelenzKuehn_osmoadaptation_hog1D"),
    "model_hog1_att": SourceAsset(
        "model_hog1_att.xml", "10.1371/journal.pcbi.1003663.s002", "sbml",
        "Corrected membrane-tethered Hog1 model identity; not an independent holdout of the published fit.",
        "PetelenzKuehn_osmoadaptation_HOG1att"),
    "model_fps1_open": SourceAsset(
        "model_fps1_open.xml", "10.1371/journal.pcbi.1003663.s003", "sbml",
        "Corrected FPS1-delta1 identity; includes an empirical time-dependent Hog1 deactivation term.",
        "PetelenzKuehn_osmoadaptation_fps1D1"),
    "raw_metabolites": SourceAsset(
        "raw_metabolites.xls", "10.1371/journal.pcbi.1003084.s001", "excel",
        "Raw and processed metabolite measurements; retain source distinctions before fitting."),
    "observations": SourceAsset(
        "observations.xls", "10.1371/journal.pcbi.1003084.s002", "excel",
        "Experimental data as used in the authors' model fitting; not an independent test of the published parameter set."),
    "western_blots": SourceAsset(
        "western_blots.xls", "10.1371/journal.pcbi.1003084.s004", "excel",
        "Western-blot measurements for strains and stress strengths; inspect normalization and fitting overlap."),
    "glycerol_assays": SourceAsset(
        "glycerol_assays.xls", "10.1371/journal.pcbi.1003084.s005", "excel",
        "Intracellular glycerol enzyme-assay measurements across strains and stresses."),
    "methods": SourceAsset(
        "methods.pdf", "10.1371/journal.pcbi.1003084.s030", "pdf",
        "Published methods, measurement processing, model equations and parameter-estimation protocol."),
}


def _download(url):
    request = urllib.request.Request(url, headers={"User-Agent": "ystwin-source-acquisition/0.1"})
    with urllib.request.urlopen(request, timeout=30) as response:
        payload = response.read(20_000_001)
        final_url = response.url
    if len(payload) > 20_000_000:
        raise ValueError("source artifact exceeds the bounded download size")
    return payload, final_url


def _validate(asset, payload):
    if not payload:
        raise ValueError(f"empty source artifact: {asset.filename}")
    metadata = {}
    if asset.kind == "sbml":
        document = libsbml.readSBMLFromString(payload.decode("utf-8-sig"))
        model = document.getModel()
        severe = [document.getError(index).getMessage() for index in range(document.getNumErrors())
                  if document.getError(index).getSeverity() >= libsbml.LIBSBML_SEV_ERROR]
        if severe:
            raise ValueError(f"invalid SBML source {asset.filename}: {severe}")
        if model is None or model.getId() != asset.model_id:
            actual = None if model is None else model.getId()
            raise ValueError(f"SBML identity mismatch for {asset.filename}: {actual!r}")
        metadata = {
            "model_id": model.getId(), "model_name": model.getName(),
            "sbml_level": document.getLevel(), "sbml_version": document.getVersion(),
            "species": model.getNumSpecies(), "reactions": model.getNumReactions(),
            "parameters": model.getNumParameters(), "rules": model.getNumRules(),
        }
    elif asset.kind == "excel":
        if not payload.startswith(bytes.fromhex("d0cf11e0a1b11ae1")) and not payload.startswith(b"PK"):
            raise ValueError(f"source is not an Excel workbook: {asset.filename}")
    elif asset.kind == "pdf" and not payload.startswith(b"%PDF"):
        raise ValueError(f"source is not a PDF: {asset.filename}")
    return metadata


def fetch_assets(output_dir, names, downloader=_download):
    output = pathlib.Path(output_dir).resolve()
    requested = tuple(dict.fromkeys(names))
    unknown = set(requested) - ASSETS.keys()
    if unknown:
        raise ValueError(f"unknown source assets: {sorted(unknown)}")
    manifest_path = output / "sources.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("publication", {}).get("doi") != "10.1371/journal.pcbi.1003084":
            raise ValueError("existing source manifest is for a different publication")
    else:
        manifest = {
            "schema_version": 1,
            "publication": {
                "title": "Quantitative analysis of glycerol accumulation, glycolysis and growth under hyper osmotic stress",
                "authors": "Petelenz-Kurdziel et al.", "year": 2013,
                "doi": "10.1371/journal.pcbi.1003084", "pmid": "23762021", "pmcid": "PMC3677637",
            },
            "correction": {
                "doi": "10.1371/journal.pcbi.1003663", "pmcid": "PMC4036733",
                "issue": "Original captions for mutant SBML supplements S3/S4/S6 were mismatched; use corrected artifacts and internal model identities.",
            },
            "license": {
                "id": "CC-BY-4.0", "url": "https://creativecommons.org/licenses/by/4.0/",
                "evidence": "License element in the primary article's Europe PMC fullTextXML, verified before acquisition.",
            },
            "interpretation": [
                "Published model parameters were fitted to experimental data; agreement with that fitting data is not independent validation.",
                "Original artifacts are retained without numerical edits. Derived observation tables must retain assay, normalization, time and genotype provenance.",
                "Intracellular glycerol retention and extracellular secretion are different observables.",
            ],
            "assets": {},
        }
    manifest["source_assumptions"] = {
        "evidence": "Primary Text S1, sections 2, 3 and 4",
        "published_fit": "A joint parameter set was estimated using wild type and multiple mutants; those mutants are not held out from the published parameters.",
        "optical_density": "The source model supplies a time-dependent OD fit and derives cell density from it; this is an experimental input, not a predicted growth law.",
        "processed_hog1": "Processed Hog1PP values assume a wild-type phosphorylation peak of 90% of total Hog1; absolute activation fractions are therefore not independently measured.",
        "metabolite_normalization": "Intracellular data were processed with a constant reference cell volume; source assignment rules map dynamic concentrations back to that measurement basis.",
        "inferred_glycerol": "Intracellular glycerol for gpd1 deletion was inferred from other pools rather than directly measured; some wild-type repeats are also inferred.",
        "mutant_forcing": "FPS1-delta1 and hog1-deletion variants contain empirical time-dependent terms introduced to reproduce observations; these are not established molecular mechanisms.",
        "trehalose": "Trehalose data were excluded from the final parameter fit because the source model lacks its regulatory dynamics.",
    }
    output.mkdir(parents=True, exist_ok=True)
    for name in requested:
        asset = ASSETS[name]
        target = output / asset.filename
        previous = manifest["assets"].get(name)
        if previous is not None and previous["url"] != asset.url:
            raise ValueError(f"source URL changed for {name}; use a separate source revision")
        if target.exists():
            if previous is None:
                raise ValueError(f"refusing to adopt an unmanifested existing artifact: {target}")
            payload = target.read_bytes()
            final_url = previous["resolved_url"]
        else:
            payload, final_url = downloader(asset.url)
        digest = hashlib.sha256(payload).hexdigest()
        if previous is not None and digest != previous["sha256"]:
            raise ValueError(f"source checksum changed for {name}; existing records are not overwritten")
        details = _validate(asset, payload)
        if not target.exists():
            with target.open("xb") as handle:
                handle.write(payload)
        manifest["assets"][name] = {
            "filename": asset.filename, "url": asset.url, "resolved_url": final_url,
            "sha256": digest, "bytes": len(payload), "kind": asset.kind,
            "purpose": asset.purpose, **details,
        }
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description="Acquire identified primary-source HOG models and measured data")
    parser.add_argument("--output-dir", type=pathlib.Path, default=paths.data_dir() / "hog2013")
    parser.add_argument("--assets", nargs="+", choices=tuple(ASSETS),
                        default=["model_wt", "observations", "western_blots", "methods"])
    args = parser.parse_args(argv)
    manifest = fetch_assets(args.output_dir, args.assets)
    for name in args.assets:
        record = manifest["assets"][name]
        print(f"{name}: {record['bytes']} bytes; sha256 {record['sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
