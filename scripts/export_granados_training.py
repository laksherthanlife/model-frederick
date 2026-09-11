#!/usr/bin/env python3
"""Granados custody tooling with a separate gated development-only workflow.

No command prints response values or response statistics. Mixed source JSON files
are custodial inputs, never learner artifacts. The original strong-protocol
workflow remains non-exporting. Development projections require their own exact
protocol, operational decision and phased release. Network calls have a 30 s curl
deadline. Existing dependencies are used; nothing is installed.
"""
from __future__ import annotations

import argparse
import ast
import base64
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import stat
import re
import subprocess
import sys
from urllib.parse import quote, urljoin
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/native_law_v2/granados"
MANIFEST = ROOT / "data/native_law_v2/granados_manifest.json"
CANDIDATE = ROOT / "data/validation_candidates/granados2018.json"
PROTOCOL = ROOT / "data/native_law_v2/protocol.json"
DEPTHS = {
    "1b38366b-7473-4ddf-9cf9-26e625533686": 1,
    "32cf23e6-ff7c-4b41-8304-f60afca3f4a6": 2,
    "721ef5f2-b858-4257-8cc8-5eea9a8599bb": 3,
    "5b7ea616-fdd5-4978-951e-8cfc3f722eb1": 3,
    "0dcc11fc-dcd1-4983-930e-26265b48a0e9": 1,
    "f0ff4386-a4a5-40a4-b501-42b3e222e698": 1,
    "5546d389-8424-40f3-8620-6be68ca9b16d": 2,
}
TIME_FIELDS = {"times", "origin"}


class CustodyError(Exception):
    """Only explicit, response-free error text belongs in this exception."""


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def relative(path):
    return str(path.relative_to(ROOT))


def load(path):
    return json.loads(path.read_bytes())


def save(path, content, sealed=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    if sealed:
        (BASE / "sealed").chmod(0o700)
        path.parent.chmod(0o700)
    encoded = (json.dumps(content, indent=2, sort_keys=False, allow_nan=False) + "\n").encode()
    path.write_bytes(encoded)
    if sealed:
        path.chmod(0o600)
    return digest(encoded)


def fetch(url, cookie=None):
    if not url.startswith("https://"):
        raise CustodyError("Only explicit HTTPS source links are permitted.")
    command = ["curl", "--fail", "--silent", "--show-error", "--location", "--max-time", "30"]
    if cookie is not None:
        command.extend(["--cookie", cookie])
    completed = subprocess.run(
        command + [url], capture_output=True, check=False, timeout=32,
    )
    if completed.returncode:
        # Never relay remote bodies or untrusted diagnostics through a data path.
        raise CustodyError(f"HTTP acquisition failed (curl exit {completed.returncode}); source URL is recorded in metadata.")
    return completed.stdout


def validate():
    import jsonschema
    manifest = load(MANIFEST)
    manifest_schema = load(BASE / "manifest.schema.json")
    jsonschema.Draft202012Validator.check_schema(manifest_schema)
    validator = jsonschema.Draft202012Validator(manifest_schema, format_checker=jsonschema.FormatChecker())
    validator.validate(manifest)
    if manifest["custody"]["raw_files_allowed_for_training"]:
        raise CustodyError("Raw mixed files may not be allowlisted.")
    # This tool is metadata-only. A valid-looking boolean or arbitrary pathname
    # cannot authorize an export or turn a source file into training data.
    if manifest["custody"]["training_allowlist"] or manifest["custody"]["export_authorized"]:
        raise CustodyError("No response-bearing export is implemented or authorized by this metadata audit.")
    checked_groups = 0
    inventory_path = BASE / "schema_inventory.json"
    if inventory_path.exists():
        inventory = load(inventory_path)
        schema = load(BASE / "inventory.schema.json")
        jsonschema.Draft202012Validator.check_schema(schema)
        inventory_validator = jsonschema.Draft202012Validator(schema, format_checker=jsonschema.FormatChecker())
        inventory_validator.validate(inventory)
        group_ids = [g["group_id"] for g in inventory["groups"]]
        if len(group_ids) != len(set(group_ids)):
            raise CustodyError("Duplicate source group IDs in the metadata inventory.")
        checked_groups = len(group_ids)
        # Regression: response-array payloads cannot be inserted into the public
        # shape objects, even if all other fields are otherwise valid.
        altered = json.loads(json.dumps(inventory))
        altered["groups"][0]["channels"]["GFP"]["fields"]["nucLoc"]["values"] = []
        if inventory_validator.is_valid(altered):
            raise CustodyError("Inventory schema failed the response-payload rejection check.")
        component_path = BASE / "dependency_components.json"
        if component_path.exists():
            members = [gid for c in load(component_path)["components"] for gid in c["group_ids"]]
            if len(members) != len(set(members)) or set(members) != set(group_ids):
                raise CustodyError("Dependency components do not cover each source group exactly once.")
        partition = manifest["proposed_partition"]
        if partition["status"] == "ineligible_no_export":
            if partition["training_group_ids"] or partition["development_group_ids"] or partition["excluded_group_ids"]:
                raise CustodyError("An ineligible candidate must not have outcome-access groups.")
            if set(partition["reserved_group_ids"]) != set(group_ids):
                raise CustodyError("Not all source groups are reserved.")
    checked_hashes = 0
    receipt_path = BASE / "acquisition.json"
    if receipt_path.exists():
        for asset in load(receipt_path)["files"]:
            path = ROOT / asset["path"]
            if digest(path.read_bytes()) != asset["sha256"] or asset["training_allowed"]:
                raise CustodyError("Acquisition integrity or custody assertion failed.")
            if path.stat().st_mode & 0o077:
                raise CustodyError("A custodial raw file has group/other permissions.")
            checked_hashes += 1
    for source in manifest["sources"]:
        if source["local_path"] and source["sha256"]:
            if digest((ROOT / source["local_path"]).read_bytes()) != source["sha256"]:
                raise CustodyError("Primary/repository source evidence checksum changed.")
    if manifest["custody"]["protocol_sha256"] is not None:
        if not PROTOCOL.exists() or digest(PROTOCOL.read_bytes()) != manifest["custody"]["protocol_sha256"]:
            raise CustodyError("The verifier protocol changed after the metadata eligibility audit.")
    print(json.dumps({"manifest_valid": True, "metadata_groups_validated": checked_groups, "raw_source_hashes_verified": checked_hashes, "response_payload_rejection_checked": bool(checked_groups), "export_authorized": False, "protocol_exists": PROTOCOL.exists()}))


def acquire():
    candidate = load(CANDIDATE)
    manifest = load(MANIFEST)
    if manifest["custody"]["training_allowlist"]:
        raise CustodyError("Acquisition is disabled after a training allowlist is published.")
    provenance = BASE / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    catalogue_bytes = fetch(candidate["source"]["item_url"])
    bundle_bytes = fetch(candidate["source"]["bundle_url"])
    catalogue = json.loads(catalogue_bytes)
    bundle = json.loads(bundle_bytes)
    if catalogue["id"] != candidate["source"]["item_id"] or bundle["page"]["totalPages"] != 1:
        raise CustodyError("Unexpected catalogue identity or pagination.")
    (provenance / "repository_item.json").write_bytes(catalogue_bytes)
    (provenance / "repository_bundle.json").write_bytes(bundle_bytes)
    records = {b["id"]: b for b in bundle["_embedded"]["bitstreams"]}
    readme_id = "c779e1c9-0b73-4824-8dae-e5f7fc39f212"
    readme = fetch(records[readme_id]["_links"]["content"]["href"])
    if digest(readme) != candidate["readme"]["sha256"]:
        raise CustodyError("README checksum differs from the verified catalogue.")
    (provenance / "readme_source.txt").write_bytes(readme)
    evidence_paths = {
        "catalogue": (provenance / "repository_item.json", catalogue_bytes),
        "bundle": (provenance / "repository_bundle.json", bundle_bytes),
        "readme": (provenance / "readme_source.txt", readme),
    }
    for source in manifest["sources"]:
        if source["source_id"] in evidence_paths:
            path, content = evidence_paths[source["source_id"]]
            source["local_path"] = relative(path)
            source["sha256"] = digest(content)
    save(MANIFEST, manifest)
    sealed = BASE / "sealed"
    raw = sealed / "raw_mixed"
    raw.mkdir(parents=True, exist_ok=True)
    sealed.chmod(0o700)
    raw.chmod(0o700)
    receipt_path = BASE / "acquisition.json"
    receipt = {"schema_version": 1, "retrieved_at_utc": now(), "response_values_disclosed": False, "files": []}
    for item in candidate["assets"]:
        record = records[item["bitstream_id"]]
        if record["name"] != item["name"] or record["sizeBytes"] != item["bytes"]:
            raise CustodyError("Live bitstream metadata differs from the verified candidate.")
        if record["_links"]["content"]["href"] != item["content_url"]:
            raise CustodyError("Content URL differs from the authoritative API link.")
        path = raw / record["name"]
        content = path.read_bytes() if path.exists() else fetch(record["_links"]["content"]["href"])
        if len(content) != record["sizeBytes"]:
            raise CustodyError("Downloaded bitstream byte count does not match the API.")
        checksum = record["checkSum"]
        if checksum["checkSumAlgorithm"] != "MD5" or hashlib.md5(content).hexdigest() != checksum["value"]:
            raise CustodyError("Downloaded bitstream checksum does not match the API.")
        if not path.exists():
            path.write_bytes(content)
        path.chmod(0o600)
        receipt["files"].append({
            "asset_id": record["id"], "name": record["name"], "source_url": record["_links"]["content"]["href"],
            "path": relative(path), "bytes": len(content), "sha256": digest(content),
            "repository_md5": checksum["value"], "training_allowed": False,
        })
        save(receipt_path, receipt)
        print(json.dumps({"acquired_asset_id": record["id"], "bytes_verified": len(content)}), flush=True)


def pointer(parts):
    return "/" + "/".join(str(p).replace("~", "~0").replace("/", "~1") for p in parts)


def groups(node, depth, path=()):
    if depth == 0:
        if not isinstance(node, dict):
            raise CustodyError("An experiment node is not an object.")
        yield path, node
    else:
        if not isinstance(node, dict):
            raise CustodyError("A grouping level is not an object.")
        for key in sorted(node):
            yield from groups(node[key], depth - 1, path + (key,))


def shape(value):
    # Numeric response values are never used in this operation.
    if value is None:
        return {"kind": "null"}
    if isinstance(value, dict):
        return {"kind": "object", "keys": sorted(value)}
    if isinstance(value, list):
        result = {"kind": "array", "outer_length": len(value)}
        nested = [len(row) for row in value if isinstance(row, list)]
        if nested:
            result["row_lengths"] = sorted(set(nested))
            result["nested_rows"] = len(nested)
        return result
    return {"kind": "scalar", "type": type(value).__name__}


def time_metadata(value):
    # This function is called only for explicitly whitelisted input/time fields.
    def flatten(x):
        if isinstance(x, list):
            for entry in x:
                yield from flatten(entry)
        elif isinstance(x, (int, float)) and not isinstance(x, bool):
            yield x
    numbers = list(flatten(value))
    result = {"shape": shape(value)}
    if numbers:
        result["minimum"] = min(numbers)
        result["maximum"] = max(numbers)
        result["distinct_value_count"] = len(set(numbers))
        if not isinstance(value, list):
            result["value"] = value
        rows = value if isinstance(value, list) and value and isinstance(value[0], list) else [value]
        positive_steps = []
        for row in rows:
            if isinstance(row, list):
                present = [x for x in row if isinstance(x, (int, float)) and not isinstance(x, bool)]
                positive_steps.extend(b - a for a, b in zip(present, present[1:]) if b > a)
        if positive_steps:
            result["positive_step_minimum"] = min(positive_steps)
            result["positive_step_maximum"] = max(positive_steps)
    return result


def metadata_ids(asset_id, path):
    # Semantics of grouping levels come from the repository README, not filename
    # heuristics. Physical dose values are deliberately not inferred from names.
    if asset_id == "1b38366b-7473-4ddf-9cf9-26e625533686":
        return {"factor_id": "sfp1", "stress_id": "gluc", "replicate_id": path[0], "condition_id": None}
    if asset_id == "32cf23e6-ff7c-4b41-8304-f60afca3f4a6":
        return {"factor_id": path[0], "stress_id": path[1], "replicate_id": "rep1", "condition_id": None}
    if asset_id == "721ef5f2-b858-4257-8cc8-5eea9a8599bb":
        return {"factor_id": path[0], "stress_id": path[1], "replicate_id": path[2], "condition_id": None}
    if asset_id == "5b7ea616-fdd5-4978-951e-8cfc3f722eb1":
        return {"factor_id": path[1], "stress_id": path[0], "replicate_id": None, "condition_id": path[2]}
    if asset_id == "0dcc11fc-dcd1-4983-930e-26265b48a0e9":
        return {"factor_id": path[0], "stress_id": "gluc", "replicate_id": None, "condition_id": "nuclear_marker_0.1_percent_glucose"}
    if asset_id == "f0ff4386-a4a5-40a4-b501-42b3e222e698":
        return {"factor_id": path[0], "stress_id": "rich_to_rich", "replicate_id": None, "condition_id": "rich_to_rich"}
    if asset_id == "5546d389-8424-40f3-8620-6be68ca9b16d":
        return {"factor_id": path[0], "stress_id": path[1], "replicate_id": None, "condition_id": None}
    raise CustodyError("No authoritative grouping contract for this asset ID.")


def inspect():
    receipt = load(BASE / "acquisition.json")
    if len(receipt["files"]) != len(DEPTHS):
        raise CustodyError("Acquisition is incomplete; inspect only a complete declared bundle.")
    inventory = {"schema_version": 1, "inspection": "Keys, shapes and input/time metadata only; no response summaries or values.", "physical_dose_values_inferred_from_names": False, "groups": []}
    raw_signatures = defaultdict(list)
    for asset in receipt["files"]:
        source = ROOT / asset["path"]
        content = source.read_bytes()
        if digest(content) != asset["sha256"]:
            raise CustodyError("Custodial source checksum changed.")
        data = json.loads(content)
        for path, experiment in groups(data, DEPTHS[asset["asset_id"]]):
            group_id = asset["asset_id"] + ":" + pointer(path)
            entry = {
                "group_id": group_id, "asset_id": asset["asset_id"], "json_pointer": pointer(path),
                "metadata": metadata_ids(asset["asset_id"], path), "channels": {},
                "biological_unit": "Whole source experiment; within-experiment rows are cells, not independent biological replicates.",
            }
            for channel, measures in sorted(experiment.items()):
                if not isinstance(measures, dict):
                    raise CustodyError("Unexpected channel schema; metadata inspection stopped.")
                entry["channels"][channel] = {
                    "fields": {field: shape(value) for field, value in sorted(measures.items())},
                    "time_metadata": {field: time_metadata(value) for field, value in sorted(measures.items()) if field in TIME_FIELDS},
                }
                # Equality-only duplicate audit. No response moments, amplitude,
                # rankings, individual values, or fitted quantities are computed.
                raw_fields = ("median", "max5", "nucLoc")
                if all(field in measures for field in raw_fields):
                    encoded = json.dumps([measures[field] for field in raw_fields], separators=(",", ":"), allow_nan=False).encode()
                    raw_signatures[(channel, digest(encoded))].append(group_id)
            inventory["groups"].append(entry)
        del data, content
    inventory_path = BASE / "schema_inventory.json"
    save(inventory_path, inventory)
    equality_sets = [sorted(set(ids)) for ids in raw_signatures.values() if len(set(ids)) > 1]
    equality_sets = sorted({tuple(ids) for ids in equality_sets})
    known_aliases = []
    ids = {entry["group_id"] for entry in inventory["groups"]}
    for entry in inventory["groups"]:
        if entry["asset_id"] == "32cf23e6-ff7c-4b41-8304-f60afca3f4a6":
            alias = "721ef5f2-b858-4257-8cc8-5eea9a8599bb:" + entry["json_pointer"] + "/rep1"
            if alias not in ids:
                raise CustodyError("Documented replicate-1 alias is absent from the actual inventory.")
            known_aliases.append([entry["group_id"], alias])
    overlaps = {
        "schema_version": 1,
        "documented_same_experiment_aliases": known_aliases,
        "exact_raw_triplet_equality_groups": [list(ids) for ids in equality_sets],
        "equality_method": "SHA256 of canonicalized [median,max5,nucLoc] arrays in the same channel; equality-only provenance audit, no numerical response summaries. Equal triplets must remain in one partition; nonmatching hashes do not prove independent experiments.",
        "remaining_independence_limit": "No assumed independence from different filenames. Dates/movie IDs and cropped/reprocessed reuse may remain unresolved. Use README replication declarations, primary methods and conservative alias components.",
    }
    save(BASE / "overlap_audit.json", overlaps)
    manifest = load(MANIFEST)
    manifest["status"] = "metadata_schema_verified_eligibility_gated"
    manifest["artifacts"]["schema_inventory"] = relative(inventory_path)
    manifest["artifacts"]["overlap_audit"] = relative(BASE / "overlap_audit.json")
    for group_set in manifest["available_group_sets"]:
        group_set["actual_group_count"] = sum(entry["asset_id"] == group_set["asset_id"] for entry in inventory["groups"])
        group_set["inventory_path"] = relative(inventory_path)
    save(MANIFEST, manifest)
    print(json.dumps({"groups_inspected": len(inventory["groups"]), "inventory_path": relative(inventory_path), "documented_alias_pairs": len(known_aliases), "exact_raw_triplet_equality_sets": len(equality_sets), "response_values_disclosed": False}))


def summarize():
    inventory = load(BASE / "schema_inventory.json")
    summary = {"schema_version": 1, "response_values_or_statistics": False, "assets": []}
    for asset in load(BASE / "acquisition.json")["files"]:
        entries = [g for g in inventory["groups"] if g["asset_id"] == asset["asset_id"]]
        summary["assets"].append({
            "asset_id": asset["asset_id"], "name": asset["name"], "group_count": len(entries),
            "group_paths": [g["json_pointer"] for g in entries],
            "channel_field_variants": sorted({json.dumps({c: sorted(m["fields"]) for c, m in g["channels"].items()}, sort_keys=True) for g in entries}),
        })
    summary["paired_channel_shape_checks"] = [{
        "group_id": g["group_id"],
        "gfp_mcherry_general_time_shapes_equal": g["channels"]["GFP"]["fields"]["nucLoc"] == g["channels"]["mCherry"]["fields"]["nucLoc"] == g["channels"]["general"]["fields"]["times"],
    } for g in inventory["groups"] if g["asset_id"] == "5546d389-8424-40f3-8620-6be68ca9b16d"]
    path = BASE / "inventory_summary.json"
    save(path, summary)
    print(json.dumps({"metadata_summary": relative(path), "response_values_disclosed": False}))


def primary():
    from bs4 import BeautifulSoup
    from pypdf import PdfReader
    provenance = BASE / "provenance"
    private = BASE / "sealed/primary_sources"
    private.mkdir(parents=True, exist_ok=True)
    private.chmod(0o700)
    lookup_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI%3A10.1073%2Fpnas.1716659115&format=json&resultType=core"
    lookup_path = provenance / "europe_pmc_lookup.json"
    lookup_bytes = lookup_path.read_bytes() if lookup_path.exists() else fetch(lookup_url)
    lookup = json.loads(lookup_bytes)
    lookup_path.write_bytes(lookup_bytes)
    records = lookup["resultList"]["result"]
    if len(records) != 1 or records[0]["doi"] != "10.1073/pnas.1716659115":
        raise CustodyError("Primary-article DOI lookup was not unique.")
    # Resolve the exact structured PMCID, never a guessed article/asset name.
    article_url = "https://pmc.ncbi.nlm.nih.gov/articles/" + records[0]["pmcid"] + "/"
    article_path = private / "primary_article.html"
    article_bytes = article_path.read_bytes() if article_path.exists() else fetch(article_url)
    article_path.write_bytes(article_bytes)
    article_path.chmod(0o600)
    soup = BeautifulSoup(article_bytes, "html.parser")
    heading = next((h for h in soup.find_all(["h2", "h3"]) if h.get_text(" ", strip=True) == "Materials and Methods"), None)
    if heading is None:
        raise CustodyError("The primary article lacks the expected Methods heading.")
    section = heading.find_parent("section")
    if section is None:
        raise CustodyError("Methods section could not be isolated safely.")
    excerpt = section.get_text("\n", strip=True)
    methods_path = provenance / "primary_methods.txt"
    methods_path.write_text(excerpt + "\n")
    # This anchor text and URL were first observed in the primary PMC article.
    supplement_anchor = next((a for a in soup.find_all("a", href=True) if a.get_text(strip=True) == "pnas.1716659115.sapp.pdf"), None)
    if supplement_anchor is None:
        raise CustodyError("The article's observed supplement link is absent.")
    supplement_url = urljoin(article_url, supplement_anchor["href"])
    supplement_path = private / "primary_supplement.pdf"
    supplement_bytes = supplement_path.read_bytes() if supplement_path.exists() else fetch(supplement_url)
    if not supplement_bytes.startswith(b"%PDF"):
        # The anonymous PMC download page explicitly supplies this SHA256 proof
        # of work to every browser. Its linked pow/vendor JS was inspected to
        # verify the algorithm; no external JavaScript is executed locally.
        challenge_page = supplement_bytes.decode("utf-8")
        constants = {key: re.search(r'const ' + key + r' = "([^"\n]+)"', challenge_page) for key in ("POW_CHALLENGE", "POW_DIFFICULTY", "POW_COOKIE_NAME")}
        if not all(constants.values()):
            raise CustodyError("The primary supplement link returned neither PDF nor its documented anonymous download challenge.")
        challenge = constants["POW_CHALLENGE"].group(1)
        difficulty = int(constants["POW_DIFFICULTY"].group(1))
        if not 1 <= difficulty <= 4:
            raise CustodyError("Unexpected PMC proof-of-work difficulty.")
        nonce = 0
        while not digest((challenge + str(nonce)).encode()).startswith("0" * difficulty):
            nonce += 1
            if nonce > 2_000_000:
                raise CustodyError("PMC anonymous download proof-of-work budget exceeded.")
        cookie = constants["POW_COOKIE_NAME"].group(1) + "=" + quote(challenge + "," + str(nonce), safe=":")
        supplement_bytes = fetch(supplement_url, cookie=cookie)
    if not supplement_bytes.startswith(b"%PDF"):
        raise CustodyError("The primary supplement download handshake did not return a PDF.")
    supplement_path.write_bytes(supplement_bytes)
    supplement_path.chmod(0o600)
    pdf = PdfReader(io.BytesIO(supplement_bytes))
    terms = ("strain", "temperature", "30", "normaliz", "normalis", "max5", "Nhp6", "replicate", "Table S1", "Materials", "Methods")
    index = {
        "schema_version": 1, "article_url": article_url, "supplement_url": supplement_url,
        "supplement_pages": len(pdf.pages),
        "page_term_locators_only": [{"page": i + 1, "terms": [term for term in terms if term.casefold() in page.extract_text().casefold()]} for i, page in enumerate(pdf.pages)],
        "response_values_disclosed": False,
    }
    save(provenance / "primary_evidence_index.json", index)
    manifest = load(MANIFEST)
    additions = [
        {"source_id": "europe_pmc_lookup", "url": lookup_url, "local_path": relative(lookup_path), "sha256": digest(lookup_bytes), "kind": "lookup_metadata", "verification": "Exact DOI lookup returned PMCID " + records[0]["pmcid"] + "; used structured PMCID to locate the primary archived article."},
        {"source_id": "primary_article", "url": article_url, "local_path": relative(article_path), "sha256": digest(article_bytes), "kind": "primary_article", "verification": "Primary PMC article; only Materials and Methods and supplement links disclosed. Full source retained in sealed primary_sources."},
        {"source_id": "primary_supplement", "url": supplement_url, "local_path": relative(supplement_path), "sha256": digest(supplement_bytes), "kind": "primary_supplement", "verification": "Exact PDF link in the primary article. Source retained in sealed/primary_sources; metadata-only page-term locators recorded. Consult the separate methods excerpt and custody manifest for audited scientific interpretation; measurement-figure values are not extracted."},
    ]
    manifest["sources"] = [s for s in manifest["sources"] if s["source_id"] not in {a["source_id"] for a in additions}] + additions
    audited_excerpt = provenance / "supplement_methods_excerpt.txt"
    manifest["artifacts"]["methods_evidence"] = relative(audited_excerpt if audited_excerpt.exists() else methods_path)
    save(MANIFEST, manifest)
    print(json.dumps({"primary_methods": relative(methods_path), "supplement_index": relative(provenance / "primary_evidence_index.json"), "supplement_pages": len(pdf.pages), "response_values_disclosed": False}))


def method_excerpts():
    from pypdf import PdfReader
    pdf = PdfReader(BASE / "sealed/primary_sources/primary_supplement.pdf")
    # Page locations were selected from the primary supplement's metadata-only
    # index. These are the contents/methods pages, not the measurement figures.
    pages = [1, 3, 4, 5, 6]
    path = BASE / "provenance/supplement_methods_excerpt.txt"
    selected = {page: pdf.pages[page - 1].extract_text() for page in pages}
    start = selected[4].index("The ratio of nuclear and cytoplasmic")
    end = selected[4].index("Therefore, we performed", start)
    selected[4] = selected[4][:start] + "[Published calibration outcome omitted from the public methods excerpt.]\n" + selected[4][end:]
    selected[5] = selected[5].split("2 Estimating ﬁtness penalties")[0]
    selected[6] = selected[6].split("3.2 Estimation of mutual information through decoding")[0]
    text = "\n\n".join(f"SUPPLEMENT PDF PAGE {page}\n" + selected[page] for page in pages)
    path.write_text(text + "\n")
    print(json.dumps({"methods_excerpt": relative(path), "pdf_pages": pages, "response_values_disclosed": False}))


def audit_protocol():
    if not PROTOCOL.exists():
        raise CustodyError("A protocol is required for candidate eligibility auditing.")
    protocol_bytes = PROTOCOL.read_bytes()
    protocol = json.loads(protocol_bytes)
    candidate_id = "granados2018_sfp1_carbon_localization"
    if protocol["protocol_id"] != "native-law-v2-preregistered-01" or candidate_id not in protocol["candidate_ids"]:
        raise CustodyError("This metadata audit is specific to the inspected preregistered Sfp1/carbon candidate.")
    inventory = load(BASE / "schema_inventory.json")
    overlaps = load(BASE / "overlap_audit.json")
    entries = {g["group_id"]: g for g in inventory["groups"]}
    parents = {group_id: group_id for group_id in entries}

    def root(group_id):
        while parents[group_id] != group_id:
            parents[group_id] = parents[parents[group_id]]
            group_id = parents[group_id]
        return group_id

    def union(group_ids):
        roots = sorted({root(group_id) for group_id in group_ids})
        for group_id in roots[1:]:
            parents[group_id] = roots[0]

    for aliases in overlaps["documented_same_experiment_aliases"] + overlaps["exact_raw_triplet_equality_groups"]:
        union(aliases)
    components = defaultdict(list)
    for group_id in sorted(entries):
        components[root(group_id)].append(group_id)
    component_records = [{
        "canonical_component_id": "granados2018:" + digest(json.dumps(members, separators=(",", ":")).encode()),
        "group_ids": members,
        "disposition": "reserved_unallocated",
        "independence_verified": False,
    } for members in components.values()]
    component_path = BASE / "dependency_components.json"
    save(component_path, {
        "schema_version": 1,
        "meaning": "Connected components of documented aliases and exact raw-array equality only; an upper bound on distinct experiments, NOT a claim of independent biological groups.",
        "unresolved_dependencies": "Primary Methods permit five factor strains in one shared-input microfluidic assay. No day, culture, movie, chamber or absolute acquisition-clock identifiers are exported. Reprocessed/cropped reuse and shared batches remain unresolved.",
        "source_group_count": len(entries), "known_alias_component_count": len(component_records),
        "components": component_records,
    })
    main_assets = {"1b38366b-7473-4ddf-9cf9-26e625533686", "32cf23e6-ff7c-4b41-8304-f60afca3f4a6", "5b7ea616-fdd5-4978-951e-8cfc3f722eb1"}
    raw_groups = [g for g in entries.values() if g["asset_id"] in main_assets and g["metadata"]["factor_id"] == "sfp1" and g["metadata"]["stress_id"] == "gluc"]
    support = [g for g in entries.values() if g["asset_id"] == "0dcc11fc-dcd1-4983-930e-26265b48a0e9" and g["metadata"]["factor_id"] == "sfp1"]
    sham = [g for g in entries.values() if g["asset_id"] == "f0ff4386-a4a5-40a4-b501-42b3e222e698" and g["metadata"]["factor_id"] == "sfp1"]
    normalized = [g for g in entries.values() if g["asset_id"] == "721ef5f2-b858-4257-8cc8-5eea9a8599bb" and g["metadata"]["factor_id"] == "sfp1" and g["metadata"]["stress_id"] == "gluc"]
    paired_metadata = load(BASE / "condition_metadata.json")["two_colour_pair_channels"]
    paired_sfp1 = [g for g in entries.values() if g["asset_id"] == "5546d389-8424-40f3-8620-6be68ca9b16d" and "sfp1" in (paired_metadata[g["metadata"]["factor_id"]]["GFP"], paired_metadata[g["metadata"]["factor_id"]]["mCherry"])]
    main_upper_bound = len({root(g["group_id"]) for g in raw_groups})
    minimum_main = sum(protocol["selection"][field] for field in ("minimum_train_independent_groups", "minimum_development_independent_groups", "minimum_final_independent_groups"))
    audit = {
        "schema_version": 1, "candidate_id": candidate_id,
        "protocol_path": relative(PROTOCOL), "protocol_sha256": digest(protocol_bytes),
        "result": "ineligible_under_current_protocol",
        "selection_used_measurement_values": False,
        "raw_matched_main_group_ids": [g["group_id"] for g in raw_groups],
        "raw_main_group_count_before_alias_merge": len(raw_groups),
        "raw_main_distinct_experiment_upper_bound": main_upper_bound,
        "required_train_plus_development_plus_final_minimum": minimum_main,
        "main_count_gate_satisfied_even_at_upper_bound": main_upper_bound >= minimum_main,
        "matched_marker_group_ids": [g["group_id"] for g in support],
        "matched_marker_group_upper_bound": len(support),
        "required_independent_support_groups": protocol["selection"]["minimum_support_independent_groups"],
        "support_count_gate_satisfied_even_at_upper_bound": len(support) >= protocol["selection"]["minimum_support_independent_groups"],
        "matched_sham_group_ids": [g["group_id"] for g in sham],
        "matched_sham_group_upper_bound": len(sham),
        "required_independent_sham_groups": 2,
        "sham_requirement_locator": "protocol.json controls.negative_biology (preregistered-01): at least two independent sham groups",
        "sham_count_gate_satisfied_even_at_upper_bound": len(sham) >= 2,
        "matched_two_colour_group_count": len(paired_sfp1),
        "two_colour_basis": "Actual pairs are Dot6/Msn2, Hog1/Msn2, Mig1/Msn2 and Yap1/Msn2; no Sfp1 pair is deposited.",
        "normalized_groups_not_raw_compatible": [g["group_id"] for g in normalized],
        "normalization_exclusion": "SI section 3.1 uses each experiment's own pre-stress 20-frame population mean and SD. Published nucLocNorm is not raw nucLoc and cannot be inverted from the deposited normalized-only replicate file. That file also lacks actual times. No raw-compatible additional Sfp1/carbon replicate is established from it.",
        "other_unresolved_gates": [
            "Distinct cultures/days for the six fig1 replicates are documented, but complete shared-movie/day/culture dependence components across all files are not identifiable from exported fields.",
            "No exported absolute acquisition clocks establish the protocol's required real crossed input/clock strata.",
            "No per-observation standard-uncertainty or independent calibration-uncertainty arrays are present in these JSON schemas.",
            "Marker localization uses an independent nuclear segmentation channel but GFP intensities are not an independent absolute detector/activity calibration; only one Sfp1 marker experiment is available and its added Nhp6A tag changes the strain context.",
            "External registration, role-separated authority, and a complete ancestral-exposure audit have not been established by this local metadata work."
        ],
        "decision": "Do not form a train/development/final partition, fit, select, calibrate, or export response-bearing projections under this protocol. All source groups remain reserved. Do not replace Sfp1 with another factor or relax criteria without a new independently registered protocol before response access.",
    }
    if all(audit[key] for key in ("main_count_gate_satisfied_even_at_upper_bound", "support_count_gate_satisfied_even_at_upper_bound", "sham_count_gate_satisfied_even_at_upper_bound")):
        raise CustodyError("Count gates changed; a new manual metadata/independence audit is required, not automatic export.")
    audit_path = BASE / "eligibility_audit.json"
    save(audit_path, audit)
    manifest = load(MANIFEST)
    manifest["status"] = "metadata_schema_verified_ineligible_current_protocol"
    manifest["custody"]["protocol_sha256"] = digest(protocol_bytes)
    manifest["custody"]["export_authorized"] = False
    manifest["custody"]["training_allowlist"] = []
    partition = manifest["proposed_partition"]
    partition["status"] = "ineligible_no_export"
    partition["training_group_ids"] = []
    partition["development_group_ids"] = []
    partition["reserved_group_ids"] = sorted(entries)
    partition["excluded_group_ids"] = []
    partition["proposal"] = audit["decision"]
    manifest["artifacts"]["eligibility_audit"] = relative(audit_path)
    manifest["artifacts"]["dependency_components"] = relative(component_path)
    save(MANIFEST, manifest)
    print(json.dumps({"candidate_id": candidate_id, "eligibility": audit["result"], "raw_main_upper_bound": main_upper_bound, "minimum_main_groups": minimum_main, "marker_group_upper_bound": len(support), "sham_group_upper_bound": len(sham), "training_allowlist": [], "response_values_disclosed": False}))


def source_links():
    """Expose only literal links in primary sources, never article results."""
    from bs4 import BeautifulSoup
    from pypdf import PdfReader
    article = BeautifulSoup((BASE / "sealed/primary_sources/primary_article.html").read_bytes(), "html.parser")
    article_links = sorted({a["href"] for a in article.find_all("a", href=True) if a["href"].startswith(("http://", "https://"))})
    pdf = PdfReader(BASE / "sealed/primary_sources/primary_supplement.pdf")
    supplement_links = []
    for number, page in enumerate(pdf.pages, start=1):
        urls = set(re.findall(r'https?://[^\s<>]+', page.extract_text()))
        for annotation in page.get("/Annots", []):
            action = annotation.get_object().get("/A", {})
            if action.get("/URI"):
                urls.add(str(action["/URI"]))
        if urls:
            supplement_links.append({"page": number, "urls": sorted(urls)})
    path = BASE / "provenance/source_links.json"
    save(path, {"schema_version": 1, "article_links": article_links, "supplement_links": supplement_links, "source_result_values_disclosed": False})
    print(json.dumps({"source_links": relative(path), "response_values_disclosed": False}))


def source_code():
    """Acquire observed, immutable GitHub blobs via gh; never execute them."""
    # URLs below were observed in GitHub API trees reached from the supplement
    # page-14 mi-by-decoding link and README's DISCO/.gitmodules chain.
    files = [
        ("swainlab/mi-by-decoding", "matlab/DataProcessingFunctions/normMeanOffset.m", "ed14a79a2fd7f7a35a9b1d3b69f0674ea2ca1e33", "normMeanOffset.m"),
        ("swainlab/mi-by-decoding", "matlab/DataProcessingFunctions/filterCellsByNucLoc.m", "7da9364e393bdda1aa8fdc4cd9908cc4a6796435", "filterCellsByNucLoc.m"),
        ("swainlab/mi-by-decoding", "matlab/ExampleScript/example_script.m", "953cf2063ab3044a4fb1a515c0221138ab9bfb0b", "example_script.m"),
        ("swainlab/mi-by-decoding", "matlab/DataProcessingFunctions/calculateInfoList.m", "5a002b5298825d36e738914385f3b0df834df63a", "calculateInfoList.m"),
        ("swainlab/mi-by-decoding", "matlab/DataProcessingFunctions/visualizeExper.m", "dcc0528b443a1d9800bccf5bd365f89cdcef2a84", "visualizeExper.m"),
        ("swainlab/mi-by-decoding", "matlab/README.md", "c07703adc20ad45ae5768f438209ccb45d450b2a", "mi_matlab_readme_source.txt"),
        ("swainlab/mi-by-decoding", "R/MIdecoding/man/YeastStressTypeResponse.Rd", "6e2473fea13fcde0bebb989dce400b14f746ef34", "YeastStressTypeResponse.Rd"),
        ("pswain/GeneralMatlabFunctions", "formatCellInf.m", "6d2cfa0bc3cc5eacb262e0c72e3a386085dd188e", "formatCellInf.m"),
        ("pswain/GeneralMatlabFunctions", "exportJSONCellInf.m", "2ebd31988c039a74e9a039b25a84b04863849056", "exportJSONCellInf.m"),
        ("swainlab/mi-by-decoding", "matlab/external/kakearney-boundedline-pkg-8179f9a/Inpaint_nans/inpaint_nans.m", "2460b512c670354cf221aef327b1bbae1da99757", "inpaint_nans.m"),
        ("swainlab/mi-by-decoding", "matlab/external/jsonlab-1.5.mltbx", "b03c48af1018b45f06a93ff9c0761ca39e536909", "jsonlab-1.5.mltbx"),
    ]
    directory = BASE / "provenance/source_code"
    directory.mkdir(parents=True, exist_ok=True)
    receipt = {"schema_version": 1, "executed_upstream_code": False, "measurement_files_downloaded_from_github": False, "files": []}
    for repository, upstream_path, blob_sha, local_name in files:
        endpoint = f"repos/{repository}/git/blobs/{blob_sha}"
        path = directory / local_name
        if path.exists():
            content = path.read_bytes()
        else:
            result = subprocess.run(["gh", "api", endpoint], capture_output=True, timeout=30, check=False)
            if result.returncode:
                raise CustodyError("GitHub source-code acquisition failed; no remote response body was logged.")
            response = json.loads(result.stdout)
            if response["sha"] != blob_sha or response["encoding"] != "base64":
                raise CustodyError("Unexpected GitHub source blob identity or encoding.")
            content = base64.b64decode(response["content"])
        if hashlib.sha1(b"blob " + str(len(content)).encode() + b"\0" + content).hexdigest() != blob_sha:
            raise CustodyError("Source-code Git blob integrity check failed.")
        path.write_bytes(content)
        receipt["files"].append({"repository": repository, "upstream_path": upstream_path, "api_url": "https://api.github.com/" + endpoint, "git_blob_sha": blob_sha, "sha256": digest(content), "local_path": relative(path)})
    archive = zipfile.ZipFile(directory / "jsonlab-1.5.mltbx")
    savejson_members = [name for name in archive.namelist() if name.endswith("/savejson.m") or name == "savejson.m"]
    if len(savejson_members) != 1:
        raise CustodyError("Could not identify a unique savejson source in the linked archive.")
    savejson_source = archive.read(savejson_members[0])
    savejson_path = directory / "jsonlab_savejson.m"
    savejson_path.write_bytes(savejson_source)
    receipt["jsonlab_installed"] = False
    receipt["extracted_static_source"] = {"archive_git_blob_sha": "b03c48af1018b45f06a93ff9c0761ca39e536909", "member": savejson_members[0], "local_path": relative(savejson_path), "sha256": digest(savejson_source)}
    lf_path = BASE / "provenance/source_code_lf/inpaint_nans.m"
    lf_path.parent.mkdir(parents=True, exist_ok=True)
    lf_source = (directory / "inpaint_nans.m").read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    lf_path.write_bytes(lf_source)
    receipt["newline_normalized_source"] = {"original_git_blob_sha": "2460b512c670354cf221aef327b1bbae1da99757", "change": "CR line endings to LF only", "local_path": relative(lf_path), "sha256": digest(lf_source)}
    save(BASE / "provenance/source_code_acquisition.json", receipt)
    print(json.dumps({"upstream_source_files_acquired": len(files), "executed_upstream_code": False, "installed_packages": False, "response_values_disclosed": False}))


def identity_checks():
    """Run only the fixed custodial source-normalization identity plan."""
    import jsonschema
    output = BASE / "identity_checks.json"
    if output.exists():
        raise CustodyError("Identity replay already recorded; inspect the existing record rather than searching for a passing representation.")
    plan_path = BASE / "identity_check_plan.json"
    plan = load(plan_path)
    manifest = load(MANIFEST)
    if manifest["status"] != "metadata_schema_verified_ineligible_current_protocol" or manifest["custody"]["training_allowlist"]:
        raise CustodyError("Identity replay requires the original ineligible, fully reserved custody state.")
    if manifest["custody"]["protocol_sha256"] != plan["original_protocol_sha256"]:
        raise CustodyError("Identity plan is not bound to the original protocol audit.")
    source_receipt = BASE / "provenance/source_code_acquisition.json"
    code = load(source_receipt)
    for record in code["files"]:
        if digest((ROOT / record["local_path"]).read_bytes()) != record["sha256"]:
            raise CustodyError("Pinned upstream source changed before identity replay.")
    normalized = code["newline_normalized_source"]
    if digest((ROOT / normalized["local_path"]).read_bytes()) != normalized["sha256"]:
        raise CustodyError("Newline-only source copy changed before identity replay.")
    raw_hashes = {}
    for record in load(BASE / "acquisition.json")["files"]:
        if record["asset_id"] in ("1b38366b-7473-4ddf-9cf9-26e625533686", "721ef5f2-b858-4257-8cc8-5eea9a8599bb"):
            if digest((ROOT / record["path"]).read_bytes()) != record["sha256"]:
                raise CustodyError("A source bitstream changed before identity replay.")
            raw_hashes[record["asset_id"]] = record["sha256"]
    # Explicitly record this newly authorized custodial computation. These
    # constants are never exposed to the learner or used for fitting/scoring.
    manifest["custody"]["measurement_statistics_computed"] = True
    manifest["custody"]["measurement_statistics_scope"] = "Baseline means/SDs are evaluated privately only inside the source-defined identity replays requested in the metadata follow-up. No response moments, constants or intermediate projections are released; no models are fitted/scored."
    manifest["custody"]["measurement_statistics_disclosed"] = False
    save(MANIFEST, manifest)
    access_path = BASE / "access_audit.json"
    access = load(access_path)
    access["measurement_bitstream_access"]["response_aggregates_amplitudes_or_fits_computed"] = True
    access["measurement_bitstream_access"]["aggregate_computation_scope"] = manifest["custody"]["measurement_statistics_scope"]
    access["measurement_bitstream_access"]["identity_replay_plan_sha256"] = digest(plan_path.read_bytes())
    save(access_path, access)
    base_escaped = str(BASE).replace("'", "''")
    expression = f"addpath('{base_escaped}'); identity_replay('{base_escaped}');"
    result = subprocess.run(
        ["octave", "--no-init-file", "--no-site-file", "--no-history", "--no-gui", "--quiet", "--eval", expression],
        capture_output=True, timeout=180, check=False,
        env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
    )
    if result.returncode:
        raise CustodyError("Custodial identity replay failed; no Octave output or response values were released.")
    report = json.loads(result.stdout)
    schema = load(BASE / "identity_checks.schema.json")
    validator = jsonschema.Draft202012Validator(schema)
    validator.validate(report)
    observed_cases = {(c["source_group_id"], c["target_group_id"], c["variant_id"]) for c in report["cases"]}
    expected_cases = {(source, target, variant["id"]) for source in plan["source_groups"] for target in plan["target_groups"] for variant in plan["fixed_replays"]}
    if observed_cases != expected_cases:
        raise CustodyError("Identity replay did not cover exactly the fixed source/target/representation plan.")
    report["plan_path"] = relative(plan_path)
    report["plan_sha256"] = digest(plan_path.read_bytes())
    report["source_code_acquisition_sha256"] = digest(source_receipt.read_bytes())
    report["runner_sha256"] = digest((BASE / "identity_replay.m").read_bytes())
    report["raw_source_sha256"] = raw_hashes

    def matched(case):
        return case["exact_double_equality"] or case["exact_source_decimal_token_equality"]

    report["known_duplicate_control_reproduced_by"] = sorted({c["variant_id"] for c in report["cases"] if c["known_duplicate_control"] and matched(c)})
    report["rep2_exact_alias_established"] = any(matched(c) for c in report["cases"] if c["target_group_id"] == plan["target_groups"][1])
    report["conclusion"] = "exact_rep2_reuse_found_keep_inseparable" if report["rep2_exact_alias_established"] else "identity_unresolved_nonmatch_does_not_establish_independence"
    report["original_candidate_remains_ineligible"] = True
    validator.validate(report)
    save(output, report)
    manifest["artifacts"]["identity_checks"] = relative(output)
    save(MANIFEST, manifest)
    print(json.dumps({"identity_metadata": relative(output), "rep2_exact_alias_established": report["rep2_exact_alias_established"], "known_duplicate_control_reproduced_by": report["known_duplicate_control_reproduced_by"], "response_values_or_statistics_released": False}))


def timing_review():
    """Inspect measured TIME fields only, not fluorescence responses."""
    from statistics import mean
    receipt = load(BASE / "acquisition.json")
    report = {
        "schema_version": 1,
        "scope": "Time/input metadata only. Nominal axes are not measured timestamps.",
        "origin_source_evidence": "example_script.m:51-70 uses deposited origin directly as a 1-based MATLAB column index. formatCellInf.m:100-117 indexes times(:,origin), subtracts its mean and remaps the origin index after time filtering.",
        "raw_time_checks": [],
        "normalized_sfp1_nominal_axes": [],
        "marker_same_experiment_structure": None,
    }
    for record in receipt["files"]:
        if record["asset_id"] == "721ef5f2-b858-4257-8cc8-5eea9a8599bb":
            # Only shapes/origin from the existing metadata inventory are needed.
            continue
        content = (ROOT / record["path"]).read_bytes()
        if digest(content) != record["sha256"]:
            raise CustodyError("Source checksum changed before timing metadata review.")
        data = json.loads(content)
        for path, experiment in groups(data, DEPTHS[record["asset_id"]]):
            channel = "general" if "general" in experiment else "GFP"
            timing = experiment[channel]
            times, origin = timing["times"], timing["origin"]
            columns = len(times[0])
            column_means = [mean(row[j] for row in times if row[j] is not None) for j in range(columns)]
            nearest = min(range(columns), key=lambda j: abs(column_means[j])) + 1
            report["raw_time_checks"].append({
                "group_id": record["asset_id"] + ":" + pointer(path),
                "time_channel": channel, "origin_stored": origin,
                "columns": columns, "nearest_zero_mean_time_column_1based": nearest,
                "origin_equals_that_column": origin == nearest,
                "mean_measured_relative_minutes_at_origin": column_means[origin - 1],
            })
    for entry in load(BASE / "schema_inventory.json")["groups"]:
        if entry["asset_id"] == "721ef5f2-b858-4257-8cc8-5eea9a8599bb" and entry["metadata"]["factor_id"] == "sfp1" and entry["metadata"]["stress_id"] == "gluc":
            origin = entry["channels"]["GFP"]["time_metadata"]["origin"]["value"]
            lengths = entry["channels"]["GFP"]["fields"]["nucLocNorm"]["row_lengths"]
            if len(lengths) != 1:
                raise CustodyError("Cannot annotate a single nominal grid for ragged normalized data.")
            report["normalized_sfp1_nominal_axes"].append({
                "group_id": entry["group_id"], "origin_frame_1based": origin, "column_count": lengths[0],
                "formula_minutes": "t_nominal(j) = 2.5 * (j - origin), j is 1-based",
                "first_nominal_minute": 2.5 * (1 - origin), "last_nominal_minute": 2.5 * (lengths[0] - origin),
                "measured_timestamps_reconstructed": False,
                "limitations": "Assumes consecutively sampled columns retaining the source origin convention. Actual per-position timing, jitter, dropped/cropped frames and exact resampling provenance are not restored by a nominal grid.",
            })
        if entry["group_id"] == "0dcc11fc-dcd1-4983-930e-26265b48a0e9:/sfp1":
            fields = entry["channels"]["GFP"]["fields"]
            selected = {name: fields[name] for name in ("nucLoc", "median", "max5", "nucInnerMedian", "cytInnerMedian", "nucCytRatio", "times")}
            report["marker_same_experiment_structure"] = {
                "group_id": entry["group_id"], "channel": "GFP", "fields_and_shapes": selected,
                "all_selected_shapes_equal": all(value == selected["times"] for value in selected.values()),
                "separate_raw_mcherry_intensity_channel_exported": "mCherry" in entry["channels"],
                "pairing_basis": "README nuclear-marker description explicitly adds mask-derived measures to this same experiment's GFP measurements; all arrays use the common cell rows and time columns. This is within-experiment pairing, not a join to marker-free fig1/fig2 cells.",
                "biological_source_group_count": 1,
            }
    path = BASE / "timing_review.json"
    save(path, report)
    manifest = load(MANIFEST)
    manifest["artifacts"]["timing_review"] = relative(path)
    save(MANIFEST, manifest)
    print(json.dumps({"timing_metadata": relative(path), "raw_groups_checked": len(report["raw_time_checks"]), "all_origins_match_1based_zero_time_columns": all(r["origin_equals_that_column"] for r in report["raw_time_checks"]), "response_values_or_statistics_released": False}))


def validate_review():
    """Validate response-free review artifacts and attach their provenance."""
    import jsonschema
    for data_name, schema_name in (
        ("metadata_supplement.json", "metadata_supplement.schema.json"),
        ("identity_checks.json", "identity_checks.schema.json"),
        ("timing_review.json", "timing_review.schema.json"),
    ):
        schema = load(BASE / schema_name)
        jsonschema.Draft202012Validator.check_schema(schema)
        jsonschema.Draft202012Validator(schema).validate(load(BASE / data_name))
    supplement = load(BASE / "metadata_supplement.json")
    if len({q["question_id"] for q in supplement["questions"]}) != 5:
        raise CustodyError("The source-question supplement must answer all five distinct questions.")
    identity = load(BASE / "identity_checks.json")
    bindings = {
        "plan_sha256": BASE / "identity_check_plan.json",
        "source_code_acquisition_sha256": BASE / "provenance/source_code_acquisition.json",
        "runner_sha256": BASE / "identity_replay.m",
    }
    for key, path in bindings.items():
        if identity[key] != digest(path.read_bytes()):
            raise CustodyError("An identity replay provenance binding changed after execution.")
    invalid = json.loads(json.dumps(identity))
    invalid["cases"][0]["normalization_mean"] = 0
    if jsonschema.Draft202012Validator(load(BASE / "identity_checks.schema.json")).is_valid(invalid):
        raise CustodyError("Identity metadata schema admits a response-derived statistic.")
    manifest = load(MANIFEST)
    if manifest["status"] != "metadata_schema_verified_ineligible_current_protocol" or manifest["custody"]["export_authorized"] or manifest["custody"]["training_allowlist"]:
        raise CustodyError("Metadata follow-up must preserve ineligibility and the empty training allowlist.")
    if manifest["custody"]["protocol_sha256"] != supplement["original_protocol_sha256"]:
        raise CustodyError("Supplement is not bound to the original protocol audit.")
    source_path = BASE / "provenance/source_code_acquisition.json"
    source = {
        "source_id": "author_processing_code",
        "url": "https://github.com/swainlab/mi-by-decoding",
        "local_path": relative(source_path),
        "sha256": digest(source_path.read_bytes()),
        "kind": "repository_metadata",
        "verification": "Supplement page 14 links mi-by-decoding; README links DISCO and its GeneralMatlabFunctions submodule. Immutable API blobs were acquired with gh, source hashes checked, and the named normalization helpers used only in the fixed custodial identity replay. No MI/classifier/plotting or model-performance code was executed; no package was installed.",
    }
    manifest["sources"] = [s for s in manifest["sources"] if s["source_id"] != source["source_id"]] + [source]
    manifest["artifacts"]["metadata_supplement"] = relative(BASE / "metadata_supplement.json")
    save(MANIFEST, manifest)
    jsonschema.Draft202012Validator(load(BASE / "manifest.schema.json")).validate(manifest)
    print(json.dumps({"review_schemas_valid": True, "identity_provenance_bindings_valid": True, "response_statistic_injection_rejected": True, "original_candidate_remains_ineligible": True, "response_export_authorized": False}))


# Separate, non-authorizing operational development workflow. These functions
# never change the original protocol, manifests, acquisition or custody records.
DEV_ROOT = BASE / "development"
DEV_PROTOCOL = ROOT / "data/native_law_v2/development_protocol.json"
DEV_PROTOCOL_ID = "granados-sfp1-published-cohort-development-01"
DEV_PROJECTION_ROOT = BASE / "development_exports" / DEV_PROTOCOL_ID
DEV_EXECUTOR_SOURCE = ROOT / "src/ystwin/analysis/native_population_development.py"
DEV_REVIEWED_PROTOCOL_SHA = "0f53f40699ebd0db6df62f00e0a0e6bfc5b56b374ec84a6e36846a8488e6dad6"
DEV_SOURCE_ID = "1b38366b-7473-4ddf-9cf9-26e625533686"
DEV_SOURCE_SHA = "1ec5a48955492aade95386d48f15109397a19acb60ef820d3c692ca59d9efc87"
DEV_SOURCE_BYTES = 4178902
DEV_ROWS = {"rep1": 55, "rep2": 136, "rep3": 202, "rep4": 195, "rep5": 157}
DEV_TRAIN = tuple(f"{DEV_SOURCE_ID}:/rep{i}" for i in (1, 2, 3))
DEV_INPUTS = tuple(f"{DEV_SOURCE_ID}:/rep{i}" for i in (4, 5))
DEV_PREFIX = tuple(range(28, 48))
DEV_RESPONSE = tuple(range(50, 74))


def dev_strict_json(content):
    def object_pairs(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise CustodyError("Duplicate JSON key rejected; source content was not logged.")
            result[key] = value
        return result

    def finite_float(token):
        value = float(token)
        if not math.isfinite(value):
            raise CustodyError("Nonfinite JSON number rejected; no value was logged.")
        return value

    def invalid_constant(_):
        raise CustodyError("Nonfinite JSON token rejected; no value was logged.")

    return json.loads(content, object_pairs_hook=object_pairs, parse_float=finite_float, parse_constant=invalid_constant)


def dev_regular_file(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file() or not stat.S_ISREG(path.stat().st_mode):
        raise CustodyError("Development artifact must be an existing regular nonsymlink file.")
    for parent in path.parents:
        if parent.is_symlink():
            raise CustodyError("Development artifact path traverses a symlink.")
    return path


def dev_save(path, value, private=False, raw_bytes=False, resume_identical=False):
    path = Path(path)
    if path.is_symlink() or any(parent.is_symlink() for parent in path.parents):
        raise CustodyError("Development output path traverses a symlink.")
    path = path.resolve()
    if not any(path.is_relative_to(root.resolve()) for root in (DEV_ROOT, DEV_PROJECTION_ROOT)):
        raise CustodyError("Development writer is restricted to custody metadata and the protocol's exact projection root.")
    content = value if raw_bytes else (json.dumps(value, indent=2, allow_nan=False) + "\n").encode()
    if path.exists():
        dev_regular_file(path)
        if not resume_identical or path.read_bytes() != content:
            raise CustodyError("Existing development artifact differs or cannot be overwritten; retain the prior record.")
        return {"path": relative(path), "sha256": digest(content), "bytes": len(content)}
    path.parent.mkdir(parents=True, exist_ok=True)
    if private:
        path.parent.chmod(0o700)
    with path.open("xb") as stream:
        stream.write(content)
    path.chmod(0o600)
    return {"path": relative(path), "sha256": digest(content), "bytes": len(content)}


def dev_artifact(path, role):
    path = dev_regular_file(path)
    content = path.read_bytes()
    return {"path": relative(path), "sha256": digest(content), "bytes": len(content), "role": role}


def dev_contract(protocol, definition):
    import jsonschema
    embedded = protocol["export_contract"]["schemas"]
    schema = {"$schema": embedded["$schema"], "$ref": "#/$defs/" + definition, "$defs": embedded["$defs"]}
    jsonschema.Draft202012Validator.check_schema(schema)
    return jsonschema.Draft202012Validator(schema)


def dev_protocol_checked():
    content = dev_regular_file(DEV_PROTOCOL).read_bytes()
    if digest(content) != DEV_REVIEWED_PROTOCOL_SHA:
        raise CustodyError("Development protocol bytes changed after independent review; no automatic amendment or release.")
    protocol = dev_strict_json(content)
    if protocol.get("protocol_id") != DEV_PROTOCOL_ID or protocol.get("document_type") != "non_authorizing_operational_development_protocol":
        raise CustodyError("Not the independently reviewed development-only protocol.")
    source = protocol["source"]
    if (source["asset_id"], source["source_sha256"], source["source_bytes"], source["custodian_only_path"]) != (
        DEV_SOURCE_ID, DEV_SOURCE_SHA, DEV_SOURCE_BYTES,
        "data/native_law_v2/granados/sealed/raw_mixed/fig1_sfp1_replicates.json",
    ):
        raise CustodyError("Development source binding differs from the reviewed pinned source.")
    if source["field_allowlist"] != ["GFP/max5", "GFP/median", "general/times", "general/origin"]:
        raise CustodyError("Development source-field contract changed.")
    partition = protocol["partition"]
    if tuple(g["group_id"] for g in partition["train"]) != DEV_TRAIN or tuple(g["group_id"] for g in partition["development"]) != DEV_INPUTS:
        raise CustodyError("Development groups changed; no substitute groups are permitted.")
    for group in partition["train"] + partition["development"]:
        if group["source_pointer"] != group["group_id"].split(":", 1)[1] or DEV_ROWS[group["source_pointer"][1:]] != group["source_rows"]:
            raise CustodyError("Development row/group metadata is inconsistent.")
    if partition["source_frame_count"] != 97 or partition["source_origin_1based"] != 49 or partition["final_groups"]:
        raise CustodyError("Unexpected source frame/origin or a final-test claim.")
    for name, offsets, count, columns in (("prefix", (-21, -2), 20, [28, 47]), ("response_coordinates", (1, 24), 24, [50, 73])):
        section = protocol[name]
        if (section["relative_column_start"], section["relative_column_end"]) != offsets or section["frame_count"] != count or section["expected_absolute_columns_1based"] != columns:
            raise CustodyError("Development window contract changed.")
    if protocol["prefix"]["minimum_prefix_cells"] != 2:
        raise CustodyError("Development prefix-coverage criterion changed.")
    return protocol, content


def dev_selected_source(protocol):
    source_path = dev_regular_file(ROOT / protocol["source"]["custodian_only_path"])
    content = source_path.read_bytes()
    if len(content) != DEV_SOURCE_BYTES or digest(content) != DEV_SOURCE_SHA:
        raise CustodyError("Pinned development source checksum/length mismatch.")
    container = dev_strict_json(content)
    selected = {}
    # Never inspect response arrays for rep6 or any excluded field. The mixed
    # container is integrity-read by the custodian, never passed to a learner.
    for rep, rows in DEV_ROWS.items():
        experiment = container[rep]
        group = {
            "max5": experiment["GFP"]["max5"],
            "median": experiment["GFP"]["median"],
            "times": experiment["general"]["times"],
            "origin": experiment["general"]["origin"],
        }
        if type(group["origin"]) is not int or group["origin"] != 49:
            raise CustodyError("Approved source origin differs from the reviewed frame index.")
        for field in ("max5", "median", "times"):
            matrix = group[field]
            if not isinstance(matrix, list) or len(matrix) != rows or any(not isinstance(row, list) or len(row) != 97 for row in matrix):
                raise CustodyError("Approved source rows/frames are not aligned as declared.")
            for row in matrix:
                for frame in DEV_PREFIX + DEV_RESPONSE:
                    value = row[frame - 1]
                    if value is not None and (isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value)):
                        raise CustodyError("An approved source scalar is neither a finite number nor null.")
        selected[rep] = group
    return selected


def dev_entry_reasons(max5, median, time, prefix):
    reasons = []
    for field, value in (("max5", max5), ("median", median), ("time", time)):
        if value is None:
            reasons.append("source_null_" + field)
        elif isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            reasons.append("invalid_nonfinite_" + field)

    def numeric(value):
        return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value)

    if numeric(max5) and max5 < 0:
        reasons.append("negative_max5")
    if numeric(median) and median <= 0:
        reasons.append("nonpositive_median")
    if numeric(time) and ((prefix and time >= 0) or (not prefix and time <= 0)):
        reasons.append("nonnegative_prefix_time" if prefix else "nonpositive_response_time")
    return reasons


def dev_mean(values):
    return math.fsum(values) / len(values) if values else None


def dev_sd(values):
    if len(values) < 2:
        return None
    centre = dev_mean(values)
    return math.sqrt(math.fsum((x - centre) ** 2 for x in values) / (len(values) - 1))


def dev_aggregate_group(rep, group):
    cohort, exclusions = [], []
    for row in range(DEV_ROWS[rep]):
        failures = []
        for frame in DEV_PREFIX:
            reasons = dev_entry_reasons(group["max5"][row][frame - 1], group["median"][row][frame - 1], group["times"][row][frame - 1], True)
            if reasons:
                failures.append({"frame_1based": frame, "reasons": reasons})
        if failures:
            exclusions.append({"cell_index": row, "failures": failures})
        else:
            cohort.append(row)
    mu, dispersion = None, None
    if len(cohort) >= 2:
        by_frame = [[group["max5"][row][frame - 1] / group["median"][row][frame - 1] for row in cohort] for frame in DEV_PREFIX]
        mu = dev_mean([dev_mean(values) for values in by_frame])
        dispersion = dev_mean([dev_sd(values) for values in by_frame])
        if not math.isfinite(mu) or not math.isfinite(dispersion):
            raise CustodyError("An approved prefix ratio/summary cannot be represented finitely.")
    frames = []
    for frame in DEV_RESPONSE:
        # Coordinates depend on the fixed prefix cohort and TIME ONLY, never
        # this frame's fluorescence validity or future valid-cell count.
        times = [group["times"][row][frame - 1] for row in cohort]
        positive_times = [t for t in times if t is not None and math.isfinite(t) and t > 0]
        coordinate = dev_mean(positive_times)
        ratios, missing, invalid = [], [], []
        for row in cohort:
            x, m, t = (group[field][row][frame - 1] for field in ("max5", "median", "times"))
            reasons = dev_entry_reasons(x, m, t, False)
            if reasons:
                record = {"cell_index": row, "reasons": reasons}
                (missing if any(value is None for value in (x, m, t)) else invalid).append(record)
            else:
                ratio = x / m
                if not math.isfinite(ratio):
                    raise CustodyError("An approved response ratio cannot be represented finitely.")
                ratios.append(ratio)
        frames.append({
            "frame_1based": frame, "time_min": coordinate, "observed_mean": dev_mean(ratios),
            "prefix_cells": len(cohort), "valid_cells": len(ratios), "source_missing_cells": len(missing),
            "invalid_cells": len(invalid), "cell_sample_sd": dev_sd(ratios),
            "source_missing_records": missing, "invalid_records": invalid,
            "sampling_under_resolved": len(ratios) < 2,
        })
    return {"group_id": f"{DEV_SOURCE_ID}:/{rep}", "origin_1based": 49, "source_rows": DEV_ROWS[rep], "eligible_cell_indices": cohort, "prefix_exclusions": exclusions, "mu": mu, "s": dispersion, "frames": frames}


def dev_source_audit():
    protocol, protocol_bytes = dev_protocol_checked()
    if (DEV_ROOT / "source_binding_audit.json").exists():
        raise CustodyError("A development source audit already exists; validate that record rather than overwriting it.")
    selected = dev_selected_source(protocol)
    aggregates = {rep: dev_aggregate_group(rep, group) for rep, group in selected.items()}
    public_groups, private_groups = [], []
    for rep, aggregate in aggregates.items():
        prefix_ok = len(aggregate["eligible_cell_indices"]) >= 2
        coordinates_ok = all(frame["time_min"] is not None for frame in aggregate["frames"])
        means_ok = all(frame["observed_mean"] is not None for frame in aggregate["frames"])
        public_groups.append({
            "group_id": aggregate["group_id"], "source_rows": aggregate["source_rows"], "source_frames": 97,
            "origin_1based": 49, "prefix_cells": len(aggregate["eligible_cell_indices"]),
            "excluded_prefix_cells": len(aggregate["prefix_exclusions"]),
            "prefix_minimum_satisfied": prefix_ok, "all_response_coordinates_available": coordinates_ok,
            "planned_response_coverage_available_custodian_checked": means_ok,
            "response_values_released": False,
        })
        # Retain all missingness/attrition reasons, but no derived response mean
        # or intensity is written before the operational release prerequisites.
        private_groups.append({
            "group_id": aggregate["group_id"], "eligible_cell_indices": aggregate["eligible_cell_indices"],
            "prefix_exclusions": aggregate["prefix_exclusions"],
            "response_coverage": [{key: value for key, value in frame.items() if key not in ("observed_mean", "cell_sample_sd", "time_min")} for frame in aggregate["frames"]],
        })
    snapshots = DEV_ROOT / "contracts"
    protocol_ref = dev_save(snapshots / "development_protocol.json", protocol_bytes, raw_bytes=True)
    protocol_ref["role"] = "protocol_snapshot"
    embedded = protocol["export_contract"]["schemas"]
    dev_save(snapshots / "protocol_schemas.json", embedded)
    for definition in ("DevelopmentReleaseRequest", "DevelopmentCustodyDecision", "Projection", "PrefixSummary", "DevelopmentPredictionFreeze", "DevelopmentBlocker"):
        dev_save(snapshots / (definition + ".schema.json"), {"$schema": embedded["$schema"], "$ref": "#/$defs/" + definition, "$defs": embedded["$defs"]})
    audit = {
        "document_type": "development_source_binding_audit_not_biological_validation",
        "protocol_id": DEV_PROTOCOL_ID, "protocol_sha256": digest(protocol_bytes),
        "source_asset_id": DEV_SOURCE_ID, "source_sha256": DEV_SOURCE_SHA, "source_bytes": DEV_SOURCE_BYTES,
        "field_allowlist": protocol["source"]["field_allowlist"],
        "prefix_frames_1based": list(DEV_PREFIX), "response_frames_1based": list(DEV_RESPONSE),
        "group_checks": public_groups,
        "source_cohort_limitation": protocol["estimand"]["sampling_qualification"],
        "ratio_rule": "Mean of individual max5/median ratios; no supplied nucLoc, Z scoring, interpolation, smoothing or post-response cohort selection.",
        "timing_rule": "Response coordinates use finite positive source times of the fixed complete-prefix cohort without consulting response intensities.",
        "groups_outside_five_ids_evaluated_or_exported": False,
        "learner_value_release_performed": False, "strong_claim_authorized": False,
    }
    audit_ref = dev_save(DEV_ROOT / "source_binding_audit.json", audit)
    audit_ref["role"] = "source_binding_audit"
    dev_save(DEV_ROOT / "private/source_coverage.json", {"document_type": "custodian_only_development_coverage_no_response_values", "protocol_id": DEV_PROTOCOL_ID, "groups": private_groups}, private=True)
    # Source eligibility is not itself permission. Metadata and projection
    # roots are distinct, and both are explicitly constrained by dev_save.
    operational_issues = [
        "A matching typed release request, real executor/settings freeze and independent operational custody decision are required before the separate development-input export. Source eligibility alone grants no value access.",
    ]
    source_failures = [g["group_id"] for g in public_groups if not g["prefix_minimum_satisfied"] or not g["all_response_coordinates_available"] or not g["planned_response_coverage_available_custodian_checked"]]
    blocker = {
        "document_type": "development_blocker_no_biological_verdict", "protocol_id": DEV_PROTOCOL_ID,
        "phase": "operational_approval", "reason_code": "custody_not_approved" if not source_failures else "timing_or_coverage_unusable",
        "reason": " ".join(operational_issues) + (" Selected-source coverage also failed; no replacement is permitted." if source_failures else " The exact selected source groups satisfy the reviewed prefix/timing/response-availability checks; this is not a biological validation result."),
        "affected_group_ids": source_failures or list(DEV_TRAIN + DEV_INPUTS),
        "available_artifacts": [protocol_ref, audit_ref],
        "operational_event_reference": "granados-development-source-review-before-value-release",
        "strong_claim_authorized": False, "further_exports_permitted": False,
    }
    dev_contract(protocol, "DevelopmentBlocker").validate(blocker)
    dev_save(DEV_ROOT / "operational_blocker.json", blocker)
    dev_save(DEV_ROOT / "operational_events.json", {"document_type": "local_development_custody_events_not_external_registry", "events": [{"event_id": blocker["operational_event_reference"], "recorded_at_utc": now(), "operator": "Granados observation/data custodian subagent; distinct from any model fitting/selection operator", "action": "Source audit and typed-contract extraction only; no response projection or operational export approval issued.", "protocol_sha256": digest(protocol_bytes), "source_sha256": DEV_SOURCE_SHA, "source_binding_audit_sha256": audit_ref["sha256"]}]})
    print(json.dumps({"development_protocol_sha256": digest(protocol_bytes), "source_binding_audit": audit_ref, "group_coverage": public_groups, "operational_approval_pending": True, "learner_allowlist": [], "response_values_disclosed": False}))


def dev_verify_ref(reference):
    path = dev_regular_file(ROOT / reference["path"])
    if path.resolve() != (ROOT / reference["path"]).absolute():
        raise CustodyError("Artifact reference is not a canonical nonsymlink path.")
    content = path.read_bytes()
    if len(content) != reference["bytes"] or digest(content) != reference["sha256"]:
        raise CustodyError("Development artifact digest/length mismatch.")
    return content


def dev_executor_literals(content):
    """Read static executor declarations through AST, not by running models."""
    values = {}
    permitted_calls = {"tuple": tuple, "list": list, "range": range}

    def literal(node):
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name) and node.id in values:
            return values[node.id]
        if isinstance(node, ast.Dict):
            return {literal(k): literal(v) for k, v in zip(node.keys, node.values)}
        if isinstance(node, (ast.Tuple, ast.List)):
            items = [literal(item) for item in node.elts]
            return tuple(items) if isinstance(node, ast.Tuple) else items
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in permitted_calls and not node.keywords:
            return permitted_calls[node.func.id](*(literal(arg) for arg in node.args))
        raise CustodyError("Executor declaration is not statically auditable; request its explicit pre-value freeze.")

    wanted = {"PROTOCOL_ID", "ASSET_ID", "SOURCE_SHA256", "QUANTITY", "PARAMETERS", "MODEL_IDS", "SEEDS", "PREFIX_FRAMES", "RESPONSE_FRAMES", "SETTINGS"}
    for node in ast.parse(content).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id in wanted:
            values[node.targets[0].id] = literal(node.value)
    if set(values) != wanted:
        raise CustodyError("Executor is missing required frozen declarations.")
    return values


def dev_freeze_executor(protocol):
    import platform
    import numpy as np
    import scipy
    code_path = DEV_ROOT / "phase_b/frozen_executor/native_population_development.py"
    environment_path = DEV_ROOT / "phase_b/executor_environment.json"
    if environment_path.exists():
        environment = dev_strict_json(dev_regular_file(environment_path).read_bytes())
        content = dev_verify_ref(environment["frozen_implementation"])
        if environment["versions"] != {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}:
            raise CustodyError("Executor environment changed before release; preserve the previous freeze.")
    else:
        content = dev_regular_file(DEV_EXECUTOR_SOURCE).read_bytes()
        code_ref = dev_save(code_path, content, raw_bytes=True, resume_identical=True)
        code_ref["role"] = "execution_code"
        declarations = dev_executor_literals(content)
        environment = {
            "document_type": "pre_value_executor_implementation_freeze_not_a_fit",
            "protocol_id": DEV_PROTOCOL_ID, "protocol_sha256": DEV_REVIEWED_PROTOCOL_SHA,
            "original_implementation_path": relative(DEV_EXECUTOR_SOURCE),
            "original_implementation_sha256_at_freeze": digest(content),
            "frozen_implementation": code_ref,
            "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
            "python_executable": sys.executable,
            "optimizer_settings": declarations["SETTINGS"],
            "canonical_model_ids": list(declarations["MODEL_IDS"]),
            "model_code_executed": False, "model_fitting_or_scoring_performed": False,
            "use_policy": "Execute this byte-bound implementation with these settings/versions. The editable original path is provenance only; it is not permission to change the frozen recipe after data release.",
        }
    declarations = dev_executor_literals(content)
    if (declarations["PROTOCOL_ID"], declarations["ASSET_ID"], declarations["SOURCE_SHA256"], declarations["QUANTITY"]) != (DEV_PROTOCOL_ID, DEV_SOURCE_ID, DEV_SOURCE_SHA, protocol["estimand"]["quantity"]):
        raise CustodyError("Executor source/estimand declarations differ from the reviewed protocol.")
    models = protocol["model_recipe"]["canonical_families"]
    if tuple(declarations["MODEL_IDS"]) != tuple(m["model_id"] for m in models) or any(tuple(declarations["PARAMETERS"][m["model_id"]]) != tuple(m["parameters"]) for m in models):
        raise CustodyError("Executor canonical family/parameter declarations differ from the protocol.")
    if tuple(declarations["PREFIX_FRAMES"]) != DEV_PREFIX or tuple(declarations["RESPONSE_FRAMES"]) != DEV_RESPONSE:
        raise CustodyError("Executor has different source windows.")
    budget = protocol["model_recipe"]["budget"]
    settings = declarations["SETTINGS"]
    if list(declarations["SEEDS"]) != budget["seeds"] or settings["max_nfev"] != budget["maximum_loss_evaluations_per_start"] or settings["starts"] != budget["maximum_starts_per_parameterized_family"]:
        raise CustodyError("Executor budget/seed declarations differ from the protocol.")
    if environment["optimizer_settings"] != settings:
        raise CustodyError("Executor settings disagree with the pre-value freeze.")
    environment_ref = dev_save(environment_path, environment, resume_identical=True)
    environment_ref["role"] = "execution_code"
    return [dev_artifact(code_path, "execution_code"), environment_ref]


def dev_source_cell(rep, group, split, window, cell, frame):
    parent = f"/{rep}"
    index = frame - 1
    return {
        "group_id": f"{DEV_SOURCE_ID}:{parent}", "split": split, "window": window,
        "cell_index": cell, "frame_1based": frame,
        "relative_time_min": group["times"][cell][index],
        "max5": group["max5"][cell][index], "median": group["median"][cell][index],
        "time_pointer": f"{parent}/general/times/{cell}/{index}",
        "max5_pointer": f"{parent}/GFP/max5/{cell}/{index}",
        "median_pointer": f"{parent}/GFP/median/{cell}/{index}",
    }


def dev_source_time(rep, group, cell, frame):
    return {"group_id": f"{DEV_SOURCE_ID}:/{rep}", "cell_index": cell, "frame_1based": frame,
            "relative_time_min": group["times"][cell][frame - 1],
            "time_pointer": f"/{rep}/general/times/{cell}/{frame - 1}"}


def dev_projection_payload(selected, phase):
    if phase == "training_values":
        ids, split, windows = DEV_TRAIN, "train", ("prefix", "response")
    elif phase == "development_inputs":
        ids, split, windows = DEV_INPUTS, "development", ("prefix",)
    else:
        raise CustodyError("No development-response projection is permitted in this release phase.")
    cells, times = [], []
    for group_id in ids:
        rep = group_id.split(":", 1)[1][1:]
        group = selected[rep]
        for window in windows:
            frames = DEV_PREFIX if window == "prefix" else DEV_RESPONSE
            for cell in range(DEV_ROWS[rep]):
                for frame in frames:
                    cells.append(dev_source_cell(rep, group, split, window, cell, frame))
        if phase == "development_inputs":
            for cell in range(DEV_ROWS[rep]):
                for frame in DEV_RESPONSE:
                    times.append(dev_source_time(rep, group, cell, frame))
    return {"group_ids": list(ids), "source_cells": cells, "source_times": times}


def dev_same_scalar(left, right):
    if left is None or right is None:
        return left is None and right is None
    return type(left) is not bool and type(right) is not bool and isinstance(left, (int, float)) and isinstance(right, (int, float)) and math.isfinite(left) and math.isfinite(right) and left == right


def dev_verify_payload(payload, selected, phase):
    """Direct scalar/locator replay; no summaries or fitted values are trusted."""
    if phase not in ("training_values", "development_inputs"):
        raise CustodyError("Only phase-b payloads may be verified/released by this entry point.")
    ids = DEV_TRAIN if phase == "training_values" else DEV_INPUTS
    if payload["group_ids"] != list(ids):
        raise CustodyError("Projection groups differ from the fixed phase partition.")
    # Build only the prescribed INDEX inventory here, independently of the
    # payload generator. Resolve each value directly from its source matrix.
    cell_inventory, time_inventory = [], []
    split = "train" if phase == "training_values" else "development"
    for group_id in ids:
        rep = group_id.split(":", 1)[1][1:]
        windows = ("prefix", "response") if split == "train" else ("prefix",)
        for window in windows:
            for cell in range(DEV_ROWS[rep]):
                for frame in (DEV_PREFIX if window == "prefix" else DEV_RESPONSE):
                    cell_inventory.append((group_id, rep, split, window, cell, frame))
        if split == "development":
            for cell in range(DEV_ROWS[rep]):
                for frame in DEV_RESPONSE:
                    time_inventory.append((group_id, rep, split, "response", cell, frame))
    scalar_checks = 0
    null_counts = {"max5": 0, "median": 0, "relative_time_min": 0}
    cell_keys = {"group_id", "split", "window", "cell_index", "frame_1based", "relative_time_min", "max5", "median", "time_pointer", "max5_pointer", "median_pointer"}
    time_keys = {"group_id", "cell_index", "frame_1based", "relative_time_min", "time_pointer"}
    for key, inventory in (("source_cells", cell_inventory), ("source_times", time_inventory)):
        if len(payload[key]) != len(inventory):
            raise CustodyError("Projection does not retain every prescribed source row/frame.")
        for actual, (group_id, rep, split, window, cell, frame) in zip(payload[key], inventory):
            if set(actual) != (cell_keys if key == "source_cells" else time_keys):
                raise CustodyError("Projection contains missing or unapproved fields.")
            if actual["group_id"] != group_id or type(actual["cell_index"]) is not int or actual["cell_index"] != cell or type(actual["frame_1based"]) is not int or actual["frame_1based"] != frame:
                raise CustodyError("Projection source row/frame/group is reordered or inconsistent.")
            if key == "source_cells" and (actual["split"] != split or actual["window"] != window):
                raise CustodyError("Projection split/window does not match its source indices.")
            fields = (("relative_time_min", "times", "general", "time_pointer"),)
            if key == "source_cells":
                fields += (("max5", "max5", "GFP", "max5_pointer"), ("median", "median", "GFP", "median_pointer"))
            for output_field, source_field, channel, pointer_field in fields:
                locator = f"/{rep}/{channel}/{source_field}/{cell}/{frame - 1}"
                if actual[pointer_field] != locator:
                    raise CustodyError("Projection scalar has an incorrect authoritative JSON pointer.")
                source_value = selected[rep][source_field][cell][frame - 1]
                if not dev_same_scalar(actual[output_field], source_value):
                    raise CustodyError("Projection scalar differs from its authoritative source location.")
                scalar_checks += 1
                null_counts[output_field] += source_value is None
    return {"source_cell_records": len(cell_inventory), "source_time_records": len(time_inventory), "source_scalar_checks": scalar_checks, "preserved_null_counts": null_counts}


def dev_saved_audit_checked(protocol):
    blocker = dev_strict_json(dev_regular_file(DEV_ROOT / "operational_blocker.json").read_bytes())
    dev_contract(protocol, "DevelopmentBlocker").validate(blocker)
    audit_reference = next(ref for ref in blocker["available_artifacts"] if ref["role"] == "source_binding_audit")
    audit = dev_strict_json(dev_verify_ref(audit_reference))
    if audit["protocol_sha256"] != DEV_REVIEWED_PROTOCOL_SHA or audit["source_sha256"] != DEV_SOURCE_SHA:
        raise CustodyError("Saved development source audit is not bound to the reviewed source/protocol.")
    if [group["group_id"] for group in audit["group_checks"]] != list(DEV_TRAIN + DEV_INPUTS):
        raise CustodyError("Saved audit groups differ from the exact release request.")
    if any(not group["prefix_minimum_satisfied"] or not group["all_response_coordinates_available"] or not group["planned_response_coverage_available_custodian_checked"] for group in audit["group_checks"]):
        raise CustodyError("A saved source/coverage failure remains blocking; no source replacement is permitted.")
    return audit_reference, audit


def dev_phase_b_export():
    """Issue an actual limited custody decision and the two phase-b projections."""
    protocol, protocol_bytes = dev_protocol_checked()
    if (DEV_ROOT / "release_manifest.json").exists():
        return dev_validate_release()
    original_paths = [PROTOCOL, MANIFEST, BASE / "acquisition.json", BASE / "access_audit.json", BASE / "eligibility_audit.json"]
    original_hashes = {relative(path): digest(path.read_bytes()) for path in original_paths}
    audit_ref, prior_audit = dev_saved_audit_checked(protocol)
    executor_refs = dev_freeze_executor(protocol)
    exporter_ref = dev_artifact(Path(__file__), "execution_code")
    request = {
        "document_type": "development_release_request", "protocol_id": DEV_PROTOCOL_ID,
        "protocol_sha256": digest(protocol_bytes), "source_asset_id": DEV_SOURCE_ID, "source_sha256": DEV_SOURCE_SHA,
        "train_ids": list(DEV_TRAIN), "development_ids": list(DEV_INPUTS),
        "exporter_code": [exporter_ref], "executor_code": executor_refs,
        "requested_phases": ["training_and_conditioning_export", "development_response_export"],
        "requester": "Parent Devin task/user instruction to the Granados custodian; encoded here against the actual resumed executor, not a fictitious external requester.",
    }
    dev_contract(protocol, "DevelopmentReleaseRequest").validate(request)
    request_ref = dev_save(DEV_ROOT / "phase_b/release_request.json", request, resume_identical=True)
    decision = {
        "document_type": "development_custody_decision_not_external_attestation", "protocol_id": DEV_PROTOCOL_ID,
        "request_sha256": request_ref["sha256"], "decision": "approve_development_exports",
        "custodian_identity": "Granados observation/data custodian subagent in this Devin session, distinct from the model fitting/selection agent; no external attestation claimed.",
        "operational_event_reference": relative(DEV_ROOT / "phase_b/approval_event.json") + "#phase-b-approved",
        "reviewed_source_asset_id": DEV_SOURCE_ID, "reviewed_source_sha256": DEV_SOURCE_SHA,
        "permitted_group_ids": list(DEV_TRAIN + DEV_INPUTS),
        "conditions": [
            "Current release is only train rep1/2/3 prefix+response, and development rep4/5 prefix plus response-frame times. Every original row and null is retained with its exact source locator.",
            "Projection files must be new regular nonsymlink files in the protocol's exact development_exports/granados-sfp1-published-cohort-development-01 root; custody metadata is separately stored in development/. No protocol or scientific criterion is changed by this writer-root correction.",
            "Use the real byte-frozen executor implementation and runtime/PRNG/optimizer settings recorded in the request. No model was fitted or scored by the custodian.",
            "Development response intensities/means remain withheld until all six canonical families, all fitting attempts and all 288 development frame predictions have been frozen, their digests received and independently checked. This command cannot release those responses.",
            "Only exact phase-b projection paths enter the learner data allowlist. Mixed raw files, rep6, every other group/field/window, and known aliases remain excluded.",
            "Failure or incompleteness is retained without source/group/cohort/window replacement. No biological ready grade, final validation, prospective-forecast certification or strong-claim authorization follows.",
        ],
        "limitations": [protocol["estimand"]["sampling_qualification"], protocol["context"]["intervention_identifiability"], protocol["partition"]["unknown_overlap"], protocol["authorization_boundary"]["security"]],
        "strong_claim_authorized": False,
    }
    dev_contract(protocol, "DevelopmentCustodyDecision").validate(decision)
    approval_ref = dev_save(DEV_ROOT / "phase_b/custody_decision.json", decision, resume_identical=True)
    approval_ref["role"] = "operational_approval"
    event_path = DEV_ROOT / "phase_b/approval_event.json"
    event = dev_strict_json(event_path.read_bytes()) if event_path.exists() else {
        "document_type": "local_operational_custody_event_not_external_registry", "event_id": "phase-b-approved",
        "recorded_at_utc": now(), "request_sha256": request_ref["sha256"], "approval_sha256": approval_ref["sha256"],
        "source_binding_audit_sha256": audit_ref["sha256"], "superseded_blocker_sha256": digest((DEV_ROOT / "operational_blocker.json").read_bytes()),
        "resolution": "Saved source coverage passes. A real resumed executor now exists and was frozen. The writer now distinguishes custody metadata root from the exact protocol-required projection root; no source, cohort, window, target or threshold was amended.",
        "values_released_before_this_approval": False, "strong_claim_authorized": False,
    }
    if event["approval_sha256"] != approval_ref["sha256"] or event["request_sha256"] != request_ref["sha256"]:
        raise CustodyError("Existing approval event is bound to different operational bytes.")
    dev_save(event_path, event, resume_identical=True)
    for reference in request["exporter_code"] + request["executor_code"]:
        dev_verify_ref(reference)
    selected = dev_selected_source(protocol)
    artifacts, bindings = [], []
    for phase, filename, role in (("training_values", "training_projection.json", "training_projection"), ("development_inputs", "development_inputs.json", "development_inputs")):
        payload = dev_projection_payload(selected, phase)
        checks = dev_verify_payload(payload, selected, phase)
        binding = {
            "document_type": "source_scalar_replay_binding_not_biological_validation", "protocol_id": DEV_PROTOCOL_ID,
            "protocol_sha256": DEV_REVIEWED_PROTOCOL_SHA, "phase": phase,
            "approval_sha256": approval_ref["sha256"], "source_asset_id": DEV_SOURCE_ID, "source_sha256": DEV_SOURCE_SHA,
            "group_ids": payload["group_ids"], "source_binding_audit": audit_ref,
            "ordered_payload_sha256": digest((json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()),
            "verification": checks,
            "development_response_intensities_in_payload": False,
            "scalar_policy": "Direct equality to deposited max5/median/times at the declared pointers, including exact null retention. No response transformations, smoothing, interpolation or row renumbering.",
            "strong_claim_authorized": False,
        }
        binding_ref = dev_save(DEV_ROOT / "phase_b" / (phase + "_binding.json"), binding, resume_identical=True)
        binding_ref["role"] = "source_binding_audit"
        projection = {
            "document_type": "source_bound_development_projection", "protocol_id": DEV_PROTOCOL_ID, "phase": phase,
            "approval_sha256": approval_ref["sha256"], "source_asset_id": DEV_SOURCE_ID, "source_sha256": DEV_SOURCE_SHA,
            **payload, "source_binding_audit": binding_ref,
        }
        dev_contract(protocol, "Projection").validate(projection)
        output_ref = dev_save(DEV_PROJECTION_ROOT / filename, projection, resume_identical=True)
        output_ref["role"] = role
        artifacts.append(output_ref)
        bindings.append(binding_ref)
    # Re-read the source container for independent scalar-to-pointer replay of
    # the serialized outputs; hashes or generator pass flags alone do not pass.
    replay_source = dev_selected_source(protocol)
    for artifact in artifacts:
        projection = dev_strict_json(dev_verify_ref(artifact))
        dev_contract(protocol, "Projection").validate(projection)
        dev_verify_payload(projection, replay_source, projection["phase"])
    manifest = {
        "document_type": "approved_development_phase_b_data_allowlist_not_biological_authorization",
        "protocol_id": DEV_PROTOCOL_ID, "protocol_sha256": DEV_REVIEWED_PROTOCOL_SHA,
        "operational_phase": "training_and_conditioning_export_complete",
        "release_request": request_ref, "custody_decision": approval_ref,
        "source_asset_id": DEV_SOURCE_ID, "source_sha256": DEV_SOURCE_SHA,
        "learner_allowlist": artifacts, "source_binding_artifacts": bindings,
        "frozen_executor_artifacts": executor_refs, "exporter_code_artifacts": [exporter_ref],
        "protocol_snapshot": dev_artifact(DEV_ROOT / "contracts/development_protocol.json", "protocol_snapshot"),
        "approved_prefix_cells_by_group": {g["group_id"]: g["prefix_cells"] for g in prior_audit["group_checks"]},
        "all_planned_coordinates_and_targets_available_custodian_checked": True,
        "coverage_details": "All approved source rows, nulls and invalid entries remain in projections. Prefix exclusions and response attrition were retained in private/source_coverage.json; development response counts/values are not a learner input at this phase.",
        "development_response_projection": None, "development_prediction_freeze_received": False,
        "development_responses_released": False,
        "further_release_requires": "Exact complete DevelopmentPredictionFreeze, matching code/data/approval digests, all six canonical classes and 288 predictions including failures; independent source/test review before any separate response release.",
        "exposed_biological_group_ids": list(DEV_TRAIN + DEV_INPUTS),
        "unexported_scope": "Rep6, all known rep6 aliases, all other source groups/fields/windows and normalized-only rep2. Unexported portions of rep1-rep5 are not independent holdouts; unknown cross-file aliases must inherit development exposure if later identified.",
        "original_strong_artifact_sha256_unchanged": original_hashes,
        "strong_claim_authorized": False, "final_test_defined": False,
        "limitations": decision["limitations"],
    }
    for path in original_paths:
        if digest(path.read_bytes()) != original_hashes[relative(path)]:
            raise CustodyError("An original strong/custody artifact changed during the development operation.")
    manifest_ref = dev_save(DEV_ROOT / "release_manifest.json", manifest, resume_identical=True)
    dev_save(DEV_ROOT / "custody_handoff.json", {
        "document_type": "custodian_handoff_to_resumed_executor_source_verifier_and_tester",
        "current_phase": manifest["operational_phase"], "release_manifest": manifest_ref,
        "learner_allowlist": artifacts, "frozen_executor_artifacts": executor_refs,
        "request_and_approval": {"request": request_ref, "decision": approval_ref},
        "next_actions": ["Independently replay approved projections against the pinned source and exact protocol locators without propagating mixed raw data to the learner.", "Use the frozen implementation/settings and only the two allowlisted phase-b data files for training and conditioning.", "Freeze every canonical family/attempt and all 288 development predictions, retaining failures, and provide the digest-bound DevelopmentPredictionFreeze to custody.", "Do not access or export development responses until the later independent custody check; do not issue a biological ready grade."],
        "development_responses_released": False, "strong_claim_authorized": False,
    }, resume_identical=True)
    return dev_validate_release()


def dev_validate_release():
    protocol, _ = dev_protocol_checked()
    manifest_path = DEV_ROOT / "release_manifest.json"
    manifest = dev_strict_json(dev_regular_file(manifest_path).read_bytes())
    if manifest["protocol_sha256"] != DEV_REVIEWED_PROTOCOL_SHA or manifest["development_responses_released"] or manifest["development_response_projection"] is not None or manifest["strong_claim_authorized"]:
        raise CustodyError("Release manifest is not the reviewed phase-b-only, non-authorizing workflow.")
    request = dev_strict_json(dev_verify_ref(manifest["release_request"]))
    decision = dev_strict_json(dev_verify_ref(manifest["custody_decision"]))
    dev_contract(protocol, "DevelopmentReleaseRequest").validate(request)
    dev_contract(protocol, "DevelopmentCustodyDecision").validate(decision)
    if decision["decision"] != "approve_development_exports" or decision["request_sha256"] != manifest["release_request"]["sha256"]:
        raise CustodyError("No matching independent operational custody approval.")
    for ref in request["executor_code"] + request["exporter_code"]:
        dev_verify_ref(ref)
    if len(manifest["learner_allowlist"]) != 2:
        raise CustodyError("Phase-b allowlist must contain exactly two data projections.")
    source = dev_selected_source(protocol)
    phases = set()
    checks = []
    for ref in manifest["learner_allowlist"]:
        path = ROOT / ref["path"]
        if path.parent != DEV_PROJECTION_ROOT or path.name not in ("training_projection.json", "development_inputs.json"):
            raise CustodyError("Learner allowlist includes a non-protocol-root or raw/unapproved path.")
        projection = dev_strict_json(dev_verify_ref(ref))
        dev_contract(protocol, "Projection").validate(projection)
        if projection["approval_sha256"] != manifest["custody_decision"]["sha256"]:
            raise CustodyError("Projection approval binding mismatch.")
        phase = projection["phase"]
        phases.add(phase)
        checks.append(dev_verify_payload(projection, source, phase))
        binding = dev_strict_json(dev_verify_ref(projection["source_binding_audit"]))
        payload = {key: projection[key] for key in ("group_ids", "source_cells", "source_times")}
        payload_hash = digest((json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode())
        if binding["ordered_payload_sha256"] != payload_hash or binding["approval_sha256"] != projection["approval_sha256"] or binding["phase"] != phase:
            raise CustodyError("Projection source-binding record does not match the actual payload.")
    if phases != {"training_values", "development_inputs"}:
        raise CustodyError("Phase-b data is incomplete or includes development responses.")
    for path, expected in manifest["original_strong_artifact_sha256_unchanged"].items():
        if digest((ROOT / path).read_bytes()) != expected:
            raise CustodyError("An original strong/custody artifact no longer matches the retained development boundary.")
    if not (DEV_ROOT / "custody_handoff.json").exists():
        manifest_ref = dev_artifact(manifest_path, "operational_event_log")
        manifest_ref.pop("role")
        dev_save(DEV_ROOT / "custody_handoff.json", {
            "document_type": "custodian_handoff_to_resumed_executor_source_verifier_and_tester",
            "current_phase": manifest["operational_phase"], "release_manifest": manifest_ref,
            "learner_allowlist": manifest["learner_allowlist"], "frozen_executor_artifacts": manifest["frozen_executor_artifacts"],
            "request_and_approval": {"request": manifest["release_request"], "decision": manifest["custody_decision"]},
            "next_actions": ["Independently verify the source-bound phase-b projections and use the exact frozen executor/settings.", "Freeze all six canonical families, every attempt and all 288 development predictions before requesting the still-withheld development responses."],
            "development_responses_released": False, "strong_claim_authorized": False,
        })
    print(json.dumps({"approved_release_manifest": dev_artifact(manifest_path, "operational_event_log"), "learner_allowlist": manifest["learner_allowlist"], "source_scalar_replay": checks, "development_responses_released": False, "strong_claim_authorized": False}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("acquire", "inspect", "summarize", "primary", "method-excerpts", "audit-protocol", "source-links", "source-code", "identity-checks", "timing-review", "validate-review", "development-audit", "development-export-inputs", "development-validate", "development-export-responses", "validate", "export"))
    args = parser.parse_args()
    if args.command == "development-export-responses":
        raise CustodyError("Development responses remain withheld: a complete digest-bound six-family/288-prediction freeze and separate custodian/tester review have not been received for release.")
    if args.command == "export":
        if not PROTOCOL.exists():
            raise CustodyError("Export refused: verifier protocol.json does not exist.")
        if load(MANIFEST)["status"] == "metadata_schema_verified_ineligible_current_protocol":
            raise CustodyError("Export refused: the registered Sfp1/carbon candidate is metadata-ineligible; every source group remains reserved and the training allowlist is empty.")
        raise CustodyError("Export refused: protocol-file existence alone is not approval. Explicit eligible whole-group IDs, an approved observation mapping and independent custody authorization are required.")
    {"acquire": acquire, "inspect": inspect, "summarize": summarize, "primary": primary, "method-excerpts": method_excerpts, "audit-protocol": audit_protocol, "source-links": source_links, "source-code": source_code, "identity-checks": identity_checks, "timing-review": timing_review, "validate-review": validate_review, "development-audit": dev_source_audit, "development-export-inputs": dev_phase_b_export, "development-validate": dev_validate_release, "validate": validate}[args.command]()


if __name__ == "__main__":
    os.umask(0o077)
    try:
        main()
    except CustodyError as error:
        sys.exit(str(error))
    except Exception as error:
        # Do not expose data-dependent exception payloads from sealed files.
        sys.exit(f"Custody operation stopped: {type(error).__name__}; no source values were logged.")
