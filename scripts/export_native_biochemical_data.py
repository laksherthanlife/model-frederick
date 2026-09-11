#!/usr/bin/env python3
"""Metadata-only inspection of pinned public native-biochemistry sources.

No measurement exporter is authorized at present. Hackett's intracellular flux
labels are model-derived, and the current protocol cannot use them. This script
therefore has no measurement-output or raw-file-output option. It retrieves two
pinned author-repository files in memory, verifies Git blob and SHA256 digests,
and prints schemas only. It neither imports the author's analysis code nor runs
R, a model, a fit, or a scoring procedure. No dependencies are installed.

The small R XDR-v2 reader is deliberately limited to the passive object types
present in the source. Unknown types, versions, references, lengths and trailing
bytes fail closed. Numeric vector payloads are never included in output. The
only numeric payload interpreted by schema() is an authoritative R `dim`
attribute (array shape), not a measurement column. Character measurement
columns are withheld. Only explicitly named specimen/condition/batch identifier
columns from chemostatInfo are exposed separately as experiment metadata.

Usage:
    python3 scripts/export_native_biochemical_data.py --self-test
    python3 scripts/export_native_biochemical_data.py --inspect-hackett

Write reviewed metadata artifacts using the agent's write/edit tools, not this
program. Any future quantitative export requires a separate protocol-approved
projection implementation, fixed whole-experiment IDs, source-cell checks,
and an independently registered declaration before measurement access.
"""

from __future__ import annotations

import argparse
import base64
import csv
from dataclasses import dataclass, field
import gzip
import hashlib
import io
import json
import struct
import subprocess
from typing import Any


REPOSITORY = "shackett/simmer"
COMMIT = "7138a45fb9653ce50d8556676a9fd1cc00bf0395"
SOURCES = {
    "protein_abundance": {
        "path": "inst/extdata/flux_input_data/proteinAbundance.tsv",
        "git_blob_sha": "b951dc2ec4bcf86534c58f6efdc331676b78c406",
        "sha256": "5227daf2a6358f188ec58b2dbfebadeac27c76d7eef03a0a82e0ad86ac9e3a9d",
        "size": 553249,
    },
    "boundary_inputs": {
        "path": "inst/extdata/flux_input_data/boundaryFluxes.Rdata",
        "git_blob_sha": "08be5778cc7dd77242bae0535bb436eb199e2522",
        "sha256": "15e9bd86d5667671f58e911fa77374bf33085b05052f10113cc7c685af8338d1",
        "size": 70639,
    },
}

TYPE_NAMES = {
    0: "NILSXP",
    1: "SYMSXP",
    2: "LISTSXP",
    6: "LANGSXP",
    9: "CHARSXP",
    10: "LGLSXP",
    13: "INTSXP",
    14: "REALSXP",
    15: "CPLXSXP",
    16: "STRSXP",
    19: "VECSXP",
    20: "EXPRSXP",
    22: "EXTPTRSXP",
    24: "RAWSXP",
    254: "NILVALUE_SXP",
}


@dataclass
class Node:
    kind: int
    length: int | None = None
    value: Any = None
    attrs: Node | None = None
    tag: Node | None = None
    car: Node | None = None
    cdr: Node | None = None
    numeric_bytes: bytes = field(default=b"", repr=False)


class RSchemaReader:
    """Passive, bounds-checked R XDR serialization-v2 subset reader."""

    def __init__(self, payload: bytes):
        if payload.startswith(b"\x1f\x8b"):
            payload = gzip.decompress(payload)
        if len(payload) > 4_000_000:
            raise ValueError("R payload exceeds the pinned-source safety limit")
        if not payload.startswith(b"RDX2\nX\n"):
            raise ValueError("Expected an RDX2 XDR serialization")
        self.payload = payload
        self.offset = 7
        self.references: list[Node] = []
        self.version = self.integer()
        self.writer_version = self.integer()
        self.minimum_reader_version = self.integer()
        if self.version != 2:
            raise ValueError("Only serialization version 2 is supported")

    def take(self, size: int) -> bytes:
        if size < 0 or self.offset + size > len(self.payload):
            raise ValueError("Invalid R object length or truncated payload")
        chunk = self.payload[self.offset : self.offset + size]
        self.offset += size
        return chunk

    def integer(self) -> int:
        return struct.unpack(">i", self.take(4))[0]

    def vector_length(self) -> int:
        length = self.integer()
        if not 0 <= length <= 1_000_000:
            raise ValueError("Unsupported R vector length")
        return length

    def node(self, depth: int = 0) -> Node:
        if depth > 300:
            raise ValueError("R object nesting exceeds safety limit")
        flags = self.integer() & 0xFFFFFFFF
        kind = flags & 255
        has_attrs = bool(flags & (1 << 9))
        has_tag = bool(flags & (1 << 10))
        if kind == 255:
            reference = flags >> 8
            if not reference:
                reference = self.integer()
            if not 1 <= reference <= len(self.references):
                raise ValueError("Invalid R object reference")
            return self.references[reference - 1]
        if kind not in TYPE_NAMES:
            raise ValueError(f"Unsupported R serialization type {kind}")
        out = Node(kind)
        if kind in (0, 254):
            return out
        if kind == 1:
            out.value = self.node(depth + 1)
            self.references.append(out)
            return out
        if kind in (2, 6):
            if has_attrs:
                out.attrs = self.node(depth + 1)
            if has_tag:
                out.tag = self.node(depth + 1)
            out.car = self.node(depth + 1)
            out.cdr = self.node(depth + 1)
            return out
        if kind == 22:
            self.references.append(out)
            out.value = [self.node(depth + 1), self.node(depth + 1)]
        elif kind == 9:
            out.length = self.integer()
            if out.length == -1:
                out.value = None
            elif 0 <= out.length <= 1_000_000:
                # Retain source strings internally; schema output uses only
                # authoritative names/class attributes, never column values.
                out.value = self.take(out.length).decode("utf-8", errors="strict")
            else:
                raise ValueError("Invalid R character length")
        elif kind in (10, 13, 14, 15, 24):
            out.length = self.vector_length()
            width = {10: 4, 13: 4, 14: 8, 15: 16, 24: 1}[kind]
            payload = self.take(out.length * width)
            # Only integer bytes are needed for named shape attributes.
            if kind == 13:
                out.numeric_bytes = payload
        elif kind in (16, 19, 20):
            out.length = self.vector_length()
            out.value = [self.node(depth + 1) for _ in range(out.length)]
        else:
            raise ValueError("Unsupported passive R object")
        if has_attrs:
            out.attrs = self.node(depth + 1)
        return out

    def parse(self) -> Node:
        root = self.node()
        if self.offset != len(self.payload):
            raise ValueError("Unconsumed bytes: schema parser did not verify source")
        return root


def strings(node: Node | None) -> list[str | None]:
    if node is None:
        return []
    if node.kind == 1:
        return strings(node.value)
    if node.kind == 9:
        return [node.value]
    if node.kind == 16:
        return [s for item in node.value for s in strings(item)]
    raise ValueError("Expected an R string metadata attribute")


def pairlist(node: Node | None) -> dict[str, Node]:
    pairs: dict[str, Node] = {}
    while node is not None and node.kind not in (0, 254):
        if node.kind != 2 or node.car is None:
            raise ValueError("Expected a tagged R pairlist")
        names = strings(node.tag)
        if len(names) != 1 or names[0] is None or names[0] in pairs:
            raise ValueError("Missing or duplicate R pairlist name")
        pairs[names[0]] = node.car
        node = node.cdr
    return pairs


def schema(node: Node, depth: int = 0) -> dict[str, Any]:
    attributes = pairlist(node.attrs)
    result: dict[str, Any] = {"r_type": TYPE_NAMES[node.kind]}
    if node.length is not None:
        result["length"] = node.length
    if "class" in attributes:
        result["class"] = strings(attributes["class"])
    if "dim" in attributes:
        dims = attributes["dim"]
        if dims.kind != 13 or dims.length is None:
            raise ValueError("Noninteger dimension attribute")
        shape = list(struct.unpack(">" + "i" * dims.length, dims.numeric_bytes))
        if any(n < 0 for n in shape):
            raise ValueError("Invalid negative array dimension")
        result["dim"] = shape
    if node.kind == 19 and depth < 4:
        names = strings(attributes.get("names"))
        if names and len(names) != len(node.value):
            raise ValueError("R names length does not match list length")
        result["fields"] = [
            {"name": names[i] if names else f"unnamed_{i + 1}", **schema(item, depth + 1)}
            for i, item in enumerate(node.value)
        ]
    return result


def experiment_metadata(node: Node) -> list[dict[str, str | None]]:
    """Return only the exact primary-source identity and acquisition fields."""
    if node.kind != 19:
        raise ValueError("Expected chemostatInfo to be an R data frame")
    names = strings(pairlist(node.attrs).get("names"))
    columns = dict(zip(names, node.value, strict=True))
    identifiers = (
        "ChemostatID", "ChemostatCond", "Sheet", "Limitation",
        "composition_data", "NMR_name", "NMR_date",
    )
    if len(columns) != len(names):
        raise ValueError("Duplicate chemostatInfo column name")
    selected = {}
    for name in identifiers:
        if columns[name].kind != 16:
            raise ValueError("An experiment identifier is not an R string vector")
        values = strings(columns[name])
        if len(values) != 25:
            raise ValueError("Unexpected number of chemostat experiment identities")
        selected[name] = values
    return [{name: values[i] for name, values in selected.items()} for i in range(25)]


def fetch(source: dict[str, Any]) -> bytes:
    endpoint = f"repos/{REPOSITORY}/contents/{source['path']}?ref={COMMIT}"
    response = subprocess.run(
        ["gh", "api", endpoint],
        capture_output=True,
        check=False,
        timeout=28,
    )
    if response.returncode:
        # Do not echo an untrusted response body that might contain data.
        raise RuntimeError(f"gh metadata retrieval failed: exit {response.returncode}")
    description = json.loads(response.stdout)
    payload = base64.b64decode(description["content"], validate=False)
    blob = hashlib.sha1(b"blob " + str(len(payload)).encode() + b"\0" + payload).hexdigest()
    if not (
        description["sha"] == source["git_blob_sha"] == blob
        and hashlib.sha256(payload).hexdigest() == source["sha256"]
        and len(payload) == source["size"]
    ):
        raise ValueError("Pinned source digest or size mismatch")
    return payload


def inspect_hackett() -> dict[str, Any]:
    output: dict[str, Any] = {
        "mode": "schema_only_no_quantitative_export",
        "repository": REPOSITORY,
        "commit": COMMIT,
        "learner_measurement_allowlist": [],
        "raw_sources_persisted": False,
        "files": {},
    }
    for name, source in SOURCES.items():
        payload = fetch(source)
        record: dict[str, Any] = dict(source)
        record["source_url"] = f"https://github.com/{REPOSITORY}/blob/{COMMIT}/{source['path']}"
        record["source_digest_verified"] = True
        if name == "protein_abundance":
            rows = csv.reader(io.StringIO(payload.decode("utf-8")), delimiter="\t")
            header = next(rows)
            expected = ["Gene"] + [
                nutrient + rate
                for nutrient in ("c", "L", "n", "p", "u")
                for rate in ("0.05", "0.11", "0.16", "0.22", "0.30")
            ]
            if header != expected:
                raise ValueError("Unexpected protein TSV header")
            count = 0
            for row in rows:
                if len(row) != len(header):
                    raise ValueError("Inconsistent protein TSV row width")
                count += 1
            record["schema"] = {"header": header, "data_rows": count, "columns": len(header)}
        else:
            reader = RSchemaReader(payload)
            objects = pairlist(reader.parse())
            record["serialization"] = {
                "format": "RDX2 XDR",
                "version": reader.version,
                "complete_byte_consumption": True,
            }
            record["schema"] = {key: schema(value) for key, value in objects.items()}
            record["experiment_metadata"] = experiment_metadata(objects["chemostatInfo"])
        output["files"][name] = record
    return output


def self_test() -> dict[str, Any]:
    def integer(n):
        return struct.pack(">i", n)

    def char(s):
        return integer(9) + integer(len(s.encode())) + s.encode()

    def symbol(s):
        return integer(1) + char(s)

    def str_vector(vals):
        return integer(16) + integer(len(vals)) + b"".join(char(v) for v in vals)

    nil = integer(254)
    names = integer(2 | (1 << 10)) + symbol("names") + str_vector(["sample_id", "rate"]) + nil
    frame_class = integer(2 | (1 << 10)) + symbol("class") + str_vector(["data.frame"]) + names
    # Synthetic engineering-only values, never biological observations.
    frame = (
        integer(19 | (1 << 9))
        + integer(2)
        + str_vector(["synthetic_a", "synthetic_b"])
        + integer(14)
        + integer(2)
        + struct.pack(">dd", 123456.25, -98765.5)
        + frame_class
    )
    encoded = b"RDX2\nX\n" + integer(2) + integer(0x30202) + integer(0x20300)
    encoded += integer(2 | (1 << 10)) + symbol("synthetic_frame") + frame + nil
    parsed = pairlist(RSchemaReader(gzip.compress(encoded)).parse())
    out = schema(parsed["synthetic_frame"])
    expected = {
        "r_type": "VECSXP",
        "length": 2,
        "class": ["data.frame"],
        "fields": [
            {"name": "sample_id", "r_type": "STRSXP", "length": 2},
            {"name": "rate", "r_type": "REALSXP", "length": 2},
        ],
    }
    if out != expected:
        raise AssertionError("Synthetic R schema mismatch")
    text = json.dumps(out)
    if any(s in text for s in ("123456", "98765", "synthetic_a", "synthetic_b")):
        raise AssertionError("Schema exposed a synthetic data value")
    for invalid in (encoded[:-1], encoded + b"unexpected", b"RDX3\nX\n"):
        try:
            RSchemaReader(invalid).parse()
        except (ValueError, struct.error):
            continue
        raise AssertionError("Malformed input was not refused")
    return {
        "self_test": "passed",
        "synthetic_only": True,
        "checks": ["schema_structure", "no_numeric_or_character_column_values", "truncation_refusal", "trailing_bytes_refusal", "unsupported_version_refusal"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--inspect-hackett", action="store_true")
    modes.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    output = self_test() if args.self_test else inspect_hackett()
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
