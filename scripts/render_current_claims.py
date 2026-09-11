from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ystwin.analysis.claims import audit_claims


_FILENAMES = ("claim_verdicts.json", "claim_inventory.csv", "claim_results.csv", "claim_report.md")


def _text(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list, bool)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False)
    return str(value)


def _cell(value) -> str:
    return _text(value).replace("|", "\\|").replace("\r", " ").replace("\n", " ")


def _inventory(rows):
    keys = ("id", "record_type", "family", "statement", "role", "method", "conditions", "model",
            "configuration", "evidence", "evaluation", "expected", "uncertainty", "coverage",
            "coverage_observed", "independent_unit", "provenance", "status", "verdict", "supported",
            "blocking", "comparison_status", "detail", "location", "selector")
    return [dict({key: _text(row.get(key)) for key in keys},
                 supported_result=_text(row.get("actual")) if row.get("supported") else "",
                 recorded_observation=_text(row.get("observed")),
                 observation_status="VERIFIED" if row.get("actual") is not None else "UNVERIFIED")
            for row in rows]


def _results(rows):
    result = []
    for row in rows:
        for index, observation in enumerate(row.get("result_rows", []), start=1):
            result.append({
                "claim_id": row["id"], "record_type": row["record_type"], "role": row.get("role", ""),
                "verdict": row["verdict"], "supported": row["supported"], "result_index": index,
                "evidence_status": observation["evidence_status"], "source": observation["source"],
                "row_key": _text(observation["row_key"]), "unit": observation["unit"],
                "input_unit": observation.get("input_unit", ""), "context": _text(observation.get("context")),
                "predicted_or_observed": _text(observation.get("observed")),
                "measured_reference": _text(observation.get("reference")),
                "value": _text(observation.get("value")), "low": _text(observation.get("low")),
                "high": _text(observation.get("high")), "recorded_result": _text(observation),
            })
    return result


def _report(rows) -> str:
    counts = Counter(row["verdict"] for row in rows)
    lines = [
        "# Structured scientific claim audit", "",
        "Only declared scientific claim families and explicit audit:value selectors are inventoried.",
        "Citation years, arbitrary decimals and historical test counts are not scientific claim selectors.", "",
        "A BLOCKED observation is cached, unverified evidence, not a supported current result. "
        "A matching copied number does not verify its method, provenance or interpretation. "
        "Refusal confirmation is support for a refusal only, never a numerical prediction.", "",
        "Rendering performs no fitting, regeneration, evidence updates or provenance restamping.", "",
        "## Verdicts", "", "| Verdict | Records |", "| --- | ---: |",
        *[f"| {key} | {value} |" for key, value in sorted(counts.items())], "",
        "## Scientific claim inventory", "",
    ]
    for row in rows:
        if row["record_type"] != "claim":
            continue
        lines += [f"### {_cell(row['id'])}", "", _cell(row["statement"]), "",
                  "| Field | Declared value or audit result |", "| --- | --- |"]
        fields = [("Role", row["role"]), ("Verdict", row["verdict"]), ("Method", row["method"]),
                  ("Conditions", row["conditions"]), ("Model", row["model"]),
                  ("Configuration", row["configuration"]), ("Evidence selectors", row["evidence"]),
                  ("Evaluation", row["evaluation"]), ("Expected / acceptance", row["expected"]),
                  ("Independent unit", row["independent_unit"]), ("Coverage", row["coverage"]),
                  ("Observed coverage", row.get("coverage_observed")), ("Uncertainty", row["uncertainty"]),
                  ("Supported current result", row["actual"] if row["supported"] else "NOT AUTHORIZED"),
                  ("Recorded observation (may be unverified)", row["observed"]),
                  ("Audit detail", row["detail"])]
        lines += [f"| {_cell(key)} | {_cell(value)} |" for key, value in fields]
        lines.append("")
    lines += ["## Explicit marker consistency", "",
              "These are copies of selected table cells, not additional independent experiments.", "",
              "| Claim id | Location | Selector | Literal expectation | Recorded observation | Cell comparison | Scientific verdict |",
              "| --- | --- | --- | --- | --- | --- | --- |"]
    for row in rows:
        if row["record_type"] == "marked_value":
            fields = [row["id"], row.get("location"), row.get("selector"), row.get("expected"),
                      row.get("observed"), row["comparison_status"], row["verdict"]]
            lines.append("| " + " | ".join(_cell(value) for value in fields) + " |")
    lines += ["", "## Registry and surface checks", "",
              "| Check | Status | Actual | Detail |", "| --- | --- | --- | --- |"]
    for row in rows:
        if row["record_type"] in {"scope", "registry"}:
            lines.append("| " + " | ".join(_cell(row.get(key)) for key in ("id", "status", "actual", "detail")) + " |")
    return "\n".join(lines) + "\n"


def _csv(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def render(rows: list[dict], output_dir: Path, *, root: Path) -> None:
    root = root.resolve()
    output_dir = output_dir.resolve()
    protected_trees = [root / part for part in ("data", "docs", "src", "scripts", "tests", ".git")]
    if output_dir in (root, root / "outputs") or any(output_dir.is_relative_to(path) for path in protected_trees):
        raise ValueError("use an isolated output directory, not a canonical repository surface")
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise ValueError("output directory must be new or empty; existing outputs are never overwritten")
    for relative in _FILENAMES:
        if (output_dir / relative).exists():
            raise ValueError(f"refusing to overwrite {relative}")
    inventory = _inventory(rows)
    results = _results(rows)
    report = _report(rows)
    serialized = json.dumps(rows, indent=2, ensure_ascii=False, sort_keys=True, allow_nan=False) + "\n"
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / _FILENAMES[0]).open("x", encoding="utf-8") as handle:
        handle.write(serialized)
    inventory_fields = tuple(inventory[0]) if inventory else ("id", "status", "verdict", "detail")
    _csv(output_dir / _FILENAMES[1], inventory, inventory_fields)
    result_fields = ("claim_id", "record_type", "role", "verdict", "supported", "result_index",
                     "evidence_status", "source", "row_key", "unit", "input_unit", "context",
                     "predicted_or_observed", "measured_reference", "value", "low", "high", "recorded_result")
    _csv(output_dir / _FILENAMES[2], results, result_fields)
    with (output_dir / _FILENAMES[3]).open("x", encoding="utf-8") as handle:
        handle.write(report)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Render declared scientific claims without modifying evidence or stamps.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--registry", type=Path, default=None, help="registry path relative to --root")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="new or empty isolated destination; relative paths are resolved from --root")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    destination = args.output_dir if args.output_dir.is_absolute() else root / args.output_dir
    rows = audit_claims(root, args.registry)
    try:
        render(rows, destination, root=root)
    except (OSError, ValueError) as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2
    blocking = sum(row["blocking"] for row in rows)
    print(f"{len(rows)} structured records; {blocking} blocking. Report: {destination / 'claim_report.md'}")
    return 1 if blocking else 0


if __name__ == "__main__":
    raise SystemExit(main())
