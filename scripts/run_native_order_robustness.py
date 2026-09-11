"""Conditional order/magnitude modelling of existing native development readouts.

No GEM, parameter fit, new assay, independent validation, or biological probability.
Only the verified native_conditions() adapter supplies measurement/precision data.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import asdict, replace
from decimal import Decimal, InvalidOperation
import io
from itertools import combinations
import json
import math
from numbers import Real
from pathlib import Path
import platform
import sys

import numpy as np
import scipy

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ystwin.analysis.partial_orders import (
    Bound, Calibration, Conversion, MagnitudeDomain, OrderClaim, Scope, Source,
    affine_extrema, compile_order, evaluate_forward, sample_uniform_box,
)
from ystwin.fba.native_reconciliation import native_conditions


# User-requested contrasts, fixed before reading data or solving; not a claim of preregistration.
PREDECLARED_CONTRASTS = {"GlucoseUptake": (0.4, 0.1), "O2uptake": (0.4, 0.28)}
YIELD_GROWTH_RATE_PER_H = 0.1
GLUCOSE_MOLAR_MASS_G_PER_MOL = 180.156
MOLAR_MASS_PROVENANCE = {
    "formula": "C6H12O6", "value": GLUCOSE_MOLAR_MASS_G_PER_MOL, "unit": "g/mol",
    "reference": "src/ystwin/mech/population.py::GLUCOSE_G_PER_MOL",
    "method": "Same conventional atomic-weight arithmetic: 6(12.011) + 12(1.008) + 6(15.999). Fixed unit conversion, not an inferred biological parameter.",
}
EXCHANGE_UNIT = "mmol/gDW/h"
DEVELOPMENT_SCOPES = {"native_development_training", "native_development_interpolation"}
REPORT_NAME = "native_order_robustness.json"
SUMMARY_NAME = "native_order_summary.csv"
OBSERVATIONS_NAME = "native_order_observations.csv"
GAP_SOURCE = Source(
    "assumed", "scripts/run_native_order_robustness.py:closed_order_relaxation",
    "gap=0: closed weak order-preserving relaxation; no positive separation or magnitude inferred from rank",
)
PROPOSAL_SOURCE = Source(
    "assumed", "scripts/run_native_order_robustness.py:bounded_uniform_proposal",
    "Independent uniform explicit nearest-rounding intervals across conditions, conditioned on numeric orders; a computational proposal, not a biological posterior or a model of measurement errors",
)
ASSUMPTIONS = [
    "Primary Table 1 tokens have closed +/- half-last-printed-decimal bands ONLY under nearest rounding; these are not statistical CIs, SDs, or bounds on biological/calibration error.",
    "Within each source/table/readout, published positive uptake magnitudes in mmol/gDW/h are comparable across the reported steady-state growth conditions. Glucose and oxygen are separate quantities/scopes.",
    "Calibration gain='known' denotes use of the published numerical coordinate, conditional on that coordinate convention; the source does not independently establish error-free instrument calibration.",
    "An edge requires both explicit bands and upper(lower_item) < lower(upper_item), compared in decimal arithmetic. Overlap, touching, or missing bounds imply incomparability, never a tie.",
    "Orders and bounded domains come from the SAME measured tokens plus the SAME rounding assumption; they are not independent evidence. Disjoint-band orders are redundant once those bounds are imposed.",
    "Order-only discards all magnitude bounds, anchors and positivity constraints. Its gap=0 LP is a closed weak relaxation of the strict qualitative order, not a physically calibrated uptake domain.",
    "Source zeros mean below detection limit; the source provides no numeric detection limit. They supply neither an exact zero nor a nearest-rounding bound/edge.",
    "Growth is an imposed reported dilution-rate setpoint, mu=D at chemostat steady state; no growth uncertainty band is inferred.",
    "Uniform proposals and a chosen shared rescaling are explicitly assumed computations, not biological probabilities or extra measurements.",
]


def _finite(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite real number")
    return float(value)


def nearest_rounding_band(record):
    """Use explicit verified precision, never pandas/float formatting or a zero token alone."""
    if record.get("rounding_is_confidence_interval") is not False:
        raise ValueError("nearest rounding must be explicitly labelled not a confidence interval")
    status = record["observation_status"]
    reported = _finite(record["reported_value"], "reported value")
    if status == "below_detection_limit":
        if reported != 0 or any(record[name] is not None for name in
                                ("value", "rounding_half_width", "detection_limit")):
            raise ValueError("native censored zeros require missing value, precision and detection limit")
        return None
    if status != "quantified":
        raise ValueError("unknown native observation status")
    value = None if record["value"] is None else _finite(record["value"], "quantified value")
    if value is not None and (value <= 0 or value != reported):
        raise ValueError("quantified value must equal its positive reported value")
    if value is None or record["rounding_half_width"] is None:
        return None  # Do not recover missing precision from the token or another condition.
    half = _finite(record["rounding_half_width"], "rounding half width")
    token = record["primary_printed_token"]
    if not isinstance(token, str):
        raise ValueError("verified primary printed token must be a string")
    try:
        center = Decimal(token)
    except InvalidOperation as exc:
        raise ValueError("invalid primary printed token") from exc
    if not center.is_finite() or float(center) != value:
        raise ValueError("primary token disagrees with quantified value")
    width = Decimal(str(half))
    if half <= 0 or width != Decimal(1).scaleb(center.as_tuple().exponent) / 2:
        raise ValueError("explicit rounding precision disagrees with verified primary token")
    return center - width, center + width


def _development_conditions(conditions):
    conditions = tuple(conditions)
    if not conditions:
        raise ValueError("native development conditions must be nonempty")
    if any(row["evidence_scope"] not in DEVELOPMENT_SCOPES for row in conditions):
        raise ValueError("only already-development native evidence scopes are allowed")
    if any(row["source"] != conditions[0]["source"] or
           row["precision_provenance"] != conditions[0]["precision_provenance"] for row in conditions):
        raise ValueError("mixed source/context/precision provenance cannot define one native scope")
    growths = [_finite(row["growth_rate_per_h"], "growth setpoint") for row in conditions]
    ids = [row["condition_id"] for row in conditions]
    if len(set(ids)) != len(ids) or len(set(growths)) != len(growths):
        raise ValueError("native condition identifiers and growth setpoints must be unique")
    return tuple(sorted(conditions, key=lambda row: row["growth_rate_per_h"]))


def build_uptake_domains(conditions, observable_id, *, tolerance=1e-8):
    """Build separate compatible public-API domains from explicit measured-readout bands."""
    if observable_id not in PREDECLARED_CONTRASTS:
        raise ValueError("this demonstration only models glucose and oxygen uptake")
    conditions = _development_conditions(conditions)
    source = conditions[0]["source"]
    if any(row["observations"][observable_id]["unit"] != EXCHANGE_UNIT for row in conditions):
        raise ValueError("native uptake unit must be mmol/gDW/h; no implicit unit conversion")
    context = json.dumps({
        "source_doi": source["doi"], "table": source["table"], **source["conditions"],
        "comparison_axis": "reported steady-state growth_rate_per_h conditions",
    }, sort_keys=True)
    scope = Scope(
        "magnitude", context, observable_id, EXCHANGE_UNIT,
        Calibration(f"{source['doi']}:{source['table']}:{observable_id}:published_coordinate", "known"),
    )
    reference = f"doi:{source['doi']}; {source['table']}, p. {source['table_page']}; {observable_id}"
    rows, bands, bounds = [], {}, []
    for condition in conditions:
        item = condition["condition_id"]
        record = condition["observations"][observable_id]
        band = nearest_rounding_band(record)
        bands[item] = band
        row = {
            "condition_id": item, "growth_rate_per_h": condition["growth_rate_per_h"],
            "evidence_scope": condition["evidence_scope"], "source_column": observable_id,
            "evidence_id": f"{item}:{observable_id}", "source_doi": source["doi"],
            "source_table_sha256": source["table_sha256"],
            "source_table_relative_path": source["table_relative_path"], **record,
            "rounding_band": None if band is None else {
                "lower": float(band[0]), "upper": float(band[1]),
                "decimal_lower": str(band[0]), "decimal_upper": str(band[1]),
            },
            "band_status": "explicit_nearest_rounding" if band is not None else
                           "unavailable_censored_or_missing_precision",
        }
        rows.append(row)
        if band is not None:
            bounds.append(Bound(item, float(band[0]), float(band[1]), scope, Source(
                "assumed", f"{reference}; {source['table_relative_path']}:{record['source_line']}",
                "Conditional closed nearest-rounding band from the verified primary token; not a statistical CI or independent evidence beyond the order",
            )))
    claims = []
    for left, right in combinations(rows, 2):
        a, b = bands[left["condition_id"]], bands[right["condition_id"]]
        if a is None or b is None:
            continue
        if a[1] < b[0]:
            lower, upper = left["condition_id"], right["condition_id"]
        elif b[1] < a[0]:
            lower, upper = right["condition_id"], left["condition_id"]
        else:
            continue
        claims.append(OrderClaim(lower, upper, scope, Source(
            "assumed", f"{reference}; TSV lines {left['source_line']}, {right['source_line']}",
            "Conditional magnitude order: explicit nearest-rounding bands strictly separated in decimal arithmetic; same measured evidence as the bounds, not enzyme importance or causal/temporal order",
        )))
    order = compile_order(bands, scope, claims)
    return {
        "rows": rows, "order": order,
        "measurement_source": Source("measured", reference, source["exchange_convention"]),
        "order_only": MagnitudeDomain(order, gap=0, gap_source=GAP_SOURCE, tolerance=tolerance),
        "bounded": MagnitudeDomain(order, bounds, gap=0, gap_source=GAP_SOURCE, tolerance=tolerance),
    }


def rank_ranges(order):
    """Possible 1-based ranks in linear extensions, not arbitrary tie-breaking or magnitudes."""
    closure = set(order.closure)
    representatives = [group[0] for group in order.classes]
    return {item: {
        "minimum_rank": 1 + sum((other, item) in closure for other in representatives),
        "maximum_rank": len(representatives) - sum((item, other) in closure for other in representatives),
    } for item in order.items}


def cover_edges(order):
    closure = set(order.closure)
    return [[a, b] for a, b in order.closure if not any(
        (a, middle) in closure and (middle, b) in closure for middle in order.items
    )]


def rescale_domain(domain, gain):
    """A declared coordinate-sensitivity check, NOT an estimated biological calibration."""
    gain = _finite(gain, "shared scale")
    if gain <= 0:
        raise ValueError("shared scale must be positive")
    source = Source(
        "assumed", "scripts/run_native_order_robustness.py:shared_scale_diagnostic",
        f"Chosen common positive coordinate multiplier {gain:g}; not a measured gain or new range evidence",
    )
    scope = replace(
        domain.order.scope, unit=f"diagnostic_coordinate({gain:g} * {domain.order.scope.unit})",
        calibration=Calibration(f"{domain.order.scope.calibration.family}:shared_scale_diagnostic", "known"),
    )
    conversion = Conversion(domain.order.scope, scope, gain, 0.0, source)
    order = compile_order(domain.order.items, scope, domain.order.claims,
                          ties=domain.order.ties, conversions=[conversion])
    return MagnitudeDomain(order, domain.bounds, domain.anchors,
                           gap=gain * domain.gap, gap_source=source, tolerance=domain.tolerance)


def apparent_specific_yield(mu_per_h, q_glucose_mmol_per_gdw_h):
    mu = _finite(mu_per_h, "growth setpoint")
    q = _finite(q_glucose_mmol_per_gdw_h, "glucose uptake")
    if mu < 0 or q <= 0:
        raise ValueError("apparent yield requires nonnegative growth and positive glucose uptake")
    return _finite(mu / (q * GLUCOSE_MOLAR_MASS_G_PER_MOL / 1000), "apparent specific yield")


def _domain_reference(result, domain_id):
    # Store every domain once, while preserving all core solver/proposal certificate fields.
    result.pop("domain")
    result["domain_id"] = domain_id
    return result


def _affine(domain, coefficients, estimand, domain_id):
    result = _domain_reference(affine_extrema(
        domain, coefficients, estimand=estimand, unit=domain.order.scope.unit,
    ), domain_id)
    return {**result, "conditional": True, "biological_probability": None}


def build_demo(*, samples=256, seed=1729, budget=4096, tolerance=1e-8, shared_scale=7.0):
    conditions = _development_conditions(native_conditions())
    indexed = {row["growth_rate_per_h"]: row for row in conditions}
    scopes, batches = {}, {}
    for key, (a, b) in PREDECLARED_CONTRASTS.items():
        built = build_uptake_domains(conditions, key, tolerance=tolerance)
        domains = {name: built[name] for name in ("order_only", "bounded")}
        domains["scaled"] = rescale_domain(domains["bounded"], shared_scale)
        items = (indexed[a]["condition_id"], indexed[b]["condition_id"])
        coefficients = {items[0]: 1.0, items[1]: -1.0}
        estimand = f"{key}({items[0]}) - {key}({items[1]})"
        refs = {name: f"#/scopes/{key}/domains/{name}" for name in domains}
        affine = {
            "order_only_absolute": _affine(domains["order_only"], {items[0]: 1.0},
                                            f"{key}({items[0]}) without magnitude bounds", refs["order_only"]),
            "order_only_contrast": _affine(domains["order_only"], coefficients, estimand, refs["order_only"]),
            "bounded_contrast": _affine(domains["bounded"], coefficients, estimand, refs["bounded"]),
            "scaled_contrast": _affine(domains["scaled"], coefficients, estimand, refs["scaled"]),
        }
        errors = {}
        for side in ("minimum", "maximum"):
            original, scaled = affine["bounded_contrast"][side], affine["scaled_contrast"][side]
            errors[side] = (abs(scaled["value"] - shared_scale * original["value"])
                            if original["validated"] and scaled["validated"] else None)
        verified = None if any(error is None for error in errors.values()) else all(
            error <= tolerance * max(1.0, abs(affine["scaled_contrast"][side]["value"]))
            for side, error in errors.items()
        )
        sampling = {}
        for name in ("order_only", "bounded"):
            batch = sample_uniform_box(domains[name], n=samples, seed=seed, budget=budget,
                                       source=PROPOSAL_SOURCE)
            sampling[name] = _domain_reference(batch.report(), refs[name])
            if name == "bounded":
                batches[key] = batch
        order = built["order"]
        scopes[key] = {
            "measurement_source": asdict(built["measurement_source"]), "rows": built["rows"],
            "predeclared_contrast": {"minuend_condition_id": items[0], "subtrahend_condition_id": items[1],
                                     "coefficients": coefficients, "offset": 0.0},
            "hierarchy": {
                "n_items": len(order.items), "n_ordered_pairs": len(order.closure),
                "n_incomparable_pairs": len(order.incomparable_pairs),
                "incomparable_pairs": [list(pair) for pair in order.incomparable_pairs],
                "cover_edges": cover_edges(order), "rank_ranges": rank_ranges(order),
                "meaning": "Magnitude hierarchy of native measured exchange readouts across reported growth conditions; ranks are not enzyme importance, parameter values, or causal/temporal precedence",
                "incomparability": "Unknown, not tied. No incomparable pairs under these explicit bands; none are manufactured." if not order.incomparable_pairs else
                                   "Unknown, not tied. Overlap, touching, or missing bands do not establish an order.",
            },
            "domains": {name: domain.report() for name, domain in domains.items()},
            "affine": affine, "sampling": sampling,
            "scaling": {
                "factor": shared_scale, "factor_is_calibration_estimate": False,
                "order_preserved": order.closure == domains["scaled"].order.closure,
                "incomparables_preserved": order.incomparable_pairs == domains["scaled"].order.incomparable_pairs,
                "endpoint_errors": errors, "verified": verified,
                "meaning": "Positive common scaling preserves order and scales the zero-offset affine contrast; it does not identify enzyme abundance E, kcat, or absolute activity. Unknown gain supplies no absolute uptake range. At fixed mu, interpreting q -> gain*q as a physical change would give apparent yield/gain.",
            },
        }
    condition = indexed[YIELD_GROWTH_RATE_PER_H]
    fixed_item = condition["condition_id"]
    forward = evaluate_forward(
        batches["GlucoseUptake"],
        lambda point: apparent_specific_yield(condition["growth_rate_per_h"], point[fixed_item]),
        estimand="apparent specific yield mu / (q_glucose * glucose_molar_mass / 1000) at D=0.10",
        unit="gDW/g_glucose", sign_tolerance=tolerance,
    )
    forward["sampling"] = _domain_reference(forward["sampling"], "#/scopes/GlucoseUptake/domains/bounded")
    for record in forward["records"]:
        record.setdefault("error", None)
        record.setdefault("error_type", None)
    forward.update(
        condition_id=fixed_item, evidence_scope=condition["evidence_scope"],
        growth_rate_per_h=condition["growth_rate_per_h"], growth_is_fixed_setpoint=True,
        molar_mass=MOLAR_MASS_PROVENANCE, is_confidence_interval=False,
        envelope_label="sampled envelope, NOT a global certificate or statistical CI",
        input_role="Only the fixed condition's glucose coordinate enters the callback; the other nine explicitly bounded proposal coordinates are nuisance readouts, not cells/replicates or extra evidence.",
        interpretation="Apparent yield is a nonlinear unit-consistent transformation of mu and measured uptake, not an independently predicted or scored biomass yield; no GEM forward solve.",
    )
    censored = [
        {"condition_id": row["condition_id"], "evidence_scope": row["evidence_scope"],
         "source_column": key, **record, "rounding_band": None, "used_for_order_or_bound": False}
        for row in conditions for key, record in row["observations"].items()
        if record["observation_status"] == "below_detection_limit"
    ]
    n_readouts = sum(len(row["observations"]) for row in conditions)
    complete = all(
        scoped["affine"]["bounded_contrast"]["status"] == "bounded" and
        scoped["affine"]["scaled_contrast"]["status"] == "bounded" and
        all(scoped["affine"][name]["status"] == "unbounded" for name in
            ("order_only_absolute", "order_only_contrast")) and
        scoped["sampling"]["bounded"]["status"] == "complete" and
        scoped["scaling"]["verified"] is True for scoped in scopes.values()
    ) and forward["failed"] == 0 and forward["successful"] == samples
    return {
        "schema_version": 1, "real_native_data": True, "law_ready": False,
        "run_status": "complete" if complete else "incomplete_or_failed",
        "biological_probability": None, "independent_validation": False,
        "source": conditions[0]["source"], "precision_provenance": conditions[0]["precision_provenance"],
        "input_entrypoint": "ystwin.fba.native_reconciliation.native_conditions()",
        "source_audit": {
            "n_conditions": len(conditions), "n_readouts": n_readouts,
            "n_quantified": n_readouts - len(censored), "n_censored": len(censored),
            "n_uptake_readouts_modelled": sum(len(scoped["rows"]) for scoped in scopes.values()),
            "evidence_scope_counts": dict(Counter(row["evidence_scope"] for row in conditions)),
            "censored_observations": censored,
            "scope": "Only already-development native data; interpolation labels retained. Pooling these conditions is not an independent validation split or new law fit.",
        },
        "predeclared_queries": {
            "contrasts_growth_rates_per_h": PREDECLARED_CONTRASTS,
            "yield_growth_rate_per_h": YIELD_GROWTH_RATE_PER_H,
            "selection": "User-requested contrasts fixed in code before data/LP calls, not post-selected signs or a claim of prospective preregistration",
        },
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "solver": "scipy.optimize.linprog / HiGHS"},
        "assumptions": list(ASSUMPTIONS), "scopes": scopes, "nonlinear": forward,
        "gem": {"loaded": False, "bounds_modified": False, "forward_probes": 0},
        "limitations": [
            "These are native measured exchange readouts, not enzyme importance, gene hierarchies, E, kcat, E*kcat, or absolute enzyme activity.",
            "Sharp affine LP ranges/signs are conditional on the declared closed numeric domains. Core certificates are residual-validated numerical SciPy/HiGHS results, not rational exact proofs or biological probabilities.",
            "The order-only closed relaxation supplies neither finite absolute ranges nor a uniform probability distribution; zero contrast is allowed by its weak embedding.",
            "Explicit nearest-rounding boxes add back magnitude assumptions derived from the same evidence. They are not measured uncertainty intervals and do not address systematic assay error.",
            "Glucose increases with growth across these observed conditions; oxygen's magnitude hierarchy need not follow growth. No claim extends outside these conditions or to interventions.",
            "The seeded independent uniform box is only a declared computational proposal. Its definedness is separate from finite-budget execution success, and all rejections/failures remain reported.",
            "The apparent-yield envelope covers successful samples only, not all feasible values or a statistical CI. Full computational digits in JSON are not extra experimental precision; CSV display uses at most six significant digits.",
            "No reserved external/product outcomes, kinetic source data, latent teacher, new experimental calibration, independent prediction score, or law-ready/biology-learning credit is used.",
        ],
    }


def _csv(rows, columns):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _display(value, digits=6):
    return None if value is None else format(value, f".{digits}g")


def summary_csv(report):
    rows = []
    for key, scoped in report["scopes"].items():
        for name, result in scoped["affine"].items():
            domain_name = result["domain_id"].rsplit("/", 1)[-1]
            sampling = scoped["sampling"].get(domain_name, {})
            scope = scoped["domains"][domain_name]["order"]["scope"]
            endpoints = [result[side] for side in ("minimum", "maximum", "feasibility")]
            residual = max(error for endpoint in endpoints for error in endpoint["residuals"].values()) if all(
                endpoint["validated"] for endpoint in endpoints
            ) else None
            rows.append({
                "query": name, "observable_id": key, "method": "conditional affine LP",
                "estimand": result["estimand"], "coefficients": json.dumps(result["coefficients"], sort_keys=True),
                "offset": result["offset"], "fixed_condition_id": None,
                "status": result["status"], "unit": result["unit"],
                "minimum": _display(result["minimum"]["value"]), "maximum": _display(result["maximum"]["value"]),
                "minimum_status": result["minimum"]["status"], "maximum_status": result["maximum"]["status"],
                "minimum_validated": result["minimum"]["validated"], "maximum_validated": result["maximum"]["validated"],
                "maximum_certificate_residual": residual,
                "sign": result["sign"], "tolerance": result["tolerance"],
                "certificate": result["certificate"], "feasibility_status": result["feasibility"]["status"],
                "domain_id": result["domain_id"], "context": scope["context"],
                "calibration_family": scope["calibration"]["family"],
                "n_incomparable_pairs": scoped["hierarchy"]["n_incomparable_pairs"],
                "incomparable_pairs": json.dumps(scoped["hierarchy"]["incomparable_pairs"]),
                "range_assumptions": "none; gap=0 weak relaxation" if domain_name == "order_only" else
                                     "explicit nearest-rounding bands; chosen common rescaling only in scaled domain",
                "sampling_status": sampling.get("status", "not_requested"),
                "requested_proposal": sampling.get("requested_proposal"),
                "proposal_source": sampling.get("proposal_source", {}).get("reference"),
                "proposal_status": sampling.get("proposal_status", "not_requested"),
                "conditioning_status": sampling.get("conditioning_status", "not_requested"),
                "distribution": sampling.get("distribution"), "seed": sampling.get("seed"),
                "budget": sampling.get("budget"), "requested_samples": sampling.get("requested", 0),
                "attempted_samples": sampling.get("attempted", 0), "accepted_samples": sampling.get("accepted", 0),
                "rejected_samples": sampling.get("rejected", 0),
                "successful_samples": 0, "failed_samples": 0, "forward_evaluated": False,
                "guarantee": "Sharp conditional affine range on the declared domain, subject to numerical certificate/tolerance; not a biological probability. Sample counts describe proposals; no callback for this LP row.",
            })
    forward = report["nonlinear"]
    sampling = forward["sampling"]
    scoped = report["scopes"]["GlucoseUptake"]
    scope = scoped["domains"]["bounded"]["order"]["scope"]
    envelope_status = "sampled_only" if forward["successful"] else "unavailable_no_successful_samples"
    rows.append({
        "query": "apparent_specific_yield", "observable_id": "GlucoseUptake", "method": "sampled nonlinear callback",
        "estimand": forward["estimand"], "coefficients": None, "offset": None,
        "fixed_condition_id": forward["condition_id"],
        "status": forward["sign_stability"], "unit": forward["unit"],
        "minimum": _display(forward["sampled_minimum"], 4), "maximum": _display(forward["sampled_maximum"], 4),
        "minimum_status": envelope_status, "maximum_status": envelope_status,
        "minimum_validated": False, "maximum_validated": False, "maximum_certificate_residual": None,
        "sign": forward["successful_sample_sign"], "tolerance": forward["sign_tolerance"],
        "certificate": "none for nonlinear envelope", "feasibility_status": sampling["feasibility"]["status"],
        "domain_id": sampling["domain_id"], "context": scope["context"],
        "calibration_family": scope["calibration"]["family"],
        "n_incomparable_pairs": scoped["hierarchy"]["n_incomparable_pairs"],
        "incomparable_pairs": json.dumps(scoped["hierarchy"]["incomparable_pairs"]),
        "range_assumptions": "explicit nearest-rounding glucose bands; fixed mu=0.10/h; conventional glucose mass 180.156 g/mol",
        "sampling_status": sampling["status"], "requested_proposal": sampling["requested_proposal"],
        "proposal_source": sampling["proposal_source"]["reference"], "proposal_status": sampling["proposal_status"],
        "conditioning_status": sampling["conditioning_status"], "distribution": sampling["distribution"],
        "seed": sampling["seed"], "budget": sampling["budget"], "requested_samples": sampling["requested"],
        "attempted_samples": sampling["attempted"], "accepted_samples": forward["accepted"],
        "rejected_samples": sampling["rejected"],
        "successful_samples": forward["successful"], "failed_samples": forward["failed"], "forward_evaluated": True,
        "guarantee": forward["envelope_label"] + "; failed samples retained in counts/records",
    })
    for row in rows:
        row.update(source_doi=report["source"]["doi"], source_table_sha256=report["source"]["table_sha256"],
                   assumptions=" | ".join(report["assumptions"]), details_file=REPORT_NAME)
    return _csv(rows, list(rows[0]))


def observations_csv(report):
    rows = []
    for key, scoped in report["scopes"].items():
        scope = scoped["domains"]["bounded"]["order"]["scope"]
        for record in scoped["rows"]:
            band = record["rounding_band"] or {}
            rows.append({
                **{name: record[name] for name in (
                    "condition_id", "growth_rate_per_h", "evidence_scope", "evidence_id", "source_column",
                    "unit", "primary_printed_token", "reported_value", "value", "observation_status",
                    "detection_limit", "rounding_half_width", "rounding_is_confidence_interval", "band_status",
                    "source_doi", "source_table_sha256", "source_table_relative_path", "source_line",
                )},
                "band_lower": band.get("decimal_lower"), "band_upper": band.get("decimal_upper"),
                **scoped["hierarchy"]["rank_ranges"][record["condition_id"]],
                "n_incomparable_pairs": scoped["hierarchy"]["n_incomparable_pairs"],
                "incomparable_with": json.dumps([
                    other for pair in scoped["hierarchy"]["incomparable_pairs"]
                    if record["condition_id"] in pair for other in pair if other != record["condition_id"]
                ]),
                "context": scope["context"], "calibration_family": scope["calibration"]["family"],
                "assumptions": "Nearest rounding only; not a CI. Strictly disjoint explicit bands only; same evidence as bounds, not independent support. Ranks are not magnitudes.",
                "result_kind": "source/rank audit, not a sampled or LP result; certificates/proposals/failures in details_file",
                "details_file": REPORT_NAME,
            })
    return _csv(rows, list(rows[0]))


def main(argv=None):
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=root / "outputs")
    parser.add_argument("--samples", type=int, default=256)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--budget", type=int, default=4096)
    parser.add_argument("--tolerance", type=float, default=1e-8)
    parser.add_argument("--shared-scale", type=float, default=7.0)
    args = parser.parse_args(argv)
    output = args.output_dir.resolve()
    if not output.is_dir() or not output.is_relative_to((root / "outputs").resolve()):
        parser.error("output directory must already exist within this project's outputs directory")
    paths = [output / name for name in (REPORT_NAME, SUMMARY_NAME, OBSERVATIONS_NAME)]
    if any(path.exists() or path.is_symlink() for path in paths):
        parser.error("refusing to overwrite existing native_order_* results; use a fresh outputs subdirectory")
    try:
        report = build_demo(samples=args.samples, seed=args.seed, budget=args.budget,
                            tolerance=args.tolerance, shared_scale=args.shared_scale)
        payloads = [json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
                    summary_csv(report), observations_csv(report)]
    except (ValueError, TypeError) as exc:
        parser.error(str(exc))
    for path, payload in zip(paths, payloads):
        try:
            with path.open("x", encoding="utf-8", newline="") as stream:
                stream.write(payload)
        except FileExistsError:
            parser.error(f"refusing to overwrite existing result: {path}")
    print("REAL native development readouts; conditional demonstration only, not a law-ready claim.")
    print("LP certificates: numerical SciPy/HiGHS, conditional on explicit nearest-rounding bands, "
          "NOT statistical CIs or biological probabilities. Orders and bounds are the same evidence.")
    print(f"Source audit: {report['source_audit']['n_censored']} censored zeros, no numeric detection limits; "
          "none used as zero anchors or order bounds. No GEM loaded or modified.")
    for key, scoped in report["scopes"].items():
        result, batch = scoped["affine"]["bounded_contrast"], scoped["sampling"]["bounded"]
        print(f"{result['estimand']}: [{_display(result['minimum']['value'])}, {_display(result['maximum']['value'])}] "
              f"{result['unit']}, {result['status']}/{result['sign']}; "
              f"incomparables={scoped['hierarchy']['n_incomparable_pairs']}; "
              f"proposal={batch['proposal_status']}/{batch['conditioning_status']}, execution={batch['status']}, "
              f"{batch['accepted']}/{batch['requested']} accepted, {batch['rejected']} rejected.")
        order_only = scoped["affine"]["order_only_contrast"]
        absolute = scoped["affine"]["order_only_absolute"]
        print(f"  Order-only: absolute range {absolute['status']}; contrast {order_only['status']}/{order_only['sign']}; "
              f"proposal={scoped['sampling']['order_only']['proposal_status']}. "
              f"Shared-scale x{scoped['scaling']['factor']:g}: covariance verified={scoped['scaling']['verified']}.")
    forward = report["nonlinear"]
    print(f"Apparent yield at {forward['condition_id']}: sampled envelope "
          f"[{_display(forward['sampled_minimum'], 4)}, {_display(forward['sampled_maximum'], 4)}] {forward['unit']}; "
          f"{forward['successful']} successful, {forward['failed']} failed; {forward['sign_stability']}; "
          f"seed={args.seed}, budget={args.budget}; NOT a global certificate/statistical CI.")
    for path in paths:
        print(path)
    return 0 if report["run_status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
