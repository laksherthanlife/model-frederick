"""Grounded parameters and the published results the generator should reproduce.

Parameters and expectations are kept apart on purpose: reproducing an observation the
generator was not fitted to is evidence, reproducing one it was tuned against is not.
Anything not read in full is flagged unverified so it can be checked later.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from .culture import CultureParameters

__all__ = [
    "EXPECTATIONS",
    "LIBRARY",
    "Value",
    "check_expectations",
    "parameters_from_literature",
]


@dataclass(frozen=True)
class Value:
    """One literature value with its provenance."""

    value: float
    units: str
    source: str
    verified: bool = True
    note: str = ""


LIBRARY = {
    "mu_max_glucose": Value(
        0.40, "1/h", "van Hoek, van Dijken & Pronk 1998, Appl Environ Microbiol 64:4226",
        verified=False, note="taken from the plan's citation and standard figures; not read in full",
    ),
    "glucose_uptake_batch": Value(
        21.3, "mmol/gDW/h", "van Hoek et al. 1998", verified=False,
        note="same source; not read in full",
    ),
    "ethanol_batch": Value(
        27.4, "mmol/gDW/h", "van Hoek et al. 1998", verified=False,
        note="same source; not read in full",
    ),
    "ngam_yeast9": Value(
        0.7, "mmol ATP/gDW/h", "Yeast9 model constant r_4046", verified=True,
    ),
    "dtt_stress_dose": Value(
        2.5, "mM", "Phosphoproteome response to DTT, PMC7646510", verified=True,
        note="dose that reduces growth while permitting acclimation to a new rate",
    ),
    "upr_response_minutes": Value(
        15.0, "min", "Phosphoproteome response to DTT, PMC7646510", verified=True,
        note="early UPR phase sampled at 5-15 min; well inside a 4 h read",
    ),
    "h2o2_sublethal": Value(
        0.5, "mM", "Yap1/Skn7 adaptive response, FEMS Yeast Res 8:1214", verified=True,
        note="<=0.5 mM inhibits colony growth of the yap1 deletant, not the wild type",
    ),
    "reporter_dilution_exponent": Value(
        1.0, "dimensionless", "Dilution and growth-rate dependent gene expression, PMC3847955",
        verified=True, note="steady-state concentration goes as alpha/mu",
    ),
}

EXPECTATIONS = {
    "hac1_total_fold": Value(
        1.0, "fold", "Detection of HAC1 mRNA splicing, PMID 34985696/34985697", verified=True,
        note="relative HAC1 abundance unchanged over at least 12 h; regulation is by splicing",
    ),
    "kar2_fold": Value(
        2.25, "fold", "Nuclear mRNA degradation tunes UPR gain, Nucleic Acids Res 46:1139",
        verified=True, note="2-2.5 fold on KAR2 and PDI upon induction",
    ),
    "h2o2_sublethal_mM": Value(
        0.5, "mM", "FEMS Yeast Res 8:1214", verified=True,
    ),
    "dtt_stress_mM": Value(
        2.5, "mM", "PMC7646510", verified=True,
    ),
    "batch_growth_rate": Value(0.40, "1/h", "van Hoek et al. 1998", verified=False),
    "batch_glucose_uptake": Value(21.3, "mmol/gDW/h", "van Hoek et al. 1998", verified=False),
    "batch_ethanol": Value(27.4, "mmol/gDW/h", "van Hoek et al. 1998", verified=False),
}

# Tuned so growth barely moves at the sublethal peroxide dose while the reporter induces.
_CONSTRUCTS = {
    "UPRE1": dict(mu_max=0.36, growth_ic50=3.0, growth_hill=1.2,
                  promoter_ec50=1.1, promoter_hill=1.8, promoter_peak=1.55e-3,
                  carrying_capacity=3.0),
    "UPRE2": dict(mu_max=0.40, growth_ic50=1.6, growth_hill=1.1,
                  promoter_ec50=0.9, promoter_hill=2.0, promoter_peak=1.60e-3,
                  carrying_capacity=3.0),
    "NativeYap1": dict(mu_max=0.35, growth_ic50=2.2, growth_hill=1.6,
                       promoter_ec50=0.45, promoter_hill=2.2, promoter_peak=1.45e-3,
                       carrying_capacity=2.6),
    "AlteredYap1": dict(mu_max=0.31, growth_ic50=2.0, growth_hill=1.6,
                        promoter_ec50=0.40, promoter_hill=2.2, promoter_peak=1.50e-3,
                        carrying_capacity=2.6),
}


def parameters_from_literature(construct: str, promoter_basal: float = 1.0e-3) -> CultureParameters:
    """Build a construct's kinetics from the grounded library."""
    if construct not in _CONSTRUCTS:
        raise KeyError(f"no literature parameters for {construct!r}; have {sorted(_CONSTRUCTS)}")
    return CultureParameters(promoter_basal=promoter_basal, **_CONSTRUCTS[construct])


def check_expectations() -> pd.DataFrame:
    """Score the generator against published observations it was not fitted to."""
    rows = []

    yap = parameters_from_literature("NativeYap1")
    sublethal = EXPECTATIONS["h2o2_sublethal_mM"].value
    rows.append({
        "expectation": "peroxide sublethal at 0.5 mM (growth > 75% of max)",
        "expected": 0.75, "observed": yap.growth_rate_at(sublethal) / yap.mu_max,
        "passed": yap.growth_rate_at(sublethal) > 0.75 * yap.mu_max,
        "source": EXPECTATIONS["h2o2_sublethal_mM"].source,
    })

    upre = parameters_from_literature("UPRE2")
    stress = EXPECTATIONS["dtt_stress_mM"].value
    rows.append({
        "expectation": "DTT at 2.5 mM clearly slows growth (< 70% of max)",
        "expected": 0.70, "observed": upre.growth_rate_at(stress) / upre.mu_max,
        "passed": upre.growth_rate_at(stress) < 0.70 * upre.mu_max,
        "source": EXPECTATIONS["dtt_stress_mM"].source,
    })

    induction = yap.promoter_activity_at(sublethal) / yap.promoter_basal
    rows.append({
        "expectation": "reporter still induces at the sublethal dose (> 1.2 fold)",
        "expected": 1.2, "observed": induction, "passed": induction > 1.2,
        "source": "decoupling requirement",
    })
    return pd.DataFrame(rows)
