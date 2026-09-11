"""Fitted constants for the shipped pathways, each beside the data that produced it.

Kept out of the TOML specs on purpose. A spec is chemistry, which is a property of the
pathway; these are fits, which are properties of a dataset -- and a fitted number is only
meaningful next to what it was fitted on and how it scored held out.

The fixed-gene flux calibration is reproduced by ``scripts/fit_pathway_flux.py`` and
checked in ``tests/test_pathway_flux.py``. Branch kinetics are a view of ELIZONDO2025,
whose fit belongs to ``scripts/fit_carotenoid_kinetics.py``. Neither shipped full-data
fit is used to select genes or fit parameters in nested product validation.
"""

from __future__ import annotations

from ..kinetic.carotenoid import ELIZONDO2025
from .flux import FluxCalibration
from .solve import CeilingRefutation, ContentCeiling, NodeKinetics

__all__ = ["ARHAR2024_BETA_CAROTENE", "BETA_CAROTENE_FLUX", "BETA_CAROTENE_KINETICS"]

# Relative CrtE expression -> total pathway flux, from Elizondo 2025's three strains at two
# dilution rates. Parameter fitting is LEAVE-ONE-STRAIN-OUT, but the gene was selected on
# this same dataset. This is a fixed-model comparison, not a prespecified nested score.
#
# Against training-only geometric means the fourteen-gene exploratory comparison gives
# CrtE +0.745, CrtYB +0.708, CrtI +0.359, and ERG9 +0.003. The remaining native genes
# range from -0.003 to -0.268. The gap does not establish a causal cassette-dosage law or
# exclude expression confounding. scripts/predict_product.py selects the gene AND model
# inside each outer training fold and reports three strain groups, six predictions.
BETA_CAROTENE_FLUX = FluxCalibration(
    alpha=1.008474e-3,
    entry_enzyme="CrtE",
    loso_rmse_log=0.2002,
    loso_skill=0.7447,
    n_states=6,
    source="PMID 40891387 (doi 10.1021/acssynbio.5c00256), Elizondo 2025; fitted by "
           "scripts/fit_pathway_flux.py",
    # The relative CrtE expressions actually seen across the six states: b-car4 at D=0.101
    # is the normalisation reference at 1.000, and b-car2 at D=0.254 is the lowest at 0.245.
    # A 4.08-fold span. The extrapolation warning fires against this rather than against a
    # literal window -- the previous hard-coded ceiling of 5.0 sat 5x above anything ever
    # measured, so no realistic input could reach it.
    expression_range=(0.245160473, 1.0),
)

ARHAR2024_BETA_CAROTENE = CeilingRefutation(
    measured_mg_per_gdcw=79.0,
    analyte="beta_carotene",
    assay="HPLC on gravimetric dry weight, lycopene below detection",
    source="PMID 39215465 (Arhar 2024, J Appl Microbiol 135:lxae224); row 39215465 of "
           "data/carotenoid_batch/published_batch_titres.tsv",
)
"""The measurement that refutes the ceiling ELIZONDO2025 implies, about 63-fold.

THE ANALYTE AND THE ASSAY ARE THE LOAD-BEARING FIELDS, not the number. Lopez 2019 reports
21 mg/gDCW, 17x the same bound, and refutes nothing: the survey's `measurand` column calls
it ``total_carotenoid`` because the paper reads absorbance at 453 nm against a beta-carotene
standard, which sums every carotenoid in the extract. Arhar's row is ``beta_carotene``, by
separation, with lycopene below detection -- so it is the smaller claim and the only one that
lands. `tests/test_carotenoid_ceiling_scope.py` cross-checks all four values against the
vendored survey row rather than trusting this copy."""

# The lycopene cyclase. vmax_per_growth and km are the Elizondo six-state fit, unchanged --
# `kinetic/carotenoid.py::ELIZONDO2025` holds the same two numbers with their intervals and
# leave-one-out score.
#
# vmax scaling with mu is BOUNDED AND NOT PINNED: the 95% profile interval on the exponent
# is [0.53, 1.78] against a point estimate of 1.04. Zero and three are excluded, so the law
# is not arbitrary; 0.6 and 1.7 are not, so it is not established either. See
# tests/test_growth_exponent_identifiability.py.
BETA_CAROTENE_KINETICS: dict[str, NodeKinetics] = {
    "lycopene": NodeKinetics(
        # The value is ELIZONDO2025's, bit for bit. What the wrapper adds is the refutation,
        # which `content_ceiling` hands on to whoever reads the bound this constant sets.
        vmax_per_growth=ContentCeiling(ELIZONDO2025.capacity_mmol_per_gdcw,
                                       refuted_by=ARHAR2024_BETA_CAROTENE),
        km=ELIZONDO2025.km_mmol_per_gdcw,
        growth_rate_range=ELIZONDO2025.growth_rate_range,
    ),
}
"""Read from ELIZONDO2025 rather than copied out of it, and LEGACY-ONLY.

The literals were duplicated here for about an hour, which is exactly long enough for a
refit to update one copy and not the other. There is one fit and it lives in
`kinetic/carotenoid.py` with its intervals, its leave-one-out score and its citation; this
is a view of it in the generic solver's vocabulary.

LEGACY-ONLY, AND THAT IS NOW A PROPERTY OF THE OBJECT RATHER THAN OF THIS PARAGRAPH. On the
step feeding a DILUTED terminal node, `vmax_per_growth` is not only a rate constant: `mu`
cancels out of `content <= vmax_per_growth`, so the same number is a hard asymptote on
beta-carotene content at every growth rate for every genotype -- 1.2483 mg/gDCW. Arhar 2024
measured 79. A content solved against this calibration is therefore a statement about the
fit and not about the strain, whatever mode produced it, and `content_ceiling(spec, kinetics)`
now returns a value whose `legacy_only` says so and whose `describe()` spells it out.

WHY IT IS NOT REFITTED INSTEAD. The repair is structural, not numerical. The capacity is a
property of the crtYB cyclase; the only genotype input the package accepts is CrtE
expression, which sets the flux INTO the chain, and raising it a hundredfold moves the answer
by 0.06%. The six calibration states differ in expression, not in cassette dosage, so the
only refit they support is a capacity scaled by a measured relative mRNA -- and on CrtE,
CrtYB and CrtI alike that fits WORSE than the flat capacity and still lands within 1.4x of
the same refuted bound. `tests/test_carotenoid_ceiling_scope.py` runs all three and fails if
any of them ever beats flat. What would settle it is content on three or more strains
differing in crtYB dosage alone; `pathway/capacity.py` refuses until then, and
`pathway/published_cassettes.py` computes from the survey why the literature cannot stand in
for that experiment."""
