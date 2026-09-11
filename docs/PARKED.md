# Parked: the product layer — mostly un-parked, August 2026

**This page's premise no longer holds and the page is kept to record that.** It said "no
carotenoid pathway is in any strain yet, so nothing below can be fitted, tested or falsified
against data". That is still true of *our* strains, and it was never true of the question:
Elizondo 2025 published three strains at two dilution rates with the pathway in them, and
those six states fit, test and falsify most of what is listed below.

What was parked was the layer. It is now assembled and the whole chain runs forward —
genotype and environment in, product out — scored leave-one-strain-out at 14.2% median error
with no product measurement of the held-out strain. `fba/carotenoid.py` is reached from the
live chain through `fba/audit.py`; `kinetic/carotenoid.py`'s fit lives on in
`pathway/calibrations.py`. See `docs/CLAIM_BOUNDARY.md`.

The rows below stay because each one settled something, and two of them settled it by being
wrong in a way worth keeping visible.

| Parked | What it settled |
| --- | --- |
| `fba/carotenoid.py` | ~~parked~~ **un-parked**: crtE / crtYB(PSY) / crtI / crtYB(LCY) installed into Yeast9 or GECKO, mass balanced, named for the real construct genes. Now reached from the live chain by `fba/audit.py`. |
| `fba/fva.py` | D1: the feasible product flux is `[0, ceiling]` at every growth level — FBA bounds the product and never predicts it. A capacity bound moves the ceiling and never lifts the floor off zero. |
| `fba/physiology.py` | ~~Plain Yeast9 overpredicts growth by 73%~~ **corrected**: it *under*predicts by 13.5% against van Hoek's own Table 1, and both models pass at 35%. The 73% and the 65% oxygen miss were artefacts of an unsourced reference and of `cap_uptake` being a no-op on the split ec model. See FINDINGS and `docs/superseded/reference-phenotype.md`. |
| `fba/surrogate.py` | LP emulator, 161x faster per call and 148,000x batched, median error 0.000%. |
| `kinetic/carotenoid.py` | ~~parked~~ **un-parked**: fitted to six measured chemostat steady states (Elizondo & Saa 2025, PMID 40891387). docs/research/KINETIC_FIT.md |
| `observation.py` → inner filter | Beta-carotene absorbs at 450 nm, on mCitrine's excitation. Inert today (no pigment) and refuses to correct without a measured coefficient. |
| `scripts/parked/run_d1.py` | Reproduces the D1 numbers. |

Tests for all of it stay in the suite so it does not rot. Nothing here is imported by
the sensor-characterisation path, which remains a separate argument on separate data.

**Still parked:** `fba/surrogate.py` (nothing calls the emulator yet), `scripts/parked/run_d1.py`,
and the inner-filter correction in `observation.py`, which is inert until a strain of ours
carries a pigment.

## Unparking checklist

1. Measure the inner-filter coefficient with a purified-carotenoid spike-in on a
   non-producing strain (`docs/PROTOCOLS.md` P3), then set `ReporterOptics.inner_filter_coeff`.
2. Add phytoene (285 nm) and lycopene (470 nm) to the beta-carotene HPLC run — same
   extraction, same injection. The intermediate ratios are what make the kinetic branch
   identifiable, and they cost almost nothing to collect.
3. Re-run `scripts/parked/run_d1.py` against the actual strain's construct manifest.
