# Generalisation: what it would take to stop being one team's project

**Question asked.** The owner's instruction was *stop calling it a beta-carotene twin — it needs
to be generalizable.* That is a scope decision, so this note works out what generalising
actually requires: which code is already host- and product-agnostic, which is not, what seam
would carry the difference, what it would have to be demonstrated on, and what breaks.

Compiled 2026-08-26 against `src/ystwin` at that date. Two files moved under it while it was
being written — `pyproject.toml`'s `description` is now *"Growth-corrected promoter activity
from reporter and optical-density time series"* and `src/ystwin/__init__.py` carries a
matching module docstring. Both changes are in the direction argued for here. Line numbers are
pointers, not identifiers.

---

## 1. Verdict

The measurement core is **already general and nobody has said so out loud**. `reporter.py`,
`growth.py`, `analysis/uncertainty.py`, `analysis/nulls.py` and the gate *shapes* contain no
yeast, no carotenoid and no instrument: they invert `dR/dt = k_synth − (μ + k_deg)·R` on any
(time, biomass, reporter) triple. Roughly 30 % of the package is host-specific, and it is
concentrated in three files — `generator/stress_panel.py`, `plate/synergy.py`,
`fba/physiology.py` — plus a scattering of frozen defaults (`gdcw_per_od = 0.42`,
`od_linear_max = 1.0`, `_MU_MAX = 0.40`, `k_deg = 0`).

The single change that generalises the most is **not a rename and not a `HostProfile`**. It is
lifting `scripts/run_sensor_characterisation.py`'s middle 60 lines into
`src/ystwin/pipeline.py::recover_activity(aligned, biomass_channel, reporter_channel, ...)`.
Today the validated core is only reachable through one exporter, one plate map and four
construct names. Make it callable on a tidy frame and every other item on this list becomes
cheap.

---

## 2. The classification table

**A finding before the table.** The brief's four categories have no slot for the largest
group. `plate/layout.py`, `qpcr.py`, `viz/figures.py` and `analysis/splits.py`'s validation
hooks are not yeast-specific and not instrument-specific — they are specific to *this
experiment*: four constructs (`UPRE1`, `UPRE2`, `NativeYap1`, `AlteredYap1`), two agents, one
plate design, one logbook keyed by date. That is a fifth category, marked **S** below, and it
matters because it is cheap to generalise (parameterise a dict) whereas **H** is expensive
(re-transcribe a literature model) and **I** is medium (write a second parser).

Categories: **G** general now · **H** host-specific, replaceable by a per-host profile ·
**P** product-specific (β-carotene) · **I** instrument-specific · **S** study-specific (this
lab's constructs, agents, plate design).

### 2.1 Core — the part that already generalises

| Module | Cat | What would have to change |
| --- | --- | --- |
| `reporter.py` | **G** | Nothing structural. `ReporterKinetics(k_deg=0.0, k_mat=None)` defaults encode "stable FP, instantaneous maturation", which is an *observation about mCitrine in yeast*, not a default. Make both required, or make the profile supply them. `default_activity_window_h`'s error-vs-window table (0.5 h → 14.1 %, 4 h → 1.1 %) was measured at 10-min sampling over 24 h with yeast-scale μ; the *shape* is general, the calibration is not. |
| `growth.py` | **G** | Argument names and one error message. `specific_growth_rate(optical_density=…)` only needs a strictly-positive biomass series; rename to `biomass` and the module works on cell counts or confluence unchanged. `_DEFAULT_DETECTION_OD = 0.02` and `_DEFAULT_RATE_QUANTILE = 0.9` are the two host/instrument constants; `window_h = 1.5` is a yeast-doubling-time choice (see §5.1). |
| `analysis/uncertainty.py` | **G** | Nothing. `MIN_PLATES_FOR_INTERVAL = 3`, the cluster bootstrap and the TOST equivalence are statistics. The one binding assumption is that the cluster is *the plate*, which holds for any plate-reader assay and not for a chemostat or a microfluidic array. |
| `analysis/nulls.py` | **G** | Nothing. Phase randomisation, shuffled loadings, rotated subspaces, `skill_score` — pure method. |
| `analysis/latent.py` | **G** | Nothing. Masked SVD, held-out dimension selection. The Ledermann ceiling it is bound by is arithmetic. |
| `analysis/transfer.py`, `recovery.py`, `matched_density.py`, `archive/code/deconvolve.py` | **G** | Nothing structural; all operate on generic reading matrices. `matched_density.py` compares at matched OD, which needs a biomass proxy, not OD specifically. |
| `estimator.py` | **G** (one **P** field) | `Observation.carotenoid` (default `0.0`) flows into `observe_rfu` → `ReporterOptics.attenuation`, so the likelihood carries a permanently inert product term. Remove the field; the filter is otherwise host-agnostic. Also needs the per-OD ↔ per-gDCW reconciliation (`ARCHITECTURE.md` §6.3) before it can run on real data at all. |
| `calib/od.py`, `calib/od_linearity.py` | **G** shape, **I** content | `ODCalibration(saturation_k, top_true_od, blank, gdcw_per_od, residual_std)` *is* the reader-response profile. Nothing to change except that `gdcw_per_od` is per-host and per-medium, and today no measured value exists for anything. Two functions named `fit_od_linearity` in two modules with different return types should be resolved before this becomes a public seam. |
| `calib/kdeg.py` | **G** | Nothing. It is the machinery that would fill `ReporterProfile.k_deg`, and it has never been run because the chase experiment has not been done. |
| `gates/g1_optical.py` | **G** shape, **I** defaults | `OpticalQualityGate` is the right pattern — frozen thresholds carried on the result. Every default is instrument- and organism-bound: `od_linear_max = 1.0` (placeholder, docstring admits it), `min_reporter_above_background = 50.0` (raw RFU at one PMT gain — not portable between two gains on the *same* reader), `min_od_above_blank = 0.02`, `min_od_fold = 1.15`, `max_decline_fraction = 0.15`. The fix is to delete the defaults, not to retune them (§5.5). |
| `diagnostics/dilution_confound.py` | **G** | `_MIN_GROWTH_FOLD = 1.15` refuses a non-growing culture; `_MODEL_INADEQUATE_FRACTION = 0.25`; the QSS test `loss > 1/(duration/2)`. All three are batch-yeast-scale (§5.1). The docstring's "physiological maximum for yeast" is a host fact used as a sanity bound and should become `HostProfile.mu_max`. |
| `analysis/power.py` | **G** with a live bug | `_one_trial` inverts dilution as `recovered = signal * growth`, i.e. `k = R·μ` — correct **only when `k_deg == 0`**, which `parameters_from_literature` defaults to. Point this at a degron reporter and it is silently wrong. `DEFAULT_WELL_CV = 0.052` is a measured yeast-plate number; `_DURATION_H = 60.0`, `_N_TIMEPOINTS = 240` are this run design. |
| `analysis/design.py`, `sensor_selection.py`, `experiment_design.py` | **G** machinery, **H** inputs | The Fisher-information and D-optimality maths is general; every input (`reporter_loadings`, `MODULES`, `REFERENCE_GROWTH_RATE = 0.22`, `_GROWTH_CV` derived from `MEASURED_GROWTH_RATE_SE / 0.22`) is yeast. `Reporter.demonstrated_in_yeast: bool` is the field that decides `RECOMMENDED_BUILD` and it is a boolean where a taxon set belongs. |
| `analysis/splits.py` | **G** machinery, **S** hooks | The `SplitKind` machinery (interpolation / extrapolation / construct- and stressor-holdout, group-aware, hashed, manifested) is domain-general and genuinely good. `_check_construct_stressor_agreement` validates against `qpcr.STRESSOR_FOR_CONSTRUCT`; that hook should take a mapping, not import one. |
| `paths.py` | **G** mechanism, **S** resolvers | The env-var → in-repo → sibling-dir contract is general and worth keeping verbatim. All nine resolvers name this lab's assets, and the two `Path.home()/"Desktop"` fallbacks are one machine's. Nine `YSTWIN_*` variables are part of the user-facing contract. |

### 2.2 Host-specific, replaceable by a per-host reference

| Module | Cat | What would have to change |
| --- | --- | --- |
| `generator/stress_panel.py` | **H** | The biggest single item. 24 `MODULES` with `driven_by` weights and PMIDs; `_CASCADE_ORDER` (24 entries, drivers-first, unchecked); 25 `STRESSORS` with `targets` and `ec50`; 24 `REPORTERS` with `also_reads` loadings. Plus scalars: `_MU_MAX = 0.40`, `_LETHAL_MULTIPLE = 3.0`, `_LETHAL_HILL = 2.5`, `_POOL_POTENCY = 0.35`, `_GENERAL_POTENCY = 2.5`, `_DEFAULT_BASAL = 0.9` (set from UPRE1/UPRE2 measurements), `CYTOSOLIC_PH = (6.4, 7.4)`. Porting this to a second host is re-transcribing a literature model, not editing a config — see §5.4. |
| `generator/context.py` | **H** | `_CARBON` (glucose/galactose/ethanol μ and derepression), `_PHASE`, `_REFERENCE_TEMPERATURE = 30.0`, `_HEAT_PER_DEGREE = 0.06`, `_HYPOXIA_THRESHOLD = 0.21`. **`CultureContext.strain: str = "BY4741"` is declared and never read anywhere** — it is the natural attachment point for a `HostProfile` and costs nothing to repurpose. |
| `fba/physiology.py` | **H** | Yeast9 reaction ids `r_2111` / `r_1714` / `r_1992` / `r_1761` / `r_1672` (each redefined in three modules), and `REFERENCE_AEROBIC_BATCH` (van Hoek 1998). A second host means a second GSMM, second id set, second measured phenotype. |
| `bridge/thermodynamic.py`, `generator/redox.py` | **H** | `STANDARD_TEMPERATURE_K = 303.15` / `_TEMPERATURE_K = 303.15` — 30 °C is a *yeast growth temperature*; E. coli and mammalian work is 37 °C (310.15 K), which moves every ΔG. `CYTOSOLIC_PH = 7.5` vs `_REFERENCE_PH = 7.0` vs `stress_panel.CYTOSOLIC_PH = (6.4, 7.4)` — one name, three meanings, all host facts. |
| `bridge/physiology_bridge.py`, `latent_bridge.py`, `tmfa.py`, `equilibrator.py`, `archive/code/module_admission.py` | **H** | Downstream of the GSMM. `_MAINTENANCE_PER_ACTIVITY = 300.0` (overridden to `1200.0` in `run_scenarios.py`) is Tier 0 in either host. |
| `generator/culture.py` | **G** shape, **H** parameters | The mechanistic ODE is general. Its parameters arrive from `literature.py` and are per-construct. |
| `generator/panel_experiment.py` | **H** | Holds the three measured yeast constants — `OBSERVED_ACTIVITY_CV = 0.146`, `MEASURED_GROWTH_RATE_SE = 0.0117`, `MEASURED_ACTIVITY_CV = 0.14`. These are *this lab's plates*; a second host needs its own three, and `sensor_selection` reads one of them at import time to compute `RECOMMENDED_BUILD`. |
| `generator/kinetics.py`, `plate_layout.py`, `unmixing.py`, `panel_gates.py`, `design.py` | **G** | Time-to-steady-state, edge-evaporation multiplier, spectral unmixing, sign/diversity gates, decoupling grid. All host-neutral. |
| `generator/literature.py` | **S** on **H** foundations | `LIBRARY` and `EXPECTATIONS` are yeast physiology values; `_CONSTRUCTS` is the four constructs. **The `Value(value, units, source, verified, note)` dataclass is the right shape for every profile field in §3** and should be reused rather than reinvented. |

### 2.3 Product-specific (β-carotene) — all currently parked

| Module | Cat | Note |
| --- | --- | --- |
| `fba/carotenoid.py` | **used** | crtE / crtYB(PSY) / crtI / crtYB(LCY) into Yeast9 or GECKO. Supplies the pathway the kinetic layer is fitted against. |
| `fba/fva.py` | **P** | D1 — feasible product flux `[0, ceiling]`, relative width 1.000. |
| `fba/surrogate.py` | **G** code, **P** justification | A generic grid interpolator that refuses to extrapolate. Its docstring justifies it by a particle-filter/FBA coupling that does not exist. |
| `kinetic/carotenoid.py` | **P** | GGPP → phytoene → lycopene → β-carotene, shared bifunctional crtYB pool. |
| `observation.py` → inner filter | **P** | `ReporterOptics.inner_filter_coeff`, `.attenuation()`, `correct_inner_filter()`, `observe_rfu(carotenoid=)`. The rest of `observation.py` is **G**. |
| `scripts/parked/run_d1.py` | **P** | Reproduces the D1 numbers. |

Total genuinely β-carotene-specific surface: **three modules, one dataclass field, four
functions, one script.** This is the smallest of the five categories by a wide margin. The
project is far less carotenoid-shaped than its own name implies — the naming problem is
almost entirely a *labelling* problem, not a code problem.

### 2.4 Instrument-specific

| Module | Cat | What would have to change |
| --- | --- | --- |
| `plate/synergy.py` | **I** | `WELL_RE = ^([A-H])(\d{1,2})$` — 96-well only (384-well is A–P, 1–24), and `\d{1,2}` accepts "13"–"99" so there is no format check either. `_parse_channel_header` parses `T° mCitrine:480,530[2]`, a BioTek Gen5 string. `_SHEET_CHANNEL_ALIASES` is nine hard-coded fluorophore spellings. `_to_hours` has five branches including *bare float = Excel serial day fraction × 24* and `"30:00"` = MM:SS = 0.5 h — both correct here, both silently 60× or 24× wrong for another exporter. `blank_subtracted = nanmin(rows) < 0` is a heuristic standing in for a fact the exporter knows (§5.6). |
| `generator/synergy.py` | **I** | A *second* Synergy reader with different assumptions, used only for the ATP-sensor plates. Its existence in-tree is a free test of the reader seam (§7). |
| `plate/dose_response.py` | **I** + **S** | The only place time arrives in **minutes**; paired with `generator/export.py`'s `time_h * 60.0`, which must stay paired. |
| `observation.py` → `ReporterOptics` | **I** | `gain`, `background`, `detector_max` are reader+gain; `autofluorescence` is host+strain. One dataclass, two profiles' worth of fields. |

### 2.5 Study-specific (the missing fifth category)

| Module | Cat | What would have to change |
| --- | --- | --- |
| `plate/layout.py` | **S** | `NEWPROTOCOL_LAYOUT`, `_STANDARD` (a second identical copy nothing keeps in sync), `_DOSES`, `RECORDED_PLATES` keyed by export date. The `recover_layout` *algorithm* is general and rather good; the constants are one plate design. |
| `qpcr.py` | **S** | `TARGET_FOR_CONSTRUCT` (HAC1, TRX2), `STRESSOR_FOR_CONSTRUCT`, `_CONSTRUCT_ALIASES`, `QPCR_PLATE`. UBC4 as a DTT-invariant reference gene is a yeast claim. The ddCq maths, the RT⁻ QC and `gdna_share` are **G**. |
| `gates/g4_anchor.py` | **G** | The anchor-agreement logic is general: "does this reporter track an independent measurement, or track growth?" Nothing in it names yeast. |
| `generator/plate.py`, `calibrate.py`, `export.py` | **S** | `DEFAULT_PANEL`, `DEFAULT_DOSES`, `PlateConditions.gdcw_per_od = 0.42`, and imports of `STRESSOR_FOR_CONSTRUCT` / `NEWPROTOCOL_LAYOUT`. |
| `viz/figures.py` | **S** | The five figures of this study. |
| `analysis/stress_model.py` | **H** | Trained artefact mapping readings to *named yeast modules*. |
| `scripts/*` | **S** | Eleven of twelve are wired to this lab's paths, constructs and agents. |

**Hardcoding census.** The four construct names appear on 40 lines across 15 files; `DTT` /
`H2O2` on 47 lines across 18. Nothing enforces agreement between the **eight** places one new
real construct must be registered (`ARCHITECTURE.md` §6.1).

---

## 3. The proposed seam

Three profiles and one pipeline function. Two of the three profiles **already exist in
fragments** — this is mostly regrouping, not new machinery. Every field follows the house
rule: *unmeasured is `None`, and the consumer refuses rather than defaults.* That rule is
already implemented in `ODCalibration.gdcw_per_od` and `ReporterOptics.inner_filter_coeff`;
the profiles extend it rather than inventing a convention.

### 3.1 `HostProfile` — genuinely new

New module src/ystwin/profiles/host.py (proposed).

```python
@dataclass(frozen=True)
class HostProfile:
    name: str                                  # "S. cerevisiae BY4741"
    taxon_id: int                              # 559292 -- an identifier, not free text
    mu_max: float                              # 1/h unstressed ceiling; a sanity bound on any recovered mu
    doubling_time_bounds_h: tuple[float, float]
    growth_temperature_c: float
    biomass_proxy: Literal["od600", "od700", "cell_count", "confluence"]
    gdcw_per_biomass_unit: float | None        # None => to_dcw() refuses, as it already does
    autofluorescence_per_biomass: float | None # None => refuse to subtract; measure on an isogenic reporter-free strain
    detection_floor: float                     # replaces growth._DEFAULT_DETECTION_OD
    cytosolic_ph_resting: float
    cytosolic_ph_range: tuple[float, float]
    regulons: Mapping[str, Module] | None       # replaces stress_panel.MODULES
    cascade_order: tuple[str, ...] | None       # replaces stress_panel._CASCADE_ORDER; validated against regulons
    gsmm: GsmmProfile | None                    # sbml path + exchange ids + measured phenotype
    source: str
```

Replaces, at these locations:

| Constant | Location |
| --- | --- |
| `_MU_MAX = 0.40` | `generator/stress_panel.py:267` |
| `_DEFAULT_DETECTION_OD = 0.02` | `growth.py:79` |
| `gdcw_per_od: float = 0.42` | `generator/plate.py:84` (`PlateConditions`) |
| `_REFERENCE_TEMPERATURE = 30.0`, `_HEAT_PER_DEGREE = 0.06`, `_HYPOXIA_THRESHOLD = 0.21`, `_CARBON`, `_PHASE` | `generator/context.py:50–52, 22, 42` |
| `STANDARD_TEMPERATURE_K = 303.15`, `CYTOSOLIC_PH = 7.5` | `bridge/thermodynamic.py:90, 104` |
| `_TEMPERATURE_K = 303.15`, `_REFERENCE_PH = 7.0` | `generator/redox.py:27–28` |
| `CYTOSOLIC_PH = (6.4, 7.4)` | `generator/stress_panel.py:39` |
| `MODULES`, `_CASCADE_ORDER` | `generator/stress_panel.py:184, 270` |
| `GLUCOSE_EXCHANGE` … `BIOMASS_REACTION`, `REFERENCE_AEROBIC_BATCH` | `fba/physiology.py:22–26, 43` (and two duplicate copies in `fba/fva.py`, `bridge/physiology_bridge.py`) |
| the "physiological maximum for yeast" sanity bound | `diagnostics/dilution_confound.py:233` (docstring) |

Attachment point: **`CultureContext.strain: str = "BY4741"` (`generator/context.py:63`) is
declared and read by nothing.** Change it to `host: HostProfile` and no behaviour moves.

### 3.2 `ReporterProfile` — composes two existing dataclasses

New module src/ystwin/profiles/reporter.py (proposed). It must **compose** `reporter.ReporterKinetics`
rather than duplicate it — `ReporterKinetics(k_deg, k_mat, k_deg_immature)` is already the
correct kinetic half and is used by `simulate_reporter` and `promoter_activity`.

```python
@dataclass(frozen=True)
class ReporterProfile:
    name: str                              # "mCitrine"
    kinetics: ReporterKinetics             # REUSE: k_deg, k_mat, k_deg_immature -- all 1/h
    kind: Kind                             # REUSE stress_panel.Kind: TRANSCRIPTIONAL | RATIOMETRIC
    excitation_nm: float
    emission_nm: float
    channel_aliases: tuple[str, ...]       # replaces plate/synergy._SHEET_CHANNEL_ALIASES
    ratiometric_partner: str | None        # the second channel; growth cancels exactly in the ratio
    validated_ph: tuple[float, float] | None
    demonstrated_in: frozenset[int]        # taxon ids -- replaces demonstrated_in_yeast: bool
    basal: float
    source: str
```

Replaces: `stress_panel.Reporter.{spectral_slot, demonstrated_in_yeast, validated_ph, basal}`
(`generator/stress_panel.py:165–168`), `plate/synergy._SHEET_CHANNEL_ALIASES:24`, and the
strain half of `observation.ReporterOptics`.

**Two fields have real values available today that the code does not use.**

- `mCitrine` matures with a half-time of **10.4 min** in budding yeast — Guerra et al., *ACS
  Synth Biol* 11(3):1129–1141 (2022), PMID 35180343, doi:10.1021/acssynbio.1c00387, which
  characterises 12 FPs in vivo. That is `k_mat ≈ 4.0 /h` against `k_mat = None`
  (instantaneous) everywhere in the repo. At μ ≈ 0.3 /h the approximation costs ~7 %, which
  is why nothing has broken. The same paper gives mKate2 at 129 min (`k_mat ≈ 0.32 /h`), i.e.
  *equal to μ* — for that reporter "instantaneous" is wrong by roughly a factor of two.
- Excitation/emission and spectra are machine-readable from FPbase — Lambert, *Nat Methods*
  16:277–278 (2019), PMID 30886412, doi:10.1038/s41592-019-0352-8 — so
  `excitation_nm`/`emission_nm` need not be hand-transcribed.

`Reporter.demonstrated_in_yeast: bool → demonstrated_in: frozenset[int]` is the single
highest-leverage one-line generalisation in the tree: it is the field
`analysis/sensor_selection` filters on to compute `RECOMMENDED_BUILD` at import time.

### 3.3 `InstrumentProfile` — mostly already there, under another name

`gates/g1_optical.OpticalQualityGate` and `calib/od.ODCalibration` are already the instrument
profile for the optical channel, and `ODCalibration.optical_gate()` is already the bridge
between them. What is missing is the *reader dialect* and a refusal.

```python
@dataclass(frozen=True)
class PlateFormat:
    rows: str          # "ABCDEFGH" or "ABCDEFGHIJKLMNOP"
    columns: int       # 12 or 24
    def well_pattern(self) -> re.Pattern: ...   # replaces plate/synergy.WELL_RE

@dataclass(frozen=True)
class InstrumentProfile:
    name: str                                  # "BioTek Synergy H1 / Gen5"
    plate_format: PlateFormat
    read: Callable[[Path], AlignedRun]         # registered, not imported
    time_unit: Literal["excel_serial_day", "hh_mm_ss", "minutes", "hours"]
    blank_convention: Literal["raw", "subtracted"]   # DECLARED, with the sniffer kept as a cross-check
    optical: ODCalibration | None              # REUSE; None => od_linear_max has no value
    detector_max: float | None
    min_reporter_above_background: float | None
    fluorescence_units: Literal["rfu", "mefl"] # RFU is reader-local; MEFL is comparable
    source: str
```

Replaces: `plate/synergy.WELL_RE:20`, `_SHEET_CHANNEL_ALIASES:24`, `_parse_channel_header:54`,
`_to_hours:76`, `KineticBlock.blank_subtracted` (as a declaration rather than a sniff), and
every default on `OpticalQualityGate:38–42`.

**The behavioural change that matters:** `OpticalQualityGate()` should stop having an
`od_linear_max` default at all. Today the placeholder `1.0` is used on every real plate, and
the docstring says it is "the single most important number to measure for the platform".
Removing the default converts a silent assumption into a `TypeError` that names the
experiment — which is precisely the house pattern (`correct_inner_filter`,
`latent_constraints`, `ODCalibration.to_dcw` all already do this).

### 3.4 `recover_activity` — the pipeline function that does not exist

New module src/ystwin/pipeline.py (proposed), ~120 lines, **no new science**.

```python
def recover_activity(
    aligned: pd.DataFrame,          # time_h x MultiIndex(channel, well) -- SynergyRun.aligned's shape
    biomass_channel: str,
    reporter_channel: str,
    host: HostProfile,
    reporter: ReporterProfile,
    instrument: InstrumentProfile,
    blank_wells: Sequence[str] | None = None,   # None => detect_blank_wells(od, rfu=...) two-channel test
) -> ActivityTable:                 # per well: mu, naive RFU/OD, activity, gate verdict, interval
```

It calls, in order: `detect_blank_wells(od, rfu=…)` (the two-channel form, which
`plate_dilution_report` and `run_gates.py` currently skip), `assess_plate` with the
instrument's gate, `growth.specific_growth_rate`, `reporter.naive_specific_fluorescence`,
`reporter.promoter_activity`, and — this is the part that is missing from the current script —
`analysis/uncertainty.fold_change` for the interval.

This is the load-bearing change. Today the only path to the validated core runs
`paths.biosensor_plates()` → `read_dose_response` → `recover_layout` → `NEWPROTOCOL_LAYOUT`
fallback → four hardcoded construct names, in a script. Nothing about a second organism or a
second instrument is hard once this function exists; everything is hard while it does not.

### 3.5 Does the seam already exist? Module by module

| Candidate | Verdict |
| --- | --- |
| `generator/literature.py` | **No.** `_CONSTRUCTS` is the four constructs' *simulator* parameters; `LIBRARY`/`EXPECTATIONS` are yeast physiology. But its `Value(value, units, source, verified, note)` dataclass is exactly the right shape for a profile field and should be lifted, not re-invented. |
| `generator/context.py` | **Closest.** `CultureContext` already separates "what the culture is in" from "what the stressor does", and its `strain` field is a dead hook waiting for a `HostProfile`. `_CARBON`/`_PHASE` are already per-host lookup tables. |
| `reporter.ReporterKinetics` | **Yes, half of one.** The kinetic half of `ReporterProfile` already exists and is correct. |
| `observation.ReporterOptics` | **Yes, but split wrong.** It mixes reader fields (`gain`, `background`, `detector_max`) with a strain field (`autofluorescence`) and a product field (`inner_filter_coeff`). Three profiles in one dataclass. |
| `gates/g1_optical.OpticalQualityGate` | **Yes.** Frozen thresholds, carried on the result, auditable without the call site. This is the pattern the other profiles should copy verbatim. |
| `calib/od.ODCalibration` | **Yes**, and `optical_gate()` is already the bridge from calibration to gate. |
| `paths.py` | **Yes, for datasets.** The env-var → in-repo → sibling contract, with `None` for absent and `require()` naming the variable, is exactly the right shape for "here is a second dataset". |

**Three of four profiles are regroupings of code that already exists. Only `HostProfile` is
new.** That is the strongest argument that this generalisation is affordable.

---

## 4. Demonstration targets, ranked

Ranked by *evidence available ÷ effort*, not by scientific interest. For "generalizable" to be
a claim rather than an aspiration, exactly one of these has to actually run.

### Rank 0 (free, do it first): the ATP-sensor plates already in this repo

`generator/synergy.py` is a **second Synergy reader with different assumptions**, used only by
`scripts/analyse_atp_sensor.py`. Running `recover_activity` over both readers is a same-day
test of the instrument seam using data already on disk. It licenses nothing scientifically and
it de-risks everything below. Not on the brief's list; it should be.

### Rank 1: *E. coli* with the Alon/Zaslaver promoter library

**By far the best ratio, and the only target reachable in weeks.**

- **Zaslaver et al., *Nat Methods* 3:623–628 (2006).** PMID 16862137, doi:10.1038/nmeth895.
  ~2,000 promoter–GFP transcriptional fusions on low-copy plasmids covering most K-12
  promoters. Supplementary data hosted by the Alon lab
  (`weizmann.ac.il/mcb/alon/EcoliTranscLibrary`).
- **Zaslaver et al., *PLoS Comput Biol* 5(10):e1000545 (2009).** PMID 19851443,
  doi:10.1371/journal.pcbi.1000545. Genome-scale promoter activity at high temporal
  resolution **as cells pass from exponential into stationary phase in different media**.
  This is the growth-dilution confound run as a designed experiment, in exactly this tool's
  data shape: 96-well, OD + GFP, dense time grid.
- **Kaplan et al., *Mol Cell* 29(6):786–792 (2008).** PMID 18374652, PMC2366073,
  doi:10.1016/j.molcel.2008.01.021. Nineteen sugar genes mapped as two-dimensional input
  functions — dose grids far denser than the four NewProtocol plates.

**What must be re-measured:** ideally nothing. **Verify before planning around it:** whether
the archived files carry *per-well OD and GFP* or only already-derived promoter activity. I
could not confirm the file contents (the Weizmann page did not respond). If only derived
activity is posted, μ cannot be re-estimated and the demonstration is much weaker. The
fallback is one plate: the library is a distributed physical resource.

**What it licenses — and this is the point.** The Alon convention *defines* promoter activity
as `(dGFP/dt)/OD`, i.e. the incumbent, widely used estimator omits the `(μ + k_deg)·R`
dilution term entirely. So this is not "our method also runs on E. coli"; it is a head-to-head
against the field-standard estimator, on the field's own data, scored with
`analysis/nulls.skill_score`, in the regime (exponential → stationary) where the two must
diverge. In the repo's own tier vocabulary that is a Tier 2 claim: agreement — or disagreement
— with a measurement the model never saw.

### Rank 2: iGEM Interlab — for units, not biology

Ranked second for a different reason from the others: it fixes the thing that stops any output
of this tool from being comparable between labs.

- **Beal et al., *PLoS ONE* 13(6):e0199432 (2018).** PMID 29928012,
  doi:10.1371/journal.pone.0199432. Bacterial fluorescence in **MEFL** (Molecules of
  Equivalent FLuorescein) via fluorescein and SpheroTech rainbow calibration particles.
- **Beal et al., *Commun Biol* 3:512 (2020).** PMID 32943734, doi:10.1038/s42003-020-01127-5
  (correction PMID 33110148). OD → cell count across **244 laboratories**; serial dilution of
  silica microspheres beats both LUDOX and CFU.
- **Beal et al., *PLoS ONE* 16(6):e0252263 (2021).** PMID 34097703. Comparative analysis
  across three interlab studies.

**What must be re-measured:** two calibration plates, run here — one fluorescein dilution
series, one microsphere series, per reader and gain. Cheapest item in this entire note.

**What it licenses:** absolute units. Every number this repo reports is currently
`RFU · OD⁻¹ · h⁻¹` on one PMT at one gain — a quantity no other lab can reproduce or compare.
A tool that reports in reader-local RFU is not general in the sense that matters, however
general its equations are. It also supplies an *inter-laboratory* variance figure to sit next
to `OBSERVED_ACTIVITY_CV = 0.146`, which is between-plate within one lab.

**Caveat:** I found no single consolidated public archive of per-team raw plate exports; the
analysed data accompanies the papers. Confirm before planning a re-analysis around it.

### Rank 3: other yeasts

**Correct the brief first: *Pichia pastoris* and *Komagataella phaffii* are the same
organism.** Kurtzman, *J Ind Microbiol Biotechnol* 36(11):1435–1438 (2009), PMID 19760441,
doi:10.1007/s10295-009-0638-4, showed that the biotechnological strains called *P. pastoris*
(GS115, X-33, CBS 7435) are *K. phaffii*. That is one target on the list, not two.

- **Vogl et al., *ACS Synth Biol* 5(2):172–186 (2016).** PMID 26592304. A functionally
  verified promoter toolbox with eGFP reporters — but the published readouts are largely
  endpoint and flow-cytometric, not the dense kinetic OD + RFU traces this tool consumes.

**What must be re-measured:** essentially the whole time series. That is the honest answer,
and it is why the most *relevant* target for bioproduction ranks third.

**What it licenses:** "runs in the production host people actually use" — persuasive to a
bioprocess audience, weak as evidence of generality, because two Saccharomycetes are close
relatives. Its real value is as a **test of the seam**: nothing structural changes between
*S. cerevisiae* and *K. phaffii* — different `mu_max`, `gdcw_per_od`, carbon axes
(methanol/glycerol), growth temperature. If adding it is a twenty-line `HostProfile` and no
edits elsewhere, the seam is real. If it is not, the seam is wrong.

### Rank 4: mammalian reporter assays — highest value, highest cost, put last deliberately

This is the target that tests the dilution correction hardest, and it is a quarter's work.

- **Eden et al., *Science* 331(6018):764–768 (2011).** PMID 21233346,
  doi:10.1126/science.1199784. Bleach-chase on 100 proteins in living human cells,
  **separating degradation from dilution explicitly**. Half-lives 45 min to 22.5 h against a
  ~24 h cell cycle — i.e. `k_deg` and `μ` are *the same size*. That is exactly the regime
  where this repo's universal `k_deg = 0` fails and where an estimator carrying both terms
  earns its keep.
- **Dénervaud et al., *PNAS* (2013), doi:10.1073/pnas.1308265110** is worth noting as the
  yeast counter-example: a microchemostat array over the yeast-GFP library under stress. In a
  chemostat μ is *set by the dilution rate*, so the confound is removed by design — a useful
  negative control for the claim, and a data shape (single-cell microscopy, no OD) this tool
  cannot currently read.

**What must be re-measured:** a full in-house time course, or a collaboration. There is no OD:
biomass comes from confluence, nuclear counts, or a constitutive reference channel. That is a
`HostProfile.biomass_proxy` field, not a rewrite — `growth.specific_growth_rate` only needs a
strictly-positive series; only its argument name and error message say "optical density".

**Realistic cost: 3–6 months.** Do not put it on a roadmap without a collaborator.

**What it licenses:** the strongest available claim — the dilution correction is not a yeast
trick — and the only test that genuinely stresses the `k_deg` / `μ` identifiability.

---

## 5. What breaks under generalisation

Adversarially, and several of these are live bugs rather than future risks.

### 5.1 The dilution correction does not assume exponential growth. Everything around it does.

The inversion `k = dR/dt + (μ + k_deg)·R` estimates μ pointwise and is agnostic. What is not:

- **`growth.specific_growth_rate(window_h=1.5, polyorder=2)`** assumes `ln OD` is locally
  quadratic over 1.5 h. For a mammalian line doubling in 24 h, 1.5 h is ~6 % of a doubling and
  the log-slope is under the reader noise.
- **In stationary phase with a stable FP the correction degenerates.** As μ → 0 and
  `k_deg = 0`, `k → dR/dt`: every bit of information now comes from a numerical derivative of
  a noisy trace. Worse, `growth_rate_uncertainty` returns an **absolute** SE — the repo
  measured 0.0117 /h and its own docstring notes this is "3 % of a healthy growth rate and
  12 % of one slowed to 0.10". At μ = 0.03 /h (mammalian) it is ~40 %; in stationary phase it
  is unbounded.
- **`dilution_confound_report` refuses outright** when the robust OD fold is below
  `_MIN_GROWTH_FOLD = 1.15`. Correct behaviour, and it means the tool currently cannot answer
  the question for stationary-phase or slow-growing systems at all — rather than answering it
  with a wide interval.
- **D2's headline statistic silently becomes undefined.** The QSS test is
  `loss > 1/max_relaxation_h` with `max_relaxation_h = duration/2`. At μ = 0.03 /h over a 48 h
  run, `loss × max_relaxation = 0.72 < 1`, so no timepoint is quasi-steady, `usable.sum() <
  _MIN_REGRESSION_POINTS`, and `dilution_explained_r2` returns NaN. **The repository's
  most-cited real-data finding — "median R² of log(RFU/OD) on growth = 0.82" — does not exist
  for slow-growing systems.** That is a correct refusal and it should be stated as a scope
  limit, not discovered by a user.

### 5.2 `k_deg = 0` is assumed universally, and one place gets it wrong rather than refusing

`ReporterKinetics.k_deg = 0.0`, `CultureParameters.k_deg = 0.0`,
`run_gates.py` passes `ReporterKinetics(k_deg=0)` explicitly,
`analysis/power.parameters_from_literature` defaults to it — and
**`analysis/power._one_trial:195` inverts dilution as `recovered = signal * growth`, i.e.
`k = R·μ`, which is valid only when `k_deg == 0`.** It is right today by coincidence of the
default and silently wrong the first time anyone passes a degron reporter.

Degron reporters are not exotic — ssrA-tagged GFP variants (Andersen et al., *Appl Environ
Microbiol* 64(6):2240–2246, 1998, doi:10.1128/aem.64.6.2240-2246.1998) have half-lives of tens
of minutes to a couple of hours, i.e. `k_deg` of roughly 0.4–1.0 /h, comparable to or larger
than μ. They are common in the iGEM registry. Aim this tool at that registry and the power
analysis is wrong on contact.

Compounding it: `calib/kdeg.py` exists precisely to measure `k_deg`, has 218 lines and a test
suite, and **has never been run** because the chase experiment has not been done. So the
profile field would have no measurement to fill it in this lab either.

### 5.3 Instantaneous maturation is a yeast-growth-rate artefact, not a reporter property

`k_mat = None` everywhere. For mCitrine in yeast that is defensible: 10.4 min half-time
(Guerra 2022) gives `k_mat ≈ 4.0 /h` against μ ≈ 0.3 /h. But:

- In *E. coli* μ reaches ~2 /h, so even a 10-minute maturation is only 2× μ. Balleza et al.,
  *Nat Methods* 15(1):47–51 (2018), PMID 29320486, doi:10.1038/nmeth.4509, measured maturation
  for **50 FPs at two temperatures in E. coli** precisely because it limits the measurement of
  rapid dynamics.
- Slow-maturing red FPs (mKate2, 129 min → `k_mat ≈ 0.32 /h`) put maturation *at* μ even in
  yeast.

`promoter_activity` already implements the two-stage inversion correctly when `k_mat` is
supplied; it differentiates twice, which roughly doubles the noise amplification. So the fix
is not code — it is that the second derivative's error budget has never been characterised the
way `default_activity_window_h` characterised the first.

### 5.4 The "seven regulons" understatement

The Ledermann (1937) ceiling `m ≤ (2p + 1 − √(8p + 1))/2` is arithmetic and travels
unchanged — at p = 4 channels, one fitted factor, in any organism. What does not travel is the
*content*: 24 `MODULES`, `_CASCADE_ORDER`, the `driven_by` weights and their PMIDs,
`reporter_loadings`, 25 `STRESSORS` and their `targets`.

Two honest points the brief's framing misses.

**The count is not stable inside yeast.** The panel went from 7 modules to 24 during this
project; two module docstrings and the README still say "seven". So "the seven regulons is
yeast-specific" understates it — the number is a modelling choice, not a fact about
*S. cerevisiae*.

**The analogue elsewhere exists and is better curated.** *E. coli* has seven sigma factors
(σ70/RpoD, σ38/RpoS, σ32/RpoH, σ24/RpoE, σ54/RpoN, σ28/FliA, σ19/FecI) plus global regulators
(CRP, ArcA, Fnr, OxyR, SoxRS, Lrp), and RegulonDB curates the topology that
`stress_panel.MODULES` had to be hand-transcribed for. In mammalian cells the canonical stress
programs are also roughly six or seven — HSF1, the three UPR arms (PERK/ATF4, IRE1/XBP1s,
ATF6), NRF2, p53, HIF1α, NF-κB — and they are also non-orthogonal (ATF4 is shared across the
integrated stress response). **The Ledermann squeeze has the same shape in all three hosts,
which is a point in favour of the finding's generality even though the panel encoding it is
not.**

### 5.5 The gates encode one reader's dynamic range, and one of them is not portable at all

- `od_linear_max = 1.0` is calibrated against nothing. `calib/od_linearity.py` turns a dilution
  series into it; the series (protocol P4) has never been run.
- `min_reporter_above_background = 50.0` is raw RFU at one PMT gain. It is not portable between
  **two gains on the same instrument** — and `plate/synergy.SynergyRun.channel` already knows
  this, refusing a bare fluorophore name when a run holds two reads "because they differ in
  gain". The gate contradicts the reader's own stated caution.
- `min_od_fold = 1.15` and `max_decline_fraction = 0.15` encode "a batch culture that grows
  several fold over the run" — false for a chemostat, false for a mammalian plate read over one
  doubling.
- `min_od_above_blank = 0.02` and `growth._DEFAULT_DETECTION_OD = 0.02` are 96-well figures
  (~0.3–0.5 cm path); a 384-well plate and a cuvette are wrong in opposite directions.

The fix is not new numbers. It is that `OpticalQualityGate` should have **no defaults**, so a
missing calibration is a refusal that names the experiment.

### 5.6 `blank_subtracted` is a heuristic where an authoritative signal exists

`KineticBlock.blank_subtracted = nanmin(rows) < 0` infers "already blanked" from the presence
of negative values. That is a property of *this exporter's habit*, not of the data. An
instrument that clips at zero, or a blank-subtracted export that happens to stay positive, is
reported as raw and `raw_channel` hands it back for a second subtraction. `ARCHITECTURE.md`
already flags this `[uncertain]`; under generalisation it stops being uncertain and becomes
certain. The instrument profile should **declare** the convention, with the sniffer retained as
a cross-check that warns on disagreement.

### 5.7 The time axis is the most dangerous single line

`_to_hours` treats a bare float as an **Excel serial day fraction × 24** and `"30:00"` as
MM:SS. Both are right for Gen5. An exporter writing elapsed *minutes* as a bare float produces
traces 24× too long, which yields a growth rate 24× too small, which yields an activity that
is wrong by 24× and looks entirely plausible. This is the failure mode most likely to produce
a confidently wrong published number.

### 5.8 The headline finding may not survive contact with a strong induction

"Most of the sensor response was dilution" rests on n = 2 plates and is INCONCLUSIVE by the
repo's own `fold_change`. It is a statement about **weak inductions near the noise floor**: a
1.33× naive fold correcting to 0.96×. In the Zaslaver diauxic-shift data, promoter activity
moves by orders of magnitude, and the dilution correction will be a small proportional
adjustment on top of a large real signal. Generalising may therefore make the correction look
*less* important, not more. Say that first, before someone else does.

### 5.9 Absolute units

Every reported number is `RFU · OD⁻¹ · h⁻¹` on one reader at one gain. Two labs cannot compare
outputs; this lab cannot compare its own outputs across a gain change. This is arguably a
larger obstacle to being a general tool than any code change, and it is the cheapest to fix
(§4, rank 2).

---

## 6. Naming and framing

### 6.1 The one sentence

> **Given a fluorescent reporter and a biomass time series from one well, recover promoter
> activity with growth dilution removed — with a calibrated interval, and an explicit refusal
> where the optics do not support the question.**

No host, no product, no instrument. It is the sentence `docs/CONTRACT.md` §"One sentence"
already writes, minus "OD" and "mCitrine". Note that `pyproject.toml` has *already* been moved
to a variant of this, so the divergence between the package **name** and the package
**description** is now live rather than hypothetical.

### 6.2 The argument for renaming

- `ystwin` expands to "yeast stress twin". Both halves are now wrong: the general claim is not
  about yeast, and `docs/CONTRACT.md` places the project at DNV level 3 at best, "honestly at
  level 2, diagnostic" — it is an estimator, not a twin.
- The name is the first thing a reader meets and it contradicts the leading sentence. A name
  that misdescribes scope is exactly the class of defect this repository spends
  `CLAIM_BOUNDARY.md` and `PARKED.md` policing.
- It gets more expensive later, not less: the moment anyone outside this team imports it, a
  rename needs a deprecation shim.

### 6.3 The argument against renaming now

- The mechanical surface is **444 import lines across 105 files**, plus **nine `YSTWIN_*`
  environment variables** that are documented in `REPRODUCING.md` and already set in people's
  shells. The imports are a `sed`; the environment variables break other people's setups.
- `scripts/audit_claims.py` checks README↔code agreement and `audit_reproducibility.py` checks
  that a stranger can rerun it. A rename means re-running and re-verifying all of that, plus
  3,241 tests, plus every path in `outputs/` provenance.
- The churn buys **zero scientific content**. Nobody is misled about scope by an `import`
  statement; they are misled by the README's first line, which is free to fix.
- The name is not yet load-bearing on any external claim, because nothing external consumes it.

### 6.4 Recommendation

**Rename the claim now. Defer the package rename until a second host actually runs, then do it
in one commit with a shim.**

Concretely:

1. Today: README first line, `pyproject.toml` description (done), `__init__.py` docstring
   (done), and `docs/CONTRACT.md`'s one sentence all say the §6.1 sentence. State the yeast
   panel, the Yeast9 validation, the Synergy reader and the carotenoid layer as **instances**,
   which is what `__init__.py` now already does.
2. Move the carotenoid layer from "parked" to "worked example". Under generalisation it is more
   parked, not less — a product-specific inner-filter correction is by definition not general —
   and calling it an example rather than a deferral is both truer and less alarming.
3. Rename only after the E. coli demonstration runs, so the new name is earned by a second
   host rather than announced ahead of one. Ship `ystwin` as a shim that re-exports and emits
   `DeprecationWarning`, and keep the `YSTWIN_*` variables as accepted aliases for one release.
4. If renamed, the name should describe the **operation**, not the organism. Lead candidate
   **`undilute`** — it names what the tool does and nothing else. Conservative alternative
   **`promact`** (promoter activity). Check availability before committing to either; I have
   not.

---

## 7. A staged plan

### Week one — about three days, no new bench data, no new science

1. **`src/ystwin/pipeline.py::recover_activity`** (§3.4). Lift the middle of
   `run_sensor_characterisation.py` into a library function on a tidy frame. The script becomes
   a thin caller. *This is the item that makes everything below cheap; do not reorder it.*
2. **`HostProfile` + a `SACCHAROMYCES` instance** populated from constants already in the
   tree, each carrying its existing `source` string. Repoint `CultureContext.strain: str` at it
   — the field is dead, so nothing moves.
3. **`ReporterProfile` composing `ReporterKinetics`.** Set `mCitrine.k_mat = 4.0 /h` from
   Guerra 2022 and run the suite. **If nothing moves, that is a result** — it says instantaneous
   maturation was safe *here* and names the FPs for which it would not be.
4. **`Reporter.demonstrated_in_yeast: bool → demonstrated_in: frozenset[int]`.** One line; it
   is the field that decides `RECOMMENDED_BUILD`.
5. **Run the whole real-data path through `generator/synergy.py`'s ATP-sensor plates** as a
   second instrument dialect. Free generality test on data already on disk.
6. **Fix `analysis/power._one_trial`** to carry `k_deg` (§5.2). It is wrong today for any
   reporter that is not a stable FP.

### Month one — about four weeks plus two calibration plates

7. **`InstrumentProfile` + a reader registry.** Parameterise `WELL_RE` by `PlateFormat`; make
   `blank_convention` declared with the sniffer kept as a warning cross-check; make the time
   convention explicit.
8. **Delete `OpticalQualityGate`'s defaults** and run protocol P4 (one dilution series, one
   plate) so `od_linear_max` becomes a measurement.
9. **Run one fluorescein/MEFL plate and one silica-microsphere plate** (Beal 2018 / 2020) so
   outputs stop being in reader-local RFU.
10. **Wire `analysis/uncertainty.fold_change` into the pipeline.** `ARCHITECTURE.md` §5.5
    already names this the highest-value change on the real-data path, and the demonstration
    below needs intervals to be worth anything.
11. **The E. coli demonstration.** Acquire the Zaslaver/Alon traces; confirm they carry per-well
    OD and GFP; run `recover_activity` over the exponential→stationary and diauxic-shift series;
    score against the incumbent `(dGFP/dt)/OD` estimator with `analysis/nulls.skill_score`
    against a baseline ladder. Write it as its own note in `docs/research/`.
12. **Add a `KOMAGATAELLA` `HostProfile`** even without data, as a seam test: if it is twenty
    lines and touches nothing else, the design in §3 is right.

### Explicitly out of scope — say no to these

- **Mammalian**, without a collaborator and 3–6 months. It is the best test and it is not
  fundable from here. Put it in "future work" and mean it.
- **Generalising the stress panel.** `MODULES`, `_CASCADE_ORDER`, `driven_by`,
  `reporter_loadings` and `STRESSORS` are a yeast literature transcription with PMIDs and one
  correction already recorded. Porting to *E. coli* means re-transcribing from RegulonDB — a
  project of its own — and **it is not needed for the core claim.** State plainly: *the general
  tool is the estimator; the panel is one host's instance.*
- **Generalising `fba/` and `bridge/`.** Yeast9 reaction ids, van Hoek phenotype, a pinned
  `pytfa`. Parked; leave parked.
- **Renaming the package** before a second host runs (§6.4).
- **Unparking the carotenoid layer.** Reclassify it as a worked example, do not revive it.

### Honest cost

Items 1–6 are **about three days** and cost nothing but review time; they make the claim in
§6.1 defensible on the code as it stands. Items 7–12 are **four to six weeks plus two or three
calibration plates**. Doing generalisation *fully* — mammalian, a second stress panel, the
prescriptive reward layer `CONTRACT.md` names as missing — is **six months and wet-lab capacity
this team does not have.**

The useful subset is the first week, plus item 11. Item 11 is the only thing on this list that
converts "generalizable" from an assertion into a measurement, and it is the one that should be
protected if anything gets cut.

---

## References

All verified against PubMed or the publisher at compile time.

1. Zaslaver A, Bren A, Ronen M, Itzkovitz S, Kikoin I, Shavit S, Liebermeister W, Surette MG,
   Alon U. A comprehensive library of fluorescent transcriptional reporters for *Escherichia
   coli*. *Nat Methods* 3(8):623–628 (2006). PMID 16862137. doi:10.1038/nmeth895
2. Zaslaver A, Bren A, Ronen M, et al. Invariant distribution of promoter activities in
   *Escherichia coli*. *PLoS Comput Biol* 5(10):e1000545 (2009). PMID 19851443.
   doi:10.1371/journal.pcbi.1000545
3. Kaplan S, Bren A, Zaslaver A, Dekel E, Alon U. Diverse two-dimensional input functions
   control bacterial sugar genes. *Mol Cell* 29(6):786–792 (2008). PMID 18374652, PMC2366073.
   doi:10.1016/j.molcel.2008.01.021
4. Beal J, Haddock-Angelli T, Baldwin G, Gershater M, Dwijayanti A, Storch M, de Mora K,
   Lizarazo M, Rettberg R; iGEM Interlab Study Contributors. Quantification of bacterial
   fluorescence using independent calibrants. *PLoS ONE* 13(6):e0199432 (2018). PMID 29928012.
   doi:10.1371/journal.pone.0199432
5. Beal J, Farny NG, Haddock-Angelli T, et al.; iGEM Interlab Study Contributors. Robust
   estimation of bacterial cell count from optical density. *Commun Biol* 3(1):512 (2020).
   PMID 32943734. doi:10.1038/s42003-020-01127-5 (author correction PMID 33110148)
6. Beal J, et al. Comparative analysis of three studies measuring fluorescence from engineered
   bacterial genetic constructs. *PLoS ONE* 16(6):e0252263 (2021). PMID 34097703
7. Balleza E, Kim JM, Cluzel P. Systematic characterization of maturation time of fluorescent
   proteins in living cells. *Nat Methods* 15(1):47–51 (2018). PMID 29320486.
   doi:10.1038/nmeth.4509
8. Guerra P, Vuillemenot L-A, Rae B, Ladyhina V, Milias-Argeitis A. Systematic in vivo
   characterization of fluorescent protein maturation in budding yeast. *ACS Synth Biol*
   11(3):1129–1141 (2022). PMID 35180343. doi:10.1021/acssynbio.1c00387
9. Eden E, Geva-Zatorsky N, Issaeva I, Cohen A, Dekel E, Danon T, Cohen L, Mayo A, Alon U.
   Proteome half-life dynamics in living human cells. *Science* 331(6018):764–768 (2011).
   PMID 21233346. doi:10.1126/science.1199784
10. Kurtzman CP. Biotechnological strains of *Komagataella* (*Pichia*) *pastoris* are
    *Komagataella phaffii* as determined from multigene sequence analysis. *J Ind Microbiol
    Biotechnol* 36(11):1435–1438 (2009). PMID 19760441. doi:10.1007/s10295-009-0638-4
11. Vogl T, Sturmberger L, Kickenweiz T, et al. A toolbox of diverse promoters related to
    methanol utilization: functionally verified parts for heterologous pathway expression in
    *Pichia pastoris*. *ACS Synth Biol* 5(2):172–186 (2016). PMID 26592304
12. Lambert TJ. FPbase: a community-editable fluorescent protein database. *Nat Methods*
    16:277–278 (2019). PMID 30886412. doi:10.1038/s41592-019-0352-8
13. Andersen JB, Sternberg C, Poulsen LK, Bjørn SP, Givskov M, Molin S. New unstable variants
    of green fluorescent protein for studies of transient gene expression in bacteria. *Appl
    Environ Microbiol* 64(6):2240–2246 (1998). doi:10.1128/aem.64.6.2240-2246.1998
    *(PMID not independently confirmed at compile time; DOI verified.)*
14. Dénervaud N, et al. A chemostat array enables the spatio-temporal analysis of the yeast
    proteome. *PNAS* (2013). doi:10.1073/pnas.1308265110 *(volume/pages not independently
    confirmed at compile time.)*

Ledermann's factor-analysis bound and the Kapteyn et al. digital-twin vocabulary are cited
where they already live, in [`IDENTIFIABILITY.md`](IDENTIFIABILITY.md) and
[`../CONTRACT.md`](../CONTRACT.md), and are not restated here.

## Provenance

Written 2026-08-26 against `src/ystwin` (53 code modules, ~13,110 LOC), `docs/CONTRACT.md`,
`docs/ARCHITECTURE.md`, `docs/CLAIM_BOUNDARY.md`, `docs/PARKED.md` and `README.md`. Module
classifications in §2 were made by reading each module, not by directory. Counts (444 import
lines / 105 files; 40 construct-name lines / 15 files; 47 agent-name lines / 18 files) were
measured, not estimated. Nothing in this note was executed against real plate data; it is a
design argument, and it is Tier 0 in this repository's own vocabulary until item 11 of §7 runs.
