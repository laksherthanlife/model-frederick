"""The yeast stress landscape: modules, their crosstalk, stressors and reporters.

Literature-informed topology with declared phenomenological response weights. A cited
interaction does not establish a measured magnitude for every coefficient, and catalogue
membership is not a complete kinetic model. Programmes overlap on shared promoters and
regulatory outputs: Rpn4 integrates Yap1/Hsf1 input, and Hog1 acts through several factors,
including Hot1/Sko1. Signaling onto a factor is not the same as transcribing its gene.
A panel assuming orthogonal sensors would over-count the states it can identify.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import numpy as np

__all__ = [
    "MODULES",
    "Kind",
    "REPORTERS",
    "STRESSORS",
    "Module",
    "Reporter",
    "Stressor",
    "combination_response",
    "demonstrated_in_yeast",
    "growth_rate",
    "healthy_ladder",
    "integration_hours",
    "module_ec50",
    "module_response",
    "ratiometric_reporters",
    "viability",
    "reporter_loadings",
    "transcriptional_reporters",
]


_DEFAULT_BASAL = 0.9
"""Constitutive floor a promoter sits on, in the same units as the loadings.

Set so a fully induced reporter reads about 1.4 times its own undosed well, which is what
UPRE1 and UPRE2 do across the measured DTT ladder. Starting from zero instead overstates
the contrast by roughly an order of magnitude, and the posterior-predictive gate refuses
a generator that does.
"""

CYTOSOLIC_PH = (6.4, 7.4)
"""Where the yeast cytosol actually goes, in pH units.

Resting is about 7.2, and glucose removal drops it to 6.4 within minutes as Pma1 loses the
ATP to pump against it (Dechant 2010, PMID 20581803). That is the same event an ATP channel
exists to report, so a sensor validated only above 7.3 is untested exactly where it matters.
"""

_MU_MAX = 0.40
"""Specific growth rate of an unstressed culture, in 1/h.

Bracketed rather than measured here. BioNumbers 108255 puts a laboratory haploid at ~90 min
in YPD (mu = 0.46) and ~140 min in synthetic medium (mu = 0.30), and 106359 gives 0.37 on
minimal medium at 30 C; 0.40 sits inside that band. Which end is right depends on the plate
medium, which `docs/DATA_INVENTORY.md` does not record, so the band is asserted in
`tests/test_constants_against_bionumbers.py` and the value is left where it was rather than
moved on a guess about which medium was used.
"""

_STABLE_FP_K_DEG = 0.0
"""Loss of mature reporter beyond dilution, in 1/h. Zero for a stable fluorescent protein,
which is what every promoter fusion in the panel uses."""


def integration_hours(mu: float = _MU_MAX, k_deg: float = _STABLE_FP_K_DEG) -> float:
    """Window a promoter fusion averages its promoter's activity over, in hours.

    Derived, not asserted, and that is the whole point: a promoter fusion obeys
    ``dR/dt = k - (mu + k_deg) R``, so its time constant is ``1/(mu + k_deg)`` and nothing
    else. Written down as a literal it silently encodes one host's growth rate -- 3.0 h is
    mu = 0.33 -- and stops being true the moment the culture, the medium or the organism
    changes. Computed, it transfers to any host for the price of that host's mu, which is
    the one number a new lab will have.
    """
    relaxation = float(mu) + float(k_deg)
    if relaxation <= 0.0:
        raise ValueError(f"a non-growing culture never relaxes; got mu={mu}, k_deg={k_deg}")
    return 1.0 / relaxation


_INTEGRATION_HOURS = integration_hours()


# QUEEN-2m's yeast demonstration is Takaine (PMID 30858198), not Yaginuma 2014
# (PMID 25283467) -- Yaginuma is real and correctly cited here, but it is an E. coli paper.
# This entry feeds RECOMMENDED_BUILD and THREE_SENSOR_BUILD, so a wrong identifier here
# reaches a build recommendation.
#
# The dermoscopy identifier this entry used to carry is retired. It now appears nowhere in
# src/ or scripts/, so refresh_citations.py has pruned its row from
# data/citations/pmid_titles.csv; the record of the error is the regression pin in
# tests/test_citations.py and docs/SOURCE_AUDIT.md 3.2.


class Kind(Enum):
    """How a sensor turns cell state into a number.

    The distinction is not cosmetic. A promoter fusion must be transcribed, translated and
    matured, accumulates over hours and is diluted by growth, so its reading is a
    convolution of past activity with the growth rate. A ratiometric biochemical sensor
    reports an equilibrium between two of its own spectral forms, so it settles in seconds,
    reverses, and its ratio cancels concentration -- growth does not enter at all.
    """

    TRANSCRIPTIONAL = "transcriptional"
    RATIOMETRIC = "ratiometric"


@dataclass(frozen=True)
class Module:
    """One transcriptional stress regulon, and the regulons that drive its gene.

    ``driven_by`` is a cascade: an upstream factor transcribing the gene for this module's
    factor, so this regulon genuinely rises. It is deliberately sparse -- most of the
    crosstalk in the yeast stress network is parallel action on shared promoters or
    signalling onto an existing factor, and neither of those belongs here.
    """

    transcription_factor: str
    response_element: str
    source: str
    driven_by: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class Stressor:
    """An agent, the modules it activates, and the dose at which it stops the cell.

    ``lethal_dose`` is where half of transcriptional capacity is gone. It defaults to a
    multiple of the inducing EC50 because an agent that induces at one concentration is
    generally lethal a few fold above it; where plates say otherwise, the measured value
    replaces the default.

    ``target_ec50`` lets one agent reach two modules at two concentrations. Left empty the
    potencies follow a default ordering: a metabolite pool moves at the concentration of
    the chemistry itself, a specific regulon at the agent's own EC50, and the general stress
    response later, once the insult is substantial.
    """

    agent: str
    units: str
    ec50: float
    targets: dict[str, float]
    source: str
    lethal_dose: float = 0.0
    target_ec50: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.lethal_dose:
            object.__setattr__(self, "lethal_dose", _LETHAL_MULTIPLE * self.ec50)


@dataclass(frozen=True)
class Reporter:
    """A promoter-fusion reporter, and every module its promoter's elements read.

    ``also_reads`` is promoter content, not regulation: a promoter carrying a second
    response element reports a module its own regulon has nothing to do with. Keeping it
    apart from the cascade stops one weight standing in for two different mechanisms.

    Its weights are relative to the reporter's own module, which `reporter_loadings` pins at
    1.0. So a weight above 1.0 is not an error: it says the literature measured the second
    element as the dominant one at this promoter, which ENA1 is. A module the reporter's own
    regulon genuinely dominates is the case, not the rule.

    ``basal`` is the constitutive floor the promoter sits on, in the same units as the
    loadings. It is not decoration: the measured induction across the whole DTT ladder is
    1.38 to 1.47 fold, so the floor is more than twice the induced component, and a
    generator starting from zero overstates the contrast by roughly an order of magnitude.

    ``ph_sensitive`` marks a sensor whose reading moves with compartment pH. It is not a
    footnote: one pH unit shifts a glutathione probe by 60 mV against a genuine peroxide
    response of 40 to 50, so such a sensor needs a pH measurement beside it or its readings
    cannot be told apart from an acidification.

    ``demonstrated_in_yeast`` is whether anyone has got it working in this organism. A
    promoter fusion always has, because any promoter can be put in front of any fluorophore.
    A biochemical sensor is one fixed protein that has to be ported, and several widely
    cited ones never have been -- recommending a build around those is recommending a
    research project rather than an experiment.

    ``spectral_slot`` is the channel a sensor occupies exclusively. A biochemical sensor is
    one fixed protein, so its spectrum is fixed too and two sensors sharing a slot cannot
    both be built. A promoter fusion has no slot: the promoter is separable from the
    protein, so it can be built onto whichever fluorophore is still free.
    """

    module: str
    element: str
    marker_gene: str
    also_reads: dict[str, float] = field(default_factory=dict)
    source: str = ""
    kind: Kind = Kind.TRANSCRIPTIONAL
    integrates_hours: float = _INTEGRATION_HOURS
    spectral_slot: str | None = None
    basal: float = _DEFAULT_BASAL
    demonstrated_in_yeast: bool = True
    validated_ph: tuple[float, float] | None = CYTOSOLIC_PH

    @property
    def ph_sensitive(self) -> bool:
        """Whether a reading cannot be told apart from an acidification on its own."""
        return self.validated_ph is None

    @property
    def ph_uncharacterised(self) -> bool:
        """Whether the cytosol goes somewhere the sensor was never shown to be flat over."""
        if self.validated_ph is None:
            return False
        low, high = self.validated_ph
        return low > CYTOSOLIC_PH[0] or high < CYTOSOLIC_PH[1]


MODULES = {
    "ESR": Module(
        "Msn2/Msn4", "STRE", "Martinez-Pastor 1996, PMID 8641288; Gasch 2000 environmental "
        "stress response"),
    "UPR": Module(
        "Hac1", "UPRE", "Cox & Walter, PMID 8898193 -- Ire1 splices HAC1 mRNA, so UPR "
        "activation is post-transcriptional and HAC1 carries no STRE"),
    "oxidative": Module(
        "Yap1/Skn7", "YRE", "Yap1/Skn7 adaptive response, FEMS Yeast Res 8:1214; no "
        "documented Msn2/4 transcription of YAP1"),
    "heat": Module(
        "Hsf1", "HSE", "Boy-Marcotte 1999, PMID 10411744 -- Hsf1 and Msn2/4 regulons "
        "overlap only slightly; Ciccarelli 2023, PMID 37467033 -- coupling is compensatory"),
    "osmotic": Module(
        "Hog1 -> Sko1/Hot1", "CRE/STRE", "Hog1 acts on Msn2/4 activity, not on MSN2/MSN4 "
        "transcription; Gorner 1998, PMID 9472026 -- nuclear entry survives hog1 deletion"),
    "proteasome": Module(
        "Rpn4", "PACE", "Owsianik 2002, PMID 11918814 and Hahn 2006, PMID 16556235 -- YRE "
        "and HSE in the RPN4 promoter, both confirmed by mutagenesis; RPN4 has no STRE",
        driven_by={"oxidative": 0.45, "heat": 0.40}),
    "iron": Module(
        "Aft1/Aft2", "FeRE", "Yamaguchi-Iwai 1996, PMID 8670839; iron sensing is Fe-S "
        "signalling and independent of the ESR. Salin 2008, PMID 18627600 puts a weak YRE "
        "upstream of AFT2, single-source and never replicated, so the weight is small",
        driven_by={"oxidative": 0.10}),
    "cell_wall": Module(
        "Slt2/Mpk1 -> Rlm1", "RLM1 box", "Jung & Levin 1999, PMID 10594829 -- Rlm1 carries "
        "most cell wall integrity transcription downstream of Slt2"),
    "dna_damage": Module(
        "Mec1/Rad53 -> Rfx1/Crt1", "X-box", "Huang 1998, PMID 9741624 -- Crt1 represses "
        "RNR genes and is released on checkpoint activation"),
    "calcium": Module(
        "Calcineurin -> Crz1", "CDRE", "Stathopoulos & Cyert 1997, PMID 9407035 -- Crz1 "
        "drives CDRE-dependent transcription under calcineurin control"),
    "carbon": Module(
        "Snf1 -> Adr1/Cat8", "CSRE", "Young 2003, PMID 12676948 -- Snf1 relieves Mig1 "
        "repression and activates Adr1 and Cat8 on glucose withdrawal"),
    "nitrogen": Module(
        "TORC1 -| Gln3/Gat1", "GATAA", "Beck & Hall 1999, PMID 10604478 -- TOR inhibition "
        "releases Gln3 into the nucleus and derepresses NCR genes"),
    "hypoxia": Module(
        "Rox1/Upc2/Hap1", "AR1/SRE", "Kwast 1998, PMID 9510529 -- heme-dependent Rox1 "
        "repression lifts under low oxygen"),
    "copper": Module(
        "Ace1/Mac1", "MRE/CuRE", "Thiele 1988, PMID 3043194 -- Ace1 activates CUP1 on "
        "copper excess; Mac1 answers copper starvation"),
    "zinc": Module(
        "Zap1", "ZRE", "Zhao & Eide 1997, PMID 9271382 -- Zap1 autoregulates and drives "
        "ZRT1 and ZRT2 under zinc limitation"),
    "sulfur": Module(
        "Met4/Met31/Met32", "Met box", "Thomas & Surdin-Kerjan 1997, PMID 9409150 -- Met4 "
        "activates the sulfur assimilation regulon"),
    "xenobiotic": Module(
        "Pdr1/Pdr3", "PDRE", "Balzi & Goffeau pleiotropic drug resistance network; Salin "
        "2008, PMID 18627600 -- Rpn4 and Pdr1 form a positive loop",
        driven_by={"proteasome": 0.20}),
    "retrograde": Module(
        "Rtg1/Rtg3", "R box", "Liao & Butow 1993, PMID 8422683 -- mitochondrial dysfunction "
        "signals to CIT2 through Rtg1 and Rtg3"),
    "alkaline_ph": Module(
        "Rim101", "Rim101 site",
        "Xu & Mitchell 2001, PMID 11698381 -- Rim101p is activated by cleavage of its "
        "C-terminal PEST-like region, in S. cerevisiae. Lamb 2001, PMID 11050096 supplies "
        "the pH half: it identifies the alkaline response genes and shows the "
        "Rim101p-dependent ones also need Rim13p, the protease that processes Rim101p. "
        "This cited Lamb & Mitchell 2003, PMID 12509465, which is a real and on-topic "
        "paper reporting something else -- Rim101p repressing NRG1 and SMP1 -- and so "
        "supported the transcription factor but not the cleavage claim beside it"),
    "redox": Module(
        "glutathione pool", "none", "Gutscher 2008, PMID 18469822 -- roGFP2-Grx1 reports "
        "the glutathione redox potential as an equilibrium ratio, not a regulon"),
    "peroxide": Module(
        "cytosolic H2O2", "none", "Pak 2020, PMID 32130885 -- HyPer7 reports the peroxide "
        "pool directly and reversibly"),
    "atp": Module(
        "cytosolic ATP", "none", "Yaginuma 2014, PMID 25283467 -- QUEEN sensor design; "
        "Takaine 2019, PMID 30858198 -- the yeast demonstration"),
    "ph": Module(
        "cytosolic pH", "none", "Miesenbock 1998, PMID 9671304 -- ratiometric pHluorin "
        "reports cytosolic pH"),
    "nadh": Module(
        "NADH:NAD ratio", "none", "Hung 2011, PMID 21982714 -- Peredox reports the cytosolic "
        "NADH:NAD ratio"),
}

_LETHAL_MULTIPLE = 3.0
_GENERAL_POTENCY = 2.5
METABOLITE_POOLS = frozenset({"redox", "peroxide", "atp", "ph", "nadh"})
"""Modules that are metabolite pools rather than transcriptional programmes.

Public because a caller reaching for E-Flux has to know them. They carry activity and have
no transcription factor, so no regulon exists to scale a reaction bound with, and
`bridge/regulation.py::gene_scales` REFUSES them rather than dropping them -- silently
excluding them would make the layer look more specific than it is. A caller excludes them
deliberately, by name, and says so on the result.
"""

_POOLS = METABOLITE_POOLS

# A pool leads the regulon only where the agent reaches it without an enzyme. A weak acid
# entering the cytosol or a blocked respiratory chain draining ATP is physical chemistry
# and moves at the concentration of the insult. Tier 0: the mechanism is sound, the factor
# is asserted.
_DIRECT_POOLS = frozenset({"atp", "ph"})
_DIRECT_POOL_LEAD = 0.35

# These have no default, because one agent's two pool arms can go opposite ways. H2O2
# reaches the cytosolic H2O2 pool by diffusion and LEADS it -- HyPer7 detects 20 uM against
# a 150 uM regulon EC50 (Kritsiligkou 2021, PMID 34118234) -- but reaches glutathione only
# through a peroxidase and LAGS it by more than tenfold, because k is 18-26 /M/s for the
# thiolate at ~2% of GSH (Winterbourn 1999, PMID 10468205) against ~10^7 /M/s for Tsa1
# (Ogusucu 2007, PMID 17210445). Same agent, same file, opposite orderings. So module_ec50
# refuses an arm here that does not carry its own value and its own source.
_POOLS_NEEDING_EVIDENCE = _POOLS - _DIRECT_POOLS
_LETHAL_HILL = 2.5

_CASCADE_ORDER = (
    "ESR", "UPR", "oxidative", "heat", "osmotic", "proteasome", "iron", "cell_wall",
    "dna_damage", "calcium", "carbon", "nitrogen", "hypoxia", "copper", "zinc", "sulfur",
    "xenobiotic", "retrograde", "alkaline_ph", "redox", "peroxide", "atp", "ph", "nadh",
)

STRESSORS = {
    "DTT": Stressor("dithiothreitol", "mM", 1.0,
                    {"UPR": 1.0, "ESR": 0.45, "redox": -0.90, "proteasome": 0.35},
                    "EC50 1.0 mM is UNSOURCED and the number to distrust in this file. It "
                    "was attributed to PMC7646510 = MacGilvray 2020, PMID 32597660, which "
                    "contains no 1.0 mM value: it doses 2.5 mM in BY4741 in YPD and the "
                    "culture acclimates to a new growth rate. The plates cannot supply it "
                    "either -- induction is still climbing at 1 mM, the last dose where the "
                    "culture is healthy, so no maximum was reached and no half of one can "
                    "be placed. The only published DTT dose-response is Pincus 2010, PMID "
                    "20625545: half-maximal UPR at 2.2 mM in W303a, under 10% of maximum at "
                    "or below 1.5 mM. Three measurements therefore disagree by more than "
                    "twofold, and 1.0 mM is kept because moving it to 2.2 would put "
                    "half-induction above the measured lethal dose while the plates read "
                    "1.38-1.47 fold induction well below it. What separates them is "
                    "testable and untested: medium (YPD supplies thiol precursors, "
                    "synthetic does not, and DATA_INVENTORY.md does not record which was "
                    "used), genotype (BY4741 is met15 delta 0, Brachmann 1998, PMID "
                    "9483801), and readout (a promoter fusion over hours against HAC1 "
                    "splicing). Lethal dose 1.55 mM IS measured -- growth falls to half of "
                    "control at 1.45 and 1.65 mM on UPRE1 and UPRE2. A window that narrow "
                    "is why the ladder cannot characterise the promoter",
                    lethal_dose=1.55,
                    target_ec50={"redox": 1.0}),
    "tunicamycin": Stressor("tunicamycin", "ug/mL", 1.0,
                            {"UPR": 1.0, "ESR": 0.35, "proteasome": 0.30},
                            "classical UPR inducer; KAR2/PDI 2-2.5 fold, NAR 46:1139"),
    "H2O2": Stressor("hydrogen peroxide", "mM", 0.15,
                     {"oxidative": 1.0, "ESR": 0.50, "peroxide": 1.0, "redox": 0.70},
                     "Oxidative EC50 0.15 mM from Goulev 2017, PMID 28418333: Yap1-GFP "
                     "nuclear entry is partial at 0.1 mM and saturates above 0.2 mM, in "
                     "S288C under sustained microfluidic delivery. It was 0.5 mM, which is "
                     "where Goulev sees growth arrest, not half-induction, and which put "
                     "the derived ESR EC50 at 1.25 mM -- above this stressor's own measured "
                     "lethal dose. The pool arms were pinned at 0.175 on the grounds that "
                     "Goulev measured Yap1, not chemistry. That was wrong twice: 0.175 was "
                     "not a measurement but 0.5 x the old blanket pool factor, the product "
                     "of the very EC50 being corrected and a default that had the ordering "
                     "gave both pool arms one ordering. They go OPPOSITE ways, which is "
                     "why no multiplier serves both. Peroxide (cytosolic H2O2, HyPer7) "
                     "LEADS at 0.05 mM: Kritsiligkou 2021, PMID 34118234 finds the lowest "
                     "detectable HyPer7 response at ~20 uM in BY4742, in S. cerevisiae with "
                     "this module's own sensor, against 100 uM for the first detectable "
                     "oxidised Yap1 (Delaunay 2000, PMID 11013218, with 25 and 50 uM tested "
                     "and below detection). 0.05 is a BRACKET between that 20 uM floor and "
                     "this agent's 0.15 mM regulon EC50, not a fitted EC50 -- no HyPer7 "
                     "half-max exists for yeast. Glutathione (redox) LAGS at 2.0 mM, which "
                     "IS a half-activation: Ayer 2013, PMID 23762325 sees no cytosolic "
                     "E_GSH shift at 0.2 mM over 60 min, its first significant shift at 1 mM "
                     "(40-50 mV, partly a pH artefact by its own admission) and half the "
                     "roGFP2 probe oxidised at 2 mM. Morgan 2013, PMID 23242256 tested "
                     "0.5-20 mM and went no lower, under a supplementary figure its authors "
                     "titled 'The cytosolic glutathione pool is robustly resistant to "
                     "perturbation'. Note 2.0 mM is TWICE this agent's lethal dose, so the "
                     "glutathione arm is unreachable inside healthy_ladder and roGFP2-Grx1 "
                     "reads flat across every well the panel can spend -- the same trap the "
                     "ESR arm fell into at 1.25 mM, and a fact about the panel rather than a "
                     "number to shrink. Lethal dose 1.0 mM IS measured -- "
                     "growth halves at 0.79 and 1.24 mM on NativeYap1 and AlteredYap1; "
                     "Goulev's 0.6 mM full arrest is a perfusion, and Tomalin 2016, PMID "
                     "26944189 explains the gap: a bolus is consumed by thiol buffering. "
                     "Its proteasome arm is not listed because Yap1 acts on RPN4, which the "
                     "cascade carries",
                     lethal_dose=1.0,
                     target_ec50={"peroxide": 0.05, "redox": 2.0}),
    "menadione": Stressor("menadione", "uM", 100.0,
                          {"oxidative": 0.90, "ESR": 0.55, "peroxide": 0.60,
                           "redox": 0.75, "nadh": -0.40},
                          "superoxide-generating redox cycler; drains reducing equivalents. "
                          "Pool arms assert no separation from the regulon because none is "
                          "measured, and the nearest evidence argues against a lead: Ayer "
                          "2013, PMID 23762325 finds 2 mM paraquat -- same redox-cycling "
                          "class -- gives no change in cytosolic E_GSH while still inducing "
                          "SOD2. Azevedo 2003, PMID 14556853 has menadione activating Yap1 "
                          "partly by direct cysteine modification, a route that never "
                          "touches the pool. The 100 uM EC50 is itself unsourced; Gasch "
                          "used 1 mM menadione bisulfite",
                          target_ec50={"peroxide": 100.0, "redox": 100.0, "nadh": 100.0}),
    "diamide": Stressor("diamide", "mM", 1.5,
                        {"oxidative": 0.80, "redox": 1.0, "ESR": 0.60},
                        "Gasch 2000: thiol-specific oxidant, strong ESR inducer at 1.5 mM. "
                        "Unlike peroxide, diamide oxidises glutathione directly, so a pool "
                        "lead is mechanistically arguable -- but no dose-response for "
                        "cytosolic E_GSH under diamide could be found in S. cerevisiae. "
                        "Every yeast paper uses it as the SATURATING calibration point for "
                        "roGFP2 and Grx1-roGFP2: 20 mM in Kojer 2012, PMID 22705944, 10 mM "
                        "in Ayer 2013, PMID 23762325. That is 7-13x this EC50, so the pool "
                        "is far from oxidised where the regulon fires. No separation is "
                        "asserted because the sign of it is not established",
                        target_ec50={"redox": 1.5}),
    "heat": Stressor("temperature shift", "degC above 30", 6.0,
                     {"heat": 1.0, "ESR": 0.70},
                     "Boy-Marcotte 1999: 25->38 C, half of induced proteins Msn2/4-dependent. "
                     "Its proteasome arm runs through Hsf1 on RPN4, carried by the cascade"),
    "NaCl": Stressor("sodium chloride", "M", 0.5,
                     {"osmotic": 1.0, "ESR": 0.65, "calcium": 0.25},
                     "Hog1 osmotic response; sodium also engages calcineurin through ENA1"),
    "sorbitol": Stressor("sorbitol", "M", 1.0, {"osmotic": 0.95, "ESR": 0.45},
                         "Gasch 2000 hyperosmotic shock without an ionic component"),
    "glucose_starvation": Stressor("glucose withdrawal", "fraction removed", 0.6,
                                   {"ESR": 1.0, "carbon": 1.0, "atp": -0.60},
                                   "carbon starvation is a core ESR trigger and relieves "
                                   "Mig1 repression through Snf1"),
    "antimycin_A": Stressor("antimycin A", "uM", 5.0,
                            {"atp": -0.85, "retrograde": 0.90, "ESR": 0.35, "nadh": 0.50},
                            "complex III inhibitor; mitochondrial dysfunction is the "
                            "canonical Rtg1/Rtg3 trigger. The NADH arm asserts no separation: "
                            "the nearest measurement is Calabrese 2019, PMID 31389622, where "
                            "10 uM -- twice this EC50 -- moves the peroxide sensor strongly "
                            "but leaves E_GSH with only a small deflection. ATP keeps the "
                            "direct-pool lead: blocking the chain drains the adenylate pool "
                            "without an enzyme in between",
                            target_ec50={"nadh": 5.0}),
    "BPS": Stressor("bathophenanthroline disulfonate", "uM", 50.0,
                    {"iron": 1.0, "ESR": 0.20}, "iron chelation activates the Aft1 regulon"),
    "MG132": Stressor("proteasome inhibitor", "uM", 40.0,
                      {"proteasome": 1.0, "ESR": 0.45, "heat": 0.30},
                      "proteasome inhibition raises Rpn4 and cytosolic misfolding"),
    "MMS": Stressor("methyl methanesulfonate", "percent", 0.02,
                    {"dna_damage": 1.0, "ESR": 0.50, "proteasome": 0.25},
                    "alkylating agent; Hahn 2006 reports RPN4 induction at 0.02 percent"),
    "hydroxyurea": Stressor("hydroxyurea", "mM", 100.0,
                            {"dna_damage": 0.95, "ESR": 0.30},
                            "ribonucleotide reductase inhibition stalls replication forks"),
    "congo_red": Stressor("congo red", "ug/mL", 50.0,
                          {"cell_wall": 1.0, "ESR": 0.30},
                          "binds nascent chitin and activates the Slt2 pathway"),
    "caffeine": Stressor("caffeine", "mM", 10.0,
                         {"cell_wall": 0.70, "nitrogen": 0.55, "ESR": 0.35},
                         "activates cell wall integrity signalling and inhibits TORC1"),
    "calcium_chloride": Stressor("calcium chloride", "mM", 200.0,
                                 {"calcium": 1.0, "ESR": 0.25},
                                 "calcineurin/Crz1 activation through cytosolic calcium"),
    "rapamycin": Stressor("rapamycin", "nM", 20.0,
                          {"nitrogen": 1.0, "ESR": 0.60, "carbon": 0.30},
                          "TORC1 inhibition releases Gln3 and induces the ESR"),
    "cobalt_chloride": Stressor("cobalt chloride", "mM", 1.0,
                                {"hypoxia": 0.90, "iron": 0.45, "oxidative": 0.30,
                                 "ESR": 0.25},
                                "hypoxia mimetic; also perturbs iron handling"),
    "copper_sulfate": Stressor("copper sulfate", "mM", 1.0,
                               {"copper": 1.0, "oxidative": 0.45, "ESR": 0.30},
                               "Ace1 activation of CUP1 with a Fenton-driven oxidative arm"),
    "TPEN": Stressor("TPEN zinc chelator", "uM", 10.0,
                     {"zinc": 1.0, "ESR": 0.25}, "zinc limitation activates Zap1"),
    "methionine_starvation": Stressor("methionine withdrawal", "fraction removed", 0.8,
                                      {"sulfur": 1.0, "ESR": 0.30},
                                      "Met4 activation on sulfur amino acid limitation"),
    "fluconazole": Stressor("fluconazole", "ug/mL", 30.0,
                            {"xenobiotic": 1.0, "hypoxia": 0.50, "cell_wall": 0.25,
                             "ESR": 0.20},
                            "azole efflux runs through Pdr5; ergosterol block engages Upc2"),
    "acetic_acid": Stressor("acetic acid", "mM", 60.0,
                            {"ph": -1.0, "ESR": 0.60, "cell_wall": 0.30},
                            "weak acid entry acidifies the cytosol at pH 4.5"),
    "sodium_hydroxide": Stressor("alkaline shift", "pH units above 5", 3.0,
                                 {"alkaline_ph": 1.0, "ph": 0.70, "iron": 0.40,
                                  "calcium": 0.35, "ESR": 0.35},
                                 "alkaline pH activates Rim101 and limits iron uptake"),
}

REPORTERS = {
    "STRE-general": Reporter(
        "ESR", "STRE", "HSP12/CTT1",
        source="STRE only; induction is lost in msn2 msn4 and retained in yap1"),
    "UPRE-ER": Reporter(
        "UPR", "UPRE", "KAR2", also_reads={"heat": 0.25},
        source="KAR2 carries both a UPRE and an HSE, so it is heat-inducible"),
    "TRX2-oxidative": Reporter(
        "oxidative", "YRE", "TRX2",
        source="Kuge & Jones 1994, PMID 8313910 -- Yap1 drives TRX2 on hydroperoxide "
               "challenge; Morgan 1997, PMID 9118942 -- the Yap1-independent arm at this "
               "promoter is Skn7, which binds TRX2 directly and co-operates with Yap1 on "
               "it. Skn7 is inside this reporter's own module, so the induction TRX2 keeps "
               "in yap1 is on-diagonal and carries no crosstalk weight. It used to carry "
               "one to the ESR, which attributed a Skn7 effect to Msn2/4"),
    "HSE-heat": Reporter(
        "heat", "HSE", "HSP104", also_reads={"ESR": 0.35},
        source="HSP104 carries both HSE and STRE"),
    "STRE-osmotic": Reporter(
        "osmotic", "STRE", "GPD1", also_reads={"ESR": 0.40},
        source="GPD1 takes Hog1 output through Sko1/Hot1 and through Msn2/4"),
    "PACE-proteasome": Reporter(
        "proteasome", "PACE", "RPT1", source="PACE only; RPN4 supplies all upstream input"),
    "FeRE-iron": Reporter(
        "iron", "FeRE", "FIT3", source="FeRE only"),
    "RLM1-cellwall": Reporter(
        "cell_wall", "RLM1 box", "PIR3/CRH1", source="Rlm1-dependent wall remodelling genes"),
    "Xbox-dna": Reporter(
        "dna_damage", "X-box", "RNR3", source="Crt1-repressed, derepressed on checkpoint"),
    "CDRE-calcium": Reporter(
        "calcium", "CDRE", "PMC1/CMK2", source="the standard calcineurin reporter"),
    "CSRE-carbon": Reporter(
        "carbon", "CSRE", "ADH2/ICL1", also_reads={"ESR": 0.25},
        source="Adr1 and Cat8 sites; glucose withdrawal also triggers the ESR"),
    "GATA-nitrogen": Reporter(
        "nitrogen", "GATAA", "DAL5/MEP2", source="nitrogen catabolite repression reporter"),
    "AR1-hypoxia": Reporter(
        "hypoxia", "AR1", "ANB1/DAN1", source="Rox1-repressed hypoxic genes"),
    "CuRE-copper": Reporter(
        "copper", "MRE", "CUP1", source="Ace1 metal response elements"),
    "ZRE-zinc": Reporter(
        "zinc", "ZRE", "ZRT1", source="Zap1 zinc responsive elements"),
    "METbox-sulfur": Reporter(
        "sulfur", "Met box", "MET17", source="Met4/Met31/Met32 sites"),
    "PDRE-xenobiotic": Reporter(
        "xenobiotic", "PDRE", "PDR5", source="Pdr1/Pdr3 drug response elements"),
    "Rbox-retrograde": Reporter(
        "retrograde", "R box", "CIT2", source="Rtg1/Rtg3 sites in the CIT2 promoter"),
    "RIM101-alkaline": Reporter(
        "alkaline_ph", "Rim101 site", "ENA1", also_reads={"calcium": 1.5},
        source="Petrezselyova 2016, PMID 27362362 -- calcineurin/Crz1 supplies about 60% of "
               "the early ENA1 alkaline response through two stress-responsive Crz1 sites, "
               "so calcium is the dominant arm at 60:40 and the weight is 1.5 against a "
               "self-loading of 1.0. It was 0.30, which had the dominant arm as the minor "
               "one. Platara 2006, PMID 17023428 makes 1.5 a floor rather than a point "
               "estimate: the residual 40% is Rim101 and Snf1 together, and only Rim101 is "
               "this reporter's own module. Renormalising to 0.4/0.6 instead would mean "
               "unpinning the diagonal for all 24 reporters to encode one measurement"),
    "roGFP2-Grx1": Reporter(
        "redox", "none", "roGFP2-Grx1", kind=Kind.RATIOMETRIC, integrates_hours=0.0,
        spectral_slot="green-ratio", validated_ph=None,
        source="Gutscher 2008, PMID 18469822. Established in yeast: Kojer 2012, PMID "
               "22705944 puts cytosolic E_GSH at -306 mV and matrix at -301 mV; Ayer 2013, "
               "PMID 23762325 reports -350 mV with a measured rather than assumed cytosolic "
               "pH. The two literatures differ by about 40 mV on that alone. Responds within "
               "5 min of a peroxide bolus and recovers within 4 min of washout"),
    "HyPer7": Reporter(
        "peroxide", "none", "HyPer7", kind=Kind.RATIOMETRIC, integrates_hours=0.0,
        spectral_slot="green-ratio",
        source="Pak 2020, PMID 32130885 -- engineered pH-insensitive, which is much of the "
               "point of it over earlier HyPer variants"),
    "QUEEN-2m": Reporter(
        "atp", "none", "QUEEN-2m", kind=Kind.RATIOMETRIC, integrates_hours=0.0,
        spectral_slot="green-ratio", validated_ph=(7.3, 8.8),
        source="Yaginuma 2014, PMID 25283467 -- ratiometric ATP sensor, characterised in "
               "E. coli, NOT in yeast: that paper is 'Diversity in ATP concentrations in "
               "a single bacterial cell population'. The yeast demonstration is Takaine "
               "2019, PMID 30858198, with the Bio-protocol imaging method in PMID "
               "33654827 -- also 2019 in PubMed, though its PMID looks 2021 because "
               "Bio-protocol indexed it late. "
               "Kd 4.5 mM ATP at 25 C, flat against Mg2+ over 1-2 mM, and flat against pH "
               "over 7.3-8.8, which stops above where a starved cytosol goes. "
               "Note Larsson 1997, PMID 9393686: in yeast chemostats glycolytic flux "
               "correlates NEGATIVELY with intracellular ATP and not at all with the "
               "ATP/ADP ratio, so this reads a pool and not a rate"),
    "pHluorin": Reporter(
        "ph", "none", "ratiometric pHluorin", kind=Kind.RATIOMETRIC, integrates_hours=0.0,
        spectral_slot="green-ratio", validated_ph=None,
        source="Miesenbock 1998, PMID 9671304 -- reports pH, so it is the partner a "
               "pH-sensitive redox sensor needs rather than a rival for the channel"),
    "Peredox": Reporter(
        "nadh", "none", "Peredox", kind=Kind.RATIOMETRIC, integrates_hours=0.0,
        spectral_slot="green-ratio", demonstrated_in_yeast=False,
        source="Hung 2011, PMID 21982714 -- cytosolic NADH:NAD ratio, but never applied in "
               "yeast. Neither has SoNar (Zhao 2015, PMID 25955212). The only genetically "
               "encoded NADH readout in S. cerevisiae is a GPD2-promoter fusion (Knudsen "
               "2014, PMID 25401080), which reports a transcriptional proxy and not a ratio"),
}


def demonstrated_in_yeast() -> list[str]:
    """Reporters someone has actually got working in S. cerevisiae.

    A build restricted to these is one the wet lab can start on; the rest would each need
    a sensor ported first, which is a different project with its own failure modes.
    """
    return [n for n, r in REPORTERS.items() if r.demonstrated_in_yeast]


def transcriptional_reporters() -> list[str]:
    """Promoter fusions, which integrate over hours and are diluted by growth."""
    return [n for n, r in REPORTERS.items() if r.kind is Kind.TRANSCRIPTIONAL]


def ratiometric_reporters() -> list[str]:
    """Biochemical sensors, whose ratio settles in seconds and ignores growth entirely."""
    return [n for n, r in REPORTERS.items() if r.kind is Kind.RATIOMETRIC]


def viability(stressor: str, dose: float) -> float:
    """Fraction of transcriptional capacity left at a dose, between one and zero.

    A dose high enough to induce strongly is also high enough to start stopping the cell,
    so the response is biphasic and a saturating Hill alone cannot produce it. This is what
    turns the curve back down past its peak.
    """
    if stressor not in STRESSORS:
        raise KeyError(f"no stressor {stressor!r}; have {sorted(STRESSORS)}")
    lethal = STRESSORS[stressor].lethal_dose
    d = max(float(dose), 0.0)
    if d <= 0.0 or lethal <= 0.0:
        return 1.0
    return float(1.0 / (1.0 + (d / lethal) ** _LETHAL_HILL))


def module_ec50(stressor: str, module: str) -> float | None:
    """Concentration at which one agent half-activates one module.

    A single saturation for every module a stressor touches makes its whole dose series
    trace one ray, so two regulons driven by the same agent can never be separated however
    finely it is dosed.

    Whether a pool leads or lags depends on how the agent reaches it, and one agent can
    do both. Peroxide reaching the cytosolic H2O2 pool is diffusion and leads; the same
    peroxide reaching glutathione needs a peroxidase and lags, because Yap1 is wired to a
    faster catalytic sensor -- which is the point of the architecture. Acid reaching
    cytosolic pH needs no enzyme either and keeps a default. Everything else carries a
    per-agent value and this function refuses to guess one. The general stress response trails because it answers substantial
    damage rather than a specific signal.

    Raises:
        KeyError: on an unknown stressor.
        ValueError: on a pool arm needing evidence with no value in ``target_ec50``.
            Inheriting a blanket multiplier here is what previously gave one agent's two
            pool arms the same ordering when the measurements disagree.
    """
    if stressor not in STRESSORS:
        raise KeyError(f"no stressor {stressor!r}; have {sorted(STRESSORS)}")
    spec = STRESSORS[stressor]
    if module not in spec.targets:
        return None
    if module in spec.target_ec50:
        return float(spec.target_ec50[module])
    if module in _POOLS_NEEDING_EVIDENCE:
        raise ValueError(
            f"{stressor!r} drives {module!r} with no EC50 for it. One agent's pool arms "
            "can lead and lag at once -- H2O2 leads cytosolic peroxide and lags "
            "glutathione -- so there is no default to fall back on; add it to target_ec50 "
            "with its source.")
    if module in _DIRECT_POOLS:
        return spec.ec50 * _DIRECT_POOL_LEAD
    if module == "ESR":
        return spec.ec50 * _GENERAL_POTENCY
    return spec.ec50


def growth_rate(stressor: str, dose: float, mu_max: float = _MU_MAX) -> float:
    """Specific growth rate under one agent at one dose.

    The lethal dose was measured as where growth falls to half of the unstressed control,
    so the same curve that gives viability gives growth. Having it here is what lets the
    generator carry the dilution confound: reporter accumulates in cells that have stopped
    diluting it, and nothing downstream can see that unless growth is part of the model.
    """
    return max(mu_max * viability(stressor, dose), 0.0)


def healthy_ladder(stressor: str, n_doses: int = 5, ceiling: float = 0.55) -> np.ndarray:
    """A dose series that stays inside the range where a culture still reports.

    Spending wells above the dose that halves growth is what the uploaded plates did, and
    those wells returned no readable activity. Spanning up to a fraction of the measured
    lethal dose keeps every rung interpretable and puts the wells where the rising limb is.

    Args:
        stressor: Agent to dose.
        n_doses: Rungs on the ladder, spaced geometrically.
        ceiling: Highest dose as a fraction of the lethal dose.
    """
    if stressor not in STRESSORS:
        raise KeyError(f"no stressor {stressor!r}; have {sorted(STRESSORS)}")
    spec = STRESSORS[stressor]
    top = ceiling * spec.lethal_dose
    return np.geomspace(top / 16.0, top, n_doses)


def _propagate(direct: dict[str, float]) -> dict[str, float]:
    """Add cascade contributions, resolving drivers before the modules they drive."""
    activity = dict(direct)
    for name in _CASCADE_ORDER:
        for driver, weight in MODULES[name].driven_by.items():
            activity[name] += weight * activity[driver]
    return activity


def module_response(stressor: str, dose: float) -> dict[str, float]:
    """Activation of every module by one stressor at one dose, saturating with dose.

    A stressor's own targets are what it acts on directly; the cascade then carries those
    into the regulons they drive, so a purely oxidative agent still raises the proteasome
    without that having to be written into the agent. Targets may be negative, because a
    reductant and an oxidant displace the same pool in opposite directions.

    Induction saturates with dose and viability falls with it, so the curve peaks and then
    turns down. Viability multiplies every module alike: a cell losing the capacity to
    transcribe cannot spare one regulon.
    """
    if stressor not in STRESSORS:
        raise KeyError(f"no stressor {stressor!r}; have {sorted(STRESSORS)}")
    spec = STRESSORS[stressor]
    d = max(float(dose), 0.0)
    alive = viability(stressor, d)
    direct = {}
    for name in MODULES:
        weight = spec.targets.get(name, 0.0)
        if weight == 0.0 or d <= 0.0:
            direct[name] = 0.0
            continue
        potency = module_ec50(stressor, name)
        direct[name] = weight * (d / (potency + d)) * alive
    return _propagate(direct)


def combination_response(doses: dict[str, float]) -> dict[str, float]:
    """Module activation under several stressors at once, additive on module activity.

    Two agents driving the same transcription factor add at that factor, so activity sums
    while each agent still saturates in its own dose. Combinations matter because a panel
    dosed one agent at a time only ever traces one ray per stressor, leaving the directions
    between them unvisited however many doses are run.
    """
    unknown = [name for name in doses if name not in STRESSORS]
    if unknown:
        raise KeyError(f"no stressor(s) {sorted(unknown)}; have {sorted(STRESSORS)}")
    total = {name: 0.0 for name in MODULES}
    for stressor, dose in doses.items():
        for module, activity in module_response(stressor, dose).items():
            total[module] += activity
    return total


def reporter_loadings(names) -> np.ndarray:
    """How each reporter loads on every module, from the elements in its promoter.

    Only promoter content appears here. A module driving another is already carried by the
    cascade in `module_response`, so reading it again at the promoter would count the same
    regulation twice.
    """
    order = list(MODULES)
    matrix = np.zeros((len(names), len(order)))
    for row, name in enumerate(names):
        if name not in REPORTERS:
            raise KeyError(f"no reporter {name!r}; have {sorted(REPORTERS)}")
        reporter = REPORTERS[name]
        matrix[row, order.index(reporter.module)] = 1.0
        for module_name, weight in reporter.also_reads.items():
            matrix[row, order.index(module_name)] += weight
    return matrix
