"""The mechanistic state vector, with every variable's units, time constant and provenance.

WHY A CONTAINER AND NOT A TUPLE. The generator this package is growing out of is algebra:
instantaneous Hill functions, no mRNA, no TF binding, no folding. Adding real ODE chains means
adding states, and a state that arrives as ``y[3]`` with no units and no measured time constant
is exactly how a model acquires parameters nobody can defend. Every :class:`StateVar` here
carries four things the integrator and the reduction audit both need:

  * its **units**, so a tolerance can be scaled and a plot can be labelled;
  * its **time constant**, in hours, as an INTERVAL -- because for two of the states that decide
    this model the literature disagrees by 9x and by 3.8x, and collapsing that to a point
    estimate is the failure this file exists to prevent;
  * its **encoding** -- whether it relaxes toward a value its driver sets, generates its own
    carrier whose FREQUENCY is the signal, or integrates without a steady state at all. The
    reduction criterion is only valid for the first;
  * what **constrains** it: a named assay, or a declared sweep axis. There is no third option.
    :class:`StateUnidentifiable` is raised at construction, not at use.

THE TIME CONSTANT MAY BE REFUSED. Seven processes the target architecture names have no
measured rate anywhere -- Ire1 clustering, Slt2 activation, Snf1 activation, Rim101 processing,
mRNA export, the yeast H2O2 permeability coefficient, and Venus chromophore oxidation.
:class:`RefusedTimescale` is what goes in the field instead of a plausible number, and it raises
on ``float()``. `pathway/capacity.py::CapacityUnmeasured` and
`pathway/environment_flux.py::EnvironmentUnmeasured` are the precedents.

SOME TIME CONSTANTS ARE NOT CONSTANTS. Dilution, plasmid loss and the generation counter all
scale as 1/mu, and mu differs between this repository's two operating windows by up to 3.6x --
0.246-0.364 /h on the plates against 0.101 /h at the fed-batch setpoint. Those states declare
``tau_growth_multiple`` instead of a fixed interval, and the window
supplies the growth rate -- so the audit computes tau where it is used rather than storing one
number that is wrong in one of the two vessels. That is the mechanism behind the non-nested
reduction sets, and it is arithmetic rather than a claim.

Every number below is MEASURED, DERIVED with the arithmetic shown, BOUNDED by a measured
completion time, or REFUSED. Provenance is the ``source`` field of each row, and the harvest it
was compiled from is `scratchpad/harvest/timescales.md`, in which every PMID was resolved
through NCBI E-utilities before it was written down.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass
from typing import Iterator, Mapping, Sequence

import numpy as np

__all__ = [
    "Encoding",
    "FEDBATCH_5D",
    "FEDBATCH_GROWTH_PER_H",
    "MechState",
    "PLATE_GROWTH_BAND_PER_H",
    "PLATE_READ_4H",
    "RefusedTimescale",
    "StateUnidentifiable",
    "StateVar",
    "TIMESCALE_CATALOGUE",
    "TimescaleUnmeasured",
    "Window",
    "catalogue_span_hours",
]

_S = 1.0 / 3600.0
_MIN = 1.0 / 60.0


class TimescaleUnmeasured(RuntimeError):
    """Raised when a caller asks for a time constant nothing has measured.

    Deliberately not a ``ValueError``: the caller did nothing wrong. The question is
    well-formed and the literature cannot answer it, which is a different thing from a bad
    argument and should not be caught by the same ``except``.
    """


class StateUnidentifiable(ValueError):
    """A state was declared with neither a named constraint nor a declared sweep axis.

    Criterion (b) of the stopping rule, enforced at construction rather than in review. A
    state that nothing measures and nobody declared as an axis is a free direction in the
    posterior wearing a mechanism's name.
    """


@dataclass(frozen=True)
class RefusedTimescale:
    """A time constant no source measures. Raises on ``float()`` rather than defaulting.

    The time-constant-shaped instance of this package's REFUSED idiom. `mech/params.py` carries
    the general one, ``Param.refused``, and these two should be folded together at integration --
    they are deliberately parallel rather than coupled, because they were written at the same
    time and importing a signature still moving would have made each file's churn the other's.

    Args:
        quantity: What would be measured, in words.
        reason: Why the number does not exist -- what the nearest source does and does not say.
        source: The nearest source, named so the refusal can be overturned by reading it.
    """

    quantity: str
    reason: str
    source: str

    def __float__(self) -> float:
        raise TimescaleUnmeasured(
            f"{self.quantity}: no measured time constant. {self.reason} (nearest: {self.source})"
        )


class Encoding(enum.Enum):
    """How a state carries information, which decides whether the ratio test applies to it.

    The reduction criterion tau/T < 0.146 is derived for a state that RELAXES toward a
    quasi-steady value its driver sets. Two other kinds of state exist here and the criterion
    is not valid for either.
    """

    LEVEL = "level"
    """Relaxes toward a value set by its driver. The only kind the ratio test licenses."""

    FREQUENCY = "frequency"
    """A self-generated carrier whose FREQUENCY is the signal, not its mean or its width.
    Cai, Dalal & Elowitz 2008, Nature 455:485-90, PMID 18818649: Crz1 bursts last ~2 min and
    calcium sets their frequency, not their duration. Hao & O'Shea 2012, Nat Struct Mol Biol
    19:31-9, PMID 22179789: Msn2 is frequency-modulated under glucose limitation and
    amplitude-modulated under oxidative stress -- the stressor picks the encoding. Replacing
    such a train by its mean is a mean-field approximation and NOT a quasi-steady-state
    reduction; Cai 2008 is a direct demonstration that the mean does not determine the output.
    """

    INTEGRATING = "integrating"
    """Accumulates; has no quasi-steady value to be reduced TO. Biomass under dx/dt = mu*x, a
    product pool, a generation counter. Its 1/mu is a doubling timescale, not a relaxation
    time, so the fast half of the criterion is meaningless for it. It can still be FROZEN when
    its total excursion over the window is below the floor."""


@dataclass(frozen=True)
class Window:
    """A declared integration window: how long the run is, and how fast the cells grow in it.

    The window is an INPUT to the model, not a runtime detail, because the irreducible set is
    different for the two vessels this project uses and the two sets are not nested. A model
    tuned to one is mis-specified for the other.

    Args:
        name: Short identifier used in the audit table.
        duration_h: Length of the run, hours.
        growth_rate_low_per_h, growth_rate_high_per_h: The band of specific growth rates the
            vessel actually runs at, MEASURED. Growth-coupled time constants are computed from
            these rather than stored, so the same catalogue serves both windows.
        source: Where the duration and the growth band come from.
    """

    name: str
    duration_h: float
    growth_rate_low_per_h: float
    growth_rate_high_per_h: float
    source: str

    def __post_init__(self) -> None:
        if self.duration_h <= 0:
            raise ValueError(f"window {self.name!r}: duration must be positive, got {self.duration_h}")
        if not 0 < self.growth_rate_low_per_h <= self.growth_rate_high_per_h:
            raise ValueError(
                f"window {self.name!r}: growth band must be positive and ordered, got "
                f"({self.growth_rate_low_per_h}, {self.growth_rate_high_per_h})"
            )


PLATE_GROWTH_BAND_PER_H = (0.24624623855022873, 0.3639682842372266)
"""Specific growth rate on the real 96-well reads, 1/h. DERIVED, and measured not assumed.

mu = ln(od_fold_change) / 4.14 h on every well that passed gate 1, pooled across:

* `outputs/g1_20260722_ERandOxidativeStress_NewProtocol_ANALYSED.csv` (67 wells),
* `outputs/g1_20260728_ERandOxidativeStress_NewProtocol_Replicate2.csv` (64 wells),
* `outputs/g1_20260803_ERandoxidativestress_Replicate3.csv` (52 wells).

The verified current gate-one tables give n = 183, median 0.3248257718, interquartile range
0.2462462386-0.3639682842, full range 0.0716175486-0.4646708144. The IQR is what is carried
here; the tails are dosed wells whose growth the stressor has moved, which is a signal
rather than a window property. This is a current-derived calibration, not a frozen
historical snapshot.

It is a chord estimate over the whole read, not the slope of a fitted window, so it is the
mean specific rate across the read -- which is exactly the quantity a window-averaged time
constant needs. `tests/test_mech_integrate.py` re-derives it from those files and fails if
they move; production keeps a literal rather than reading generated outputs at import.
"""

FEDBATCH_GROWTH_PER_H = 0.101
"""Setpoint of the exponential fed-batch, 1/h.

The lower Elizondo chemostat state, which is the one `fba/fedbatch.py` reproduces comfortably.
Its own `REFUSES_THE_UPPER_ELIZONDO_STATE` records why the upper state (mu = 0.254 /h) is
refused for this vessel: 0.254 is also that strain's measured mu_max, so a fed-batch designed
for it has zero margin and its feed no longer sets mu. A single value rather than a band,
because a fed-batch setpoint is chosen and then MEASURED back off the biomass series -- this
project does not run chemostats and does not assert mu from a pump.
"""

PLATE_READ_4H = Window(
    name="PLATE_READ_4H",
    duration_h=4.14,
    growth_rate_low_per_h=PLATE_GROWTH_BAND_PER_H[0],
    growth_rate_high_per_h=PLATE_GROWTH_BAND_PER_H[1],
    source=(
        "4.14 h is generator/plate.py::PlateConditions.duration_h, the length of the real "
        "Synergy reads; growth band is the current gate-one IQR over 183 passing wells "
        "on the 20260722, 20260728 and 20260803 plates named above"
    ),
)
"""The 96-well plate read. Reporter dilution, FP maturation and the transcript pool all live
inside it; plasmid loss and adaptation do not."""

FEDBATCH_5D = Window(
    name="FEDBATCH_5D",
    duration_h=120.0,
    growth_rate_low_per_h=FEDBATCH_GROWTH_PER_H,
    growth_rate_high_per_h=FEDBATCH_GROWTH_PER_H,
    source=(
        "120 h is the fed-batch run length in docs/research/CALIBRATION_DATA.md row 4 "
        "(Nunta/Watcharawipas 2024, PMID 38921419, 5-L semi-defined fed-batch); mu is the "
        "lower Elizondo state fba/fedbatch.py reproduces"
    ),
)
"""The vessel that makes the titre. Plasmid loss dominates here and is negligible on a plate;
everything the plate read is sensitive to has collapsed to algebra."""


@dataclass(frozen=True)
class StateVar:
    """One named variable of the mechanistic state, with its units and its time constant.

    Exactly one of ``tau_h`` and ``tau_growth_multiple`` must be supplied. The second is for
    states whose time constant is 1/mu times a constant, which the window supplies.

    Args:
        name: Identifier used in the state vector and the audit table.
        units: Physical units of the variable itself, not of its time constant.
        source: Citation for the time constant, with its provenance tag.
        tau_h: ``(low, high)`` interval in hours, equal when one number, or a
            :class:`RefusedTimescale`. An interval because two of the states that decide this
            model have conflicting measurements -- mRNA half-life spans 9x across methods and
            FP maturation 3.8x across two yeast papers -- and a point estimate would hide the
            straddle the audit exists to report.
        tau_growth_multiple: ``tau = multiple / mu``. 1.0 for dilution of a stable protein,
            ln(2) for a generation counter, 1/(-ln(1-loss_per_generation)) for plasmid loss.
        encoding: See :class:`Encoding`. Decides whether the ratio test applies at all.
        constrained_by: The named assay that pins this state on the instruments this project
            owns. Empty only when ``sweep_axis`` is set.
        sweep_axis: Declared sensitivity axis: nothing here measures it and the model reports a
            band rather than a point.
        note: One line of context the source alone does not carry -- a caveat on the
            provenance, or what the reduction would delete.
    """

    name: str
    units: str
    source: str
    tau_h: tuple[float, float] | RefusedTimescale | None = None
    tau_growth_multiple: float | None = None
    encoding: Encoding = Encoding.LEVEL
    constrained_by: str = ""
    sweep_axis: bool = False
    note: str = ""

    def __post_init__(self) -> None:
        if (self.tau_h is None) == (self.tau_growth_multiple is None):
            raise ValueError(
                f"state {self.name!r}: supply exactly one of tau_h and tau_growth_multiple"
            )
        if isinstance(self.tau_h, tuple):
            low, high = self.tau_h
            if not 0 < low <= high:
                raise ValueError(f"state {self.name!r}: tau interval must be positive and ordered")
        if self.tau_growth_multiple is not None and self.tau_growth_multiple <= 0:
            raise ValueError(f"state {self.name!r}: tau_growth_multiple must be positive")
        if not self.constrained_by and not self.sweep_axis:
            raise StateUnidentifiable(
                f"state {self.name!r} declares neither a measurement that constrains it nor a "
                f"sweep axis. Criterion (b): name the assay, or declare the axis."
            )

    @property
    def timescale_is_refused(self) -> bool:
        return isinstance(self.tau_h, RefusedTimescale)

    def taus_in(self, window: Window) -> tuple[float, float]:
        """The ``(low, high)`` time-constant interval in hours, for this window.

        Growth-coupled states are evaluated at the window's own growth band, which is why a
        catalogue with one row per state can serve two vessels whose irreducible sets differ.

        Raises:
            TimescaleUnmeasured: if the time constant is refused.
        """
        if isinstance(self.tau_h, RefusedTimescale):
            float(self.tau_h)
        if self.tau_growth_multiple is not None:
            m = self.tau_growth_multiple
            return (m / window.growth_rate_high_per_h, m / window.growth_rate_low_per_h)
        low, high = self.tau_h  # type: ignore[misc]
        return (float(low), float(high))


@dataclass(frozen=True)
class MechState:
    """An ordered, named state vector: the bridge between mechanism and ``solve_ivp``.

    Holds the variables in integration order and converts between a name-keyed mapping, which
    is how a right-hand side should be written, and the flat array the solver wants. Nothing
    here integrates; :mod:`ystwin.mech.integrate` does that.
    """

    variables: tuple[StateVar, ...]

    def __post_init__(self) -> None:
        names = [v.name for v in self.variables]
        duplicates = sorted({n for n in names if names.count(n) > 1})
        if duplicates:
            raise ValueError(f"duplicate state names: {duplicates}")

    def __len__(self) -> int:
        return len(self.variables)

    def __iter__(self) -> Iterator[StateVar]:
        return iter(self.variables)

    def __getitem__(self, name: str) -> StateVar:
        for v in self.variables:
            if v.name == name:
                return v
        raise KeyError(f"no state named {name!r}; have {list(self.names)}")

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(v.name for v in self.variables)

    @property
    def units(self) -> tuple[str, ...]:
        return tuple(v.units for v in self.variables)

    def index(self, name: str) -> int:
        try:
            return self.names.index(name)
        except ValueError:
            raise KeyError(f"no state named {name!r}; have {list(self.names)}") from None

    def pack(self, values: Mapping[str, float]) -> np.ndarray:
        """Name-keyed mapping -> flat array in integration order. Missing names are an error."""
        missing = sorted(set(self.names) - set(values))
        if missing:
            raise KeyError(f"no initial value for {missing}")
        extra = sorted(set(values) - set(self.names))
        if extra:
            raise KeyError(f"{extra} are not states of this system; have {list(self.names)}")
        return np.array([float(values[n]) for n in self.names], dtype=float)

    def unpack(self, vector: Sequence[float]) -> dict[str, float]:
        """Flat array -> name-keyed mapping. The inverse of :meth:`pack`."""
        arr = np.asarray(vector, dtype=float)
        if arr.shape[-1] != len(self):
            raise ValueError(f"expected {len(self)} states, got {arr.shape[-1]}")
        return {n: float(arr[..., i]) for i, n in enumerate(self.names)}

    def subset(self, names: Sequence[str]) -> "MechState":
        """The same system restricted to ``names``, in this system's order."""
        wanted = set(names)
        unknown = sorted(wanted - set(self.names))
        if unknown:
            raise KeyError(f"no state named {unknown}; have {list(self.names)}")
        return MechState(tuple(v for v in self.variables if v.name in wanted))


_PLASMID_LOSS_PER_GENERATION = 0.05
_LN2 = math.log(2.0)

TIMESCALE_CATALOGUE = MechState((
    # -- Transport and phosphorelay: the fast band that sets the stiffness ceiling ----------
    StateVar(
        name="ssk1_p", units="fraction phosphorylated",
        tau_h=(1.0 / 160.0 * _S, 1.0 / 160.0 * _S),
        source="MEASURED. Janiak-Spens, Cook & West 2005, Biochemistry 44:377-86, PMID 15628880: "
               "Ypd1 -> Ssk1-R2 phosphotransfer k = 160 /s, tau = 1/160 s = 6.25 ms",
        sweep_axis=True,
        note="The fastest resolved process in the architecture; it is what makes RK45 need "
             "7.7e5 steps for a 4 h read. No reverse transfer from Ssk1-R2~P was detected.",
    ),
    StateVar(
        name="ypd1_p", units="fraction phosphorylated",
        tau_h=(1.0 / 29.0 * _S, 1.0 / 29.0 * _S),
        source="MEASURED. Janiak-Spens 2005, PMID 15628880: SLN1-R1~P -> Ypd1, k = 29 /s, "
               "tau = 34.5 ms",
        sweep_axis=True,
    ),
    StateVar(
        name="acid_cytosolic", units="mM undissociated acid",
        tau_h=(10.0e-3 * _S, 10.0e-3 * _S),
        source="DERIVED from MEASURED-not-in-yeast. tau = r/(3P) = 1/(P * A/V) with "
               "P(acetic) = 6.9e-3 cm/s (BNID 110821, Walter & Gutknecht, egg-PC/decane "
               "bilayer) and A/V = 1.39e4 /cm from a 42 um^3 cell (BNID 100427): tau = 10 ms",
        sweep_axis=True,
        note="The permeability is an ARTIFICIAL BILAYER, not a yeast membrane, and ignores the "
             "cell wall. Order of magnitude only; Gabba 2020 measures the yeast membrane 192-707x "
             "less permeable than a model bilayer.",
    ),
    StateVar(
        name="cell_volume", units="fraction of isosmotic volume",
        tau_h=(3.0 * _S, 3.0 * _S),
        source="BOUNDED by MEASURED. Braz J Microbiol 52:895-903 (2021), PMID 33476034: a 60% "
               "volume loss completes in < 8 s and reverses ~16 s after return to isosmotic "
               "medium, so tau ~ 3 s",
        sweep_axis=True,
        note="A completion time, not a fitted exponential constant.",
    ),
    StateVar(
        name="atp", units="mM",
        tau_h=(0.6 * _S, 0.8 * _S),
        source="DERIVED. Pool 2.5 mM (BNID 115810) over a fermentative ATP flux of "
               "~30 mmol/gDCW/h in 1.7-2.0 mL/gDCW (BNIDs 112806, 111353) = 4.2-4.9 mM/s",
        constrained_by="the ATP biosensor plates resolved by paths.atp_sensor_plates",
        note="The q_glucose in the numerator is a textbook figure, not a BNID. The ~1 s answer "
             "survives a 3x error in it; the exact value is not claimed.",
    ),
    # -- Sensors and transcription factors --------------------------------------------------
    StateVar(
        name="hsf1_free", units="fraction free of Hsp70",
        tau_h=(1.0 / 2.783 * _MIN, 1.0 / 2.783 * _MIN),
        source="FITTED SOURCE PRIOR. Zheng et al. 2016, eLife 5:e18638, PMID 27831465; "
               "Table 1 DOI 10.7554/eLife.18638.021 and parameter-screen Methods: "
               "k2 = 2.783 /min was selected by fitting/screening a model in arbitrary "
               "abundance units, not measured as an intrinsic binding off-rate. "
               "tau = 1/k2 = 21.6 s is a conditional model timescale, not a measured "
               "Hsf1 activation time. Verified primary XML PMC5127643; locator and SHA256 "
               "are recorded in data/parameter_evidence.json, sources.zheng2016.",
        sweep_axis=True,
        note="The switch only. Zheng's phosphorylation gain peaks at 20 min and the feedback "
             "closes at 60 min -- a 167x split on one input, uncoupled by their own account.",
    ),
    StateVar(
        name="yap1_nuclear", units="fraction nuclear",
        tau_h=(120.0 * _S, 120.0 * _S),
        source="MEASURED. Goulev et al. 2017, eLife 6:e23971, PMID 28418333, Fig 1-fig supp 1C, "
               "30 s sampling: Yap1 nuclear relocation ~120 s",
        sweep_axis=True,
    ),
    StateVar(
        name="crz1_nuclear", units="fraction nuclear",
        tau_h=(2.0 * _MIN, 2.0 * _MIN),
        encoding=Encoding.FREQUENCY,
        source="MEASURED. Cai, Dalal & Elowitz 2008, Nature 455:485-90, PMID 18818649: bursts "
               "last ~2 min and calcium sets the FREQUENCY, not the duration",
        sweep_axis=True,
        note="The state the ratio test gets wrong. T/tau = 124 on a plate read and reducing it "
             "is still wrong: frequency modulation is what makes targets scale proportionally "
             "and independently of promoter characteristics.",
    ),
    StateVar(
        name="msn2_nuclear", units="fraction nuclear",
        tau_h=(5.0 * _MIN, 5.0 * _MIN),
        encoding=Encoding.FREQUENCY,
        source="MEASURED. Hao & O'Shea 2012, Nat Struct Mol Biol 19:31-9, PMID 22179789: 5 min "
               "bursts, frequency-modulated under glucose limitation",
        sweep_axis=True,
        note="Same TF, two encodings: amplitude-modulated under oxidative stress. The STRESSOR "
             "picks the encoding, so a dose ladder in H2O2 and one in glucose limitation are "
             "not measuring the same thing about Msn2.",
    ),
    StateVar(
        name="hog1_p", units="fraction dually phosphorylated",
        tau_h=(5.0 * _MIN, 15.0 * _MIN),
        source="BOUNDED by MEASURED. Mettetal, Muzzey, Gomez-Uribe & van Oudenaarden 2008, "
               "Science 319:482-4, PMID 18218902: the protein-synthesis-independent feedback "
               "starts recovering within 5 min and is finished by 15 min",
        sweep_axis=True,
        note="A second, slower feedback through gene expression appears only after large shocks, "
             "so the effective constant is dose-dependent.",
    ),
    # -- Slow chemistry, transcripts and the reporter's own two clocks -----------------------
    StateVar(
        name="permeability_h2o2", units="cm/s",
        tau_h=(15.0 * _MIN, 15.0 * _MIN),
        source="BOUNDED by MEASURED. Branco, Marinho, Cyrne & Antunes 2004, JBC 279:6501, "
               "PMID 14645222: plasma-membrane permeability to H2O2 falls 40% by 15 min and "
               "2-fold on full adaptation",
        sweep_axis=True,
        note="Transport is the fastest process everywhere else in this table. Here the "
             "COEFFICIENT is itself a slow state, so eliminating transport as quasi-steady is "
             "right for the concentration and silently deletes the adaptation.",
    ),
    StateVar(
        name="mrna", units="transcripts/cell",
        tau_h=(3.6 * _MIN / _LN2, 32.0 * _MIN / _LN2),
        source="MEASURED, two methods that disagree 9x. Chan, Mugler, Heinrich, Vallotton & Weis "
               "2018, eLife 7:e32536, PMID 30192227, 4tU non-invasive: median t_half 3.6 min "
               "(tau 5.2 min). Transcription-shutoff: typical t_half ~32 min (BNID 114181, "
               "tau 46.2 min). Chan states shutoff gives up to 26x longer half-lives",
        constrained_by="RT-qPCR fold-change (ystwin/qpcr.py)",
        note="Shalem 2008 (PMID 18854817) and Castells-Roca 2011 (PMID 21364882): induced "
             "transcripts are DESTABILISED and repressed ones stabilised, their decay-rate ratio "
             "moving 2.7x inside 10 min, so one global k_deg has the wrong sign for one class.",
    ),
    StateVar(
        name="reporter_immature", units="a.u./cell",
        tau_h=(10.4 * _MIN / _LN2, 39.0 * _MIN / _LN2),
        source="MEASURED, two yeast papers that disagree 3.8x. Guerra, Vuillemenot, Rae, "
               "Ladyhina & Milias-Argeitis 2022, ACS Synth Biol 11:1129-41, PMID 35180343: "
               "mCitrine t_half 10.4 min at 30 C by optogenetic induction (tau 15.0 min). "
               "Gordon et al. 2007, Nat Methods 4:175-81, PMID 17237792: YFP 39 +/- 7 min by "
               "cycloheximide chase (tau 56.3 min)",
        constrained_by="mCitrine RFU rise on the real plates -- the channel exists, the number "
                       "does not (reporter.ReporterKinetics.k_mat defaults to None)",
        note="The conflict is methodological: Guerra's objection to a cycloheximide chase is "
             "that translation inhibitors perturb pH, which is exactly the artefact "
             "generator/redox.py already measures as larger than the signal it contaminates.",
    ),
    StateVar(
        name="trehalose", units="g/g protein",
        tau_h=(60.0 * _MIN, 60.0 * _MIN),
        source="BOUNDED by MEASURED. Hottiger, Schmutz & Wiemken 1987, J Bacteriol 169:5518-22, "
               "PMID 2960663: 0.01 -> 1 g/g protein within 1 h on a 27 -> 40 C shift",
        sweep_axis=True,
        note="Neutral trehalase rises 3x alongside a 6x rise in Tps, so the pool is a "
             "fast-turnover buffer maintained by a futile cycle, not a slow store.",
    ),
    StateVar(
        name="glutathione", units="mM GSH",
        tau_h=(90.0 * _MIN / _LN2, 90.0 * _MIN / _LN2),
        source="MEASURED. BNID 113864: GSH pool t_half ~90 min, tau = 130 min",
        constrained_by="Grx1-roGFP2 ratio, subject to the pH confound generator/redox.py "
                       "measures as 60 mV against a 40-50 mV real signal",
    ),
    StateVar(
        name="upr_output", units="fraction HAC1 spliced",
        tau_h=(2.0, 4.0),
        source="MEASURED. Pincus et al. 2010, PLoS Biol 8(7):e1000415, PMID 20625545: HAC1 "
               "splicing deactivates within 2 h at 1.5 mM DTT and within 4 h at 2.2 mM",
        constrained_by="RT-qPCR of HAC1 on the UPRE1 strain (ystwin/qpcr.py)",
        note="A measured DEACTIVATION time, not a fitted relaxation constant, and it is "
             "dose-dependent: at 5 mM the response is maximal and sustained, so the apparent "
             "dose-response over a fixed read is a duration curve in an amplitude curve's "
             "clothes.",
    ),
    StateVar(
        name="cross_protection", units="fold protection",
        tau_h=(4.0 * 1.5, 5.0 * (140.0 / 60.0)),
        source="DERIVED from MEASURED. Guan, Haroon, Bravo, Will & Gasch 2012, Genetics "
               "192:495-505, PMID 22851651: H2O2 resistance persists 4-5 generations after the "
               "salt pre-treatment is removed, carried by long-lived Ctt1p rather than by new "
               "synthesis. At 90-140 min per generation (BNID 108255) that is 6.0-11.7 h",
        sweep_axis=True,
        note="Longer than a plate read, so it sets the INITIAL CONDITION rather than moving "
             "inside the window. A model that starts every well from a defined t = 0 is "
             "asserting the pre-culture carried no history.",
    ),
    # -- Growth-coupled: tau = multiple / mu, so the window supplies the number ---------------
    StateVar(
        name="reporter_mature", units="a.u./cell",
        tau_growth_multiple=1.0,
        source="DERIVED. Loss rate is mu + k_deg and D2 measured k_deg ~ 0 for mCitrine on the "
               "real plates, so dilution is the whole of it: tau = 1/mu. Every dynamic-SILAC "
               "median protein half-life (2.0-8.8 h, BNID 115093) is slower than dilution at "
               "any growth rate this project runs",
        constrained_by="mCitrine RFU on the real plates",
        note="This is the observable. Its time-integral over the window is what a stable-FP "
             "plate read records.",
    ),
    StateVar(
        name="biomass", units="gDCW/L",
        tau_growth_multiple=1.0,
        encoding=Encoding.INTEGRATING,
        source="MEASURED per window. mu is fitted from the biomass series, never asserted from "
               "a pump: 0.245-0.363 /h on the real plates, 0.101 /h at the fed-batch setpoint",
        constrained_by="OD600 on the real plates",
    ),
    StateVar(
        name="product", units="mmol/gDCW",
        tau_growth_multiple=1.0,
        encoding=Encoding.INTEGRATING,
        source="DEFINITION, not a literature number. A terminal intracellular pool is held down "
               "by washout alone, X = v_in/mu (pathway/solve.py), so its timescale is 1/mu",
        constrained_by="HPLC titre. A plate absorbance is NOT specific -- it sums lycopene with "
                       "beta-carotene and reads above the ceiling",
    ),
    StateVar(
        name="plasmid_bearing", units="fraction of cells",
        tau_growth_multiple=1.0 / -math.log1p(-_PLASMID_LOSS_PER_GENERATION),
        source="DERIVED from a REVIEW-LEVEL figure. Engineered 2-micron plasmids are lost at up "
               "to ~5%/generation without selection (Futcher & Cox 1984, J Bacteriol 157:283-90, "
               "PMID 6361000, established the negative copy-number/instability correlation). "
               "tau = 1/(mu * -ln(1-0.05)) = 19.50/mu",
        sweep_axis=True,
        note="The 5%/generation could not be tied to one primary table and is an order of "
             "magnitude. It is the ONLY process that flips between the two windows: 4 h is "
             "1.0-1.5 generations, 120 h at the setpoint is 17.5.",
    ),
    StateVar(
        name="generations", units="generations",
        tau_growth_multiple=_LN2,
        encoding=Encoding.INTEGRATING,
        source="DEFINITION, not a literature number. dg/dt = mu/ln2, so one generation is ln2/mu",
        constrained_by="OD600 on the real plates",
    ),
    # -- Refused: named so nobody re-derives them, and raising rather than defaulting ---------
    StateVar(
        name="ire1_cluster",
        units="fraction of Ire1 in clusters",
        tau_h=RefusedTimescale(
            quantity="Ire1 oligomerisation and cluster dissolution",
            reason="clusters are imaged at a fixed 45 min timepoint (10 mM DTT) and splicing "
                   "assayed at a fixed 2 h (2 mM DTT); no formation or dissolution rate is "
                   "reported anywhere in the paper",
            source="van Anken et al. 2014, eLife 3:e05031",
        ),
        source="REFUSED -- see the RefusedTimescale",
        sweep_axis=True,
    ),
    StateVar(
        name="slt2_p",
        units="fraction dually phosphorylated",
        tau_h=RefusedTimescale(
            quantity="Slt2/Mpk1 activation by heat",
            reason="the paper reports 'strongly activated by mild heat shock, sustained during "
                   "growth at high temperature'; the time course is in figures and no minutes "
                   "are stated in text",
            source="Kamada, Jung, Piotrowski & Levin 1995, Genes Dev 9:1559-71",
        ),
        source="REFUSED -- see the RefusedTimescale",
        sweep_axis=True,
    ),
    StateVar(
        name="snf1_p",
        units="fraction active",
        tau_h=RefusedTimescale(
            quantity="Snf1/AMPK activation",
            reason="no activation time constant found in any source. The nearest number is not "
                   "Snf1 at all: Gorner 2002's '> 2 min' is PKA dephosphorylating the Msn2 NLS",
            source="Gorner et al. 2002, EMBO J 21:135-44, PMID 11782433 (the WRONG kinase)",
        ),
        source="REFUSED -- see the RefusedTimescale",
        sweep_axis=True,
    ),
))
"""Every state the target architecture names whose time constant is measured, derived, bounded
or refused, in one table so the audit has something real to run on.

Twenty-six rows spanning 6.25 ms to 193 h -- a factor of 1.1e8. That spread is the whole reason
:mod:`ystwin.mech.integrate` pins a stiff solver and applies a window-aware reduction: an
explicit method chooses its step from the fastest row and the observable only sees the slowest.
"""


def catalogue_span_hours(states: MechState = TIMESCALE_CATALOGUE) -> tuple[float, float]:
    """Fastest and slowest time constant across both declared windows, hours.

    Computed rather than asserted, because the stiffness ratio is the argument for the whole
    integrator and an asserted 1e8 is exactly the kind of number this repository has been
    burned by. Refused rows are skipped: an unmeasured process cannot widen a measured span.
    """
    taus: list[float] = []
    for v in states:
        if v.timescale_is_refused:
            continue
        for window in (PLATE_READ_4H, FEDBATCH_5D):
            taus.extend(v.taus_in(window))
    if not taus:
        raise TimescaleUnmeasured("no state in this system has a measured time constant")
    return (min(taus), max(taus))
