"""Where the *fluorescence* detector stops responding linearly, from two gains on one plate.

Why this can be measured when `calib/od_linearity.py` cannot
------------------------------------------------------------
`od_linearity.py` needs a dilution series nobody has run: to ask whether a reader is
linear you need to know what it *should* have read, and only a dilution factor supplies
that. The fluorescence channel has the same answer for free and nobody had used it.
Seventeen of the 22 Gen5 ``.xpt`` files read the same fluorophore **twice per step at two
photomultiplier gains** -- ``mCitrine:480,530`` and ``mCitrine:480,530[2]``, from
protocols named ``...mCitrine_OD600+75100.prt``. Two gains on one well are a dilution
series in the light: the same photons, scaled by a known factor. If the reported scale is
linear the two channels differ by a constant multiple at every signal level, and where one
of them compresses the multiple moves.

What that measurement says
--------------------------
Fitting ``high = ratio * low + offset`` on the mid-range and asking how the residual
behaves at every level, over 96,459 paired readings from the ten pairs whose protocol name
declares gain 75 then 100 (``scripts/score_gain_linearity.py``):

* **The multiple is 8.15** (per-plate 8.079 to 8.217, a 0.5% spread), and it does not
  drift with signal. The median of ``measured / predicted`` stays inside **0.999 to
  1.001** from the fitting window at 2,000 RFU to 99,998 -- the largest reading in the
  entire corpus.
* **There is no roll-off.** Nothing compresses, softens or bends. The channel is linear
  right up to the last reading it reports.
* **It fails by refusing, not by lying.** Above the top of scale Gen5 stops reporting a
  number and writes :data:`RFU_OVERFLOW` in its place. 11,149 of the corpus's 303,552
  fluorescence readings are that sentinel, across 8 channels, and wherever a lower-gain
  partner survived it predicts at least 98.5% of :data:`RFU_CEILING` for the flagged
  reading -- so the flag is the top of scale and not an early surrender.

So the thing this was meant to find -- a silently wrong reading in the region where the
detector is compressing -- **does not exist for fluorescence on this instrument**. The
hazard is the opposite one: ``plate/gen5.py`` returns the sentinel as an ordinary float, so
a caller reading a saturated channel gets ``-99999.0`` where it expects RFU. Use
:func:`is_overflow` before arithmetic on any ``[2]``/``[3]`` channel.

The gain law, as a check on the above
-------------------------------------
``20260707_oxidative_stress_preliminary.xpt`` reads mCitrine three times in one endpoint
step, at 50, 75 and 100 (``endpoint_mCitrine_shake+OD600+5075100.prt``). The two
independent steps give ratios 18.902 and 8.079, and

    ln(18.902) / ln(75/50) = 7.249       ln(8.079) / ln(100/75) = 7.262

which is :data:`RFU_GAIN_EXPONENT`. Two gain steps that share no readings agree to 0.2%
that RFU goes as gain to a fixed power -- so the ratios above are the gain relation and
not a coincidence of one plate.

**This does not set** ``OpticalQualityGate.od_linear_max``
----------------------------------------------------------
It is a different detector. Fluorescence is read by a photomultiplier looking at emitted
light through the 530 nm filter; absorbance is read by a photodiode measuring transmitted
600 nm light against a reference, and its nonlinearity is stray light and pathlength, not
anode current. A number measured on one says nothing about the other, and moving
``od_linear_max`` on this evidence would be exactly the kind of laundering this repository
exists to refuse. ``od_linear_max`` stays at its placeholder 1.0 and still needs the
dilution series `od_linearity.py` describes.

The absorbance channel does carry one *negative* result, from the same files: OD600 tops
out at 1.958 over 20 plates, whose maxima take 20 distinct values, and the most readings
any one plate has at its own maximum is **2** -- against 11,149 fluorescence readings
pinned to one number. The reader is not clipping absorbance anywhere these cultures
reached, so the corpus contains no top-of-scale for OD to find. That bounds nothing about
linearity, which is the point: saturation of the *number* and departure from
*proportionality* are different questions, and only the second one matters for OD.

Caveats that travel with the ratio
----------------------------------
* The two reads are **52 seconds apart**, not simultaneous -- Gen5 sweeps the whole plate
  at one gain and then again at the next, and the offset is exactly 52 s at every timepoint
  of all three AFL-debugging plates. At mu = 0.3 /h that is 0.43% of growth between them,
  which is inside the 0.5% plate-to-plate spread of the ratio and cannot be separated from
  it here. It is a fixed multiplicative bias on the ratio and does not touch the linearity
  conclusion, which is about the *shape* of the response.
* Below the fitting window the test has no power: at 2,000 RFU on the high-gain channel
  the low-gain partner is near 245, and by 500 it is near 60, where one count of the
  instrument's integer quantisation is 1.6%. :func:`fit_gain_pair` measures that scatter
  and refuses to call it saturation.
* **The two channels are checked over different ranges, and the asymmetry is the whole
  caveat.** Gain 75 is checked against gain 100 from 40 to **12,677** RFU, above which its
  partner has overflowed and can check nothing. Gain 100 is checked against gain 75 from
  383 to **99,998**. The largest gain-75 reading in the corpus is 48,244 -- 3.8x past
  where gain 75 was directly checked -- so calling it unsaturated rests on the reported
  RFU scale having been linear at 48,244 *on the other channel of the same detector*,
  which is an inference and is weaker than the direct check.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..plate import gen5

__all__ = [
    "BIN_EDGES",
    "GainPairFit",
    "MCITRINE_GAIN_RATIO_75_TO_100",
    "RFU_CEILING",
    "RFU_GAIN_EXPONENT",
    "RFU_LINEAR_MAX_MEASURED",
    "RFU_OVERFLOW",
    "fit_gain_pair",
    "is_overflow",
]

RFU_OVERFLOW = gen5.OVERFLOW
"""What Gen5 writes in a fluorescence well whose reading passed the top of scale.

Measured, not documented by the vendor: 11,149 of the 303,552 fluorescence readings in the
22-file corpus carry exactly this value, spread over 8 channels of 7 files. **No reading
anywhere lies between 99,998 and this sentinel**, which is what makes it a flag rather than
a reading -- a detector running out of range would report values approaching its limit, and
this one reports none.

Six of those channels are the higher-gain member of a pair whose lower-gain partner
survived, and there the flag lands where it should: the lowest value any flagged reading's
partner predicts is 98,492, **98.5% of** :data:`RFU_CEILING`, so the instrument is not
giving up early. The other two are both channels of ``24h_30min-interval_mCitrine.xpt``
(``mCitrine:510,530``), which is 100% sentinel at both gains and therefore says nothing
about the detector -- only that the read was misconfigured.

Defined in `plate/gen5.py` as :data:`~ystwin.plate.gen5.OVERFLOW`, because it is a property
of the file format; the evidence for it is here because this is where it was measured.
`Gen5Channel.values` still returns it as a float -- that reader does not repair files -- but
`Gen5Channel.frame` refuses a channel containing it unless the caller says what to do, and
anything working from `values` must mask it first (see :func:`is_overflow`) or it will
average a reporter trace with -99999 in it and get a number that looks like a measurement.
"""

RFU_CEILING = 100_000.0
"""Top of the reported fluorescence scale, in RFU, at any gain.

The largest non-sentinel reading anywhere in the 22-file corpus is 99,998, and in each of
the six overflowing channels that also carry real readings the largest is 99,998 / 99,994 /
99,993 / 99,983 / 99,965 / 99,957. Readings stop at a round decimal number and jump
straight to :data:`RFU_OVERFLOW`: a hard limit on the reported value, with no compression
approaching it (:data:`RFU_LINEAR_MAX_MEASURED`).

It is stated per *channel*, not per plate: at gain 75 the same ceiling takes 8.15x more
incident light to reach than at gain 100, which is why no gain-75 reading in this corpus
comes near it -- the largest is 48,244, 48% of the ceiling.
"""

RFU_LINEAR_MAX_MEASURED = 99_998.0
"""Highest fluorescence reading shown to be on the linear response, in RFU.

Not a fitted breakpoint -- there is no breakpoint. This is the largest reading in the
corpus, and the bin containing it holds 325 pairs with a median ``measured / predicted``
of 0.9991 against a fit anchored a decade and a half lower. The value is the *edge of the
evidence*, exactly as `od_linearity.fit_linear_range` refuses to extrapolate past the
densest well actually run.

Which channel this was measured on matters, so it is stated: the gain-100 channel is
checked against gain 75 over 383 to 99,998 RFU, and the gain-75 channel against gain 100
over 40 to only 12,677, above which its partner has overflowed and can check nothing. See
the last bullet of the module docstring -- a gain-75 reading above 12,677 is covered by
inference from the other channel and not by a direct check.
"""

MCITRINE_GAIN_RATIO_75_TO_100 = 8.15
"""RFU at gain 100 divided by RFU at gain 75, mCitrine 480/530.

Per-plate affine fits ``high = ratio * low + offset`` on the pairs with the high channel
between 2,000 and 20,000 RFU, over the ten ``.xpt`` pairs whose protocol name declares
these two gains: mean 8.149, range 8.079 to 8.217, spread 0.5%. **The spread is the number
to carry, not the third digit of the mean** -- per-plate ratios and offsets are columns of
``outputs/gain_linearity.csv`` for exactly that reason.

Reproduced by ``scripts/score_gain_linearity.py``. Nothing in this repository converts
between gains; this is here so that a future analysis reaching for the ``[2]`` channel has
the factor and its uncertainty rather than deriving one from whatever plate is at hand.
"""

RFU_GAIN_EXPONENT = 7.26
"""``RFU ~ gain ** k``. Fitted k, from the one plate read at three gains.

``20260707_oxidative_stress_preliminary.xpt``, protocol
``endpoint_mCitrine_shake+OD600+5075100.prt``: the 50 -> 75 step gives 18.902 and the
75 -> 100 step gives 8.079, and ``ln(18.902)/ln(1.5) = 7.249`` against
``ln(8.079)/ln(4/3) = 7.262``. One plate, one fluorophore, two gain steps that share no
readings. Recorded because it is what makes the pair ratios interpretable as a property of
the detector rather than of a plate; it is not used to convert anything.
"""

BIN_EDGES = (0.0, 500.0, 1_000.0, 2_000.0, 5_000.0, 10_000.0, 20_000.0, 30_000.0,
             40_000.0, 50_000.0, 60_000.0, 70_000.0, 80_000.0, 90_000.0, 95_000.0,
             98_000.0, 100_000.0)
"""Signal levels the departure from proportionality is reported at, high-gain RFU.

Coarse at the bottom where the check has no power and fine near
:data:`RFU_CEILING` where a roll-off would appear if there were one. Fixed rather than
quantile-derived so that two plates' tables have the same rows and can be read side by
side.
"""

_MIN_BIN_POINTS = 30
"""Fewest pairs in a bin before its median departure is allowed to fail a channel.

A bin of five readings has a median that moves several percent on one outlier, and calling
that saturation would put a breakpoint wherever the last few wells happened to land. A bin
thinner than this is measured and reported and cannot conclude anything -- which is also
why :attr:`GainPairFit.linear_max` stops at the highest bin that reached this count, rather
than at the highest reading on the plate.
"""

_MIN_ANCHOR_POINTS = 10
"""Fewest pairs in the anchor window before the multiple may be fitted at all.

Lower than :data:`_MIN_BIN_POINTS` on purpose: the anchor is a straight line fitted to
points spanning a decade of signal, which ten points determine, while a bin's job is to
resolve a one-percent shift in a median, which ten points do not.
"""


def is_overflow(values) -> np.ndarray:
    """Which readings are :data:`RFU_OVERFLOW` rather than measurements.

    Args:
        values: Fluorescence readings, any shape.

    Returns:
        Boolean array, ``True`` where the instrument reported an overflow.
    """
    return np.asarray(values, dtype=float) <= RFU_OVERFLOW + 1.0


@dataclass(frozen=True)
class GainPairFit:
    """One channel pair on one plate: the multiple, and where it stops holding."""

    ratio: float
    """Fitted ``high / low`` slope, from the anchor window."""

    offset: float
    """Fitted intercept, RFU on the high channel. Small compared with the anchor window."""

    linear_max: float
    """Highest high-channel reading inside a bin that held, RFU.

    The edge of the *evidence*, not of the plate: a bin too thin to resolve a departure
    cannot extend it, so a plate whose top decade holds nine readings reports a linear
    range ending where the last well-populated bin ended. The same doctrine as
    `od_linearity.fit_linear_range`, which never extrapolates past the densest well
    actually run.
    """

    saturation_observed: bool
    """``True`` only if a bin below :data:`RFU_CEILING` departed by more than ``tolerance``."""

    worst_departure: float
    """Largest ``|median(measured / predicted) - 1|`` over the bins that were tested."""

    scatter: float
    """Half the 5th-to-95th-percentile spread of ``measured / predicted`` in the anchor.

    The floor under ``tolerance``: a departure smaller than this is not distinguishable
    from the reading noise of the pair.
    """

    tolerance: float
    n_pairs: int
    """Pairs where both channels reported a number."""

    n_overflow: int
    """Readings where the high channel reported :data:`RFU_OVERFLOW` instead."""

    def summary(self) -> str:
        edge = ("compression seen" if self.saturation_observed
                else f"no compression up to {self.linear_max:,.0f} RFU")
        return (f"x{self.ratio:.3f} +/- {self.scatter:.1%}, {edge} "
                f"({self.n_pairs:,} pairs, {self.n_overflow:,} overflowed)")


def fit_gain_pair(low_gain_rfu, high_gain_rfu, tolerance: float = 0.05,
                  anchor: tuple[float, float] = (2_000.0, 20_000.0),
                  bin_edges: tuple[float, ...] = BIN_EDGES) -> GainPairFit:
    """Where two simultaneous reads of one plate stop being a fixed multiple of each other.

    Fits ``high = ratio * low + offset`` inside ``anchor`` -- chosen high enough that
    integer quantisation on the low channel is negligible and far enough below
    :data:`RFU_CEILING` that no compression could be inside it -- then walks the bins
    upward from ``anchor[0]`` asking whether the median of ``measured / predicted`` is
    still 1. Bins below ``anchor[0]`` are measured and reported but cannot fail the
    channel: down there the low-gain partner is a two-digit integer and its quantisation,
    not the detector, sets the departure.

    Args:
        low_gain_rfu: Readings from the lower-gain channel, any shape.
        high_gain_rfu: The same wells and timepoints at the higher gain.
        tolerance: Largest median departure from the anchor line still called linear.
        anchor: High-channel RFU window the multiple is fitted in.
        bin_edges: Ascending high-channel RFU edges to report departures at.

    Returns:
        A :class:`GainPairFit`. Overflowed readings are excluded from every fit and
        counted separately, never treated as the number -99999.

    Raises:
        ValueError: if the two arrays differ in shape, if the anchor is not ascending, or
            if the anchor window holds fewer than ``_MIN_ANCHOR_POINTS`` pairs -- a plate
            that never reached the signal level the multiple is fitted at cannot answer
            this question, and a fit through three points would not say so.
    """
    low = np.asarray(low_gain_rfu, dtype=float).ravel()
    high = np.asarray(high_gain_rfu, dtype=float).ravel()
    if low.shape != high.shape:
        raise ValueError("the two gain channels must hold the same readings")
    if not anchor[0] < anchor[1]:
        raise ValueError(f"anchor window {anchor} is not ascending")

    overflow = is_overflow(high) | is_overflow(low)
    n_overflow = int(overflow.sum())
    low, high = low[~overflow], high[~overflow]

    inside = (high >= anchor[0]) & (high <= anchor[1])
    if int(inside.sum()) < _MIN_ANCHOR_POINTS:
        raise ValueError(
            f"only {int(inside.sum())} pairs inside the anchor window {anchor}; this plate "
            f"never reached the signal level the multiple is fitted at")
    design = np.vstack([low[inside], np.ones(int(inside.sum()))]).T
    (ratio, offset), *_ = np.linalg.lstsq(design, high[inside], rcond=None)

    predicted = ratio * low + offset
    usable = np.abs(predicted) > 1.0
    departure = np.ones_like(high)
    departure[usable] = high[usable] / predicted[usable]

    anchor_spread = departure[inside & usable]
    scatter = float(np.percentile(anchor_spread, 95) - np.percentile(anchor_spread, 5)) / 2.0

    # The linear range runs to the top of the highest bin that BOTH held a departure of
    # one and carried enough readings to have shown otherwise. It starts at the top of the
    # anchor, because that window was fitted and cannot fail its own fit.
    linear_max = float(high[inside].max())
    saturation, worst = False, 0.0
    for low_edge, high_edge in zip(bin_edges[:-1], bin_edges[1:]):
        if low_edge < anchor[0]:
            continue
        in_bin = usable & (high >= low_edge) & (high < high_edge)
        if int(in_bin.sum()) < _MIN_BIN_POINTS:
            continue
        off_by = abs(float(np.median(departure[in_bin])) - 1.0)
        worst = max(worst, off_by)
        if off_by > tolerance:
            saturation = True
            break
        linear_max = max(linear_max, float(high[in_bin].max()))

    return GainPairFit(
        ratio=float(ratio), offset=float(offset), linear_max=linear_max,
        saturation_observed=saturation, worst_departure=float(worst), scatter=scatter,
        tolerance=tolerance, n_pairs=int(low.size), n_overflow=n_overflow,
    )
