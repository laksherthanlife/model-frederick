"""The second gain nobody read: what it says about the fluorescence detector's linear range.

    YSTWIN_GEN5_XPT="<Experiments dir>" python3 scripts/score_gain_linearity.py

Writes ``outputs/gain_linearity.csv``. Needs the Gen5 ``.xpt`` files; refuses cleanly
without them.

**The channel that was already there.** Twelve of the 22 ``.xpt`` files read the same
fluorophore twice in one step at two photomultiplier gains -- ``mCitrine:480,530`` and
``mCitrine:480,530[2]``, from protocols named ``...mCitrine_OD600+75100.prt``. This whole
repository reads the first and has never touched the second. Two gains on one well are a
dilution series in the light: the same photons scaled by a known factor, which is exactly
the input `calib/od_linearity.py` says it needs and has never had.

**What it establishes, and what it refuses to.** See
``src/ystwin/calib/gain_linearity.py`` for the argument. In one line: the reported RFU
scale is linear over its whole range and then stops reporting, so the failure mode is an
overflow flag rather than a compressed number -- and it is the *fluorescence* detector, so
it does not move ``OpticalQualityGate.od_linear_max``, which is a photodiode.

**Why not only the AFL plates.** The three AFL-debugging runs are the cleanest pairs here
and they never exceed 8,498 RFU at gain 75. They can show the multiple is constant; they
cannot show where it stops, because nothing on them gets near the top of scale. Six other
plates do overflow, and they are what locates the ceiling. Which plates are the AFL runs is
a column, so the two questions stay separable in the table.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.calib.gain_linearity import (BIN_EDGES, RFU_CEILING, RFU_OVERFLOW,
                                         fit_gain_pair, is_overflow)
from ystwin.plate import gen5

#: ``.xpt`` -> the gains its protocol name declares, in read order.
#:
#: Taken from the protocol string Gen5 stores in the file, not guessed from the readings:
#: ``24h-10min_mCitrine_OD600+75100.prt`` says the fluorophore was read at 75 then 100, and
#: ``endpoint_mCitrine_shake+OD600+5075100.prt`` says 50, 75 then 100. Files whose protocol
#: name does not carry the gains are scored too, under ``GAIN_UNSTATED`` -- the ratio is
#: still measurable, only its interpretation as a gain step is not.
DECLARED_GAINS = {
    "20260705_AFL&control_1st.xpt": (75, 100),
    "20260707_oxidative_stress_preliminary.xpt": (50, 75, 100),
    "20260708_ER_stress_1st.xpt": (75, 100),
    "20260712_AFL&control_2nd_AND_oxidative_stress_1st.xpt": (75, 100),
    "20260713_AFL&Control3_ER Stress1.xpt": (75, 100),
    "20260802_AFL-debugging_1st.xpt": (75, 100),
    "20260803_AFL-debugging_2nd.xpt": (75, 100),
    "20260804_AFL-debugging_3rd.xpt": (75, 100),
    "20260805_ATP-sensor.xpt": (75, 100),
    "20260808_ATP-sensors(ACS-ICL).xpt": (75, 100),
}

#: Written in the gain columns when the protocol name does not state them. Not zero and not
#: blank: a reader must be able to tell "read at gain 0" from "the file does not say".
GAIN_UNSTATED = -1

#: The auto feedback loop plates, which are what ``scripts/score_afl_circuit.py`` scores.
#: Listed here rather than imported so that this table keeps meaning the same thing if that
#: script's series change; the question asked of them is only "does the reporter channel
#: this analysis reads ever approach the ceiling".
AFL_PLATES = ("20260705_AFL&control_1st.xpt",
              "20260802_AFL-debugging_1st.xpt",
              "20260803_AFL-debugging_2nd.xpt",
              "20260804_AFL-debugging_3rd.xpt")

#: High-channel RFU window the multiple is fitted in. Above the quantisation floor -- one
#: instrument count on the low channel is then under 0.5% -- and five times below
#: :data:`RFU_CEILING`, so nothing inside it can already be compressed.
ANCHOR = (2_000.0, 20_000.0)

#: Largest median departure from the fitted line still called linear. Five percent is the
#: same figure `calib/od_linearity.fit_od_calibration` defaults to, and it is roughly four
#: times the pair scatter, so a bin failing it is not failing on noise.
TOLERANCE = 0.05


def _pairs(read: gen5.Gen5Read, low: str, high: str) -> tuple[np.ndarray, np.ndarray]:
    """The two channels flattened onto a common grid of readings.

    A run stopped early can leave the second channel one timepoint shorter than the first
    -- the sweep at the next gain had not finished -- so the shared prefix is used and the
    difference is never padded.
    """
    a, b = read.channel(low).values, read.channel(high).values
    n = min(a.shape[0], b.shape[0])
    return a[:n].ravel(), b[:n].ravel()


def gain_pairs(read: gen5.Gen5Read) -> list[tuple[gen5.Gen5Channel, gen5.Gen5Channel]]:
    """Consecutive channels that are the same fluorophore and optics read again.

    ``mCitrine:480,530`` then ``mCitrine:480,530[2]`` is a pair; ``OD600:600`` then
    ``mCitrine:480,530`` is not. Adjacency in file order is what makes the second read the
    next gain up rather than an unrelated step, which is also how the ``[n]`` suffix is
    assigned -- see ``plate/gen5.py``.
    """
    fluorescent = [c for c in read.channels if c.fluorophore != "OD600"]
    return [(lower, upper) for lower, upper in zip(fluorescent[:-1], fluorescent[1:])
            if (lower.fluorophore, lower.optics) == (upper.fluorophore, upper.optics)]


def score_plate(path: pathlib.Path) -> pd.DataFrame:
    """Every gain pair on one plate, as departure-from-proportionality by signal level.

    Returns:
        For a pair that could be fitted: one row per RFU bin with ``status='scored'``,
        plus a row at :data:`RFU_CEILING` holding the overflowed readings, which have no
        departure because they have no value. For a pair that could not: a single row
        whose ``status`` is the reason, so a table reader sees the refusal rather than an
        absence. Empty when the file holds no channel read twice.
    """
    read = gen5.read_xpt(path)
    fluorescent = [c for c in read.channels if c.fluorophore != "OD600"]
    gains = DECLARED_GAINS.get(path.name)
    rows = []
    for lower, upper in gain_pairs(read):
        low, high = _pairs(read, lower.name, upper.name)
        identity = {
            "plate": path.stem,
            "low_channel": lower.name, "high_channel": upper.name,
            "gain_low": gains[fluorescent.index(lower)] if gains else GAIN_UNSTATED,
            "gain_high": gains[fluorescent.index(upper)] if gains else GAIN_UNSTATED,
            "is_afl_plate": path.name in AFL_PLATES,
        }
        try:
            fit = fit_gain_pair(low, high, tolerance=TOLERANCE, anchor=ANCHOR)
        except ValueError as refusal:
            rows.append({
                **identity, "status": str(refusal),
                "plate_ratio": np.nan, "plate_offset": np.nan, "plate_scatter": np.nan,
                "plate_linear_max_rfu": np.nan, "plate_saturation_observed": False,
                "plate_n_pairs": int((~(is_overflow(low) | is_overflow(high))).sum()),
                "plate_n_overflow": int((is_overflow(low) | is_overflow(high)).sum()),
                "bin_low_rfu": np.nan, "bin_high_rfu": np.nan, "n_pairs": 0,
                "departure_median": np.nan, "departure_p05": np.nan,
                "departure_p95": np.nan, "high_max_rfu": np.nan,
            })
            continue

        valid = ~(is_overflow(low) | is_overflow(high))
        low, high = low[valid], high[valid]
        predicted = fit.ratio * low + fit.offset
        usable = np.abs(predicted) > 1.0
        departure = np.where(usable, high / np.where(usable, predicted, 1.0), np.nan)

        shared = {
            **identity, "status": "scored",
            "plate_ratio": fit.ratio, "plate_offset": fit.offset,
            "plate_scatter": fit.scatter,
            "plate_linear_max_rfu": fit.linear_max,
            "plate_saturation_observed": fit.saturation_observed,
            "plate_n_pairs": fit.n_pairs, "plate_n_overflow": fit.n_overflow,
        }
        for low_edge, high_edge in zip(BIN_EDGES[:-1], BIN_EDGES[1:]):
            in_bin = usable & (high >= low_edge) & (high < high_edge)
            n = int(in_bin.sum())
            rows.append({
                **shared, "bin_low_rfu": low_edge, "bin_high_rfu": high_edge,
                "n_pairs": n,
                "departure_median": float(np.median(departure[in_bin])) if n else np.nan,
                "departure_p05": float(np.percentile(departure[in_bin], 5)) if n else np.nan,
                "departure_p95": float(np.percentile(departure[in_bin], 95)) if n else np.nan,
                "high_max_rfu": float(high[in_bin].max()) if n else np.nan,
            })
        rows.append({
            **shared, "bin_low_rfu": RFU_CEILING, "bin_high_rfu": np.inf,
            "n_pairs": fit.n_overflow, "departure_median": np.nan,
            "departure_p05": np.nan, "departure_p95": np.nan,
            "high_max_rfu": RFU_OVERFLOW if fit.n_overflow else np.nan,
        })
    return pd.DataFrame(rows)


def overflow_census(directory: pathlib.Path) -> dict[str, float]:
    """Every overflowed fluorescence reading in the corpus, counted per channel.

    The per-pair count in the main table charges a timepoint once however many of its
    channels overflowed. This counts readings, which is the number that matters to a caller
    doing arithmetic on one channel, and it reaches files the pair analysis refuses.
    """
    flagged, readings, channels, highest = 0, 0, 0, -np.inf
    for path in sorted(directory.glob("*.xpt")):
        try:
            read = gen5.read_xpt(path)
        except gen5.Gen5FormatError:
            continue
        for channel in read.channels:
            if channel.fluorophore == "OD600":
                continue
            values = channel.values.ravel()
            over = is_overflow(values)
            channels += int(over.any())
            flagged += int(over.sum())
            readings += int(values.size)
            if (~over).any():
                highest = max(highest, float(values[~over].max()))
    return {"flagged": flagged, "readings": readings,
            "channels_with_an_overflow": channels, "highest_reading": highest}


def absorbance_top_of_scale(directory: pathlib.Path) -> dict[str, float]:
    """The negative result for OD600: does the absorbance channel clip anywhere?

    A reader at the top of its scale pins many wells to one number -- the fluorescence
    channel does exactly that, 11,149 times, at :data:`RFU_OVERFLOW`. So the probe is the
    number of readings equal to the plate's own maximum. One, on a plate of 13,920, is a
    reader still resolving; a plate whose maximum is shared by hundreds of readings is one
    that has stopped. The plate maxima themselves are also counted: a common ceiling would
    put every plate at the same number.

    This says nothing about *linearity*, which is a different question and the one
    ``od_linear_max`` is about.
    """
    highest, most_pinned, plates, maxima = 0.0, 0, 0, set()
    for path in sorted(directory.glob("*.xpt")):
        try:
            read = gen5.read_xpt(path)
        except gen5.Gen5FormatError:
            continue
        for channel in read.channels:
            if channel.fluorophore != "OD600":
                continue
            values = channel.values.ravel()
            values = values[values > 0]
            if values.size < 20:
                continue
            plates += 1
            top = float(values.max())
            highest = max(highest, top)
            maxima.add(round(top, 4))
            most_pinned = max(most_pinned, int((values == top).sum()))
    return {"max_od": highest, "plates": plates,
            "most_readings_at_a_plate_max": most_pinned,
            "distinct_plate_maxima": len(maxima)}


def main() -> int:
    directory = paths.gen5_xpt_dir()
    if directory is None:
        raise SystemExit(
            "the Gen5 .xpt files are not here. Set YSTWIN_GEN5_XPT to the Experiments "
            "directory. They are the raw instrument files and are not redistributable.")

    frames = []
    for path in sorted(directory.glob("*.xpt")):
        try:
            frame = score_plate(path)
        except gen5.Gen5FormatError as bad:
            print(f"  skipped {path.name}: {bad}")
            continue
        if not frame.empty:
            frames.append(frame)
    if not frames:
        raise SystemExit(f"no .xpt in {directory} carries a channel read at two gains")
    table = pd.concat(frames, ignore_index=True)

    pairs = table.drop_duplicates(["plate", "low_channel", "high_channel"])
    scored = pairs[pairs.status == "scored"]
    declared = scored[scored.gain_low == 75]

    print("Two gains on one plate: the fluorescence detector's linear range.\n")
    print("  Every pair the corpus holds, ratio fitted at "
          f"{ANCHOR[0]:,.0f}-{ANCHOR[1]:,.0f} RFU on the high channel.\n")
    print(f"    {'plate':52} {'gains':>9} {'ratio':>7} {'+/-':>6} "
          f"{'linear to':>10} {'overflowed':>11}")
    for _, row in scored.iterrows():
        gains = ("unstated" if row.gain_low == GAIN_UNSTATED
                 else f"{row.gain_low:.0f}->{row.gain_high:.0f}")
        print(f"    {row.plate[:52]:52} {gains:>9} {row.plate_ratio:7.3f} "
              f"{row.plate_scatter:5.1%} {row.plate_linear_max_rfu:10,.0f} "
              f"{row.plate_n_overflow:11,d}")

    refused = pairs[pairs.status != "scored"]
    if not refused.empty:
        print(f"\n  {len(refused)} pair(s) the plate could not answer with, and why. A pair "
              "that never reached")
        print("  the anchor window has no multiple to check; it is named, not dropped.\n")
        for _, row in refused.iterrows():
            print(f"    {row.plate[:44]:44} {row.low_channel}->{row.high_channel}")
            print(f"      {row.status} ({row.plate_n_overflow:,} of "
                  f"{row.plate_n_overflow + row.plate_n_pairs:,} readings overflowed)")

    print(f"\n  The {len(declared)} plates whose protocol declares gain 75 -> 100 give "
          f"{declared.plate_ratio.mean():.3f} "
          f"({declared.plate_ratio.min():.3f} to {declared.plate_ratio.max():.3f}), "
          f"a {declared.plate_ratio.std() / declared.plate_ratio.mean():.1%} spread.")

    three = table[(table.gain_low == 50) & (table.status == "scored")]
    two = table[(table.plate == "20260707_oxidative_stress_preliminary")
                & (table.gain_low == 75) & (table.status == "scored")]
    if not three.empty and not two.empty:
        low_step, high_step = float(three.plate_ratio.iloc[0]), float(two.plate_ratio.iloc[0])
        print("\n  And the check that makes those ratios a property of the DETECTOR rather "
              "than of a plate:")
        print("  one endpoint plate read at 50, 75 and 100 gives two gain steps that share "
              "no readings.")
        print(f"    50 -> 75  x{low_step:6.3f}   exponent "
              f"ln(r)/ln(75/50)  = {np.log(low_step) / np.log(75 / 50):.3f}")
        print(f"    75 -> 100 x{high_step:6.3f}   exponent "
              f"ln(r)/ln(100/75) = {np.log(high_step) / np.log(100 / 75):.3f}")

    pooled = table[(table.gain_low == 75) & (table.status == "scored")
                   & np.isfinite(table.departure_median)]
    print("\n  Departure from that line, pooled over those plates, by signal level.")
    print("  The column that matters is the median: a detector compressing would drift "
          "below 1.")
    print(f"  Below {ANCHOR[0]:,.0f} the check has no power -- the gain-75 partner is then "
          "a two-digit integer,")
    print("  and its quantisation, not the detector, is what the spread is measuring.\n")
    print(f"    {'high-gain RFU':>22} {'pairs':>8} {'median':>8} {'p05':>7} {'p95':>7}")
    for (low_edge, high_edge), block in pooled.groupby(["bin_low_rfu", "bin_high_rfu"]):
        n = int(block.n_pairs.sum())
        if not n:
            continue
        weights = block.n_pairs.to_numpy(dtype=float)
        median = float(np.average(block.departure_median, weights=weights))
        p05 = float(np.average(block.departure_p05, weights=weights))
        p95 = float(np.average(block.departure_p95, weights=weights))
        print(f"    {low_edge:9,.0f}-{high_edge:<12,.0f} {n:8,d} {median:8.4f} "
              f"{p05:7.4f} {p95:7.4f}")

    census = overflow_census(directory)
    print(f"\n  It never bends. It stops: {census['flagged']:,} of the "
          f"{census['readings']:,} fluorescence readings in the corpus")
    print(f"  are the overflow flag {RFU_OVERFLOW:,.0f} rather than a number, across "
          f"{census['channels_with_an_overflow']} channels, and the largest reading")
    print(f"  anywhere is {census['highest_reading']:,.0f} against a ceiling of "
          f"{RFU_CEILING:,.0f}. Nothing in between was ever reported.")

    # And the flag is not raised early: what the surviving partner says a flagged reading
    # WOULD have been. If the detector gave up before the top of scale this figure would
    # come in well under it.
    weakest, overflowing = np.inf, 0
    for name, gains in DECLARED_GAINS.items():
        path = directory / name
        if not path.is_file() or gains[-2:] != (75, 100):
            continue
        read = gen5.read_xpt(path)
        low, high = _pairs(read, *[c.name for c in gain_pairs(read)[-1]])
        flagged = is_overflow(high) & ~is_overflow(low)
        if not flagged.any():
            continue
        overflowing += 1
        fit = fit_gain_pair(low, high, tolerance=TOLERANCE, anchor=ANCHOR)
        weakest = min(weakest, float((fit.ratio * low[flagged] + fit.offset).min()))
    if overflowing:
        print(f"  Nor is the flag raised early: over the {overflowing} declared plates that "
              f"overflow, the lowest value")
        print(f"  any flagged reading's surviving partner predicts is {weakest:,.0f} -- "
              f"{weakest / RFU_CEILING:.1%} of the ceiling.")

    print("\n  Which channel is checked over which range, because they are not the same "
          "range.\n")
    lows, highs = [], []
    for name, gains in DECLARED_GAINS.items():
        path = directory / name
        if not path.is_file() or gains[-2:] != (75, 100):
            continue
        read = gen5.read_xpt(path)
        low, high = _pairs(read, *[c.name for c in gain_pairs(read)[-1]])
        both = ~(is_overflow(low) | is_overflow(high))
        lows.append(low[both])
        highs.append(high[both])
    low, high = np.concatenate(lows), np.concatenate(highs)
    print(f"    gain 75, checked against gain 100 : {low.min():,.0f} to {low.max():,.0f} "
          f"RFU   ({low.size:,} pairs)")
    print(f"    gain 100, checked against gain 75 : {high.min():,.0f} to "
          f"{high.max():,.0f} RFU   (the same pairs)")
    print("    Above the gain-75 figure its partner has overflowed and can check nothing, "
          "so a gain-75")
    print("    reading larger than that rests on the gain-100 channel having been linear "
          "there instead.")
    unchecked = float(max(
        gen5.read_xpt(directory / name).channel("mCitrine:480,530").values.max()
        for name, gains in DECLARED_GAINS.items()
        if (directory / name).is_file() and gains[-2:] == (75, 100)))
    print(f"    The largest gain-75 reading in the corpus is {unchecked:,.0f} RFU, which is "
          f"{unchecked / low.max():.1f}x past that")
    print("    and is therefore an inference from the other channel rather than a "
          "measurement of this one.")

    print("\n  What the analyses that read gain 75 are actually sitting at.\n")
    for name in AFL_PLATES:
        path = directory / name
        if not path.is_file():
            continue
        values = gen5.read_xpt(path).channel("mCitrine:480,530").values
        print(f"    {path.stem[:52]:52} max {values.max():8,.0f} RFU  "
              f"= {values.max() / RFU_CEILING:5.1%} of the ceiling, "
              f"{int(is_overflow(values).sum())} overflowed")
    plate_dir = paths.REPO_ROOT / "data" / "plates"
    committed = sorted(plate_dir.glob("*mCitrine*.csv"))
    if committed:
        n, highest = 0, -np.inf
        for path in committed:
            block = pd.read_csv(path).drop(columns=["elapsed_hms", "temperature_c"],
                                           errors="ignore").to_numpy(dtype=float)
            n += int(np.isfinite(block).sum())
            highest = max(highest, float(np.nanmax(block)))
        print(f"\n    data/plates, {len(committed)} committed mCitrine tables: "
              f"{n:,} readings, max {highest:,.0f} RFU")
        print(f"    = {highest / RFU_CEILING:.1%} of the ceiling. None is the overflow flag.")

    absorbance = absorbance_top_of_scale(directory)
    print(f"\n  And the negative, on the OTHER detector. OD600 reaches "
          f"{absorbance['max_od']:.3f} across {absorbance['plates']} plates, whose maxima "
          f"take {absorbance['distinct_plate_maxima']}")
    print(f"  distinct values, and the most readings any one plate has AT its own maximum "
          f"is {absorbance['most_readings_at_a_plate_max']}.")
    print("  There is no top of scale for absorbance in this corpus -- but that is "
          "clipping, not linearity,")
    print("  and it is a photodiode measuring transmission rather than a photomultiplier "
          "counting emission.")
    print("  NOTHING here sets OpticalQualityGate.od_linear_max. It stays the placeholder "
          "1.0 and still")
    print("  needs the dilution series calib/od_linearity.py asks for.")

    out = paths.outputs_dir() / "gain_linearity.csv"
    table.to_csv(out, index=False)
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
