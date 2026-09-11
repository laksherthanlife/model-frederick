"""What the three ATP-sensor builds actually measured, and what they could not have.

Three constructs, all reported as not working: an ICL-UAS build (2026-08-05) and an ACS-UAS
and a second ICL-UAS build (2026-08-09), each read across four galactose:glucose ratios at
several total sugar concentrations.

One thing here is solid and it is a correction to how the plates were read. At a single
timepoint the plate spreads three- to fourfold, monotone in both axes, and that looks like a
dose response. It is growth phase: faster wells reach any point of the growth curve sooner,
and conditions sit up to eleven hours apart at the same density. Matched on density the
whole plate spreads 1.1-fold. There is no carbon-source response.

The second thing is that the plates do contain a response, and the single-timepoint
comparison is what hid it. Each condition's RFU/OD is flat through exponential growth and
then rises to a peak: 2.84-fold at 2% glucose peaking at h12.6, 2.12-fold at 1.5% peaking at
h11.5, 1.48-fold at 1% peaking at h9.0, and nothing at 0.5%. The peak time moves with how
much glucose was there, which is what exhausting it looks like, and the response needs
glucose -- the galactose-rich columns do not peak at all.

Reading every well at h10 mixes that up with itself. At h10 the 2% and 1.5% wells are still
climbing toward peaks they have not reached, while the 1% well peaked an hour earlier and
the 0.5% well never peaked. A bar chart at one timepoint then reports how far through its
own transient each condition happens to be, not how much signal it makes.

Two controls would have settled it and neither plate carried them: a promoterless strain, so
baseline fluorescence can be separated from autofluorescence, and a constitutive-promoter
strain, so the growth-rate coupling can be subtracted. Nor was the axis these promoters read
ever varied -- CSRE answers to Cat8 and Adr1, which answer to running out of fermentable
sugar. ACS1 derepresses several hundredfold on ethanol (Kratzer & Schuller 1997,
PMID 9427394); galactose is another fermentable sugar and the plates held nothing else.
"""

import numpy as np
import pytest

from ystwin.analysis.matched_density import (
    fold_range,
    specific_fluorescence_at_density,
    time_to_density,
)
from ystwin.generator.synergy import read_export
from ystwin.readings import CorrectedOD, CorrectedRFU

RATIOS = {1: "glu 10:0", 2: "gal 3:7", 3: "gal 6:4", 4: "gal 9:1"}


def wells_by_sample(layout, keep_rows=None):
    out = {}
    for well, sample in layout.items():
        if sample == "BLK" or (keep_rows and well[0] not in keep_rows):
            continue
        out.setdefault(sample, []).append(well)
    return out


def condition(wells):
    """Plate row is total sugar, column group is the galactose:glucose ratio."""
    first = sorted(wells, key=lambda w: int(w[1:]))[0]
    return first[0], RATIOS[(int(first[1:]) - 1) // 3 + 1]


def series(path, keep_rows=None):
    layout, tables = read_export(path)
    od, fl = tables["Blank OD600:600"], tables["Blank mCitrine:480,530"]
    out = {}
    for sample, wells in wells_by_sample(layout, keep_rows).items():
        out[condition(wells)] = (od.index.values, od[wells].mean(axis=1).values,
                                 fl[wells].mean(axis=1).values)
    return out


def at_time(path, hour, keep_rows=None):
    vals = []
    for t, o, f in series(path, keep_rows).values():
        i = int(np.abs(t - hour).argmin())
        vals.append(f[i] / max(o[i], 0.02))
    return vals


def at_density(path, target, keep_rows=None):
    return [specific_fluorescence_at_density(t, CorrectedOD(o), CorrectedRFU(f), target)
            for t, o, f in series(path, keep_rows).values()]


PLATES = [("icl", None), ("acs_icl", set("ABC")), ("acs_icl", set("DEF"))]


@pytest.fixture(scope="module")
def plates(atp_sensor_icl, atp_sensor_acs_icl):
    return {"icl": atp_sensor_icl, "acs_icl": atp_sensor_acs_icl}


class TestTheSignalIsRealNotNoise:
    @pytest.mark.parametrize("key,rows", PLATES)
    def test_replicates_agree_to_within_a_few_percent(self, plates, key, rows):
        """Typical condition, not the worst: one well group on each plate sits near 15%."""
        layout, tables = read_export(plates[key])
        od, fl = tables["Blank OD600:600"], tables["Blank mCitrine:480,530"]
        spread = []
        for wells in wells_by_sample(layout, rows).values():
            a = fl[wells].values / np.clip(od[wells].values, 0.02, None)
            spread.append(float(np.median(a.std(axis=1) / a.mean(axis=1))))

        assert float(np.median(spread)) < 0.08

    @pytest.mark.parametrize("key,rows", PLATES)
    def test_and_at_one_timepoint_it_looks_like_a_strong_dose_response(self, plates, key, rows):
        assert fold_range(at_time(plates[key], 12.0, rows)) > 2.5


class TestButMostOfThatIsGrowthPhase:
    @pytest.mark.parametrize("key,rows", PLATES)
    def test_matched_on_density_the_early_response_nearly_vanishes(self, plates, key, rows):
        assert fold_range(at_density(plates[key], 0.5, rows)) < 1.3

    @pytest.mark.parametrize("key,rows", PLATES)
    def test_which_is_far_less_than_the_same_plate_read_at_one_hour(self, plates, key, rows):
        at_od = fold_range(at_density(plates[key], 0.5, rows))
        at_clock = fold_range(at_time(plates[key], 12.0, rows))

        assert at_clock > 2 * at_od

    def test_because_the_conditions_are_hours_apart_at_the_same_density(self, plates):
        reached = [time_to_density(t, CorrectedOD(o), 1.0) for t, o, _ in series(plates["icl"]).values()]
        reached = [h for h in reached if np.isfinite(h)]

        assert max(reached) - min(reached) > 4.0


SUGAR = {"A": "2%", "B": "1.5%", "C": "1%", "D": "0.5%"}


def induction(path, row, ratio):
    """Peak RFU/OD against the same well's own exponential-phase baseline."""
    t, o, f = series(path)[(row, ratio)]
    ratio_series = f / np.clip(np.maximum.accumulate(o), 0.02, None)
    baseline = float(np.median(ratio_series[(t > 2) & (t < 5)]))
    return float(ratio_series.max() / baseline), float(t[int(np.argmax(ratio_series))])


class TestThereIsAResponseAndItIsOrdered:
    """What the fixed-timepoint comparison threw away."""

    def test_the_richest_glucose_well_nearly_triples(self, plates):
        fold, _ = induction(plates["icl"], "A", "glu 10:0")

        assert fold > 2.5

    def test_induction_falls_with_the_glucose_available(self, plates):
        folds = [induction(plates["icl"], row, "glu 10:0")[0] for row in "ABCD"]

        assert folds == sorted(folds, reverse=True)
        assert folds[-1] < 1.2

    def test_and_the_peak_arrives_later_the_more_glucose_there_was(self, plates):
        """A dose response in time, which is what running out of a sugar looks like."""
        peaks = [induction(plates["icl"], row, "glu 10:0")[1] for row in "ABC"]

        assert peaks == sorted(peaks, reverse=True)

    def test_the_galactose_rich_columns_never_peak(self, plates):
        for ratio in ("gal 6:4", "gal 9:1"):
            fold, when = induction(plates["icl"], "A", ratio)

            assert fold < 1.5, ratio
            assert when < 2.0, ratio


class TestWhyTheSingleTimepointHidIt:
    def test_at_ten_hours_the_conditions_are_at_different_points_of_their_own_transient(
            self, plates):
        still_rising = [induction(plates["icl"], row, "glu 10:0")[1] > 10.0 for row in "ABCD"]

        assert still_rising == [True, True, False, False]

    def test_so_the_ordering_at_one_timepoint_is_not_the_ordering_of_response(self, plates):
        """2% reads higher than 1% at h10 partly because 1% has already come back down."""
        t, o, f = series(plates["icl"])[("C", "glu 10:0")]
        ratio_series = f / np.clip(np.maximum.accumulate(o), 0.02, None)
        at_ten = ratio_series[int(np.abs(t - 10.0).argmin())]

        assert at_ten < ratio_series.max()


class TestWhatStillCannotBeConcludedFromThesePlates:
    def test_no_promoterless_control_was_on_either_plate(self, plates):
        """Mitochondrial biogenesis at the diauxic shift raises flavin autofluorescence at
        480/530, and it would scale with glucose consumed exactly as the signal does. That
        is the one alternative that fits everything here, and only this control excludes it.
        """
        for path in plates.values():
            layout, _ = read_export(path)

            assert set(layout.values()) - {"BLK"} == {
                f"SPL{i}" for i in range(1, len(set(layout.values())))}

    def test_and_raw_fluorescence_falls_after_its_peak_which_nothing_stable_does(self, plates):
        """mCitrine is not degraded and cells do not lose it, so a real decline means the
        measurement is losing signal in a way none of the kinetics above models."""
        _, tables = read_export(plates["icl"])
        layout, _ = read_export(plates["icl"])
        wells = [w for w, s in layout.items() if s == "SPL1"]
        raw = tables["Blank mCitrine:480,530"][wells].mean(axis=1).values

        assert raw[-1] < 0.85 * raw.max()


class TestTheModelAgreesWithWhatWasMeasured:
    """The panel already routes these constructs to the carbon regulon, not to ATP: the
    CSRE reporter's marker genes are ADH2 and ICL1 and its module is Snf1 -> Adr1/Cat8."""

    def test_a_CSRE_reporter_reads_carbon_and_not_atp(self):
        from ystwin.generator.stress_panel import REPORTERS

        assert REPORTERS["CSRE-carbon"].module == "carbon"
        assert "ICL1" in REPORTERS["CSRE-carbon"].marker_gene

    def test_the_only_atp_reporter_is_ratiometric_not_transcriptional(self):
        """Which is the design conclusion: no promoter fusion reads the adenylate pool."""
        from ystwin.generator.stress_panel import Kind, REPORTERS

        atp = [name for name, r in REPORTERS.items() if r.module == "atp"]

        assert atp == ["QUEEN-2m"]
        assert REPORTERS["QUEEN-2m"].kind is Kind.RATIOMETRIC

    def test_ethanol_derepresses_the_carbon_regulon_far_harder_than_galactose(self):
        """The axis the plates never varied. Not set from this data: nothing here can
        calibrate galactose without a control that says the reporter reads anything."""
        from ystwin.generator.context import CultureContext, baseline_activity

        on_galactose = baseline_activity(CultureContext(carbon_source="galactose"))["carbon"]
        on_ethanol = baseline_activity(CultureContext(carbon_source="ethanol"))["carbon"]

        assert on_ethanol > on_galactose

    def test_and_glucose_starvation_is_what_actually_drives_it(self):
        from ystwin.generator.stress_panel import STRESSORS

        assert "carbon" in STRESSORS["glucose_starvation"].targets


class TestWhatQUEENWasActuallyValidatedOver:
    """The ATP channel of the three-sensor build, checked against its own paper.

    Yaginuma 2014 (PMID 25283467) gives Kd 4.5 mM for ATP at 25 C, no response to Mg2+ at
    1-2 mM, and no response to pH over 7.3-8.8. The first two are fine: the model runs free
    magnesium at 1 mM, inside the tested band. The third stops at 7.3, and a glucose-starved
    yeast cytosol drops to 6.4 (Dechant 2010, PMID 20581803) -- the very event an ATP
    channel is for. The protocol paper that followed (Takaine 2019) does not characterise
    pH either, so the gap is in the record, not in this model.
    """

    def test_the_magnesium_the_model_assumes_is_inside_the_tested_band(self):
        from ystwin.bridge.thermodynamic import FREE_MAGNESIUM_M

        assert 1e-3 <= FREE_MAGNESIUM_M <= 2e-3

    def test_queen_is_not_pH_coupled_so_it_needs_no_pH_channel(self):
        from ystwin.analysis.sensor_selection import THREE_SENSOR_BUILD, interpretable

        assert interpretable(THREE_SENSOR_BUILD)

    def test_but_it_was_never_tested_where_a_starved_cytosol_goes(self):
        from ystwin.analysis.sensor_selection import THREE_SENSOR_BUILD, uncharacterised_ph
        from ystwin.generator.stress_panel import CYTOSOLIC_PH, REPORTERS

        low, _ = REPORTERS["QUEEN-2m"].validated_ph

        assert low > CYTOSOLIC_PH[0]
        assert uncharacterised_ph(THREE_SENSOR_BUILD) == ["QUEEN-2m"]

    def test_the_peroxide_sensor_does_not_have_this_problem(self):
        """HyPer7 was engineered pH-insensitive across the physiological range."""
        from ystwin.generator.stress_panel import REPORTERS

        assert not REPORTERS["HyPer7"].ph_uncharacterised

    def test_a_pH_coupled_probe_is_a_different_failure_and_stays_separate(self):
        from ystwin.analysis.sensor_selection import interpretable, uncharacterised_ph

        build = ["UPRE-ER", "roGFP2-Grx1", "QUEEN-2m"]

        assert not interpretable(build)
        assert "roGFP2-Grx1" not in uncharacterised_ph(build)
