"""The generator's acceptance test: synthetic plates must survive the real pipeline.

Write a synthetic plate as a Synergy export, then run it through the same reader,
the same optical-quality gate, the same layout recovery and the same dilution
analysis used on wet-lab files. If any stage behaves differently, the generator is
producing something that only looks right to itself.
"""

import numpy as np
import pytest

from ystwin.diagnostics.dilution_confound import detect_blank_wells, dilution_confound_report
from ystwin.gates.g1_optical import OpticalQualityGate, assess_plate
from ystwin.generator.export import write_synergy_export
from ystwin.generator.plate import DEFAULT_PANEL, PlateConditions, generate_plate
from ystwin.plate.synergy import read_synergy_kinetic
from ystwin.reporter import ReporterKinetics
from ystwin.readings import RawOD, RawRFU


@pytest.fixture(scope="module")
def written(tmp_path_factory):
    plate = generate_plate(DEFAULT_PANEL, PlateConditions(seed=7))
    path = write_synergy_export(plate, tmp_path_factory.mktemp("gen") / "synthetic.xlsx")
    return plate, path


def test_the_real_reader_parses_the_generated_export(written):
    _, path = written

    run = read_synergy_kinetic(path)

    assert set(run.channel_names) == {"OD600", "mCitrine"}


def test_the_values_survive_the_write_and_read(written):
    plate, path = written

    run = read_synergy_kinetic(path)
    od = run.channel("OD600").data

    assert od["A1"].to_numpy() == pytest.approx(plate.od["A1"].to_numpy(), rel=1e-6)


def test_the_timebase_survives(written):
    plate, path = written

    times = np.asarray(read_synergy_kinetic(path).channel("OD600").times_h)

    assert times == pytest.approx(plate.times_h, abs=1e-3)


def test_blank_detection_finds_the_wells_that_were_generated_as_blanks(written):
    plate, path = written
    run = read_synergy_kinetic(path)
    aligned = run.aligned("OD600")

    found = detect_blank_wells(RawOD(aligned["OD600"]), rfu=RawRFU(aligned["mCitrine"]))

    assert set(found) == set(plate.blank_wells)


def _g1(plate, path):
    aligned = read_synergy_kinetic(path).aligned("OD600")
    cultures = [w for w in aligned["OD600"].columns if w not in plate.blank_wells]
    return assess_plate(
        RawOD(aligned["OD600"][cultures]), od_blank=plate.conditions.od_blank,
        gate=OpticalQualityGate(od_linear_max=1.0),
    )


def test_the_optical_gate_gives_a_pass_rate_like_a_real_plate(written):
    """The 2026-07-22 plate passed 67 of 84 wells, 80%."""
    plate, path = written

    assert 0.7 < _g1(plate, path).passed.mean() < 0.95


def test_no_well_fails_on_linear_range(written):
    """The NewProtocol inoculation keeps every culture inside the reader's range."""
    plate, path = written
    table = _g1(plate, path)

    assert "linear_range" not in set(table[~table.passed].failures.str.split(",").explode())


def test_the_wells_that_fail_are_the_high_dose_ones_as_on_the_real_plates(written):
    """Strong inhibition leaves a near-flat trace, which reads as decline under noise."""
    plate, path = written
    table = _g1(plate, path).set_index("well")
    truth = plate.truth.set_index("well")

    failed = [w for w in table.index if not table.loc[w, "passed"] and w in truth.index]
    passed = [w for w in table.index if table.loc[w, "passed"] and w in truth.index]

    assert truth.loc[failed].dose_index.mean() > truth.loc[passed].dose_index.mean()


def test_the_failure_modes_match_the_ones_the_real_plates_show(written):
    plate, path = written
    table = _g1(plate, path)
    reasons = set(table[~table.passed].failures.str.split(",").explode())

    assert reasons <= {"sustained_decline", "dynamic_range"}


def test_the_layout_is_recovered_from_the_generated_plate(written):
    from ystwin.generator.export import derived_dose_response
    from ystwin.plate.layout import recover_layout

    plate, path = written
    aligned = read_synergy_kinetic(path).aligned("OD600")
    cultures = [w for w in aligned["OD600"].columns if w not in plate.blank_wells]
    normalised = (
        (aligned["mCitrine"][cultures] - plate.conditions.optics.background)
        / (aligned["OD600"][cultures] - plate.conditions.od_blank)
    )

    found = recover_layout(normalised, derived_dose_response(plate))

    assert found.layout.construct_columns == plate.conditions.layout.construct_columns


def _report(plate, aligned, well, subtract_autofluorescence=False):
    background = plate.conditions.optics.background
    rfu = aligned["mCitrine"][well].to_numpy()
    od = aligned["OD600"][well].to_numpy()
    if subtract_autofluorescence:
        # What a reporter-free parent (BY4741) contributes, per unit biomass.
        biomass = (od - plate.conditions.od_blank) * plate.conditions.gdcw_per_od
        rfu = rfu - plate.conditions.optics.autofluorescence * biomass
    return dilution_confound_report(
        aligned.index.to_numpy(), RawOD(od), RawRFU(rfu),
        kinetics=ReporterKinetics(k_deg=0.0),
        od_blank=plate.conditions.od_blank, rfu_background=background,
    )


def test_the_dilution_analysis_recovers_the_promoter_activity_that_was_generated(written):
    """Known input in, known input out.

    The observable is ``S = (RFU-bg)/(OD-blank) = (autofluorescence + gain*R) * gdcw``,
    so inverting it returns ``gdcw*gain*k_synth`` plus a term ``mu*gdcw*autofluorescence``
    contributed by the cells' own background emission.
    """
    plate, path = written
    aligned = read_synergy_kinetic(path).aligned("OD600")
    truth = plate.truth.set_index("well")
    conditions = plate.conditions

    errors = []
    for well in ("A1", "D1", "A4", "D4"):
        report = _report(plate, aligned, well)
        recovered = np.nanmedian(report.promoter_activity[5:-5])
        mu = np.nanmedian(report.growth_rate[5:-5])
        expected = (
            conditions.gdcw_per_od * conditions.optics.gain * truth.loc[well].promoter_activity
            + mu * conditions.gdcw_per_od * conditions.optics.autofluorescence
        )
        errors.append(abs(recovered - expected) / expected)

    assert np.median(errors) < 0.15


class TestAutofluorescenceIsItsOwnConfound:
    """Cellular background emission adds a growth-rate-proportional term.

    Inverting the raw ratio returns ``gdcw*gain*k_synth + mu*gdcw*autofluorescence``.
    The second term is not promoter activity, and because it scales with mu it moves
    with dose exactly as the dilution artefact does. Subtracting what a reporter-free
    parent emits removes it -- which is what the BY4741 control plates are for.
    """

    def test_leaving_it_in_biases_the_estimate_upward(self, written):
        plate, path = written
        aligned = read_synergy_kinetic(path).aligned("OD600")

        raw = np.nanmedian(_report(plate, aligned, "A1").promoter_activity[5:-5])
        corrected = np.nanmedian(
            _report(plate, aligned, "A1", subtract_autofluorescence=True).promoter_activity[5:-5]
        )

        assert raw > corrected

    def test_subtracting_it_recovers_the_generated_activity_directly(self, written):
        plate, path = written
        aligned = read_synergy_kinetic(path).aligned("OD600")
        truth = plate.truth.set_index("well")
        conditions = plate.conditions

        errors = []
        for well in ("A1", "D1", "A4", "D4"):
            report = _report(plate, aligned, well, subtract_autofluorescence=True)
            recovered = np.nanmedian(report.promoter_activity[5:-5])
            expected = (
                conditions.gdcw_per_od * conditions.optics.gain
                * truth.loc[well].promoter_activity
            )
            errors.append(abs(recovered - expected) / expected)

        assert np.median(errors) < 0.10


def test_a_generated_plate_shows_the_dilution_confound_the_real_ones_do(written):
    """Naive fold should exceed corrected fold, as on every real plate."""
    plate, path = written
    aligned = read_synergy_kinetic(path).aligned("OD600")

    report = dilution_confound_report(
        aligned.index.to_numpy(),
        RawOD(aligned["OD600"]["F4"].to_numpy()),      # UPRE2 at a high dose
        RawRFU(aligned["mCitrine"]["F4"].to_numpy()),
        kinetics=ReporterKinetics(k_deg=0.0),
        od_blank=plate.conditions.od_blank,
        rfu_background=plate.conditions.optics.background,
    )

    assert report.negative_activity_fraction < 0.25
