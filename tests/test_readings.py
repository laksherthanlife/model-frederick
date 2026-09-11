"""The correction-state types, and the mistakes they exist to make impossible.

`docs/ARCHITECTURE.md` §6.3 carried three tables saying which functions want a raw reading
and which want a corrected one. Every entry was true and none was enforced, so the tables
were a promise about what callers would remember. These tests are the same statements with
a failure mode.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ystwin.readings import (
    CorrectedOD,
    CorrectedRFU,
    RawOD,
    RawRFU,
    Reading,
    ReporterPerGDCW,
    SpecificFluorescence,
    require,
)


class TestTheConversionsAreTheOnlyRouteBetweenStates:
    def test_subtracting_the_blank_gives_a_corrected_density(self):
        corrected = RawOD([0.19, 0.29, 0.39]).minus_blank(0.09)

        assert isinstance(corrected, CorrectedOD)
        assert corrected.values == pytest.approx([0.10, 0.20, 0.30])

    def test_subtracting_the_background_gives_a_corrected_fluorescence(self):
        corrected = RawRFU([350.0, 450.0]).minus_background(250.0)

        assert isinstance(corrected, CorrectedRFU)
        assert corrected.values == pytest.approx([100.0, 200.0])

    def test_a_corrected_density_cannot_have_its_blank_taken_off_twice(self):
        """The state machine is one-way on purpose: there is no second blank to remove."""
        assert not hasattr(CorrectedOD([0.1, 0.2]), "minus_blank")

    def test_a_corrected_fluorescence_over_a_corrected_density_is_specific(self):
        specific = RawRFU([350.0, 450.0]).minus_background(250.0).per(
            RawOD([0.19, 0.29]).minus_blank(0.09))

        assert isinstance(specific, SpecificFluorescence)
        assert specific.values == pytest.approx([1000.0, 1000.0])

    def test_dividing_by_an_uncorrected_density_is_refused(self):
        """The composite mistake, and the one that looks most like ordinary code: a
        fluorescence that HAS been corrected over a density that has not."""
        with pytest.raises(TypeError, match="minus_blank"):
            CorrectedRFU([100.0, 200.0]).per(RawOD([0.19, 0.29]))

    def test_a_blank_that_is_not_a_number_is_refused_at_the_subtraction(self):
        """A NaN blank turns a plate into NaNs several functions downstream, where the
        cause is no longer visible."""
        with pytest.raises(ValueError, match="finite"):
            RawOD([0.19]).minus_blank(float("nan"))

        with pytest.raises(ValueError, match="finite"):
            RawRFU([350.0]).minus_background(float("inf"))


class TestTheMistakeThatActuallyHappens:
    """The wrong variable of the right shape. Both are float arrays of the same length
    off the same plate, which is why the parameter name never caught it."""

    def test_a_bare_array_is_not_a_reading(self):
        with pytest.raises(TypeError, match="bare ndarray"):
            require(np.array([0.1, 0.2]), CorrectedOD, name="optical_density")

    def test_the_message_for_a_bare_array_says_what_to_write(self):
        with pytest.raises(TypeError, match=r"CorrectedOD\(values\)"):
            require(np.array([0.1]), CorrectedOD, name="optical_density")

    def test_the_message_for_the_wrong_state_names_the_conversion(self):
        """"expected CorrectedOD, got RawOD" says what is wrong. It does not say what to
        do, and a caller who knew that would not have made the mistake."""
        with pytest.raises(TypeError, match=r"call \.minus_blank\(od_blank\) on it"):
            require(RawOD([0.19]), CorrectedOD, name="optical_density")

    def test_asking_for_a_raw_reading_and_getting_a_corrected_one_is_also_caught(self):
        with pytest.raises(TypeError, match="before the blank came off"):
            require(CorrectedOD([0.1]), RawOD, name="optical_density")

    def test_a_reading_of_the_right_kind_passes_through_unchanged(self):
        reading = CorrectedOD([0.1, 0.2])

        assert require(reading, CorrectedOD, name="optical_density") is reading

    def test_more_than_one_kind_may_be_accepted(self):
        """`detect_blank_wells` is looking FOR the blank, so it takes either state."""
        raw = RawOD([0.09, 0.19])

        assert require(raw, RawOD, CorrectedOD, name="od") is raw

    def test_wrapping_a_reading_in_another_reading_is_refused(self):
        """A slip rather than an intention, and the outer wrapper would silently claim a
        state the inner one had already contradicted."""
        with pytest.raises(TypeError, match="already a reading"):
            CorrectedOD(RawOD([0.19, 0.29]))


class TestShapeIsNotPartOfTheState:
    """A well is a trace and a plate is a frame, and either can be raw or corrected."""

    def test_a_frame_survives_the_round_trip_as_a_frame(self):
        frame = pd.DataFrame({"A1": [0.19, 0.29], "H12": [0.09, 0.09]})

        corrected = RawOD(frame).minus_blank(0.09)

        assert isinstance(corrected.values, pd.DataFrame)
        assert list(corrected.values.columns) == ["A1", "H12"]
        assert corrected.values["A1"].tolist() == pytest.approx([0.10, 0.20])

    def test_a_trace_is_coerced_to_a_float_array(self):
        """So a list from a parser and an array from a fit are the same thing here."""
        reading = CorrectedOD([1, 2, 3])

        assert isinstance(reading.values, np.ndarray)
        assert reading.values.dtype == float

    def test_array_reaches_the_numbers_whichever_container_they_came_in(self):
        frame = CorrectedOD(pd.DataFrame({"A1": [0.1], "A2": [0.2]}))
        trace = CorrectedOD([0.1, 0.2])

        assert frame.array.ravel().tolist() == pytest.approx([0.1, 0.2])
        assert trace.array.tolist() == pytest.approx([0.1, 0.2])

    def test_numpy_can_consume_a_reading_directly(self):
        """The escape hatch for arithmetic inside a consumer that has already checked the
        type. It does not weaken `require`, which asks about the type and not the buffer."""
        assert np.asarray(CorrectedOD([0.1, 0.2])).tolist() == pytest.approx([0.1, 0.2])
        assert len(CorrectedOD([0.1, 0.2])) == 2

    def test_asking_numpy_for_a_copy_gets_a_copy(self):
        """The first version of `__array__` ignored the flag and handed back the reading's
        own buffer, so writing to the copy wrote through to the reading. A type added to
        stop one silent aliasing mistake must not ship another."""
        reading = CorrectedOD([1.0, 2.0, 3.0])

        duplicate = np.array(reading, copy=True)
        duplicate[0] = 99.0

        assert reading.values[0] == 1.0

    def test_refusing_to_copy_is_honoured_and_refused_when_impossible(self):
        """`copy=False` means "a view or an error", never a silent copy."""
        reading = CorrectedOD([1.0, 2.0])

        assert np.shares_memory(np.array(reading, copy=False), reading.values)  # a view, not a copy
        with pytest.raises(ValueError, match="without copying"):
            np.array(reading, dtype=np.int64, copy=False)


class TestTheTypesAreDistinctAndNotInterchangeable:
    def test_every_state_is_its_own_type(self):
        kinds = [RawOD, CorrectedOD, RawRFU, CorrectedRFU, SpecificFluorescence]

        assert len({k.__name__ for k in kinds}) == 5
        assert all(issubclass(k, Reading) for k in kinds)

    def test_a_density_is_not_a_fluorescence_even_though_both_are_corrected(self):
        with pytest.raises(TypeError):
            require(CorrectedRFU([100.0]), CorrectedOD, name="optical_density")

    def test_a_specific_fluorescence_is_not_a_corrected_one(self):
        """§6.3's third hazard: `promoter_activity` documents a per-cell concentration and
        every caller passes RFU/OD. Those are different quantities and this says so."""
        with pytest.raises(TypeError):
            require(SpecificFluorescence([1000.0]), CorrectedRFU, name="reporter")


class TestTheConverterArchitectureCalledMissing:
    """§6.3 recorded that per-OD and per-gDCW "are never used together and there is no
    converter", and named the reconciliation as what wiring `estimator.py` needs first.

    The filter's reporter state is per gDCW; every real caller of `promoter_activity` passes
    RFU/OD. The two differ by `gdcw_per_od`, and that factor is confounded with the optics
    gain -- `generator/calibrate.py` fits their product and reports the pair as
    unidentifiable without a measured dry weight. So the converter exists and refuses to
    guess, which is the only honest shape it can have.
    """

    def test_it_converts_on_the_stated_factor(self):
        converted = SpecificFluorescence([1000.0, 2000.0]).per_gdcw(0.42)

        assert isinstance(converted, ReporterPerGDCW)
        assert converted.values == pytest.approx([1000.0 / 0.42, 2000.0 / 0.42])

    def test_the_factor_cannot_be_defaulted(self):
        """A default would invent the half of `gain x gdcw_per_od` that nobody separated."""
        with pytest.raises(TypeError):
            SpecificFluorescence([1000.0]).per_gdcw()

    @pytest.mark.parametrize("factor", [0.0, -0.42, float("nan"), float("inf")])
    def test_a_factor_that_is_not_a_positive_number_is_refused(self, factor):
        """Zero divides the trace to infinity and a negative flips its sign; both would read
        as the reporter doing something."""
        with pytest.raises(ValueError, match="positive finite"):
            SpecificFluorescence([1000.0]).per_gdcw(factor)

    def test_the_two_bases_are_different_types(self):
        """Which is the point. They were both bare arrays, so nothing stopped a per-OD trace
        reaching a per-gDCW slot -- a silent error of exactly `gdcw_per_od`, about 2.4x."""
        with pytest.raises(TypeError):
            require(SpecificFluorescence([1.0]), ReporterPerGDCW, name="reporter")
        with pytest.raises(TypeError):
            require(ReporterPerGDCW([1.0]), SpecificFluorescence, name="reporter")

    def test_the_conversion_is_one_way(self):
        """There is no inverse, deliberately. Going back would need the same unmeasured
        factor, and a round trip would make it look free."""
        assert not hasattr(ReporterPerGDCW([1.0]), "per_gdcw")
        assert not hasattr(ReporterPerGDCW([1.0]), "per_od")


@pytest.mark.parametrize("density", [0.0, -0.1, np.nan, np.inf])
def test_specific_fluorescence_requires_finite_positive_density(density):
    with pytest.raises(ValueError, match="positive|finite"):
        CorrectedRFU([1.0, 2.0]).per(CorrectedOD([0.1, density]))


@pytest.mark.parametrize("signal", [np.nan, np.inf, -np.inf])
def test_specific_fluorescence_refuses_nonfinite_signal(signal):
    with pytest.raises(ValueError, match="finite"):
        CorrectedRFU([1.0, signal]).per(CorrectedOD([0.1, 0.2]))


def test_specific_fluorescence_cannot_broadcast_unpaired_measurements():
    with pytest.raises(ValueError, match="shape"):
        CorrectedRFU([1.0, 2.0]).per(CorrectedOD([0.1]))


@pytest.mark.parametrize("axis", ["index", "columns"])
def test_specific_fluorescence_requires_matching_frame_labels(axis):
    signal = pd.DataFrame({"A1": [1.0, 2.0]}, index=[0.0, 1.0])
    density = pd.DataFrame({"A1": [0.1, 0.2]}, index=signal.index)
    setattr(density, axis, [2.0, 3.0] if axis == "index" else ["B1"])

    with pytest.raises(ValueError, match="labels"):
        CorrectedRFU(signal).per(CorrectedOD(density))


def test_no_copy_numpy_protocol_shares_memory_rather_than_only_matching_values():
    reading = CorrectedOD([1.0, 2.0])
    view = np.array(reading, copy=False)

    assert np.shares_memory(view, reading.values)
    view[0] = 3.0
    assert reading.values[0] == 3.0
