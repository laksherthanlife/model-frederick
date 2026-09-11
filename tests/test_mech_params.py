"""The parameter type, the REFUSED sentinel, and the free-scalar gate.

Two of these classes are calibrations rather than unit tests and are marked as such in
their docstrings: ``TestTheGapsDocumentsPhSliceCountsEight`` rebuilds
`ARCHITECTURE_GAPS.md` 0.3's pH slice and checks the gate reproduces its published count,
and ``TestTheNoiseFloorsAreTheRepositorysOwn`` pins the two measured floors criterion (a)
scores against so a target cannot be given an invented one.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ystwin.generator.panel_experiment import (
    MEASURED_GROWTH_RATE_SE,
    OBSERVED_ACTIVITY_CV,
)
from ystwin.mech.params import (
    FREE_TAGS,
    PINNED_TAGS,
    FreeScalarGateFailed,
    NoSingleValue,
    Param,
    ParamRegistry,
    RefusedValue,
    SweptValue,
    Tag,
    Target,
    provenance_table,
    tag_census_table,
)


def a_measured() -> Param:
    """This repository's own plate CV, which is what the reduction criterion is set from."""
    return Param.measured(
        "OBSERVED_ACTIVITY_CV", OBSERVED_ACTIVITY_CV, "dimensionless",
        "generator/panel_experiment.py -- 28 matched conditions across the 2026-07-22 and "
        "2026-08-03 plates")


def a_refusal() -> Param:
    """`ARCHITECTURE_TARGET.md` DEPTH 2's anion exporter: the shipped example of a refusal."""
    return Param.refused(
        "Vmax_ex", "mmol/gDCW/h",
        "ARCHITECTURE_TARGET.md section 2 DEPTH 2; searched for Tpo2/Tpo3/Pdr12 capacity",
        reason="no measured Tpo2, Tpo3 or Pdr12 anion-export capacity exists in "
               "S. cerevisiae.",
        missing="an export rate at a known intracellular anion content, which also fixes "
                "K_ex.")


def a_sweep() -> Param:
    """Cytosolic buffering capacity: the architecture's headline swept axis."""
    return Param.swept(
        "beta", "mM/pH",
        "Gabba 2020, fitted 160-320 mM phosphate, which converts to 87-175 "
        "mM/pH at pH 7; Gerber 2016 reports 200 in yeast",
        bounds=(90.0, 200.0))


class TestTheTagIsMandatoryAtConstruction:
    """The audit's finding was 65% ASSERTED against 14% MEASURED. A default tag would be
    that finding with a type annotation on it."""

    def test_a_param_with_no_tag_is_not_constructible(self):
        with pytest.raises(TypeError) as excinfo:
            Param(1.0, "mM")

        assert "tag" in str(excinfo.value)

    def test_a_none_tag_is_refused_rather_than_treated_as_unknown(self):
        with pytest.raises(ValueError, match="is not one of"):
            Param(1.0, "mM", None, "somewhere")

    def test_an_unknown_tag_names_all_seven_and_says_why_it_matters(self):
        with pytest.raises(ValueError) as excinfo:
            Param(1.0, "mM", "PROBABLY_FINE", "somewhere")
        message = str(excinfo.value)

        for tag in Tag.ALL:
            assert tag in message
        assert "65%" in message and "14%" in message

    def test_there_is_one_constructor_per_tag_and_each_sets_its_own(self):
        built = {
            a_measured().tag,
            Param.bounded("ec50", 0.05, "mM", "Kritsiligkou 2021, PMID 34118234 floor to "
                          "the 0.15 mM regulon EC50", bounds=(0.02, 0.15),
                          missing="a HyPer7 half-max in yeast").tag,
            Param.borrowed("g_H", 25.0, "uS/cm2", "Sanders 1981", organism="Neurospora "
                           "crassa").tag,
            Param.derived("S_over_V", 1.39e4, "1/cm", "3/r for a 42 um3 sphere, BNID "
                          "100427").tag,
            Param.asserted("_LETHAL_MULTIPLE", 3.0, "x EC50",
                           "generator/stress_panel.py -- the panel's lethal dose is 3x "
                           "its EC50 because the ladder needed a top; the hardcoding "
                           "audit ranks it #1 by damage").tag,
            a_sweep().tag,
            a_refusal().tag,
        }

        assert built == set(Tag.ALL)


class TestTheRefusedSentinel:
    """`pathway/solve.py` refuses a node one constant short by name and
    `kinetic/carotenoid.py` refuses the unidentified half of a branch. This is the same
    act, moved into the type."""

    def test_float_raises_rather_than_returning_a_default(self):
        with pytest.raises(RefusedValue):
            float(a_refusal())

    def test_the_refusal_carries_the_reason_the_missing_measurement_and_the_citation(self):
        param = a_refusal()

        with pytest.raises(RefusedValue) as excinfo:
            float(param)
        message = str(excinfo.value)

        assert param.reason in message
        assert param.missing in message
        assert param.source in message
        assert param.name in message and param.units in message

    def test_it_is_a_notimplementederror_like_calibrated_kinetics(self):
        from ystwin.kinetic.carotenoid import calibrated_kinetics

        assert issubclass(RefusedValue, NotImplementedError)
        assert issubclass(NoSingleValue, NotImplementedError)
        with pytest.raises(NotImplementedError):
            calibrated_kinetics(0.1)

    def test_numpy_cannot_coerce_it_either(self):
        """The obvious way around a __float__ that raises is to let numpy do the cast."""
        with pytest.raises(RefusedValue):
            np.array([a_refusal()], dtype=float)

    def test_a_refusal_carrying_a_value_is_refused_at_construction(self):
        with pytest.raises(ValueError, match="decorative"):
            Param(0.0, "mmol/gDCW/h", Tag.REFUSED, "somewhere", name="Vmax_ex",
                  reason="none exists.", missing="a capacity.")

    def test_a_refusal_must_say_why_and_what_would_close_it(self):
        with pytest.raises(ValueError, match="no reason"):
            Param(None, "mM", Tag.REFUSED, "searched", name="x", missing="a number.")
        with pytest.raises(ValueError, match="would close it"):
            Param(None, "mM", Tag.REFUSED, "searched", name="x", reason="none exists.")

    def test_only_a_refusal_carries_a_reason(self):
        with pytest.raises(ValueError, match="carries a refusal reason"):
            Param(1.0, "mM", Tag.ASSERTED, "somewhere", name="x", reason="none exists.")


class TestASweptParamRefusesItsOwnMidpoint:
    """GAPS 0.3: scoring on whether some point in a sweep reproduces a target is fitting
    with extra steps. Taking a midpoint silently is the first step of that."""

    def test_float_raises_and_points_at_the_two_ways_through(self):
        with pytest.raises(SweptValue) as excinfo:
            float(a_sweep())
        message = str(excinfo.value)

        assert "band, not a point" in message
        assert ".sweep_points(" in message and ".at(" in message
        assert "90" in message and "200" in message

    def test_the_grid_includes_both_endpoints_and_is_increasing(self):
        points = a_sweep().sweep_points(5)

        assert points[0] == 90.0 and points[-1] == 200.0
        assert len(points) == 5
        assert list(points) == sorted(points)

    def test_a_grid_of_one_point_is_a_sweep_that_reports_one_answer(self):
        with pytest.raises(ValueError, match="both endpoints"):
            a_sweep().sweep_points(1)

    def test_at_pins_one_point_inside_the_band_and_refuses_one_outside(self):
        param = a_sweep()

        assert param.at(120.0) == 120.0
        with pytest.raises(ValueError, match="outside it"):
            param.at(320.0)

    def test_a_point_valued_param_has_no_grid(self):
        with pytest.raises(ValueError, match="not SWEPT"):
            a_measured().sweep_points(3)
        with pytest.raises(ValueError, match="not SWEPT"):
            a_measured().at(0.1)

    def test_a_sweep_needs_bounds_and_refuses_a_value(self):
        with pytest.raises(ValueError, match="needs bounds"):
            Param(None, "mM/pH", Tag.SWEPT, "somewhere", name="beta")
        with pytest.raises(ValueError, match="decorative"):
            Param(145.0, "mM/pH", Tag.SWEPT, "somewhere", name="beta", bounds=(90.0, 200.0))

    def test_bounds_must_increase(self):
        with pytest.raises(ValueError, match="increasing"):
            Param.swept("beta", "mM/pH", "somewhere", bounds=(200.0, 90.0))


class TestBoundedNamesTheFloorTheCeilingAndTheAbsentMeasurement:
    """The audit's own definition, from `stress_panel.py`'s H2O2 target EC50."""

    def test_the_value_must_lie_inside_its_own_bracket(self):
        with pytest.raises(ValueError, match="outside its own bracket"):
            Param.bounded("ec50", 0.5, "mM", "Kritsiligkou 2021, PMID 34118234",
                          bounds=(0.02, 0.15), missing="a HyPer7 half-max in yeast")

    def test_a_bracket_with_no_way_out_of_it_is_refused(self):
        with pytest.raises(ValueError, match="two out of three"):
            Param(0.05, "mM", Tag.BOUNDED, "Kritsiligkou 2021, PMID 34118234",
                  name="ec50", bounds=(0.02, 0.15))

    def test_it_still_has_a_float_and_still_costs_a_degree_of_freedom(self):
        param = Param.bounded("ec50", 0.05, "mM", "Kritsiligkou 2021, PMID 34118234 floor",
                              bounds=(0.02, 0.15), missing="a HyPer7 half-max in yeast")

        assert float(param) == 0.05
        assert param.is_free


class TestBorrowedCarriesItsOrganismOrItIsNotConstructible:
    """MEASURED-dagger in `ARCHITECTURE_TARGET.md`: import only with the flag attached."""

    def test_the_flag_cannot_be_detached(self):
        with pytest.raises(ValueError, match="names no organism"):
            Param(25.0, "uS/cm2", Tag.BORROWED, "Sanders 1981", name="g_H")

    @pytest.mark.parametrize("organism", ["S. cerevisiae", "yeast",
                                          "Saccharomyces cerevisiae"])
    def test_borrowing_from_the_host_is_a_measured_row_graded_down_by_mistake(self, organism):
        with pytest.raises(ValueError, match="graded down by mistake"):
            Param.borrowed("g_H", 25.0, "uS/cm2", "Sanders 1981", organism=organism)

    def test_only_borrowed_carries_an_organism(self):
        with pytest.raises(ValueError, match="Only BORROWED"):
            Param(0.146, "dimensionless", Tag.MEASURED, "this repo", name="cv",
                  organism="Neurospora crassa")

    def test_a_borrowed_point_is_pinned_not_free(self):
        param = Param.borrowed("g_H", 25.0, "uS/cm2", "Sanders 1981",
                               organism="Neurospora crassa")

        assert float(param) == 25.0
        assert not param.is_free


class TestUnitsSourceAndTheNumberItself:
    def test_units_are_required_because_three_time_bases_collide_here(self):
        with pytest.raises(ValueError, match="units are required"):
            Param(1.0, "  ", Tag.ASSERTED, "somewhere", name="k")

        assert float(Param.asserted("_LETHAL_MULTIPLE", 3.0, "x EC50",
                                    "generator/stress_panel.py, no source")) == 3.0

    @pytest.mark.parametrize("tag,extra", [
        (Tag.MEASURED, {}),
        (Tag.ASSERTED, {}),
        (Tag.REFUSED, {"value": None, "reason": "none exists.", "missing": "a number."}),
    ])
    def test_source_is_required_for_every_grade_including_refused(self, tag, extra):
        kwargs = {"value": 1.0, **extra}
        with pytest.raises(ValueError, match="source is required"):
            Param(kwargs.pop("value"), "mM", tag, "   ", name="x", **kwargs)

    @pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
    def test_a_non_finite_value_is_not_a_measurement(self, value):
        with pytest.raises(ValueError, match="must be finite"):
            Param.measured("x", value, "mM", "somewhere")

    def test_a_ci_that_does_not_contain_its_value_is_refused(self):
        with pytest.raises(ValueError, match="does not contain its value"):
            Param.measured("k_ref", 0.0798, "1/h", "Branco 2004",
                           ci95=(0.09, 0.11))

    def test_the_repositorys_own_scanner_finds_the_identifier_in_a_source(self):
        param = Param.bounded(
            "h2o2_target_ec50", 0.05, "mM",
            "Kritsiligkou 2021, PMID 34118234, HyPer7 in BY4742: a 20 uM floor under this "
            "module's own sensor, bracketed against a 0.15 mM regulon EC50",
            bounds=(0.02, 0.15), missing="a HyPer7 half-max in yeast")

        assert param.identifiers == ("pmid:34118234",)
        assert a_measured().identifiers == ()


class TestTheRegistry:
    def test_add_hands_the_param_back_so_a_module_constant_reads_normally(self):
        registry = ParamRegistry("mech.demo")

        cv = registry.add(a_measured())

        assert float(cv) == OBSERVED_ACTIVITY_CV
        assert registry["OBSERVED_ACTIVITY_CV"] is cv
        assert "OBSERVED_ACTIVITY_CV" in registry
        assert len(registry) == 1

    def test_an_anonymous_param_cannot_be_registered(self):
        with pytest.raises(ValueError, match="needs a name"):
            ParamRegistry("mech.demo").add(Param(1.0, "mM", Tag.ASSERTED, "somewhere"))

    def test_registering_one_constant_twice_is_how_copies_drift_apart(self):
        registry = ParamRegistry("mech.demo")
        registry.add(a_measured())

        with pytest.raises(ValueError, match="already registered"):
            registry.add(a_measured())

    def test_a_missing_name_lists_what_is_there(self):
        registry = ParamRegistry("mech.demo")
        registry.add(a_measured())

        with pytest.raises(KeyError, match="OBSERVED_ACTIVITY_CV"):
            registry["beta"]

    def test_the_census_covers_every_tag_and_sums_to_the_registry(self):
        registry = ParamRegistry("mech.demo")
        for param in (a_measured(), a_sweep(), a_refusal()):
            registry.add(param)

        counts = registry.by_tag()

        assert set(counts) == set(Tag.ALL)
        assert sum(counts.values()) == len(registry) == 3
        assert counts[Tag.MEASURED] == counts[Tag.SWEPT] == counts[Tag.REFUSED] == 1
        assert counts[Tag.ASSERTED] == 0

    def test_uncited_finds_a_cited_grade_with_no_identifier_and_ignores_asserted(self):
        registry = ParamRegistry("mech.demo")
        registry.add(a_measured())
        registry.add(Param.bounded("ec50", 0.05, "mM", "Kritsiligkou 2021, PMID 34118234",
                                   bounds=(0.02, 0.15),
                                   missing="a HyPer7 half-max in yeast"))
        registry.add(Param.asserted("_LETHAL_MULTIPLE", 3.0, "x EC50",
                                    "generator/stress_panel.py, no source"))

        names = [p.name for p in registry.uncited()]

        assert names == ["OBSERVED_ACTIVITY_CV"]

    def test_borrowed_and_refusals_are_reportable_as_sets(self):
        registry = ParamRegistry("mech.demo")
        registry.add(Param.borrowed("g_H", 25.0, "uS/cm2", "Sanders 1981",
                                    organism="Neurospora crassa"))
        registry.add(a_refusal())

        assert [p.organism for p in registry.borrowed()] == ["Neurospora crassa"]
        assert [p.name for p in registry.refusals()] == ["Vmax_ex"]

    def test_a_registry_needs_the_name_of_its_piece(self):
        with pytest.raises(ValueError, match="needs the name"):
            ParamRegistry("  ")


def ph_slice() -> ParamRegistry:
    """`ARCHITECTURE_GAPS.md` 0.3's pH slice, transcribed. Sources are the target
    document's own DEPTH 2 table; nothing here is a new number."""
    registry = ParamRegistry("mech.ph")
    registry.add(Param.swept(
        "beta", "mM/pH",
        "Gabba 2020, 160-320 mM phosphate -> 87-175 mM/pH at pH 7; Gerber "
        "2016 reports 200 in yeast; Kahm's 200 traces to Neurospora",
        bounds=(90.0, 200.0)))
    registry.add(Param.swept(
        "Jmax_pump", "mmol ATP/gDCW/h",
        "anchored on yeast physiology: 2x Abbott's delta q_glucose at half yield to "
        "Verduyn 1992's 8x qO2 rise under benzoate",
        bounds=(12.0, 24.0)))
    registry.add(Param.swept(
        "pKa_pma1", "pH",
        "Zhao 2021 established the mechanism at two pH points; no "
        "titration exists",
        bounds=(6.0, 7.4), missing="a Pma1 activity titration across pH."))
    registry.add(Param.swept(
        "n", "dimensionless", "Hill coefficient of the same autoinhibition; unmeasured",
        bounds=(1.0, 4.0)))
    registry.add(Param.refused(
        "Vmax_ex", "mmol/gDCW/h", "searched for Tpo2/Tpo3/Pdr12 capacity",
        reason="no measured yeast anion-export capacity exists.",
        missing="an export rate at a known intracellular anion content."))
    registry.add(Param.refused(
        "K_ex", "mmol/gDCW", "searched for Tpo2/Tpo3/Pdr12 capacity",
        reason="the same absent measurement fixes both constants.",
        missing="a second export rate at a different anion content."))
    registry.add(Param.swept(
        "g_H", "uS/cm2",
        "Sanders 1981 measured 25 uS/cm2 in Neurospora crassa and no yeast value was "
        "found. GAPS 0.3 counts this axis swept but states no band, so the endpoints "
        "here are an ASSERTED decade either side of the borrowed point and exist only "
        "to make this fixture constructible",
        bounds=(2.5, 250.0), missing="a yeast proton conductance."))
    registry.add(Param.swept(
        "V_CYT_L_PER_GDCW", "L/kgDCW",
        "BNIDs 112806 and 111353; not measured for these strains",
        bounds=(1.7, 2.0)))

    registry.add_target(Target(
        name="E4_lactate",
        observable="growth rate across a lactate ladder at matched undissociated "
                   "concentration",
        assay="OD600 kinetic read; growth rate as the slope of log OD against time",
        noise_floor=MEASURED_GROWTH_RATE_SE,
        units="1/h",
        source="Gabba 2020, Table 2 measured P_lactic = 0, so "
               "Henderson-Hasselbalch fails by 9x for free"))
    registry.add_target(Target(
        name="E6_orij_slope",
        observable="growth rate against cytosolic pH",
        assay="OD600 kinetic read; growth rate as the slope of log OD against time",
        noise_floor=MEASURED_GROWTH_RATE_SE,
        units="1/h",
        source="Orij 2009"))
    return registry


class TestTheGapsDocumentsPhSliceCountsEight:
    """A calibration, not a unit test. `ARCHITECTURE_GAPS.md` 0.3 counts the pH slice at
    eight free scalars against two usable targets and concludes the block adds mechanism
    without adding refutability. If the rule in :data:`FREE_TAGS` ever stops reproducing
    that eight, the rule has drifted from the document that motivated it."""

    def test_the_free_count_is_eight_and_the_names_are_the_documents_own(self):
        gate = ph_slice().gate()

        assert len(gate.free) == 8
        assert set(gate.free) == {"beta", "Jmax_pump", "pKa_pma1", "n", "Vmax_ex", "K_ex",
                                  "g_H", "V_CYT_L_PER_GDCW"}

    def test_the_target_count_is_two(self):
        assert len(ph_slice().gate().targets) == 2

    def test_eight_against_two_does_not_pass_and_the_refusal_lists_all_eight(self):
        registry = ph_slice()

        assert registry.gate().passes is False
        with pytest.raises(FreeScalarGateFailed) as excinfo:
            registry.require_gate()
        message = str(excinfo.value)

        assert "8 free scalars against 2 independent targets" in message
        for name in registry.gate().free:
            assert name in message

    def test_a_fitted_target_does_not_pay_for_a_sweep(self):
        registry = ph_slice()
        for i in range(6):
            registry.add_target(Target(
                name=f"fitted_{i}", observable="a swept point that matched",
                assay="OD600 kinetic read", noise_floor=MEASURED_GROWTH_RATE_SE,
                units="1/h", source="scored by searching the sweep", fitted=True))

        gate = registry.gate()

        assert len(gate.targets) == 2
        assert len(gate.fitted_targets) == 6
        assert gate.passes is False
        assert "not counted, fitted" in gate.summary()

    def test_the_gate_passes_when_the_targets_outnumber_the_free_scalars(self):
        registry = ParamRegistry("mech.demo")
        registry.add(a_sweep())
        registry.add_target(Target("t1", "growth rate", "OD600 kinetic read",
                                   MEASURED_GROWTH_RATE_SE, "1/h", "Orij 2009"))

        gate = registry.require_gate()

        assert gate.passes and "PASSES" in gate.summary()

    def test_a_piece_with_no_free_scalars_passes_with_no_targets(self):
        registry = ParamRegistry("mech.params")
        registry.add(a_measured())

        assert registry.require_gate().passes

    def test_the_two_tag_sets_partition_the_seven(self):
        assert set(PINNED_TAGS) | set(FREE_TAGS) == set(Tag.ALL)
        assert not set(PINNED_TAGS) & set(FREE_TAGS)


class TestTheNoiseFloorsAreTheRepositorysOwn:
    """Criterion (a) names an assay and its floor. A target with an invented floor is the
    error `REVISED_BUILD_LIST.md` records against Route C2 -- a plate-reader activity CV
    used to score an FBA product ceiling."""

    def test_a_target_needs_an_assay_and_a_positive_floor(self):
        with pytest.raises(ValueError, match="assay is required"):
            Target("t", "growth rate", "  ", MEASURED_GROWTH_RATE_SE, "1/h", "Orij 2009")
        with pytest.raises(ValueError, match="noise_floor must be positive"):
            Target("t", "growth rate", "OD600 kinetic read", 0.0, "1/h", "Orij 2009")

    def test_the_two_measured_floors_are_the_ones_the_gate_examples_use(self):
        assert OBSERVED_ACTIVITY_CV == 0.146
        assert MEASURED_GROWTH_RATE_SE == 0.0117
        assert all(t.noise_floor == MEASURED_GROWTH_RATE_SE for t in ph_slice().targets)

    def test_one_measurement_cannot_be_counted_twice(self):
        registry = ph_slice()

        with pytest.raises(ValueError, match="already registered"):
            registry.add_target(Target("E4_lactate", "the same thing again",
                                       "OD600 kinetic read", MEASURED_GROWTH_RATE_SE,
                                       "1/h", "Gabba 2020"))


class TestTheProvenanceTable:
    def test_every_parameter_reaches_the_table_with_its_tag_and_its_freedom(self):
        table = provenance_table(ph_slice())

        assert "mech.ph" in table
        assert "8 free scalars against 2 independent targets -- REFUSED" in table
        for param in ph_slice().params:
            assert f"`{param.name}`" in table
            assert f"**{param.tag}**" in table

    def test_a_refusal_renders_as_a_dash_and_never_as_a_number(self):
        table = provenance_table(ph_slice())
        row = next(line for line in table.splitlines() if "`Vmax_ex`" in line)

        assert "| -- |" in row
        assert "no measured yeast anion-export capacity exists." in row
        assert "would close it:" in row

    def test_a_swept_axis_renders_as_its_band(self):
        table = provenance_table(ph_slice())
        row = next(line for line in table.splitlines() if "`beta`" in line)

        assert "90 - 200" in row

    def test_the_targets_carry_their_assay_and_floor_into_the_table(self):
        table = provenance_table(ph_slice())

        assert "OD600 kinetic read" in table
        assert "0.0117 1/h" in table

    def test_a_pipe_in_a_source_does_not_break_the_table(self):
        registry = ParamRegistry("mech.demo")
        registry.add(Param.asserted("x", 1.0, "mM", "a | b"))

        row = next(line for line in provenance_table(registry).splitlines()
                   if "`x`" in line)

        assert r"a \| b" in row
        assert row.replace(r"\|", "").count("|") == 7  # six columns, seven delimiters

    def test_several_registries_render_in_one_document(self):
        other = ParamRegistry("mech.oxidative")
        other.add(Param.measured("k_ref", 0.0798, "1/h",
                                 "Branco 2004, refit by Selvaggio 2018"))

        table = provenance_table(ph_slice(), other)

        assert "### mech.ph" in table and "### mech.oxidative" in table

    def test_the_table_needs_a_registry(self):
        with pytest.raises(ValueError, match="at least one registry"):
            provenance_table()


class TestTheTagCensus:
    def test_it_has_a_column_per_tag_and_reports_the_asserted_share(self):
        registry = ParamRegistry("mech.demo")
        registry.add(a_measured())
        registry.add(Param.asserted("_LETHAL_MULTIPLE", 3.0, "x EC50",
                                    "generator/stress_panel.py, no source"))
        registry.add(Param.asserted("c", 0.1, "dimensionless", "convention"))

        table = tag_census_table(registry, ph_slice())
        header, _, demo, ph = table.splitlines()

        for tag in Tag.ALL:
            assert tag in header
        assert demo.endswith("| 2 | 67% |")
        assert ph.endswith("| 8 | 0% |")

    def test_it_needs_a_registry(self):
        with pytest.raises(ValueError, match="at least one registry"):
            tag_census_table()
