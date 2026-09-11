"""The scenario script: what it refuses without a model, and the names it indexes by.

``scripts/run_scenarios.py`` writes ``scenario_predictions.csv``, whose latent rows are
refused by G4 and must never be quoted. Its loop indexes four separate tables by construct
and by reaction id, and every one of those lookups is a plain dict access -- so a name that
drifts is a KeyError minutes into a run, or worse a row that silently never gets built.

The one GSMM test skips when the model is absent, including when its variable has been
aimed somewhere the model is not.
"""

from __future__ import annotations

import importlib.util
import os
import pathlib

import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "run_scenarios.py"


@pytest.fixture(scope="module")
def outputs(tmp_path_factory):
    return tmp_path_factory.mktemp("scenario_outputs")


@pytest.fixture(scope="module")
def gsmm():
    """Yeast9, or a skip naming the variable that would find it.

    Deliberately not the shared ``yeast_gem`` fixture: that one skips on ``None`` only, and
    ``paths._resolve`` returns whatever an override points at as long as it exists -- so an
    override aimed at a directory resolves to the directory and the load raises "the SBML
    model is not valid" instead of skipping. Pointing the variable at an empty directory is
    how this suite is meant to be run without its assets, so the guard belongs here.
    """
    from ystwin import paths

    path = paths.yeast_gem()
    if path is None or not path.is_file():
        pytest.skip("yeast-GEM not present; set YSTWIN_YEAST_GEM to yeast-GEM.xml")
    import cobra

    return cobra.io.read_sbml_model(str(path))


@pytest.fixture(scope="module")
def script(outputs):
    """The script as a module, with its output directory redirected away from outputs/.

    ``OUT`` binds at import time and ``main`` writes ``scenario_predictions.csv`` into it
    unconditionally, so the redirect has to be in place before the import.
    """
    previous = os.environ.get("YSTWIN_OUTPUTS")
    os.environ["YSTWIN_OUTPUTS"] = str(outputs)
    try:
        spec = importlib.util.spec_from_file_location("run_scenarios", _SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if previous is None:
            os.environ.pop("YSTWIN_OUTPUTS", None)
        else:
            os.environ["YSTWIN_OUTPUTS"] = previous
    return module


class TestNoModelMeansNoTable:
    """``scenario_predictions.csv`` is written in one call after the whole sweep, so the
    only way a partial one appears is if the refusal comes late. It comes first."""

    def test_an_absent_gsmm_exits_naming_the_variable(self, script, monkeypatch,
                                                      tmp_path):
        monkeypatch.setenv("YSTWIN_YEAST_GEM", str(tmp_path / "no_such_model.xml"))

        with pytest.raises(SystemExit) as refusal:
            script.main()

        assert "YSTWIN_YEAST_GEM" in str(refusal.value)

    def test_the_refusal_happens_before_anything_is_written(self, script, monkeypatch,
                                                            tmp_path, outputs):
        monkeypatch.setenv("YSTWIN_YEAST_GEM", str(tmp_path / "no_such_model.xml"))

        with pytest.raises(SystemExit):
            script.main()

        assert not (outputs / "scenario_predictions.csv").exists()


class TestTheReactionsAreNamedOnceAndOnlyOnce:
    """Each row is built as ``{NAMES[k]: r[f"flux_{k}"] for k in REACTIONS}``, so the two
    constants have to agree exactly: a reaction with no name is a KeyError, and a name with
    no reaction is a column that never appears."""

    def test_every_reported_reaction_has_a_name(self, script):
        assert set(script.REACTIONS) <= set(script.NAMES)

    def test_every_name_belongs_to_a_reported_reaction(self, script):
        assert set(script.NAMES) <= set(script.REACTIONS)

    def test_no_two_reactions_share_a_column_name(self, script):
        """Two ids mapping to one name would silently drop a flux: the later key wins and
        the earlier reaction disappears from the table with no error."""
        assert len(set(script.NAMES.values())) == len(script.NAMES)

    def test_the_printed_summary_columns_are_among_the_named_ones(self, script):
        """The console table prints glucose and oxygen per branch. Both have to survive the
        rename above or the summary prints a header over nothing."""
        assert {"glucose", "oxygen"} <= set(script.NAMES.values())

    @pytest.mark.integration
    def test_every_reported_reaction_is_in_the_model(self, script, gsmm):
        """Four of the five are exchange reactions and one is biomass. An id that is not in
        Yeast9 raises inside ``predict_flux``, after the plate has been generated."""
        present = {reaction.id for reaction in gsmm.reactions}

        assert set(script.REACTIONS) <= present


class TestEveryConstructOnThePlateCanBeDosed:
    """The loop reads ``STRESSOR_FOR_CONSTRUCT[construct]`` and then
    ``DEFAULT_DOSES[stressor]``. Both are dict accesses on a panel defined elsewhere."""

    def test_every_panel_construct_has_a_stressor(self, script):
        from ystwin.generator.plate import DEFAULT_PANEL
        from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT

        assert set(DEFAULT_PANEL) <= set(STRESSOR_FOR_CONSTRUCT)

    def test_every_panel_stressor_has_a_dose_ladder(self, script):
        from ystwin.generator.plate import DEFAULT_DOSES, DEFAULT_PANEL
        from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT

        for construct in DEFAULT_PANEL:
            assert STRESSOR_FOR_CONSTRUCT[construct] in DEFAULT_DOSES

    def test_every_ladder_starts_at_an_unstressed_rung(self, script):
        """The zero dose is the control every other rung is read against. A ladder without
        one gives a scenario table with no baseline row for that construct."""
        from ystwin.generator.plate import DEFAULT_DOSES, DEFAULT_PANEL
        from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT

        for construct in DEFAULT_PANEL:
            ladder = DEFAULT_DOSES[STRESSOR_FOR_CONSTRUCT[construct]]
            assert min(ladder) == 0.0


class TestTheLatentBranchIsMarkedNotQuotable:
    """Every latent number the script prints is refused by G4. The refusal travels with the
    rows rather than living only in the banner, because the CSV outlives the console."""

    def test_the_measured_branch_is_named_the_string_the_summary_selects_on(self, script):
        """The summary selects ``sub.branch == "physiology"`` per dose and prints nothing
        at all where that selection is empty -- so a rename in the bridge would empty the
        console table and the pivot below it without raising anywhere."""
        from ystwin.bridge.physiology_bridge import MeasuredState, physiology_constraints

        constraints = physiology_constraints(MeasuredState(
            biomass_gl=1.0, growth_rate=0.2, glucose_mM=20.0, time_h=10.0))

        assert constraints.branch == "physiology"

    def test_the_latent_branch_is_not_load_bearing_by_construction(self, script):
        """``compare_branches`` marks it, and the script copies the flag into every row
        rather than deciding for itself."""
        from ystwin.bridge import latent_bridge

        assert "REFUSED by G4" in latent_bridge._G4_STATUS


class TestTheLatentBranchIsCapableOfDifferingFromThePhysiologyBranch:
    """The bug this class exists for was not a wrong number. It was a branch that could
    not produce one.

    ``latent_constraints`` normalises promoter activity by ``_ACTIVITY_FULL_SCALE``, which
    is 1305 **RFU/OD/h**. The script used to divide its recovered activity by
    ``gdcw_per_od * gain`` first, handing over a synthesis rate around 1e-3. The stress
    fraction was then about 3e-7, the latent NGAM came out 0.7000002 against a resting 0.7,
    and the script's closing line reported "glucose uptake shifts by -0.0% to +0.0%" --
    which reads as a finding about the latent layer and was arithmetic about units.

    Everything else in the file passed throughout. The name lookups were right, the
    refusal path was right, the G4 label was right. What nothing asserted was that the
    branch *moves*, and a branch that cannot move is the one failure mode a comparison
    script has no defence against: its output looks like a clean negative result.
    """

    @pytest.fixture(scope="class")
    def activities(self, script):
        """The corrected activity at every rung of one construct's ladder, RFU/OD/h."""
        from ystwin.generator.plate import (
            DEFAULT_DOSES,
            DEFAULT_PANEL,
            PlateConditions,
            generate_plate,
        )
        from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT

        conditions = PlateConditions(seed=42)
        plate = generate_plate(DEFAULT_PANEL, conditions)
        construct = next(iter(DEFAULT_PANEL))
        ladder = DEFAULT_DOSES[STRESSOR_FOR_CONSTRUCT[construct]]
        return [script._corrected_activity(plate, conditions, plate.times_h, construct, i)
                for i in range(len(ladder))]

    def test_the_recovered_activity_is_on_the_scale_the_bridge_normalises_by(self,
                                                                             activities):
        """The direct statement of the unit. ``_ACTIVITY_FULL_SCALE`` is the measured
        dynamic range of these constructs, so an activity fed to the bridge belongs within
        an order of magnitude of it -- not four below."""
        from ystwin.bridge.latent_bridge import _ACTIVITY_FULL_SCALE

        assert max(activities) > 0.1 * _ACTIVITY_FULL_SCALE

    def test_the_top_dose_raises_maintenance_well_above_resting(self, script, activities):
        from ystwin.bridge.latent_bridge import (
            _RESTING_MAINTENANCE,
            LatentState,
            latent_constraints,
        )
        from ystwin.bridge.physiology_bridge import MeasuredState, physiology_constraints

        physical = physiology_constraints(MeasuredState(
            biomass_gl=1.0, growth_rate=0.05, glucose_mM=20.0, time_h=24.0))
        top = latent_constraints(
            LatentState(max(activities), "UPRE2", 24.0), physical,
            acknowledge_unvalidated=True, basal_activity=activities[0])

        assert top.atp_maintenance > 2.0 * _RESTING_MAINTENANCE

    def test_the_unstressed_rung_leaves_maintenance_exactly_at_resting(self, script,
                                                                       activities):
        """The other half, and the reason the fix is a *baseline* rather than a rescaling:
        excess is measured from the construct's own 0 mM rung, so an undosed culture
        carries no stress maintenance at all. Without this the same change that made the
        branch move would have made it move everywhere, including where nothing happened.
        """
        from ystwin.bridge.latent_bridge import (
            _RESTING_MAINTENANCE,
            LatentState,
            latent_constraints,
        )
        from ystwin.bridge.physiology_bridge import MeasuredState, physiology_constraints

        physical = physiology_constraints(MeasuredState(
            biomass_gl=1.0, growth_rate=0.4, glucose_mM=20.0, time_h=24.0))
        control = latent_constraints(
            LatentState(activities[0], "UPRE2", 24.0), physical,
            acknowledge_unvalidated=True, basal_activity=activities[0])

        assert control.atp_maintenance == pytest.approx(_RESTING_MAINTENANCE)

    def test_maintenance_never_leaves_the_range_yeast9_can_solve(self, script, activities):
        """The constant this replaced demanded an NGAM around 3e5. Yeast9 goes infeasible
        at 19.19 and growth is already zero there, so the old value sat about 2e4 past the
        point where the cell stops growing -- and an infeasible LP is a caught exception
        and a missing row, not a loud failure."""
        from ystwin.bridge.latent_bridge import LatentState, latent_constraints
        from ystwin.bridge.physiology_bridge import MeasuredState, physiology_constraints

        physical = physiology_constraints(MeasuredState(
            biomass_gl=1.0, growth_rate=0.05, glucose_mM=20.0, time_h=24.0))

        for activity in activities:
            got = latent_constraints(
                LatentState(activity, "UPRE2", 24.0), physical,
                acknowledge_unvalidated=True, basal_activity=activities[0])

            assert 0.0 < got.atp_maintenance < 19.19
