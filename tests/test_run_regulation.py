"""The D3 regulation script: what it refuses without a model, and what it declares.

``scripts/run_regulation.py`` writes five ``regulation_*`` tables and is the source of the
D3 sections of ``docs/FINDINGS.md``. It is one long ``main`` around the GECKO model, so the
properties worth pinning are the ones that hold before any LP is solved:

* with no enzyme-constrained model it exits naming the variable, and writes nothing. A
  partial ``regulation_product_range.csv`` would be read as a result;
* the modules it declares unmappable really have no regulon, and every module its probes
  drive really has one. ``eflux_layer`` raises on any module with activity and no regulon,
  so a drift in either direction turns into a crash halfway through a long run -- after
  ``regulation_module_reach.csv`` has already been written;
* an attainable capacity is never negative, because a negative one would flip an E-Flux
  upper bound below its lower bound and make the model infeasible for arithmetic reasons.

The GSMM sections are covered by ``tests/test_regulation.py`` against ``bridge/regulation``
itself; nothing here loads a genome-scale model.
"""

from __future__ import annotations

import importlib.util
import json
import os
import pathlib

import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "run_regulation.py"


@pytest.fixture(scope="module")
def outputs(tmp_path_factory):
    """Where the script is allowed to write during these tests, and nowhere else."""
    return tmp_path_factory.mktemp("regulation_outputs")


@pytest.fixture(scope="module")
def script(outputs):
    """The script as a module, with its output directory redirected away from outputs/.

    ``OUT`` binds at import time, so redirecting after the import would not help: the
    refusal test below calls ``main``, and on a machine that does have the GECKO model a
    stray success would overwrite five tracked tables.
    """
    previous = os.environ.get("YSTWIN_OUTPUTS")
    os.environ["YSTWIN_OUTPUTS"] = str(outputs)
    try:
        spec = importlib.util.spec_from_file_location("run_regulation", _SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if previous is None:
            os.environ.pop("YSTWIN_OUTPUTS", None)
        else:
            os.environ["YSTWIN_OUTPUTS"] = previous
    return module


@pytest.fixture
def regulons():
    """The vendored SGD regulons, or a skip naming what is missing."""
    from ystwin.bridge.regulation import load_regulons

    try:
        return load_regulons()
    except FileNotFoundError as exc:
        pytest.skip(f"SGD regulation records not present: {exc}")


class TestNoModelMeansNoTable:
    """The five tables are the D3 result. Written from a run that could not load the model
    they would be worse than absent, because a stale one on disk reads as a fresh one."""

    def test_an_absent_gecko_model_exits_naming_the_variable(self, script, monkeypatch,
                                                             tmp_path):
        """``YSTWIN_EC_YEAST_GEM`` set to somewhere the model is not resolves to nothing
        rather than falling through to a default -- that is the whole point of the
        override -- so this is also the path a machine without the asset takes."""
        monkeypatch.setenv("YSTWIN_EC_YEAST_GEM", str(tmp_path / "no_such_model.xml"))

        with pytest.raises(SystemExit) as refusal:
            script.main(["--output-dir", str(tmp_path / "out")])

        assert "YSTWIN_EC_YEAST_GEM" in str(refusal.value)

    def test_the_refusal_leaves_no_partial_table_behind(self, script, monkeypatch,
                                                        tmp_path, outputs):
        monkeypatch.setenv("YSTWIN_EC_YEAST_GEM", str(tmp_path / "no_such_model.xml"))

        with pytest.raises(SystemExit):
            script.main(["--output-dir", str(tmp_path / "out")])

        assert not list(outputs.glob("regulation_*.csv"))


class TestTheUnmappableModulesAreDeclaredAndCorrect:
    """``eflux_layer`` insists that a module with activity and no regulon be declared
    rather than filtered, and ``UNMAPPED`` is that declaration. It is a literal, so it can
    go stale in both directions: too small and a long run dies partway through, too large
    and a module that has since acquired a regulon is silently excluded from the layer."""

    def test_every_declared_module_really_has_no_regulon(self, script, regulons):
        """The five metabolite pools have no transcription factor at all; calcium has Crz1
        but no target genes at the MacIsaac evidence level. Either way the declaration has
        to be true of the records on disk."""
        for module in script.UNMAPPED:
            regulon = regulons.get(module)
            assert regulon is None or not regulon.genes, f"{module} now has a regulon"

    def test_every_module_the_probes_drive_is_either_mapped_or_declared(self, script,
                                                                       regulons):
        """This is the condition ``eflux_layer`` raises on, checked against the four probe
        agents the script actually doses."""
        from ystwin.generator.stress_panel import module_response

        for stressor, dose in script.PROBES:
            for module, value in module_response(stressor, dose).items():
                if value == 0.0 or module in script.UNMAPPED:
                    continue
                regulon = regulons.get(module)
                assert regulon and regulon.genes, (
                    f"{stressor} drives {module}, which has no regulon and is not declared")

    def test_the_pools_are_a_subset_of_the_declaration(self, script):
        assert script.POOLS <= script.UNMAPPED

    def test_every_declared_module_is_a_module(self, script):
        """A name that is not a module would sit in ``allow_unmapped`` doing nothing, and
        the module it was meant to excuse would still raise."""
        from ystwin.generator.stress_panel import MODULES

        assert script.UNMAPPED <= set(MODULES)


class TestTheProbesAreRealAgentsAtRealDoses:
    """One agent per mechanism the panel covers, each dosed at its own EC50 so the
    activities are mid-scale rather than saturated."""

    def test_every_probe_agent_exists(self, script):
        from ystwin.generator.stress_panel import STRESSORS

        assert {name for name, _ in script.PROBES} <= set(STRESSORS)

    def test_every_probe_is_dosed_at_its_own_ec50(self, script):
        """A dose that drifted off the EC50 would put the layer at a saturated activity,
        where the scale is 1.0 for every gene and the layer stops distinguishing modules."""
        from ystwin.generator.stress_panel import STRESSORS

        for name, dose in script.PROBES:
            assert dose == pytest.approx(STRESSORS[name].ec50), name

    def test_every_probe_moves_at_least_one_module(self, script):
        from ystwin.generator.stress_panel import module_response

        for name, dose in script.PROBES:
            assert any(v != 0.0 for v in module_response(name, dose).values()), name


class TestTheGrowthGridIsD1s:
    """The before/after columns are only comparable to the published D1 table because the
    growth fractions are the same ones. They are copied, so they can drift."""

    def test_the_grid_runs_from_maximum_growth_downward(self, script):
        assert script.GROWTH_FRACTIONS[0] == 1.0
        assert list(script.GROWTH_FRACTIONS) == sorted(script.GROWTH_FRACTIONS,
                                                       reverse=True)

    @pytest.mark.integration
    def test_the_tracked_table_was_swept_over_this_grid(self, script):
        """The tracked ``regulation_product_range.csv`` is quoted as the D3 result. If the
        constant here and the fractions in that file disagree, the file predates the
        constant and the numbers being read are from a different sweep."""
        import pandas as pd

        from ystwin import paths

        path = paths.outputs_dir() / "regulation_product_range.csv"
        if not path.exists():
            pytest.skip(f"{path.name} not present; run scripts/run_regulation.py")
        swept = set(pd.read_csv(path).growth_fraction.round(6))

        assert swept == {round(f, 6) for f in script.GROWTH_FRACTIONS}


class TestAnAttainableCapacityIsNeverNegative:
    """The tighter of the two E-Flux references, computed rather than asserted. It is an
    FVA maximum, and an irreversible-backwards reaction has a negative one -- which as a
    capacity would put an upper bound below the lower bound and make the model infeasible
    for a reason that has nothing to do with regulation."""

    @staticmethod
    def _toy_model(script):
        import cobra

        model = cobra.Model("toy")
        a = cobra.Metabolite("a_c")
        supply = cobra.Reaction("EX_a", lower_bound=0.0, upper_bound=10.0)
        supply.add_metabolites({a: 1.0})
        biomass = cobra.Reaction(script.BIOMASS, lower_bound=0.0, upper_bound=1000.0)
        biomass.add_metabolites({a: -1.0})
        # No metabolites and a wholly negative range: always feasible, and its FVA
        # maximum is -1.
        backwards = cobra.Reaction("only_backwards", lower_bound=-5.0, upper_bound=-1.0)
        model.add_reactions([supply, biomass, backwards])
        model.objective = script.BIOMASS
        return model

    def test_a_backwards_only_reaction_reports_zero_rather_than_a_negative_capacity(
            self, script):
        model = self._toy_model(script)

        capacities = script._attainable_capacities(model, ["only_backwards"])

        assert capacities["only_backwards"] == 0.0

    def test_a_forward_reaction_reports_its_own_maximum(self, script):
        """The counterpart, so the clamp above is about the sign and not about the
        function returning zero for everything."""
        model = self._toy_model(script)

        capacities = script._attainable_capacities(model, ["EX_a"])

        assert capacities["EX_a"] > 0.0

    def test_the_capacity_is_read_at_the_stated_growth_fraction(self, script):
        """The default is 90% of maximum growth, which is what the counterfactual section
        says it computed. At a fraction that low the supply reaction is still pinned near
        its maximum, and a function that ignored the argument would return the same number
        at every fraction."""
        model = self._toy_model(script)

        loose = script._attainable_capacities(model, ["EX_a"], growth_fraction=0.0)
        tight = script._attainable_capacities(model, ["EX_a"], growth_fraction=1.0)

        assert tight["EX_a"] >= loose["EX_a"] - 1e-9
        assert tight["EX_a"] == pytest.approx(10.0)


class TestTheCitedMechanismsAreNamedByPmid:
    """Both mechanisms are printed with their PMIDs in the banners, and the banners are
    what a reader checks the method against."""

    def test_the_eflux_and_cosmic_pmids_come_from_the_module_that_implements_them(
            self, script):
        from ystwin.bridge import regulation

        assert script.COLIJN_PMID == regulation.COLIJN_PMID
        assert script.COSMIC_PMID == regulation.COSMIC_PMID

    def test_eflux_is_recorded_as_unable_to_lift_a_floor(self, script):
        """The whole D3 answer rests on this: a capacity constraint bounds |v| and never
        requires flux, so no E-Flux layer can move the product floor off zero. If this
        constant ever reads True the script's section 3 conclusion is wrong."""
        assert script.EFLUX_LIFTS_LOWER_BOUND is False


@pytest.mark.parametrize("widths, expected", [
    pytest.param([(1.0, 1.0)] * 24 + [(float("nan"), 1.0)] * 4,
                 "0 of 24 defined comparisons; 4 undefined", id="registered-availability"),
    pytest.param([(1.0, 0.5), (1.0, 1.0), (1.0, 1.5), (0.0, 0.0), (1.0, 0.0)],
                 "2 of 5 defined comparisons; 0 undefined", id="finite-including-zero"),
    pytest.param([(1.0, 1.0 - 1e-9), (1.0, 1.0 - 0.5e-9), (1.0, 1.0 - 2e-9)],
                 "1 of 3 defined comparisons; 0 undefined", id="unchanged-threshold"),
    pytest.param([(float("nan"), 1.0), (1.0, float("nan")),
                  (float("nan"), float("nan")), (float("inf"), 1.0),
                  (1.0, float("inf")), (-float("inf"), 1.0), (1.0, -float("inf"))],
                 "0 of 0 defined comparisons; 7 undefined", id="nonfinite-either-arm"),
    pytest.param([], "0 of 0 defined comparisons; 0 undefined", id="empty"),
])
def test_relative_width_summary_counts_only_finite_pairs(script, capsys, widths, expected):
    import pandas as pd

    # Width availability, not stressor names or growth fractions, defines the denominator.
    frame = pd.DataFrame(widths, columns=["relative_width_unregulated", "relative_width_regulated"],
                         dtype=float)
    before = frame.copy(deep=True)

    script._print_relative_width_summary(frame)

    assert capsys.readouterr().out == (
        f"\n  rows where regulation narrowed the relative width: {expected}\n")
    pd.testing.assert_frame_equal(frame, before)


@pytest.mark.parametrize("maximum", [0.0, -1e-14])
def test_nonpositive_ceiling_widths_remain_undefined_in_summary(script, capsys, maximum):
    import pandas as pd

    from ystwin.fba.fva import FluxRange

    width = FluxRange(reaction="product", minimum=0.0, maximum=maximum, max_growth=1.0,
                      enforced_growth=0.5, glucose_uptake=1.0, pathway_capacity=None).relative_width
    frame = pd.DataFrame({"relative_width_unregulated": [width, 1.0],
                          "relative_width_regulated": [1.0, width]})

    script._print_relative_width_summary(frame)

    assert capsys.readouterr().out == (
        "\n  rows where regulation narrowed the relative width: "
        "0 of 0 defined comparisons; 2 undefined\n")


def test_measured_runner_keeps_all_tasks_and_conditions_before_ranking(script, monkeypatch):
    import cobra

    model = cobra.Model("condition-toy")
    model.add_reactions([cobra.Reaction(script.BIOMASS)])
    calls = []

    def constrain(m, exchange, magnitude, **kwargs):
        calls.append((exchange, magnitude, kwargs))

    def refuse(m, tasks, **kwargs):
        assert {task.name for task in tasks} == {"biomass", "ethanol", "CO2"}
        assert len(calls) == 2
        assert {row[0] for row in calls} == {"r_1714", "r_1992"}
        raise script.InfeasibleRegulation("full set infeasible")

    monkeypatch.setattr(script, "constrain_uptake", constrain)
    monkeypatch.setattr(script, "task_priority_order", refuse)
    result = script.measured_task_report(model)
    assert result["status"] == "infeasible"
    assert result["tasks"] == ["biomass", "ethanol", "CO2"]
    assert result["order"] == []
    assert len(calls) == 2


def test_matched_runner_refuses_plain_models_in_both_comparison_arms(script):
    import cobra

    plain = cobra.Model("yeastGEM_v9__46__0__46__2")
    ec = cobra.Model("yeastGEM_v9__46__0__46__2")
    result = script.matched_task_comparison(plain, ec, strain_transfer="synthetic check")
    assert result["status"] == "incompatible"
    assert "enzyme" in result["refusal"]


def test_new_runner_requires_an_explicit_output_destination(script):
    with pytest.raises(SystemExit):
        script.parse_args([])


def test_task_report_uses_portable_full_model_paths_and_exact_hashes(script, tmp_path, monkeypatch):
    from ystwin.artifacts import sha256

    plain = _REPO / "data/gem/yeast-GEM.xml.gz"
    ec = _REPO / "data/gem/ecYeastGEM_yeast902.xml.gz"
    monkeypatch.setattr(script, "load_model", lambda path: (object(), "GLPK"))
    monkeypatch.setattr(script, "matched_task_comparison", lambda *args, **kwargs: {
        "status": "completed", "models": {}, "biological_validation": False})
    output = tmp_path / "report"
    script.main(["--tasks-only", "--plain-model", str(plain), "--ec-model", str(ec),
                 "--output-dir", str(output)])
    record = json.loads((output / "regulation_matched_tasks.json").read_text())
    assert record["model_paths"] == {"plain": "data/gem/yeast-GEM.xml.gz",
                                     "ec": "data/gem/ecYeastGEM_yeast902.xml.gz"}
    assert record["model_sha256"] == {"plain": sha256(plain), "ec": sha256(ec)}
