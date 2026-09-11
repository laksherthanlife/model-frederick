"""Do the audits detect what they claim to?

An audit that always passes is worse than no audit: it is a green light nobody has
checked, and it makes the next person trust the ones that do work. So every check in
``scripts/audit_reproducibility.py`` and ``scripts/audit_determinism.py`` is exercised
twice here -- once against a fixture that satisfies the invariant, once against a fixture
built to violate it -- and the violation has to be caught.

The fixtures are miniature repositories and miniature packages under ``tmp_path``. The
real audits are deliberately *not* run over the real repository: that would be slow, and
it would couple this suite to whatever happens to be in ``outputs/`` today, so a stale
table would show up as a broken test rather than as the audit finding it is meant to be.
"""
from __future__ import annotations

import importlib
import importlib.util
import json
import os
import pathlib
import subprocess
import sys

import numpy as np
import pytest

from ystwin import artifacts

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"


def _load(name: str):
    """Import an audit script by path; they live in ``scripts/`` and are not a package.

    Registered in ``sys.modules`` before execution because ``@dataclass`` resolves
    annotations by looking its own module up there, and a module absent from it fails at
    class-definition time rather than at use.
    """
    key = f"_audit_{name}"
    spec = importlib.util.spec_from_file_location(key, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[key] = module
    spec.loader.exec_module(module)
    return module


repro = _load("audit_reproducibility")
determinism = _load("audit_determinism")


def _status(rows, fragment: str) -> str:
    """The status of the one row whose check name contains ``fragment``."""
    matched = [r for r in rows if fragment in r["check"]]
    assert len(matched) == 1, f"expected one row matching {fragment!r}, got {[r['check'] for r in matched]}"
    return matched[0]["status"]


# ===========================================================================
# audit_reproducibility: where the data lives
# ===========================================================================


class TestAssetResolvers:
    """Absent data is a documented state; a path only this machine has is not."""

    def test_a_path_inside_the_repo_passes(self, tmp_path):
        record = {"name": "vendored", "env_var": "YS_X", "resolved": tmp_path / "data",
                  "defaults": []}

        rows = repro.check_asset_resolvers([record], tmp_path, environ={})

        assert _status(rows, "vendored") == "PASS"

    def test_an_absent_asset_is_not_a_failure(self, tmp_path):
        """The whole repository is built to work without the wet-lab exports."""
        record = {"name": "plates", "env_var": "YS_PLATES", "resolved": None, "defaults": []}

        rows = repro.check_asset_resolvers([record], tmp_path, environ={})

        assert _status(rows, "plates") == "PASS"

    def test_a_path_outside_the_repo_with_no_override_fails(self, tmp_path):
        record = {"name": "plates", "env_var": "YS_PLATES",
                  "resolved": pathlib.Path("/somewhere/else"), "defaults": []}

        rows = repro.check_asset_resolvers([record], tmp_path, environ={})

        assert _status(rows, "plates") == "FAIL"

    def test_the_same_path_passes_once_a_variable_points_at_it(self, tmp_path):
        """An override is what makes an outside path portable: it is written down."""
        record = {"name": "plates", "env_var": "YS_PLATES",
                  "resolved": pathlib.Path("/somewhere/else"), "defaults": []}

        rows = repro.check_asset_resolvers([record], tmp_path,
                                           environ={"YS_PLATES": "/somewhere/else"})

        assert _status(rows, "plates") == "PASS"

    def test_a_resolver_that_raises_fails(self, tmp_path):
        record = {"name": "broken", "env_var": None, "resolved": "error: OSError: nope",
                  "defaults": []}

        rows = repro.check_asset_resolvers([record], tmp_path, environ={})

        assert _status(rows, "broken") == "FAIL"


class TestDescribeResolvers:
    """The environment variable is read off the call, not off the source text."""

    def test_it_recovers_the_variable_each_resolver_consults(self):
        module = importlib.import_module("ystwin.paths")

        records = {r["name"]: r["env_var"] for r in repro.describe_resolvers(module)}

        assert records["yeast_gem"] == "YSTWIN_YEAST_GEM"
        assert records["biosensor_plates"] == "YSTWIN_PLATES"

    def test_a_resolver_consulting_no_variable_reports_none(self):
        """``data_dir`` is vendored and ``outputs_dir`` is a destination; neither overrides."""
        module = importlib.import_module("ystwin.paths")

        records = {r["name"]: r["env_var"] for r in repro.describe_resolvers(module)}

        assert records["data_dir"] is None

    def test_it_leaves_the_module_as_it_found_it(self):
        module = importlib.import_module("ystwin.paths")
        before = module._resolve

        repro.describe_resolvers(module)

        assert module._resolve is before


class TestEnvVarsDocumented:
    def test_a_documented_variable_passes(self):
        records = [{"name": "plates", "env_var": "YS_PLATES"}]

        rows = repro.check_env_vars_documented(records, "set `YS_PLATES` to the directory")

        assert _status(rows, "YS_PLATES") == "PASS"

    def test_an_undocumented_variable_fails(self):
        records = [{"name": "plates", "env_var": "YS_PLATES"}]

        rows = repro.check_env_vars_documented(records, "no mention of it anywhere")

        assert _status(rows, "YS_PLATES") == "FAIL"


# ===========================================================================
# audit_reproducibility: the environment the numbers were produced in
# ===========================================================================


PINNED = '[project]\nrequires-python = ">=3.11"\ndependencies = [\n  "numpy==2.4.1",\n]\n'


class TestDependencyPins:
    def test_a_pin_matching_the_installed_version_passes(self):
        rows = repro.check_dependency_pins(PINNED, lambda name: "2.4.1")

        assert _status(rows, "pin matches installed: numpy") == "PASS"

    def test_a_pin_that_is_not_what_is_installed_fails(self):
        """The silent-divergence case: the file still reads correct."""
        rows = repro.check_dependency_pins(PINNED, lambda name: "2.3.0")

        assert _status(rows, "pin matches installed: numpy") == "FAIL"

    def test_a_range_fails(self):
        ranged = PINNED.replace("numpy==2.4.1", "numpy>=2.0")

        rows = repro.check_dependency_pins(ranged, lambda name: "2.4.1")

        assert _status(rows, "dependency pinned: numpy") == "FAIL"

    def test_a_missing_runtime_dependency_fails(self):
        rows = repro.check_dependency_pins(PINNED, lambda name: None)

        assert _status(rows, "pin matches installed: numpy") == "FAIL"

    def test_a_missing_optional_dependency_is_skipped_not_failed(self):
        """``[thermo]`` is documented as optional and its tests skip without it."""
        text = PINNED + '\n[project.optional-dependencies]\nthermo = ["pytfa==0.9.1"]\n'

        rows = repro.check_dependency_pins(text, lambda name: None if name == "pytfa" else "2.4.1")

        assert _status(rows, "pin matches installed: pytfa") == "SKIP"

    def test_a_tool_section_is_not_read_as_a_list_of_dependencies(self):
        """The defect that made this a TOML parse instead of a regex.

        Adding `[tool.ruff]` to this repository's own pyproject.toml produced four FAIL
        rows -- ``dependency pinned: E``, ``F``, ``E501`` and ``archive`` -- because the
        group discovery matched every ``name = [...]`` in the document and skipped seven
        key names by hand. A linter rule is not a dependency and no blocklist can be
        relied on to know that; PEP 621 says which keys hold requirements, so this asks it.
        """
        text = (PINNED
                + '\n[tool.ruff]\nextend-exclude = ["archive"]\n'
                + '\n[tool.ruff.lint]\nselect = ["E", "F"]\nignore = ["E501"]\n'
                + '\n[tool.pytest.ini_options]\nmarkers = ["slow: long-running"]\n')

        rows = repro.check_dependency_pins(text, lambda name: "2.4.1")

        assert [r["check"] for r in rows] == ["pin matches installed: numpy [runtime]"]

    def test_a_pyproject_that_is_not_valid_toml_is_refused_not_ignored(self):
        """A file this cannot read must fail loudly. Returning no rows would report a
        repository with unreadable pins as one with no unpinned dependency."""
        rows = repro.check_dependency_pins("[project\ndependencies = [", lambda name: "1.0")

        assert _status(rows, "dependency pins readable") == "FAIL"


class TestPythonVersion:
    def test_a_new_enough_interpreter_passes(self):
        rows = repro.check_python_version(PINNED, version_info=(3, 13, 0))

        assert _status(rows, "requires-python") == "PASS"

    def test_too_old_an_interpreter_fails(self):
        rows = repro.check_python_version(PINNED, version_info=(3, 9, 0))

        assert _status(rows, "requires-python") == "FAIL"


class TestSolver:
    def test_glpk_passes(self):
        rows = repro.check_solver("optlang.glpk_interface", ["glpk"])

        assert _status(rows, "cobra solver") == "PASS"

    def test_a_commercial_solver_fails(self):
        """Not wrong, but not the vertex the committed FVA widths came from."""
        rows = repro.check_solver("optlang.gurobi_interface", ["glpk", "gurobi"])

        assert _status(rows, "cobra solver") == "FAIL"

    def test_no_cobra_is_a_skip(self):
        rows = repro.check_solver(None)

        assert _status(rows, "cobra solver") == "SKIP"


# ===========================================================================
# audit_reproducibility: what is in outputs/
# ===========================================================================


@pytest.fixture
def repo(tmp_path):
    """A miniature checkout: git index, outputs/, scripts/, src/ystwin/."""
    (tmp_path / "outputs").mkdir()
    (tmp_path / "scripts").mkdir()
    (tmp_path / "src" / "ystwin").mkdir(parents=True)
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True, capture_output=True)
    return tmp_path


def _stage(repo_path: pathlib.Path, relative: str) -> None:
    """Add to the index. ``git ls-files`` reads the index, so no commit is needed."""
    subprocess.run(["git", "add", relative], cwd=repo_path, check=True, capture_output=True)


class TestOutputsTracked:
    def test_a_staged_table_passes(self, repo):
        (repo / "outputs" / "result.csv").write_text("a,b\n1,2\n")
        _stage(repo, "outputs/result.csv")

        rows = repro.check_outputs_tracked(repo)

        assert _status(rows, "outputs/result.csv") == "PASS"

    def test_an_untracked_table_fails(self, repo):
        """A result nobody else can see is not a result."""
        (repo / "outputs" / "result.csv").write_text("a,b\n1,2\n")

        rows = repro.check_outputs_tracked(repo)

        assert _status(rows, "outputs/result.csv") == "FAIL"

    def test_a_deliberately_ignored_file_passes(self, repo):
        """The .gitignore states the policy; being on either side of it is fine."""
        (repo / ".gitignore").write_text("*.png\n")
        (repo / "outputs" / "figure.png").write_bytes(b"\x89PNG")

        rows = repro.check_outputs_tracked(repo)

        assert _status(rows, "outputs/figure.png") == "PASS"


class TestOutputTables:
    def test_a_table_with_rows_passes(self, repo):
        (repo / "outputs" / "result.csv").write_text("a,b\n1,2\n")

        rows = repro.check_output_tables(repo)

        assert _status(rows, "outputs/result.csv") == "PASS"

    def test_an_empty_file_fails(self, repo):
        (repo / "outputs" / "result.csv").write_text("")

        rows = repro.check_output_tables(repo)

        assert _status(rows, "outputs/result.csv") == "FAIL"

    def test_headers_with_no_data_fail(self, repo):
        """What a crashed script leaves behind, and it reads exactly like a result."""
        (repo / "outputs" / "result.csv").write_text("a,b\n")

        rows = repro.check_output_tables(repo)

        assert _status(rows, "outputs/result.csv") == "FAIL"

    def test_something_that_is_not_a_csv_fails(self, repo):
        (repo / "outputs" / "result.csv").write_text('a,b\n"unclosed,2\n1\n')

        rows = repro.check_output_tables(repo)

        assert _status(rows, "outputs/result.csv") == "FAIL"


# ===========================================================================
# audit_reproducibility: which script writes which file
# ===========================================================================


class TestWritersFromScripts:
    def test_it_finds_a_literal_destination(self, repo):
        (repo / "scripts" / "run_x.py").write_text(
            'OUT = paths.outputs_dir()\nframe.to_csv(OUT / "result.csv", index=False)\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["outputs/result.csv"]}

    def test_a_computed_name_becomes_a_pattern(self, repo):
        (repo / "scripts" / "run_x.py").write_text('g1.to_csv(OUT / f"g1_{tag}.csv")\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["g1_*.csv"]}

    def test_a_read_is_not_a_write(self, repo):
        """``run_calibration.py`` names the characterisation table; it reads it."""
        (repo / "scripts" / "run_x.py").write_text(
            'MEASURED = OUT / "sensor.csv"\nframe.to_csv(OUT / "result.csv")\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["result.csv"]}

    def test_it_finds_a_model_saved_rather_than_written(self, repo):
        (repo / "scripts" / "run_x.py").write_text('model.save(OUT / "model.npz")\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["model.npz"]}

    def test_a_figure_is_a_write(self, repo):
        """``savefig`` was not in the set, and its absence reported all twelve committed
        figures as having no writer at all -- "nobody can regenerate it" -- while
        ``make_figures.py`` regenerates them on demand. A false alarm from an integrity
        gate is worse than a quiet one: the honest response is to stop believing it."""
        (repo / "scripts" / "run_x.py").write_text('fig.savefig(OUT / "fig01.svg")\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["fig01.svg"]}

    def test_a_destination_bound_to_a_name_first_is_still_found(self, repo):
        """The pattern can be lines above the call, and a walk that only reads the call
        sees a bare ``Name`` and gives up. ``viz/figures.py`` writes this way, though its
        stem is the caller's -- see the bare-extension tests below for that half."""
        (repo / "scripts" / "run_x.py").write_text(
            'destination = OUT / "result.csv"\nframe.to_csv(destination, index=False)\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["result.csv"]}

    def test_a_write_inside_a_library_module_belongs_to_the_script_that_reaches_it(
            self, repo):
        """``make_figures.py`` writes nothing itself: every figure leaves through
        ``fig.savefig`` inside ``ystwin.viz.figures``. Attribution follows the same
        reachability walk the freshness check already runs, because a script that calls a
        module owns what that module writes."""
        (repo / "src" / "ystwin" / "drawing.py").write_text(
            'def draw(out):\n    png = out / "fig09.png"\n    fig.savefig(png)\n')
        (repo / "scripts" / "run_x.py").write_text(
            'from ystwin.drawing import draw\n\ndraw(OUT)\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["fig09.png"]}

    def test_a_bare_extension_claims_nothing(self, repo):
        """``viz/figures.py`` writes ``out_dir / f"{stem}.svg"`` and the stem is the
        caller's, so resolving the binding yields ``*.svg`` -- true, and useless. It
        matches every SVG in the repository, so every script reaching that module would
        claim every figure and "nobody can regenerate this" would never fire for a figure
        again. Attribution has to come from somewhere that knows which figure."""
        (repo / "src" / "ystwin" / "drawing.py").write_text(
            'def draw(out, stem):\n    svg = out / f"{stem}.svg"\n    fig.savefig(svg)\n')
        (repo / "scripts" / "run_x.py").write_text(
            'from ystwin.drawing import draw\n\ndraw(OUT, "anything")\n')

        assert repro.writers_from_scripts(repo) == {}

    def test_a_wildcard_with_a_literal_prefix_still_claims(self, repo):
        """The other side of the line. ``g1_*.csv`` names a family of files this script and
        no other writes; only the bare extension is too weak."""
        (repo / "scripts" / "run_x.py").write_text('g1.to_csv(OUT / f"g1_{tag}.csv")\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["g1_*.csv"]}

    def test_a_loop_over_literal_names_resolves_to_those_names(self, repo):
        """``run_training.py`` writes four tables through one call:
        ``for frame, name in ((a, "training_recovery"), ...): frame.to_csv(OUT / f"{name}.csv")``.
        Read without the loop that is ``*.csv``, which claims every table in the repository
        -- and once bare extensions stopped attributing, it claimed none. The literal names
        are two lines above the call."""
        (repo / "scripts" / "run_x.py").write_text(
            'for frame, name in ((a, "one"), (b, "two")):\n'
            '    frame.to_csv(OUT / f"{name}.csv", index=False)\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["one.csv", "two.csv"]}

    def test_a_loop_over_a_bare_list_of_names_resolves_too(self, repo):
        (repo / "scripts" / "run_x.py").write_text(
            'for name in ["alpha", "beta"]:\n'
            '    frame.to_csv(OUT / f"{name}.csv")\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["alpha.csv", "beta.csv"]}

    def test_a_computed_loop_variable_keeps_the_wildcard(self, repo):
        """And is then dropped as too weak, which is the honest outcome: the script writes
        some CSV whose name nobody can know without running it, so it cannot be said to be
        the writer of any particular one."""
        (repo / "scripts" / "run_x.py").write_text(
            'for name in discover():\n'
            '    frame.to_csv(OUT / f"{name}.csv")\n')

        assert repro.writers_from_scripts(repo) == {}

    def test_a_prefix_before_the_loop_variable_survives_the_substitution(self, repo):
        (repo / "scripts" / "run_x.py").write_text(
            'for name in ("one",):\n'
            '    frame.to_csv(OUT / f"table_{name}.csv")\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["table_one.csv"]}

    def test_the_figure_manifest_supplies_the_attribution_the_source_cannot(self, repo):
        """``make_figures.py`` records ``figure_produced_by`` beside every file it emits.
        That is the authoritative statement of what produced what, written by the thing
        that produced it -- a structured record rather than a glob over source text."""
        (repo / "scripts" / "make_figures.py").write_text("VALUE = 1\n")
        (repo / "outputs" / "figures_manifest.csv").write_text(
            "figure,svg,png,figure_produced_by\n"
            "fig01,outputs/fig01.svg,outputs/fig01.png,python3 scripts/make_figures.py\n")

        assert repro.writers_from_scripts(repo)["make_figures.py"] == [
            "outputs/fig01.png", "outputs/fig01.svg"]

    def test_a_figure_no_manifest_row_claims_still_has_no_writer(self, repo):
        """The case this whole distinction exists for: three committed ``atp_sensor_*.png``
        are produced by nothing in ``scripts/``, and a ``*.png`` glob was hiding that."""
        (repo / "scripts" / "make_figures.py").write_text("VALUE = 1\n")
        (repo / "outputs" / "figures_manifest.csv").write_text(
            "figure,svg,png,figure_produced_by\n"
            "fig01,outputs/fig01.svg,outputs/fig01.png,python3 scripts/make_figures.py\n")

        writers = repro.writers_from_scripts(repo)

        assert repro._writers_of("atp_sensor_confidence.png", writers) == []

    def test_an_absent_manifest_contributes_nothing_rather_than_raising(self, repo):
        """The audit runs on repositories that have never generated a figure."""
        (repo / "scripts" / "run_x.py").write_text('frame.to_csv(OUT / "result.csv")\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["result.csv"]}

    def test_a_manifest_missing_its_columns_contributes_nothing(self, repo):
        (repo / "scripts" / "run_x.py").write_text('frame.to_csv(OUT / "result.csv")\n')
        (repo / "outputs" / "figures_manifest.csv").write_text("figure\nfig01\n")

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["result.csv"]}

    def test_a_module_the_script_does_not_import_is_not_attributed_to_it(self, repo):
        """The other half. Attribution that followed every module in ``src`` would blame
        every script for every output and the freshness check would become noise."""
        (repo / "src" / "ystwin" / "elsewhere.py").write_text(
            'def draw(out):\n    png = out / "fig09.png"\n    fig.savefig(png)\n')
        (repo / "scripts" / "run_x.py").write_text('frame.to_csv(OUT / "result.csv")\n')

        assert repro.writers_from_scripts(repo) == {"run_x.py": ["result.csv"]}


class TestWritersFromDoc:
    def test_it_reads_the_table_in_the_guide(self):
        text = (
            "| Script | Writes into `outputs/` | Needs |\n"
            "| --- | --- | --- |\n"
            "| `run_gates.py` | `g1_<plate>.csv`, `d2_g1passed_<plate>.csv` | data |\n"
        )

        assert repro.writers_from_doc(text) == {
            "run_gates.py": ["d2_g1passed_*.csv", "g1_*.csv"]}


class TestWriterMapMatchesDoc:
    def test_agreement_passes(self):
        rows = repro.check_writer_map_matches_doc({"run_x.py": ["result.csv"]},
                                                  {"run_x.py": ["result.csv"]})

        assert _status(rows, "run_x.py") == "PASS"

    def test_a_documented_table_nothing_writes_fails(self):
        rows = repro.check_writer_map_matches_doc({"run_x.py": ["result.csv"]},
                                                  {"run_x.py": ["result.csv", "ghost.csv"]})

        assert _status(rows, "run_x.py") == "FAIL"

    def test_a_new_undocumented_table_fails(self):
        rows = repro.check_writer_map_matches_doc({"run_x.py": ["result.csv", "extra.csv"]},
                                                  {"run_x.py": ["result.csv"]})

        assert _status(rows, "run_x.py") == "FAIL"

    def test_the_guide_may_be_more_specific_than_the_source(self):
        """``f"d2_{tag}.csv"`` and ``d2_<plate>__<reporter>.csv`` are the same files."""
        rows = repro.check_writer_map_matches_doc({"run_d2.py": ["d2_*.csv"]},
                                                  {"run_d2.py": ["d2_*__*.csv"]})

        assert _status(rows, "run_d2.py") == "PASS"


class TestLibraryDependencies:
    def test_it_follows_imports_transitively_through_relative_ones(self, repo):
        src = repo / "src" / "ystwin"
        (src / "__init__.py").write_text("")
        (src / "top.py").write_text("from .middle import thing\n")
        (src / "middle.py").write_text("from ..leaf import other\n")
        (src / "leaf.py").write_text("")
        script = repo / "scripts" / "run_x.py"
        script.write_text("from ystwin.top import thing\n")

        found = {p.name for p in repro.library_dependencies(script, src)}

        assert found == {"top.py", "middle.py", "leaf.py"}

    def test_it_ignores_modules_the_script_never_reaches(self, repo):
        src = repo / "src" / "ystwin"
        (src / "__init__.py").write_text("")
        (src / "used.py").write_text("")
        (src / "unused.py").write_text("")
        script = repo / "scripts" / "run_x.py"
        script.write_text("from ystwin import used\n")

        found = {p.name for p in repro.library_dependencies(script, src)}

        assert found == {"used.py", "__init__.py"}


def _age(path: pathlib.Path, seconds: float) -> None:
    """Backdate a file, so mtime ordering can be stated rather than raced."""
    when = 1_700_000_000 - seconds
    os.utime(path, (when, when))


class TestOutputFreshness:
    def test_a_table_newer_than_its_code_passes(self, repo):
        script = repo / "scripts" / "run_x.py"
        script.write_text('frame.to_csv(OUT / "result.csv")\n')
        table = repo / "outputs" / "result.csv"
        table.write_text("a\n1\n")
        _age(script, 3600)
        _age(table, 0)

        rows = repro.check_output_freshness(repo, repro.writers_from_scripts(repo))

        assert _status(rows, "outputs/result.csv") == "PASS"

    def test_a_table_older_than_its_script_fails(self, repo):
        script = repo / "scripts" / "run_x.py"
        script.write_text('frame.to_csv(OUT / "result.csv")\n')
        table = repo / "outputs" / "result.csv"
        table.write_text("a\n1\n")
        _age(table, 3600)
        _age(script, 0)

        rows = repro.check_output_freshness(repo, repro.writers_from_scripts(repo))

        assert _status(rows, "outputs/result.csv") == "FAIL"

    def test_a_table_older_than_the_library_it_used_fails(self, repo):
        """The scripts are thin; the behaviour is under src/."""
        src = repo / "src" / "ystwin"
        (src / "__init__.py").write_text("")
        module = src / "engine.py"
        module.write_text("")
        script = repo / "scripts" / "run_x.py"
        script.write_text('from ystwin import engine\nframe.to_csv(OUT / "result.csv")\n')
        table = repo / "outputs" / "result.csv"
        table.write_text("a\n1\n")
        for path in (script, table, src / "__init__.py"):
            _age(path, 3600)
        _age(module, 0)

        rows = repro.check_output_freshness(repo, repro.writers_from_scripts(repo))

        assert _status(rows, "outputs/result.csv") == "FAIL"
        assert "engine.py" in [r for r in rows if "result.csv" in r["check"]][0]["detail"]

    def test_a_library_edit_the_script_never_imports_is_not_blamed(self, repo):
        src = repo / "src" / "ystwin"
        (src / "__init__.py").write_text("")
        (src / "unrelated.py").write_text("")
        script = repo / "scripts" / "run_x.py"
        script.write_text('frame.to_csv(OUT / "result.csv")\n')
        table = repo / "outputs" / "result.csv"
        table.write_text("a\n1\n")
        _age(script, 3600)
        _age(table, 3600)
        _age(src / "unrelated.py", 0)

        rows = repro.check_output_freshness(repo, repro.writers_from_scripts(repo))

        assert _status(rows, "outputs/result.csv") == "PASS"

    def test_a_table_no_script_writes_fails(self, repo):
        """Present, citable, and with no reproducer at all."""
        (repo / "outputs" / "orphan.csv").write_text("a\n1\n")

        rows = repro.check_output_freshness(repo, {})

        assert _status(rows, "outputs/orphan.csv") == "FAIL"

    def test_the_more_specific_pattern_claims_the_file(self, repo):
        """``d2_*.csv`` also matches ``d2_g1passed_x.csv``; the longer match wins."""
        (repo / "outputs" / "d2_g1passed_x.csv").write_text("a\n1\n")
        writers = {"run_d2.py": ["d2_*.csv"], "run_gates.py": ["d2_g1passed_*.csv"]}
        (repo / "scripts" / "run_d2.py").write_text("")
        (repo / "scripts" / "run_gates.py").write_text("")

        rows = repro.check_output_freshness(repo, writers)

        matched = [r for r in rows if "d2_g1passed_x.csv" in r["check"]]
        assert len(matched) == 1
        assert "run_gates.py" in matched[0]["expected"]


@pytest.mark.slow
class TestOfflineSuite:
    """The check that proves the repository works for someone with no plate reader."""

    def _suite(self, tmp_path, body: str) -> pathlib.Path:
        (tmp_path / "tests").mkdir()
        (tmp_path / "tests" / "test_mini.py").write_text(body)
        return tmp_path

    def test_a_suite_that_skips_without_data_passes(self, tmp_path):
        root = self._suite(tmp_path, (
            "import os, pytest\n"
            "def test_needs_data():\n"
            "    if not os.path.isdir(os.environ.get('MINI_DATA', '/nowhere')):\n"
            "        pytest.skip('no data')\n"
            "    assert True\n"))

        rows = repro.check_offline_suite(root, ["MINI_DATA"])

        assert _status(rows, "suite passes with no real data") == "PASS"

    def test_a_suite_that_errors_without_data_fails(self, tmp_path):
        """A fixture that reads a file directly instead of going through paths.py."""
        root = self._suite(tmp_path, (
            "import os\n"
            "def test_needs_data():\n"
            "    assert os.path.isdir(os.environ.get('MINI_DATA', '/nowhere'))\n"))

        rows = repro.check_offline_suite(root, ["MINI_DATA"])

        assert _status(rows, "suite passes with no real data") == "FAIL"


# ===========================================================================
# audit_determinism: comparing two results exactly
# ===========================================================================


class TestDigest:
    def test_equal_values_digest_equally(self):
        assert determinism.digest({"a": np.arange(3)}) == determinism.digest({"a": np.arange(3)})

    def test_a_difference_in_the_last_bit_is_a_difference(self):
        """Two bootstraps agreeing to fifteen decimals are not the same bootstrap."""
        one = 0.1 + 0.2
        other = 0.30000000000000004 + 1e-17

        assert determinism.digest(one) != determinism.digest(other + 1e-16)

    def test_nan_matches_nan(self):
        assert determinism.digest(np.array([np.nan])) == determinism.digest(np.array([np.nan]))

    def test_dictionary_order_does_not_matter(self):
        assert determinism.digest({"a": 1, "b": 2}) == determinism.digest({"b": 2, "a": 1})

    def test_an_array_is_not_its_list(self):
        assert determinism.digest(np.array([1.0])) != determinism.digest([1.0])

    def test_it_sees_inside_an_object(self):
        class Result:
            def __init__(self, value):
                self.value = value

        assert determinism.digest(Result(1.0)) != determinism.digest(Result(2.0))


# ===========================================================================
# audit_determinism: the two properties
# ===========================================================================


def _specimen(call, **kwargs):
    """A specimen whose target resolves to something harmless; ``call`` ignores it."""
    return determinism.Specimen("builtins::len", lambda target, seed: call(seed), **kwargs)


class TestCheckSpecimen:
    def test_a_properly_seeded_function_passes_both(self):
        rows = determinism.check_specimen(
            "good", _specimen(lambda seed: np.random.default_rng(seed).normal(size=4)))

        assert [r["status"] for r in rows] == ["PASS", "PASS"]

    def test_an_unseeded_draw_fails_repeatability(self):
        """The global generator: same seed, different answer."""
        rows = determinism.check_specimen(
            "leaky", _specimen(lambda seed: np.random.default_rng().normal(size=4)))

        assert _status(rows, "same seed, same answer") == "FAIL"

    def test_a_function_that_ignores_its_seed_fails_sensitivity(self):
        """Reproducible, and the seed is decoration. Invisible without this check."""
        rows = determinism.check_specimen("inert", _specimen(lambda seed: 42.0))

        assert _status(rows, "same seed, same answer") == "PASS"
        assert _status(rows, "seed changes the answer") == "FAIL"

    def test_a_declared_insensitive_function_passes_when_it_is(self):
        rows = determinism.check_specimen(
            "choice", _specimen(lambda seed: 2, seed_sensitive=False,
                                reason="an argmin over integer widths"))

        assert [r["status"] for r in rows] == ["PASS", "PASS"]

    def test_a_declared_insensitive_function_fails_when_it_is_not(self):
        """The declaration is asserted in both directions, so it cannot hide a finding."""
        rows = determinism.check_specimen(
            "surprise", _specimen(lambda seed: float(seed), seed_sensitive=False,
                                  reason="claimed to be an argmin"))

        assert _status(rows, "seed changes the answer") == "FAIL"

    def test_declaring_insensitivity_without_a_reason_fails(self):
        rows = determinism.check_specimen(
            "unjustified", _specimen(lambda seed: 2, seed_sensitive=False))

        assert [r["status"] for r in rows] == ["FAIL"]

    def test_a_specimen_that_raises_fails(self):
        def boom(seed):
            raise ValueError("need at least 5 timepoints")

        rows = determinism.check_specimen("broken", _specimen(boom))

        assert rows[0]["status"] == "FAIL"
        assert "ValueError" in rows[0]["actual"]

    def test_a_specimen_resolves_the_object_it_names(self):
        specimen = determinism.Specimen(
            "builtins::sorted", lambda target, seed: target([seed, 0]))

        rows = determinism.check_specimen("sorted", specimen)

        assert _status(rows, "same seed, same answer") == "PASS"


# ===========================================================================
# audit_determinism: discovery and coverage
# ===========================================================================


@pytest.fixture
def fake_package(tmp_path):
    """A two-module package with one seeded function, one without, and a re-export."""
    package = tmp_path / "faketwin"
    package.mkdir()
    (package / "__init__.py").write_text("from .core import seeded\n")
    (package / "core.py").write_text(
        "def seeded(data, seed: int = 0):\n    return data\n\n\n"
        "def unseeded(data):\n    return data\n\n\n"
        "def _private(data, seed: int = 0):\n    return data\n")
    sys.path.insert(0, str(tmp_path))
    try:
        yield importlib.import_module("faketwin")
    finally:
        sys.path.remove(str(tmp_path))
        for name in [n for n in sys.modules if n.startswith("faketwin")]:
            del sys.modules[name]


class TestDiscovery:
    def test_it_finds_a_seeded_function(self, fake_package):
        found = determinism.discover_seeded(fake_package)

        assert "faketwin.core::seeded" in found

    def test_it_ignores_a_function_with_no_seed(self, fake_package):
        found = determinism.discover_seeded(fake_package)

        assert "faketwin.core::unseeded" not in found

    def test_it_ignores_private_names(self, fake_package):
        found = determinism.discover_seeded(fake_package)

        assert not any(name.endswith("_private") for name in found)

    def test_a_re_export_is_counted_once_at_its_definition(self, fake_package):
        """``__init__`` re-exports ``seeded``; it is one entry point, not two."""
        found = determinism.discover_seeded(fake_package)

        assert [n for n in found if n.endswith("::seeded")] == ["faketwin.core::seeded"]


class TestCoverage:
    def test_a_covered_entry_point_passes(self):
        discovered = {"pkg.mod::f": None}
        specimens = {"f": determinism.Specimen("pkg.mod::f", lambda t, s: 0)}

        rows = determinism.check_coverage(discovered, specimens, {})

        assert _status(rows, "every seeded entry point is exercised") == "PASS"

    def test_a_new_seeded_function_with_no_specimen_fails(self):
        """The nag: a seeded function added today is outside the audit until covered."""
        rows = determinism.check_coverage({"pkg.mod::brand_new": None}, {}, {})

        assert _status(rows, "every seeded entry point is exercised") == "FAIL"
        assert "brand_new" in [r for r in rows if "exercised" in r["check"]][0]["detail"]

    def test_an_exempt_entry_point_passes(self):
        rows = determinism.check_coverage({"pkg.mod::needs_gsmm": None}, {},
                                          {"pkg.mod::needs_gsmm": "needs the GSMM"})

        assert _status(rows, "every seeded entry point is exercised") == "PASS"

    def test_a_specimen_for_a_deleted_function_fails(self):
        specimens = {"gone": determinism.Specimen("pkg.mod::gone", lambda t, s: 0)}

        rows = determinism.check_coverage({}, specimens, {})

        assert _status(rows, "no longer exists") == "FAIL"

    def test_an_unimportable_module_is_skipped_not_failed(self):
        """An optional dependency this machine lacks is not a determinism finding."""
        marker = f"pkg.thermo::{determinism.IMPORT_FAILED}"

        rows = determinism.check_coverage({marker: None}, {}, {})

        assert _status(rows, "module importable") == "SKIP"


class TestProvenanceCatchesWhatAnMtimeCannot:
    """The hole that let eight artefacts go stale behind a green gate.

    `check_output_freshness` compares modification times, and its own docstring defends the
    choice: *"a fresh clone writes every file at once, so mtimes there are equal and nothing
    fires"*. That was written as a reassurance. It is the vulnerability -- `git checkout`
    stamps every file with the checkout time, so on any clone the gate passes **vacuously**,
    on content that need not match the code at all. On 2026-08-31 it reported 284/284 while
    `outputs/heldout_scores.csv` held 179 training rows against the code's 287 and five
    biosensor tables moved enough to break ten pinned prose numbers.

    A digest of the producers fires in exactly the case an mtime cannot: every file the same
    age, and the content disagreeing anyway.
    """

    def _tree(self, tmp_path):
        (tmp_path / "scripts").mkdir()
        (tmp_path / "outputs").mkdir()
        (tmp_path / "src" / "ystwin").mkdir(parents=True)
        (tmp_path / "scripts" / "run_thing.py").write_text(
            "from pathlib import Path\nPath('outputs/thing.csv').write_text('a,b\\n1,2\\n')\n")
        (tmp_path / "outputs" / "thing.csv").write_text("a,b\n1,2\n")
        (tmp_path / "data").mkdir()
        run = {"id": "thing", "lifecycle": "current", "status": "investigation_pending",
               "reason": "Synthetic reproducibility fixture", "producer": {"path": "scripts/run_thing.py", "callable": "__main__"},
               "dependencies": [], "parameters": {"rows": 1}, "runtime": {"kind": "python"},
               "model_identity": {"kind": "none", "reason": "literal fixture"},
               "reproduction": {"argv": ["{python}", "-B", "scripts/run_thing.py"]},
               "artifacts": [{"path": "outputs/thing.csv", "comparison": {"format": "csv", "row_keys": ["a"], "units": {"b": "dimensionless"}}}]}
        (tmp_path / "data/artifact_registry.json").write_text(json.dumps({"schema_version": 1, "runs": [run]}))
        return tmp_path, {"run_thing.py": ["outputs/thing.csv"]}

    def _verify(self, repo):
        destination = repo.with_name(repo.name + "-reproduction")
        receipt = artifacts.reproduce_run(repo, "thing", destination)
        assert receipt["success"]
        proof = repo / "data/receipt.json"
        proof.write_text(json.dumps(receipt))
        document = json.loads((repo / "data/artifact_registry.json").read_text())
        document["runs"][0]["status"] = "verified"
        document["runs"][0]["verification"] = {"path": "data/receipt.json", "sha256": artifacts.sha256(proof)}
        (repo / "data/artifact_registry.json").write_text(json.dumps(document))

    def test_a_stamped_tree_passes(self, tmp_path):
        repo, writers = self._tree(tmp_path)
        self._verify(repo)

        rows = repro.check_output_provenance(repo, writers)

        assert [r["status"] for r in rows] == ["PASS"]

    def test_editing_the_script_fails_even_with_every_mtime_equal(self, tmp_path):
        """The decisive case. Both files are rewritten to the same instant, which is what a
        clone looks like, and the digest still fires."""
        import os
        import time

        repo, writers = self._tree(tmp_path)
        self._verify(repo)
        (repo / "scripts" / "run_thing.py").write_text("# edited\nimport pandas as pd\n")
        now = time.time()
        for path in list(repo.rglob("*")):
            if path.is_file():
                os.utime(path, (now, now))

        rows = repro.check_output_provenance(repo, writers)

        assert [r["status"] for r in rows] == ["FAIL"]
        assert "identity changed; numerical effect is untested" in rows[0]["detail"]

    def test_an_unstamped_tree_says_so_rather_than_passing_silently(self, tmp_path):
        repo, writers = self._tree(tmp_path)

        rows = repro.check_output_provenance(repo, writers)

        assert len(rows) == 1
        assert rows[0]["status"] == "FAIL"
        assert "successful reproduction and semantic comparison have not been recorded" in rows[0]["detail"]

    def test_a_new_table_is_reported_as_unstamped_not_as_stale(self, tmp_path):
        """A table nobody has stamped has not drifted -- it has never been blessed. Failing
        it would train an author to re-stamp reflexively, which is how the check dies."""
        repo, writers = self._tree(tmp_path)
        self._verify(repo)
        (repo / "outputs" / "second.csv").write_text("c\n3\n")
        writers["run_thing.py"].append("second.csv")

        rows = repro.check_output_provenance(repo, writers)
        new = [r for r in rows if "second.csv" in r["check"]]

        assert new and new[0]["status"] == "FAIL"
        assert new[0]["actual"] == "investigation_pending"
        assert "Not in the authoritative artifact registry" in new[0]["detail"]

    def test_the_digest_depends_on_content_not_on_names_alone(self, tmp_path):
        repo, _ = self._tree(tmp_path)
        script = repo / "scripts" / "run_thing.py"
        before = repro.producer_digest([script], repo)
        script.write_text(script.read_text() + "# one more line\n")

        assert repro.producer_digest([script], repo) != before
