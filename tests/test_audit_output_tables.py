"""The gate that says a committed table still comes out of its generator must itself run.

`scripts/audit_output_tables.py` closes the unguarded end of this repository's one
verification chain: `audit_claims.py` binds prose to a CELL, and until 2026-09-04 nothing
asked whether the cell still came from the code. It did not, in at least one case --
``outputs/proteomics_gain.csv`` had drifted so far from `simulate_proteomics_gain.py` that one
column changed SIGN.

A gate that nobody runs is the same as no gate, so the suite runs it over the fastest
registered pairs. Not all of them: the carotenoid fit takes fourteen seconds, which belongs on
a command line and not in a unit test. What the suite guards is that the machinery works and
that the registry has not quietly emptied.
"""

from __future__ import annotations

import importlib.util
import pathlib

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def gate():
    path = REPO / "scripts" / "audit_output_tables.py"
    spec = importlib.util.spec_from_file_location("_audit_output_tables", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestTheRegistryIsRealAndNotEmpty:
    def test_it_holds_the_table_it_was_written_for(self, gate):
        """`proteomics_gain.csv` is the proof case: it is the artefact whose drift made this
        gate necessary, so removing it from the registry has to fail."""
        assert "proteomics_gain.csv" in gate.REGISTRY["simulate_proteomics_gain"]

    def test_every_registered_script_exists(self, gate):
        missing = [name for name in gate.REGISTRY
                   if not (REPO / "scripts" / f"{name}.py").exists()]

        assert not missing, missing

    def test_every_registered_table_is_committed(self, gate):
        missing = [table for tables in gate.REGISTRY.values() for table in tables
                   if not (REPO / "outputs" / table).exists()]

        assert not missing, missing

    def test_every_exemption_states_a_reason(self, gate):
        """An unnamed omission is indistinguishable from an oversight, which is the argument
        this repository makes for `audit_determinism.py`'s exemption list too."""
        assert gate.EXEMPT
        for script, reason in gate.EXEMPT.items():
            assert len(reason) > 20, f"{script} is exempt without saying why"

    def test_nothing_is_both_registered_and_exempt(self, gate):
        assert not set(gate.REGISTRY) & set(gate.EXEMPT)

    def test_the_fast_subset_is_a_subset(self, gate):
        assert set(gate.FAST) <= set(gate.REGISTRY)


class TestTheGateActuallyRuns:
    """A smoke check over the cheap pairs, so the machinery cannot silently stop working."""

    @pytest.mark.parametrize("script", ("environment_channel", "score_phb_environment",
                                        "error_budget"))
    def test_a_registered_table_reproduces_from_its_generator(self, gate, script):
        rows = gate.check_pair(script, gate.REGISTRY[script])

        assert rows
        for row in rows:
            assert row["ok"], f"{row['table']}: {row['detail']}"


class TestItCanActuallyFail:
    """A gate that cannot report a mismatch is decoration."""

    def test_a_changed_cell_is_reported_with_its_row_and_column(self, gate, tmp_path):
        import pandas as pd

        committed = tmp_path / "committed.csv"
        regenerated = tmp_path / "regenerated.csv"
        pd.DataFrame({"a": [1.0, 2.0], "b": ["x", "y"]}).to_csv(committed, index=False)
        pd.DataFrame({"a": [1.0, 2.5], "b": ["x", "y"]}).to_csv(regenerated, index=False)

        reason = gate.compare_tables(committed, regenerated)

        assert "'a'" in reason and "row 1" in reason

    def test_a_changed_text_cell_is_reported_too(self, gate, tmp_path):
        import pandas as pd

        committed = tmp_path / "committed.csv"
        regenerated = tmp_path / "regenerated.csv"
        pd.DataFrame({"law": ["mu only", "z only"]}).to_csv(committed, index=False)
        pd.DataFrame({"law": ["mu only", "z alone"]}).to_csv(regenerated, index=False)

        assert "row 1" in gate.compare_tables(committed, regenerated)

    def test_a_new_column_is_reported_rather_than_ignored(self, gate, tmp_path):
        import pandas as pd

        committed = tmp_path / "committed.csv"
        regenerated = tmp_path / "regenerated.csv"
        pd.DataFrame({"a": [1.0]}).to_csv(committed, index=False)
        pd.DataFrame({"a": [1.0], "b": [2.0]}).to_csv(regenerated, index=False)

        assert "columns changed" in gate.compare_tables(committed, regenerated)

    def test_a_boolean_column_is_compared_as_a_verdict_and_not_as_a_number(self, gate,
                                                                          tmp_path):
        """`True - True` raises in numpy, and a boolean is a verdict with no tolerance."""
        import pandas as pd

        committed = tmp_path / "committed.csv"
        regenerated = tmp_path / "regenerated.csv"
        pd.DataFrame({"is_a_prediction": [True, False]}).to_csv(committed, index=False)
        pd.DataFrame({"is_a_prediction": [True, True]}).to_csv(regenerated, index=False)

        assert "row 1" in gate.compare_tables(committed, regenerated)

    def test_an_identical_pair_reports_nothing(self, gate, tmp_path):
        import pandas as pd

        committed = tmp_path / "committed.csv"
        regenerated = tmp_path / "regenerated.csv"
        for path in (committed, regenerated):
            pd.DataFrame({"a": [1.0, 2.0], "ok": [True, False]}).to_csv(path, index=False)

        assert gate.compare_tables(committed, regenerated) == ""


def test_the_output_override_every_generator_must_honour():
    """`YSTWIN_OUTPUTS` is how this gate redirects a table, so no script may crash on it.

    Seven scripts printed their destination with ``out.relative_to(paths.REPO_ROOT)``, which
    RAISES for a path outside the repository -- so they could not run under the very override
    `paths.outputs_dir` documents, and two of them failed this gate on its first run.
    `paths.display_path` replaced it.
    """
    from ystwin import paths

    inside = paths.REPO_ROOT / "outputs" / "x.csv"
    outside = pathlib.Path("/tmp/somewhere/x.csv")

    assert paths.display_path(inside) == "outputs/x.csv"
    assert paths.display_path(outside) == str(outside)
    assert "relative_to(paths.REPO_ROOT)" not in "".join(
        path.read_text(encoding="utf-8", errors="replace")
        for path in sorted((REPO / "scripts").glob("*.py")))
