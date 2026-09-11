from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from scripts import run_hog_learning as runner
from ystwin.fba import dynamic_rates, product_panel, solver


def test_gem_probe_uses_the_real_coupled_rates_interface(monkeypatch):
    monkeypatch.setattr(runner.paths, "ec_yeast_gem", lambda: "fixture.xml")
    monkeypatch.setattr(solver, "load_model", lambda _: (object(), solver.SolverSettings("glpk", 1e-7, "fixture")))
    reference = pd.DataFrame({"orf": ["YDL022W"], "gene": ["GPD1"]})
    monkeypatch.setattr(product_panel, "prepare_panel_model", lambda _: (object(), reference))
    monkeypatch.setattr(product_panel, "apply_relative_abundances", lambda model, inputs, ratios: (
        ratios["YDL022W"], pd.DataFrame({"applied_cap_mmol_per_gdcw": [ratios["YDL022W"] * 1e-5]})))

    def install(model, name, *, output_compartment):
        assert name == "glycerol" and output_compartment == "c"
        return model, SimpleNamespace(reaction_id="DM_glycerol", output_basis="intracellular retention")

    def build(ratio, *args, **kwargs):
        assert kwargs["product_reaction_id"] == "DM_glycerol"
        rates = dynamic_rates.CoupledRates(
            np.array([0.0, 1.0, 2.0]), np.array([[0, 0, 0], [0.1, 1, ratio], [0.1, 2, 2 * ratio]]),
            np.array([False, True, True]), 0.5, "DM_glycerol", ("infeasible", "optimal", "optimal"))
        return rates, SimpleNamespace(ngam=1.0)

    monkeypatch.setattr(product_panel, "install_product", install)
    monkeypatch.setattr(dynamic_rates, "build_network_rates", build)
    signals = pd.DataFrame({"time_s": [0, 7200], "nacl_molar": [0.4, 0.4], "gpd1_capacity_multiplier": [1.0, 2.0]})
    result = runner._gem_capacity_check(signals)
    np.testing.assert_allclose(result.intracellular_glycerol_capacity_mmol_per_gdcw_h, [2, 4])
    assert result.product_reaction_id.eq("DM_glycerol").all()


def test_report_weights_measurement_blocks_not_technical_trace_counts():
    frame = pd.DataFrame({"split": ["heldout_dose"] * 3, "experiment_id": ["a", "a", "b"],
                          "observations": [13] * 3, "fitted_shape_rmse": [0.1, 0.3, 0.4],
                          "published_shape_rmse": [0.1, 0.3, 0.4], "persistent_shape_rmse": [0.1, 0.3, 0.4]})
    result = runner._summary(frame)["heldout_dose"]
    assert result["fitted_shape_rmse"] == pytest.approx(0.3)
    assert result["measurement_blocks"] == 2 and result["trace_count"] == 3


def test_runner_preserves_existing_output_artifacts(tmp_path):
    marker = tmp_path / "existing.csv"
    marker.write_text("preserve")
    with pytest.raises(FileExistsError, match="fresh"):
        runner.main(["--output-dir", str(tmp_path)])
    assert marker.read_text() == "preserve"
