"""A second compound, and the model it falsifies.

`solve_pathway`'s degradation outlet was added on 2026-08-31 with a synthetic test and no
measurement behind it. Boender 2009 supplies one, and supplies something better first: a
falsification of the two-outlet model that needs no fit at all.
"""

from __future__ import annotations

import importlib.util
import pathlib

import pytest

_SCRIPT = (pathlib.Path(__file__).resolve().parents[1]
           / "scripts" / "calibrate_glycogen.py")


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("calibrate_glycogen", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def frame(script):
    return script.states()


class TestTheTwoOutletModelIsFalsifiedWithoutFittingAnything:
    """The part that does not depend on any parameter. Growth dilution alone says the pool
    rises 25x as mu drops from 0.025 to 0.001; it rises 2.12x. No choice of flux changes
    that ratio, because flux cancels out of it."""

    def test_the_pool_does_not_diverge_as_growth_stops(self, script, frame):
        chem = float(frame.loc["chemostat_D0025", script.CONTENT])
        ret = float(frame.loc["retentostat_22d", script.CONTENT])

        assert 1.5 < ret / chem < 3.0

    def test_growth_dilution_alone_would_predict_an_order_of_magnitude_more(self, script,
                                                                            frame):
        chem_mu = float(frame.loc["chemostat_D0025", "mu_per_h"])
        ret_mu = float(frame.loc["retentostat_22d", "mu_upper_bound_per_h"])
        chem = float(frame.loc["chemostat_D0025", script.CONTENT])
        ret = float(frame.loc["retentostat_22d", script.CONTENT])

        assert (chem_mu / ret_mu) / (ret / chem) > 8.0


class TestTheThirdOutletIsCalibratedNotAsserted:
    def test_the_fitted_turnover_is_physically_plausible(self, script, frame):
        """GPH1 and SGA1 turning over a storage pool on a timescale of tens of hours."""
        flux, k_deg = script.calibrate(frame)

        assert 0.005 < k_deg < 0.10
        assert flux > 0

    def test_the_solver_reproduces_the_states_it_was_fitted_on(self, script, frame):
        from ystwin.pathway.solve import NodeKinetics, solve_pathway
        from ystwin.pathway.spec import load_pathway

        flux, k_deg = script.calibrate(frame)
        spec = load_pathway("glycogen")
        kinetics = {"glycogen": NodeKinetics(degradation_rate_per_h=k_deg)}

        for state, mu_col in (("chemostat_D0025", "mu_per_h"),
                              ("retentostat_22d", "mu_upper_bound_per_h")):
            mu = float(frame.loc[state, mu_col])
            solved = solve_pathway(spec, flux, mu, kinetics).terminal.content_mmol_per_gdcw
            assert solved == pytest.approx(float(frame.loc[state, script.CONTENT]), rel=1e-6)

    def test_the_fit_is_saturated_and_the_docstring_says_so(self, script):
        """Two states, two unknowns: it identifies the pair and cannot fail on them. A test
        that let this read as validation would be the misleading kind."""
        prose = " ".join(script.__doc__.split())   # the phrase spans a line wrap

        assert "saturated" in prose
        assert "cannot be wrong on its own data" in prose

    def test_the_out_of_sample_point_is_a_different_limitation(self, script, frame):
        """Nitrogen starvation, not glucose limitation. Within 36%, and reported as a
        comparison rather than a held-out test."""
        flux, k_deg = script.calibrate(frame)
        measured = float(frame.loc["n_starved_shake_flask", script.CONTENT])

        assert 0.6 < (flux / k_deg) / measured < 1.0
