"""An emulator for the constraint-to-flux map, so state estimation is affordable.

A particle filter over a dynamic GSMM needs one LP per particle per interval. At a
thousand particles and a hundred intervals that is 10^5 solves per culture, inside a
parameter-inference loop. On an enzyme-constrained model with ~8000 reactions that
is hours per culture.

The map from uptake bounds to optimal fluxes is piecewise linear -- it is the
solution of an LP as its right-hand side moves -- so linear interpolation over a
sampled grid is exact wherever the active basis does not change, and degrades
gracefully where it does. It refuses to extrapolate, because outside the sampled box
the basis is unknown and a wrong answer would be indistinguishable from a right one.
"""

import time
from copy import deepcopy

import numpy as np
import pytest

from ystwin.fba.surrogate import build_surrogate

pytestmark = pytest.mark.integration

INPUTS = {"r_1714": (-20.0, -2.0), "r_1992": (-15.0, -1.0)}
OUTPUTS = ["r_2111", "r_1761"]


@pytest.fixture(scope="module")
def _surrogate_template(yeast_gem_factory):
    return build_surrogate(yeast_gem_factory(), inputs=INPUTS, outputs=OUTPUTS, resolution=6)


@pytest.fixture
def surrogate(_surrogate_template):
    return deepcopy(_surrogate_template)


def test_it_reproduces_the_lp_at_the_points_it_was_trained_on(surrogate, yeast_gem):
    point = {"r_1714": -20.0, "r_1992": -15.0}

    with yeast_gem as m:
        for rid, value in point.items():
            m.reactions.get_by_id(rid).lower_bound = value
        truth = m.slim_optimize()

    assert surrogate.predict(point)["r_2111"] == pytest.approx(truth, rel=1e-6)


def test_it_reproduces_the_lp_at_points_between_the_grid_nodes(surrogate, yeast_gem):
    rng = np.random.default_rng(0)
    errors = []
    for _ in range(12):
        point = {rid: float(rng.uniform(lo, hi)) for rid, (lo, hi) in INPUTS.items()}
        with yeast_gem as m:
            for rid, value in point.items():
                m.reactions.get_by_id(rid).lower_bound = value
            truth = m.slim_optimize()
        predicted = surrogate.predict(point)["r_2111"]
        errors.append(abs(predicted - truth) / max(abs(truth), 1e-9))

    assert np.median(errors) < 0.02


def test_a_finer_grid_reduces_the_interpolation_error(yeast_gem):
    rng = np.random.default_rng(1)
    points = [
        {rid: float(rng.uniform(lo, hi)) for rid, (lo, hi) in INPUTS.items()} for _ in range(10)
    ]
    truths = []
    for point in points:
        with yeast_gem as m:
            for rid, value in point.items():
                m.reactions.get_by_id(rid).lower_bound = value
            truths.append(m.slim_optimize())

    def median_error(resolution):
        s = build_surrogate(yeast_gem, inputs=INPUTS, outputs=["r_2111"], resolution=resolution)
        return np.median([
            abs(s.predict(p)["r_2111"] - t) / max(abs(t), 1e-9) for p, t in zip(points, truths)
        ])

    assert median_error(9) <= median_error(3)


def test_it_is_orders_of_magnitude_faster_than_solving_the_lp(surrogate, yeast_gem):
    point = {"r_1714": -11.3, "r_1992": -7.1}

    start = time.perf_counter()
    for _ in range(200):
        surrogate.predict(point)
    surrogate_s = time.perf_counter() - start

    start = time.perf_counter()
    for _ in range(5):
        with yeast_gem as m:
            for rid, value in point.items():
                m.reactions.get_by_id(rid).lower_bound = value
            m.slim_optimize()
    lp_s = (time.perf_counter() - start) / 5

    assert surrogate_s / 200 < lp_s / 20, "surrogate should be at least 20x faster per call"


def test_it_refuses_to_extrapolate_outside_the_sampled_box(surrogate):
    with pytest.raises(ValueError, match="outside the sampled range"):
        surrogate.predict({"r_1714": -50.0, "r_1992": -7.0})


def test_it_reports_which_grid_points_were_infeasible(yeast_gem):
    s = build_surrogate(
        yeast_gem, inputs={"r_1714": (-10.0, -0.001), "r_1992": (-10.0, -0.001)},
        outputs=["r_2111"], resolution=4,
    )

    assert s.n_infeasible >= 0
    assert s.n_samples == 16


def test_every_requested_output_is_predicted(surrogate):
    prediction = surrogate.predict({"r_1714": -10.0, "r_1992": -8.0})

    assert set(prediction) == set(OUTPUTS)


def test_it_can_be_evaluated_on_many_points_at_once(surrogate):
    points = {"r_1714": np.linspace(-18, -4, 32), "r_1992": np.linspace(-12, -3, 32)}

    batch = surrogate.predict_many(points)

    assert batch["r_2111"].shape == (32,)


def test_an_unknown_input_reaction_is_rejected_when_building(yeast_gem):
    with pytest.raises(KeyError, match="r_9999"):
        build_surrogate(yeast_gem, inputs={"r_9999": (-1.0, 0.0)}, outputs=["r_2111"], resolution=3)
