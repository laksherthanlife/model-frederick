"""The mechanistic seam in the data generator: what it moves, and what it must not.

`generator/panel_experiment.py::panel_dataset` is what produces the training data, and until
this seam existed it imported nothing from `mech/`. These tests pin the two halves of that
change, and the second half is the one that matters more.

**Off by default, and byte-identical when off.** Every committed table and every pinned
number in this repository was computed on the algebraic path. So the assertion here is
``np.array_equal`` on the readings, not ``approx``: with ``mechanism=None`` nothing moved,
and with a corner declared nothing moved outside
:data:`~ystwin.generator.panel_experiment.MECHANISED_CHANNELS` either. A mechanistic path
that quietly perturbs the other 23 stressors would make every downstream diff unreadable.

**The geometry result, stated as it came out.** The median dose x time ray share over the
25 stressors is 0.9986 on the algebraic path and 0.9985 with the mechanism on. It did not
move, and the reason is arithmetic rather than interesting: two of the 25 stressors have a
mech block and 23 do not, so a median over 25 cannot see them. What did move is H2O2, from
0.9951 to 0.5838-0.8560 across the swept corners, and DTT from 0.9993 to 0.9949-0.9979.
Quoting the median as evidence that the mechanism changed the dataset's geometry would be
false; quoting H2O2 without saying it is one stressor in 25 would be worse.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from ystwin.generator.panel_experiment import (
    MECHANISED_CHANNELS,
    OBSERVED_ACTIVITY_CV,
    Mechanism,
    _cascade_column,
    _COMMITTED_UPR_FIT,
    _mechanistic_direct,
    dose_time_ray_shares,
    panel_dataset,
    pool_replicates,
    with_growth_channel,
)
from ystwin.generator.stress_panel import (
    MODULES,
    STRESSORS,
    module_response,
    transcriptional_reporters,
)
from ystwin.mech import oxidative as mech_oxidative
from ystwin.mech import upr as mech_upr

GEOMETRY_KWARGS = dict(reporters=transcriptional_reporters(), seed=0, noise_cv=0.0,
                       replicates=1, read_times_h=[0.5, 1.0, 2.0, 4.0, 8.0])
"""The design the geometry claim is quoted on, exactly as the build list states it."""

COMMITTED_RAY_SHARES = (0.9986, 0.9945, 0.9998)
"""Median, minimum and maximum ray share of the algebraic generator, to 4 decimal places.

The pre-build baseline. If this moves, the algebraic path moved, and every committed table
computed on it is stale -- which is the failure this file exists to catch early.
"""


@pytest.fixture(scope="module")
def algebraic():
    return panel_dataset(**GEOMETRY_KWARGS)


@pytest.fixture(scope="module")
def corners():
    return Mechanism.corners()


# --------------------------------------------------------------------------------------
# Off by default, and inert everywhere it is not declared
# --------------------------------------------------------------------------------------

def test_the_default_is_off_and_reproduces_the_committed_baseline(algebraic):
    shares = np.array(list(dose_time_ray_shares(algebraic).values()))
    got = (float(np.median(shares)), float(shares.min()), float(shares.max()))
    assert [round(value, 4) for value in got] == list(COMMITTED_RAY_SHARES)


def test_a_declared_corner_leaves_every_other_stressor_byte_identical(algebraic, corners):
    """Not ``approx``: the untouched stressors must not move by one bit.

    The mechanistic branch is entered only where a block drives a module, and the algebraic
    activity is subtracted before the quadrature runs, so the incumbent's closed form is what
    produces every other reading rather than a quadrature that happens to agree with it.
    """
    driven = {stressor for stressor, _module in MECHANISED_CHANNELS}
    keep = ~np.isin(algebraic.labels, list(driven))
    for spec in corners:
        moved = panel_dataset(mechanism=spec, **GEOMETRY_KWARGS)
        assert np.array_equal(algebraic.readings[keep], moved.readings[keep])
        assert np.array_equal(algebraic.modules, moved.modules)
        assert np.array_equal(algebraic.labels, moved.labels)
        assert np.array_equal(algebraic.doses, moved.doses)


_ROUTE_NOISE = 1e-12
"""Largest reading difference attributable to recomputation rather than to the mechanism.

The mechanistic path rebuilds an untouched column through a different but algebraically
identical route, which costs at most 2.2e-16 here -- measured across every reporter under
H2O2. The smallest change the mechanism actually makes is 5.3e-2, on the xenobiotic arm two
cascade edges from the one it drives. The two are ten orders of magnitude apart, so an exact
comparison tests the floating-point route and not the seam.
"""


def test_the_untouched_modules_keep_the_incumbent_closed_form(corners):
    """Inside a driven stressor, a reporter reading no module of the cascade keeps its value.

    H2O2 drives the general stress response as well as Yap1, and STRE-general reads only the
    former. The mechanism must leave that column exactly where the Hill put it, or the seam
    is leaking quadrature error into channels it does not claim.
    """
    kwargs = dict(GEOMETRY_KWARGS, stressors=["H2O2"])
    base = panel_dataset(**kwargs)
    column = base.reporters.index("STRE-general")
    for spec in corners:
        moved = panel_dataset(mechanism=spec, **kwargs)
        assert np.max(np.abs(base.readings[:, column] - moved.readings[:, column])) <= _ROUTE_NOISE


def test_pooling_and_the_growth_channel_still_work_on_a_mechanistic_dataset():
    data = panel_dataset(mechanism=Mechanism.midpoint(), **GEOMETRY_KWARGS)
    assert np.all(np.isfinite(data.readings))
    assert pool_replicates(data).readings.shape[0] == len(set(zip(data.labels, data.doses)))
    assert with_growth_channel(data).readings.shape[1] == len(data.reporters) + 1


# --------------------------------------------------------------------------------------
# The geometry, both numbers
# --------------------------------------------------------------------------------------

def test_the_median_ray_share_does_not_move_and_that_is_the_finding(algebraic, corners):
    """23 of the 25 stressors have no mech block, so a median over 25 cannot move.

    This is the honest headline and it is asserted rather than described: the mechanism is
    wired, it works, and the panel-level median it was measured against is unchanged at the
    fourth decimal place. Reporting it as a success would be the worst outcome available.
    """
    before = np.array(list(dose_time_ray_shares(algebraic).values()))
    for spec in corners:
        after = np.array(list(dose_time_ray_shares(
            panel_dataset(mechanism=spec, **GEOMETRY_KWARGS)).values()))
        assert abs(float(np.median(after)) - float(np.median(before))) < 5e-4
        assert round(float(after.max()), 4) == round(float(before.max()), 4)


def test_exactly_the_two_declared_stressors_move(algebraic, corners):
    before = dose_time_ray_shares(algebraic)
    driven = {stressor for stressor, _module in MECHANISED_CHANNELS}
    for spec in corners:
        after = dose_time_ray_shares(panel_dataset(mechanism=spec, **GEOMETRY_KWARGS))
        assert {name for name in before if before[name] != after[name]} == driven


def test_the_peroxide_ray_share_collapses_and_the_dtt_one_barely_moves(algebraic, corners):
    """The two blocks are not equally consequential, and the numbers say so.

    H2O2 loses between 14% and 42% of its surface off the leading ray because the pipetted
    dose is consumed inside the read and the whole Yap1 arm of the panel then decays while
    the general stress response it is read beside does not. DTT moves too, and by a factor
    on the off-ray energy rather than a rounding error, but it stays above 0.99: the UPR
    module reaches exactly one of the nineteen promoters and drives nothing downstream.
    """
    before = dose_time_ray_shares(algebraic)
    peroxide, dithiothreitol = [], []
    for spec in corners:
        after = dose_time_ray_shares(panel_dataset(mechanism=spec, **GEOMETRY_KWARGS))
        peroxide.append(after["H2O2"])
        dithiothreitol.append(after["DTT"])
    assert 0.55 < min(peroxide) and max(peroxide) < 0.90
    assert (1.0 - min(peroxide)) / (1.0 - before["H2O2"]) > 25.0
    assert 0.99 < min(dithiothreitol) < max(dithiothreitol) < before["DTT"]
    assert (1.0 - min(dithiothreitol)) / (1.0 - before["DTT"]) > 2.5


def test_the_peroxide_reporter_peaks_and_falls_where_the_incumbent_only_rises():
    """The structural difference behind the geometry, on one trace.

    The algebraic path holds the promoter at a level the pipetted dose sets and never lets
    it go, so a stable fluorophore climbs for the whole read. The block consumes the dose in
    minutes, after which the same stable fluorophore is diluted by a growing culture. One is
    monotone over 8 h and the other is not, and no rescaling of the first produces the
    second -- which is what a ray share of 0.58 is measuring.
    """
    kwargs = dict(GEOMETRY_KWARGS, stressors=["H2O2"])
    base = panel_dataset(**kwargs)
    moved = panel_dataset(mechanism=Mechanism.midpoint(), **kwargs)
    column = base.reporters.index("TRX2-oxidative")
    top = base.doses.max()
    old = base.readings[base.doses == top][:, column]
    new = moved.readings[moved.doses == top][:, column]
    assert np.all(np.diff(old) > 0)
    assert np.argmax(new) not in (0, len(new) - 1)
    assert new[-1] < new.max()


def test_the_grid_is_fine_enough_to_be_a_tolerance_and_not_a_model_choice():
    """Doubling the quadrature grid moves nothing a plate could see."""
    reference = panel_dataset(
        mechanism=replace(Mechanism.midpoint(), grid_points=3201), **GEOMETRY_KWARGS)
    coarse = panel_dataset(mechanism=Mechanism.midpoint(), **GEOMETRY_KWARGS)
    relative = np.max(np.abs(coarse.readings - reference.readings) / reference.readings)
    assert relative < 1e-4 * OBSERVED_ACTIVITY_CV


# --------------------------------------------------------------------------------------
# The override is the block's own arithmetic, not a re-implementation
# --------------------------------------------------------------------------------------

def test_the_peroxide_arm_starts_exactly_where_the_incumbent_hill_starts():
    """At ``t = 0`` the two paths are the same number, and that is what makes them nested.

    ``yap1_nuclear_fraction`` is ``E/(E + K_ex)`` with ``K_ex`` BOUNDED at 0.15 mM from
    Goulev 2017, and `stress_panel.STRESSORS['H2O2'].ec50` is 0.15 mM from a different
    source entirely. They agree, so the only thing the mechanism adds on this channel is
    that ``E`` moves. If either constant is ever revised this test fails, which is the point:
    the nesting is a property of two independent numbers agreeing, not of the code.
    """
    grid = np.linspace(0.0, 4.14, 401)
    spec = Mechanism.midpoint()
    for dose in (0.05, 0.15, 0.55):
        module, activation = _mechanistic_direct(spec, "H2O2", dose, 0.32, grid)
        assert module == "oxidative"
        assert activation[0] == pytest.approx(module_response("H2O2", dose)["oxidative"],
                                              rel=1e-12)


def test_the_peroxide_arm_carries_the_cascade_the_hill_it_replaces_carried():
    """Overriding a direct activity must travel `stress_panel`'s cascade, not bypass it."""
    column = _cascade_column("oxidative")
    order = list(MODULES)
    assert column[order.index("oxidative")] == 1.0
    assert column[order.index("proteasome")] == pytest.approx(0.45)
    assert column[order.index("iron")] == pytest.approx(0.10)
    assert column[order.index("xenobiotic")] == pytest.approx(0.45 * 0.2)
    assert column[order.index("ESR")] == 0.0

    kwargs = dict(GEOMETRY_KWARGS, stressors=["H2O2"])
    base = panel_dataset(**kwargs)
    moved = panel_dataset(mechanism=Mechanism.midpoint(), **kwargs)
    # Real cascade effects are >= 5e-2; recomputing an untouched column costs <= 2.2e-16.
    changed = {name for index, name in enumerate(base.reporters)
               if np.max(np.abs(base.readings[:, index] - moved.readings[:, index])) > _ROUTE_NOISE}
    assert changed == {"TRX2-oxidative", "PACE-proteasome", "FeRE-iron", "PDRE-xenobiotic"}


def test_the_overridden_modules_are_the_ones_module_response_reports_directly():
    """The invariant the override rests on: nothing drives oxidative or UPR in the cascade.

    If either acquires a driver, `module_response` stops reporting its direct activity and
    subtracting it would double-count. The code raises there; this states it as a fact about
    the table so the failure arrives at the table rather than deep in a dataset.
    """
    for _stressor, module in MECHANISED_CHANNELS:
        assert not MODULES[module].driven_by


def test_the_upr_arm_is_the_blocks_own_occupancy_normalised_by_its_own_saturation():
    """Zero at zero dose, and it reaches full UPRE occupancy only as the load saturates."""
    grid = np.linspace(0.0, 4.14, 401)
    spec = Mechanism.midpoint()
    assert _mechanistic_direct(spec, "DTT", 0.0, 0.32, grid) is None
    _module, low = _mechanistic_direct(spec, "DTT", 0.05, 0.32, grid)
    _module, high = _mechanistic_direct(spec, "DTT", 0.85, 0.32, grid)
    assert low[0] == pytest.approx(0.0, abs=1e-9)
    assert high.max() > 10.0 * low.max()
    assert high.max() <= STRESSORS["DTT"].targets["UPR"]


# --------------------------------------------------------------------------------------
# Every axis this path opens is declared, and every refusal is still a refusal
# --------------------------------------------------------------------------------------

def test_the_committed_scalars_are_the_ones_mech_upr_committed():
    """The four restated scalars are pinned against `mech/upr.py`'s own record of them."""
    doc = mech_upr.PlateUpr.__doc__
    assert "0.2424" in doc and "0.9900" in doc and "8.0000" in doc and "1.8093" in doc
    assert _COMMITTED_UPR_FIT == {"basal_load": 0.2424, "top_load": 0.99, "hill_n": 8.0,
                                  "basal_share": 1.8093}
    low, high = mech_upr.IRE1_HILL_N.bounds
    assert _COMMITTED_UPR_FIT["hill_n"] == high, "hill_n sits on its own upper bound"
    assert Mechanism.midpoint().od600 == float(mech_oxidative.PLATE_INOCULUM_OD600)


def test_the_corners_are_the_corners_of_the_two_swept_axes(corners):
    cells = mech_oxidative.CELLS_PER_ML_PER_OD600.bounds
    folding = mech_upr.FOLDING_RATE_PER_H.bounds
    assert len(corners) == 4
    assert {spec.cells_per_ml_per_od600 for spec in corners} == set(cells)
    assert {spec.folding_rate_per_h for spec in corners} == set(folding)


def test_a_point_off_a_swept_axis_is_refused():
    low, high = mech_oxidative.CELLS_PER_ML_PER_OD600.bounds
    with pytest.raises(ValueError, match="swept"):
        replace(Mechanism.midpoint(), cells_per_ml_per_od600=high * 2.0)
    with pytest.raises(ValueError, match="swept"):
        replace(Mechanism.midpoint(), cells_per_ml_per_od600=low / 2.0)
    with pytest.raises(ValueError, match="Korennykh|bracket"):
        replace(Mechanism.midpoint(), hill_n=12.0)


def test_the_dtt_entry_stays_refused_and_only_a_declaration_gets_past_it():
    """The dose axis is not measured, and this seam does not quietly supply one."""
    assert mech_upr.DTT_ENTRY.tag == "REFUSED"
    with pytest.raises(Exception):
        float(mech_upr.DTT_ENTRY)
    with pytest.raises(mech_upr.EntryRefused):
        mech_upr.DoseEntry(top_load=0.5, top_dose_mM=1.0)
    assert Mechanism.midpoint().upr_entry().declared is True


def test_an_inadmissible_er_load_is_refused():
    with pytest.raises(ValueError, match="basal_load"):
        replace(Mechanism.midpoint(), basal_load=0.995)
    with pytest.raises(ValueError, match="basal_load"):
        replace(Mechanism.midpoint(), top_load=0.1, basal_load=0.2424)


def test_a_dtt_dose_above_the_measured_lethal_one_is_refused_rather_than_extrapolated():
    lethal = float(mech_upr.DTT_LETHAL_MM)
    over = lethal * 1.2 / STRESSORS["DTT"].ec50
    with pytest.raises(ValueError, match="lethal"):
        panel_dataset(stressors=["DTT"], doses=[over], replicates=1, noise_cv=0.0,
                      read_times_h=[4.14], mechanism=Mechanism.midpoint())
    assert panel_dataset(stressors=["DTT"], doses=[over], replicates=1, noise_cv=0.0,
                         read_times_h=[4.14]).readings.shape[0] == 1


def test_the_mechanism_needs_a_read_schedule():
    """A mechanistic reading is a time course; the timeless steady state does not exist."""
    with pytest.raises(ValueError, match="read_times_h"):
        panel_dataset(stressors=["H2O2"], replicates=1, noise_cv=0.0,
                      mechanism=Mechanism.midpoint())
    with pytest.raises(ValueError, match="finite and non-negative"):
        panel_dataset(stressors=["H2O2"], replicates=1, noise_cv=0.0, read_times_h=[0.0],
                      mechanism=Mechanism.midpoint())


def test_the_declared_dose_axis_brackets_the_measured_induction_rather_than_matching_it():
    """What the axis costs, MEASURED, instead of a sentence saying it is declared.

    The plates put the UPRE fold across the whole DTT ladder at 1.38-1.47
    (`stress_panel._DEFAULT_BASAL`), and the algebraic generator is calibrated to sit inside
    that. The mechanistic arm's fold at 4.14 h is a function of ``top_load``, which is the
    axis `mech/upr.py::DTT_ENTRY` refuses to supply: it runs from 1.00 at the bottom of the
    admissible range to 1.61 at the committed fit's 0.99, which is on its own box. The
    measurement is inside that span and does not pick a point in it.
    """
    kwargs = dict(reporters=transcriptional_reporters(), seed=0, noise_cv=0.0, replicates=1,
                  read_times_h=[4.14], stressors=["DTT"])
    column = transcriptional_reporters().index("UPRE-ER")
    basal = 0.9
    folds = {}
    for top_load in (0.3, 0.5, 0.7, 0.99):
        readings = panel_dataset(
            mechanism=replace(Mechanism.midpoint(), top_load=top_load), **kwargs
        ).readings[:, column]
        folds[top_load] = float(readings.max() / basal)
    assert folds[0.3] < 1.05
    assert folds[0.5] < 1.38 < folds[0.7]
    assert 1.55 < folds[0.99] < 1.70
    assert folds[0.3] < folds[0.5] < folds[0.7] < folds[0.99]
