"""What four channels actually resolve, as opposed to what the panel simulates.

The panel carries 25 stressors over 24 state variables. The build reads four modules, and
the honest characterisation is not "four stresses" but three: through these channels the 25
stressors collapse onto three axes carrying 97% of the variance, and within each axis they
are indistinguishable.

That is the rank ceiling doing what it says rather than anything going wrong, but it is
worth stating plainly. A build that cannot separate hydrogen peroxide from menadione should
not be described as measuring oxidative stress in general, and sixteen of the panel's
stressors move it barely at all.

Which nine it does see moved when selection was re-derived around the three named axes.
Trading the proteasome channel for an ATP one lost MG132 and MMS, whose signal was the PACE
element, and gained antimycin_A and copper_sulfate. They do not land where the names
suggest: copper joins the oxidative cluster, which is right for a Fenton metal, and the
adenylate axis is antimycin_A beside glucose starvation -- a respiratory chain inhibitor and
a culture with nothing to respire, which is the same state reached two ways.
"""

import numpy as np

from ystwin.analysis.sensor_selection import RECOMMENDED_BUILD
from ystwin.generator.stress_panel import MODULES, STRESSORS, module_response, reporter_loadings

DETECTABLE = ["H2O2", "menadione", "tunicamycin", "antimycin_A", "DTT", "heat", "diamide",
              "glucose_starvation", "copper_sulfate"]


def signature(stressor):
    order = list(MODULES)
    loadings = reporter_loadings(RECOMMENDED_BUILD.stress_reporters)
    activity = np.array([module_response(stressor, STRESSORS[stressor].ec50)[m] for m in order])
    return loadings @ activity


def unit(stressor):
    raw = signature(stressor)
    return raw / max(float(np.linalg.norm(raw)), 1e-12)


class TestHowMuchOfThePanelItSees:
    def test_most_stressors_barely_move_it(self):
        weak = [s for s in STRESSORS if float(np.abs(signature(s)).sum()) < 0.25]

        assert len(weak) > len(STRESSORS) / 2

    def test_the_ones_it_does_see_are_the_ones_expected(self):
        strong = {s for s in STRESSORS if float(np.abs(signature(s)).sum()) >= 0.25}

        assert strong == set(DETECTABLE)

    def test_an_iron_chelator_is_invisible_to_it(self):
        """Nothing in the build reads the iron regulon."""
        assert float(np.abs(signature("BPS")).sum()) < 0.1


class TestItResolvesThreeAxesNotFour:
    def test_three_axes_carry_most_of_the_variance(self):
        rows = np.array([unit(s) for s in DETECTABLE])
        share = np.linalg.svd(rows, compute_uv=False) ** 2
        share = share / share.sum()

        assert share[:3].sum() > 0.95

    def test_the_fourth_axis_is_nearly_empty(self):
        rows = np.array([unit(s) for s in DETECTABLE])
        share = np.linalg.svd(rows, compute_uv=False) ** 2

        assert (share / share.sum())[3] < 0.1


class TestWithinAnAxisStressorsAreIndistinguishable:
    def test_peroxide_and_menadione_look_the_same(self):
        assert float(unit("H2O2") @ unit("menadione")) > 0.95

    def test_tunicamycin_and_dtt_look_the_same(self):
        assert float(unit("tunicamycin") @ unit("DTT")) > 0.95

    def test_the_oxidative_cluster_swallows_copper(self):
        """Copper is a Fenton metal, so the oxidative channel is where it should land."""
        for other in ("H2O2", "menadione", "diamide"):
            assert float(unit("copper_sulfate") @ unit(other)) > 0.95, other

    def test_the_adenylate_axis_is_antimycin_beside_starvation(self):
        """A blocked respiratory chain and nothing to respire reach the same state."""
        assert float(unit("antimycin_A") @ unit("glucose_starvation")) > 0.85

    def test_and_that_axis_is_clear_of_the_oxidative_one(self):
        assert float(unit("antimycin_A") @ unit("H2O2")) < 0.2

    def test_but_the_axes_are_separable_from_each_other(self):
        """Which is what makes it three axes rather than one."""
        assert float(unit("H2O2") @ unit("DTT")) < 0.4
