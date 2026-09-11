"""One agent reaching two regulons at two concentrations, which it could not before.

A single saturation term drove every module a stressor targets, so peroxide was forced to
engage Yap1 and Msn2/4 at exactly the same concentration. That is not how any of it works,
and it quietly flattened the dose axis: every module a stressor touched rose in lockstep, so
a dose series traced one ray and no amount of dosing could separate two regulons the same
agent drives.

One ordering is defensible as a default, and it is not the one this file used to assert.
The general stress response is the later, coarser arm: a specific sensor fires first and
Msn2/4 follows once the insult is substantial, which is the ordering Gasch describes.

The other default -- that a metabolite pool always responds below the dose that moves any
regulon, because "peroxide oxidises glutathione stoichiometrically with no threshold to
cross" -- was false chemistry, and the interesting part is that H2O2's two pool arms go
OPPOSITE ways.

Cytosolic H2O2 leads. It is what Yap1's sensor reads, H2O2 reaches it by diffusion, and
HyPer7 detects ~20 uM against 100 uM for the first detectable oxidised Yap1.

Glutathione lags, by more than tenfold. H2O2 reacts with the glutathione thiolate at
18-26 /M/s and only ~2% of GSH is thiolate at cytosolic pH (Winterbourn 1999, PMID
10468205), against ~10^7 /M/s for the peroxiredoxin Tsa1 (Ogusucu 2007, PMID 17210445).
GSSG appears only through a peroxidase, so the pool is the sink at the end of the relay
rather than the first thing peroxide touches, and Yap1 -- wired to the fast catalytic
Orp1/Gpx3 sensor (Delaunay 2002, PMID 12437921) -- fires well before it moves.

One agent, one insult, two pool arms, two orderings. No multiplier can encode that, so
there is no default and each arm carries its own source. The retraction is in
``docs/superseded/pool-before-regulon.md``.
"""

import numpy as np
import pytest

from ystwin.generator.stress_panel import (
    MODULES,
    STRESSORS,
    module_ec50,
    module_response,
)


class TestPerModulePotency:
    def test_the_glutathione_pool_lags_the_regulon_the_agent_drives(self):
        """Yap1 fires through a catalytic relay long before the bulk pool moves.

        Ayer 2013, PMID 23762325 sees no cytosolic E_GSH shift at 0.2 mM over an hour and
        half-oxidation only at 2 mM, against Yap1 nuclear entry from 0.1 mM.
        """
        assert module_ec50("H2O2", "redox") > module_ec50("H2O2", "oxidative")

    def test_the_peroxide_pool_does_lead_the_regulon_it_drives(self):
        """And this is the same agent, going the other way, which is the whole point.

        `peroxide` is cytosolic H2O2 -- what Yap1's sensor reads -- not glutathione. H2O2
        reaches it by diffusion with no enzyme in between, and HyPer7 detects it at ~20 uM
        in BY4742 (Kritsiligkou 2021, PMID 34118234, S. cerevisiae with this module's own
        sensor) against 100 uM for the first detectable oxidised Yap1 (Delaunay 2000, PMID
        11013218).
        """
        assert module_ec50("H2O2", "peroxide") < module_ec50("H2O2", "oxidative")

    def test_so_one_agent_holds_both_orderings_at_once(self):
        """Which is what no single multiplier can express, and why the default is gone."""
        assert (module_ec50("H2O2", "peroxide")
                < module_ec50("H2O2", "oxidative")
                < module_ec50("H2O2", "redox"))

    def test_a_pool_arm_with_no_measured_ec50_is_refused_not_guessed(self):
        """The bug this replaces: a blanket multiplier gave every pool arm a number, so a
        refuted ordering was inherited silently by eight agents."""
        from ystwin.generator.stress_panel import STRESSORS as PANEL, Stressor

        agent = Stressor("test", "mM", 1.0, {"oxidative": 1.0, "redox": 0.5}, "test")
        PANEL["_probe"] = agent
        try:
            with pytest.raises(ValueError, match="no default to fall back on"):
                module_ec50("_probe", "redox")
        finally:
            del PANEL["_probe"]

    def test_a_pool_reached_without_an_enzyme_still_leads(self):
        """Acid entering the cytosol is physical chemistry, so this half of the old rule
        survives -- and it is the only half that does."""
        assert module_ec50("acetic_acid", "ph") < module_ec50("acetic_acid", "ESR")
        assert module_ec50("glucose_starvation", "atp") < STRESSORS["glucose_starvation"].ec50

    def test_the_general_response_needs_more_than_the_specific_one(self):
        assert module_ec50("H2O2", "ESR") > module_ec50("H2O2", "oxidative")
        assert module_ec50("DTT", "ESR") > module_ec50("DTT", "UPR")

    def test_the_specific_target_keeps_the_agents_own_ec50(self):
        assert module_ec50("H2O2", "oxidative") == pytest.approx(STRESSORS["H2O2"].ec50)

    def test_an_explicit_override_wins_over_the_default(self):
        from ystwin.generator.stress_panel import Stressor

        agent = Stressor("test", "mM", 1.0, {"UPR": 1.0, "ESR": 0.5}, "test",
                         target_ec50={"ESR": 0.25})

        assert agent.target_ec50["ESR"] == 0.25

    def test_a_module_the_agent_does_not_target_has_no_potency(self):
        assert module_ec50("BPS", "heat") is None

    def test_an_unknown_agent_is_refused(self):
        with pytest.raises(KeyError, match="unobtainium"):
            module_ec50("unobtainium", "ESR")


class TestItUnflattensTheDoseAxis:
    def test_two_modules_of_one_agent_no_longer_rise_in_lockstep(self):
        """The point of the change: their ratio has to move with dose."""
        low = module_response("H2O2", 0.05)
        high = module_response("H2O2", 2.0)

        assert (low["redox"] / low["oxidative"]) != pytest.approx(
            high["redox"] / high["oxidative"], rel=0.05)

    def test_the_peroxide_arm_separates_too_now_that_it_is_measured(self):
        """Both pool arms unflatten the axis, in opposite directions."""
        low = module_response("H2O2", 0.05)
        high = module_response("H2O2", 2.0)

        assert (low["peroxide"] / low["oxidative"]) != pytest.approx(
            high["peroxide"] / high["oxidative"], rel=0.05)

    def test_the_glutathione_pool_trails_at_low_dose(self):
        """The inverse of what this asserted before. At a twentieth of the lethal dose the
        regulon is well up and the pool has barely moved."""
        response = module_response("H2O2", 0.05)
        pool = response["redox"] / STRESSORS["H2O2"].targets["redox"]
        regulon = response["oxidative"] / STRESSORS["H2O2"].targets["oxidative"]

        assert pool < regulon

    def test_a_dose_series_no_longer_traces_a_single_ray(self):
        doses = np.array([0.02, 0.05, 0.1, 0.3, 0.8])
        vectors = np.array([[module_response("H2O2", d)[m] for m in MODULES] for d in doses])
        normalised = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)

        assert np.linalg.matrix_rank(normalised, tol=1e-6) > 1

    def test_zero_dose_still_moves_nothing(self):
        assert all(v == 0.0 for v in module_response("H2O2", 0.0).values())

    def test_every_module_still_saturates_upward(self):
        for name in ("oxidative", "peroxide", "ESR"):
            low = module_response("H2O2", 0.01)[name]
            high = module_response("H2O2", 0.4)[name]
            assert high > low, name
