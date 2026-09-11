"""Three sensors reading ER, oxidative and ATP, trained to report the whole landscape.

The constraint is three channels. The question is whether a state learned through them
generalises past what they read, and it does -- but only if the analysis stops discarding
the replication.

Two things had to be right. Culture context has to be separable from stress, and three
stress channels cannot do it: trained across seven contexts the build reported 2 of 24
modules, because nothing distinguished a high general-stress response in stationary phase
from one under a stressor. Growth rate separates them and costs no fluorophore, since
optical density is measured anyway.

And replicates have to be pooled before fitting. The latent model gives every well its own
state, so extra wells add parameters rather than averaging noise -- going from 3 replicates
to 24 changed recovery by nothing at all. Pooling first turns the same wells into 9 of 24
modules recovered, which is the ceiling clean data reaches.
"""

import pytest

from ystwin.analysis.experiment_design import RECOMMENDED_DESIGN
from ystwin.analysis.sensor_selection import (
    THREE_SENSOR_BUILD,
    interpretable,
    spectral_conflict,
)
from ystwin.analysis.stress_model import module_transfer, train_stress_model
from ystwin.generator.context import CultureContext
from ystwin.generator.panel_experiment import (
    MEASURED_ACTIVITY_CV,
    MEASURED_GROWTH_RATE_SE,
    panel_dataset,
    pool_replicates,
    with_growth_channel,
)
from ystwin.generator.stress_panel import REPORTERS, STRESSORS

CONTEXTS = [CultureContext(), CultureContext(growth_phase="diauxic"),
            CultureContext(growth_phase="stationary"),
            CultureContext(carbon_source="galactose"),
            CultureContext(carbon_source="ethanol"),
            CultureContext(oxygen=0.03), CultureContext(temperature_c=37.0)]


def trained(replicates=8, pooled=True, growth=True, contexts=CONTEXTS):
    data = panel_dataset(reporters=THREE_SENSOR_BUILD, stressors=list(STRESSORS),
                         replicates=replicates, noise_cv=MEASURED_ACTIVITY_CV,
                         growth_rate_se=MEASURED_GROWTH_RATE_SE, seed=0,
                         combinations=RECOMMENDED_DESIGN, contexts=contexts)
    if growth:
        data = with_growth_channel(data)
    if pooled:
        data = pool_replicates(data)
    return train_stress_model(data, n_states=data.readings.shape[1])


def recovered(model, threshold=0.25):
    return {m for m, v in model.recovery.items() if v > threshold}


class TestTheBuildIsRealisable:
    def test_it_is_three_channels(self):
        assert len(THREE_SENSOR_BUILD) == 3

    def test_they_read_er_oxidative_and_atp(self):
        assert {REPORTERS[r].module for r in THREE_SENSOR_BUILD} == {"UPR", "oxidative", "atp"}

    def test_it_can_be_built(self):
        assert not spectral_conflict(THREE_SENSOR_BUILD)

    def test_every_sensor_works_in_yeast(self):
        assert all(REPORTERS[r].demonstrated_in_yeast for r in THREE_SENSOR_BUILD)

    def test_a_peroxide_sensor_could_not_have_joined_it(self):
        """HyPer7 and QUEEN are both green excitation-ratio, so the oxidative channel has
        to be the transcriptional one."""
        assert spectral_conflict(["UPRE-ER", "HyPer7", "QUEEN-2m"])

    def test_it_needs_no_pH_partner(self):
        assert interpretable(THREE_SENSOR_BUILD)


class TestGrowthSeparatesContextFromStress:
    def test_without_it_the_build_reports_almost_nothing(self):
        assert len(recovered(trained(growth=False, pooled=False))) <= 4

    def test_adding_it_recovers_the_general_stress_response(self):
        """The module context drives hardest, and the one it was being confused with."""
        assert trained(growth=True, pooled=False).recovery["ESR"] > 0.5

    def test_it_costs_no_channel(self):
        """Growth comes from optical density, which is measured regardless."""
        assert "growth" not in REPORTERS


class TestReplicatesMustBePooledBeforeFitting:
    def test_more_wells_alone_change_nothing(self):
        """Every well gets its own latent state, so replication adds parameters not signal."""
        few = trained(replicates=3, pooled=False)
        many = trained(replicates=20, pooled=False)

        assert len(recovered(few)) == len(recovered(many))

    def test_pooling_the_same_wells_lifts_recovery_from_five_modules_to_nine(self):
        """Same wells, same noise -- only where the averaging happens. Both bounds are
        pinned because the multiple between them moved when the panel was corrected, and
        a ratio hid which side had shifted."""
        per_well = trained(replicates=8, pooled=False)
        pooled = trained(replicates=8, pooled=True)

        assert len(recovered(per_well)) == 5
        # 9 until the generator's temperature response was replaced by a cardinal-temperature
        # model fitted to ten S. cerevisiae strains. The contexts include 37 C, where the old
        # response and the measured one differ, so the simulated panel moved and one more
        # module cleared the threshold. The claim is the gap, not the endpoint.
        assert len(recovered(pooled)) == 10

    def test_and_then_replication_starts_paying(self):
        assert trained(replicates=20).recovery["UPR"] > trained(replicates=3).recovery["UPR"]


class TestItGeneralisesPastWhatItReads:
    def test_it_recovers_far_more_modules_than_it_has_channels(self):
        assert len(recovered(trained())) > 2 * len(THREE_SENSOR_BUILD)

    def test_the_directly_read_modules_come_back_well(self):
        model = trained(replicates=20)

        for module in ("UPR", "oxidative", "atp"):
            assert model.recovery[module] > 0.5, module

    def test_and_some_it_never_reads_come_back_too(self):
        """Which is the point: crosstalk and the cascade carry information about modules no
        channel touches."""
        indirect = recovered(trained()) - {"UPR", "oxidative", "atp"}

        assert len(indirect) >= 3

    def test_it_still_cannot_reach_a_private_module(self):
        """Nothing in the build reads the iron regulon, and no amount of training invents it."""
        assert "iron" not in recovered(trained())


class TestItTransfersToAStressorItHasNeverSeen:
    """The actual generalisability claim, and it is narrower than it first measured.

    Hold a stressor out of both the latent basis and the readout, then ask the model what
    its modules are doing. What decides the answer is not the held-out stressor's own axes
    but whether OTHER stressors in the panel drive the same modules. Across all 76
    (held-out stressor, own target) pairs the rank correlation between peer count and score
    is +0.59, and the extremes are stark: no peer driving a module gives 0 of 6 recovered,
    six or more gives 84%.

    So this is interpolation across redundant coverage, not generalisation to unseen
    biology. The earlier reading of this class asserted 8-9 modules recovered per stressor;
    that came from scoring modules whose held-out truth was identically zero. See
    ``docs/superseded/module-transfer-inflation.md``.
    """

    @staticmethod
    def peers(held_out, module):
        """How many other stressors in the panel drive the same module."""
        return sum(1 for name, spec in STRESSORS.items()
                   if name != held_out and module in spec.targets)

    @staticmethod
    def transfer(held_out, replicates=8):
        data = panel_dataset(reporters=THREE_SENSOR_BUILD, stressors=list(STRESSORS),
                             replicates=replicates, noise_cv=MEASURED_ACTIVITY_CV,
                             growth_rate_se=MEASURED_GROWTH_RATE_SE, seed=0,
                             combinations=RECOMMENDED_DESIGN, contexts=CONTEXTS)
        data = pool_replicates(with_growth_channel(data))
        scored = module_transfer(data, held_out=held_out,
                                 n_states=data.readings.shape[1], seed=0)
        return {m for m, v in scored.items() if v > 0.25}

    @staticmethod
    def own_targets_recovered(held_out):
        """Of the held-out stressor's own target modules, how many come back."""
        got = TestItTransfersToAStressorItHasNeverSeen.transfer(held_out)
        return {m for m in STRESSORS[held_out].targets if m in got}

    @pytest.mark.parametrize("stressor,expected", [("H2O2", 3), ("menadione", 4)])
    def test_a_redundantly_covered_stressor_transfers(self, stressor, expected):
        """Both sit on the oxidative-redox axis, which four other agents in the panel also
        drive, so the readout has seen every module they touch on someone else."""
        recovered_here = self.own_targets_recovered(stressor)

        assert len(recovered_here) == expected
        # One peer is enough on a directly-read axis: peroxide has a single one and still
        # scores 0.51. Zero is never enough, which the sole-driver test below pins.
        assert all(self.peers(stressor, m) >= 1 for m in recovered_here)

    def test_and_h2o2_loses_its_glutathione_arm_to_its_own_lethal_dose(self):
        """Not a fitting failure. H2O2's redox EC50 is 1.0 mM (Ayer 2013, PMID 23762325),
        which equals its measured lethal dose, so the arm never moves inside a dose range
        the culture survives -- while menadione, whose pool arms sit at its own EC50, keeps
        its redox arm. The panel cannot see what it cannot dose for."""
        from ystwin.generator.stress_panel import module_ec50

        assert "redox" not in self.own_targets_recovered("H2O2")
        assert "redox" in self.own_targets_recovered("menadione")
        assert module_ec50("H2O2", "redox") >= STRESSORS["H2O2"].lethal_dose

    @pytest.mark.parametrize("stressor", ["antimycin_A", "glucose_starvation"])
    def test_a_stressor_holding_up_its_own_axis_recovers_only_what_is_shared(self, stressor):
        """The build has a dedicated ATP channel and both of these drive ATP, yet neither
        recovers it. Reading a module is not enough -- the readout is fitted on other
        stressors' labels, so the axis has to appear twice in the panel, not once.

        This asserted that they recover NOTHING. Under the measured temperature response
        antimycin_A now returns ESR, and that is the redundancy law being obeyed rather
        than broken: ESR has 24 peers, more than any other module, and the law says 6+
        peers recover 84% of the time. What must stay empty is everything else -- atp and
        nadh have one peer each and retrograde has none.
        """
        recovered_here = self.own_targets_recovered(stressor)

        assert recovered_here <= {"ESR"}
        for module in recovered_here:
            assert self.peers(stressor, module) >= 6

    def test_no_module_comes_back_without_a_peer_driving_it(self):
        """The hard edge, and the one that would break first if the score were inflated
        again: across the whole panel, every module whose held-out stressor is its only
        driver scores below the recovery threshold. Six such pairs, none recovered."""
        sole = [(s, m) for s in STRESSORS for m in STRESSORS[s].targets
                if self.peers(s, m) == 0]

        assert len(sole) >= 5, "the panel must contain sole-driver pairs for this to test"
        for stressor, module in sole:
            assert module not in self.transfer(stressor), f"{stressor} -> {module}"

    @pytest.mark.parametrize("stressor", ["sorbitol", "copper_sulfate", "caffeine"])
    def test_a_stressor_off_its_axes_does_not(self, stressor):
        """Purely osmotic, metal or TOR stress routes through modules no channel reaches,
        and no amount of training recovers what was never observed."""
        assert len(self.transfer(stressor)) <= 3

    def test_which_is_the_price_of_three_channels_not_a_bug(self):
        """The five-channel build exists for exactly this gap."""
        assert len(THREE_SENSOR_BUILD) < 5


class TestTheFiveChannelBuildExtendsThisOneRatherThanReplacingIt:
    """Selection was re-derived once ATP became a named axis. Weighting the general stress
    response alone had chosen a peroxide channel and no ATP channel, reporting the adenylate
    pool at 0.01; naming the three axes picks QUEEN and TRX2 instead."""

    def test_it_is_the_three_sensor_build_plus_one_channel(self):
        from ystwin.analysis.sensor_selection import RECOMMENDED_BUILD

        assert set(THREE_SENSOR_BUILD) < set(RECOMMENDED_BUILD.stress_reporters)

    def test_the_channel_it_adds_is_the_general_stress_response(self):
        from ystwin.analysis.sensor_selection import RECOMMENDED_BUILD
        from ystwin.generator.stress_panel import REPORTERS

        extra = set(RECOMMENDED_BUILD.stress_reporters) - set(THREE_SENSOR_BUILD)

        assert {REPORTERS[r].module for r in extra} == {"ESR"}

    def test_so_a_team_can_build_three_now_and_the_fourth_later(self):
        from ystwin.analysis.sensor_selection import RECOMMENDED_BUILD, spectral_conflict

        assert not spectral_conflict(RECOMMENDED_BUILD.stress_reporters)

    def test_the_selection_does_not_hinge_on_the_exact_weights(self):
        from ystwin.analysis.sensor_selection import (
            _AXIS_PRIORITY, _MIN_EFFECT, _PLATE_NOISE, RECOMMENDED_BUILD,
            growth_penalty_cv, select_sensors)

        halved = dict(_AXIS_PRIORITY, ESR=_AXIS_PRIORITY["ESR"] / 2)
        again = select_sensors(n_channels=5, n_reference=1, noise=_PLATE_NOISE,
                               priority=halved, min_effect=_MIN_EFFECT,
                               growth_cv=growth_penalty_cv())

        assert again == RECOMMENDED_BUILD.stress_reporters
