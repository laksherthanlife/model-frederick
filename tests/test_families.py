"""A family has to be a different structure, and the baseline has to be the same one.

Two properties carry the whole cross-family argument and everything here serves one of them.

The baseline family must reproduce `panel_dataset` exactly. If it does not, every per-family
score is measured from a different origin than the numbers already reported in
`docs/NULL_RESULTS.md`, and a displacement between families stops being interpretable.

Every other family must differ from the baseline *in direction* and not merely in scale. A
configuration that rescales the readings is a units change: any subspace method is invariant
to it, so a robustness claim that survived only such families would have survived nothing.
The test that pins this is the one that fits the best scalar between a family and the baseline
and demands a residual -- it is the operational form of "structural, not a re-draw".

The declarations are tested too, which is unusual and deliberate. `perturbs` and
`plausible_because` are the scientific content of a family: a robustness claim is worth the
plausibility of the alternatives it survived, so a family that cannot defend itself makes the
claim cheaper rather than stronger.
"""

from __future__ import annotations

import copy

import numpy as np
import pytest

from ystwin.analysis.experiment_design import RECOMMENDED_DESIGN
from ystwin.generator.families import (
    FAMILY_BUILDERS,
    PanelFamily,
    baseline_family,
    build_family,
    family_dataset,
    family_names,
    module_activity,
    sample_family,
)
from ystwin.generator.panel_experiment import (
    MEASURED_ACTIVITY_CV,
    MEASURED_GROWTH_RATE_SE,
    panel_dataset,
)
from ystwin.generator.stress_panel import MODULES, REPORTERS

READERS = ["STRE-general", "UPRE-ER", "TRX2-oxidative", "HSE-heat",
           "STRE-osmotic", "PACE-proteasome", "FeRE-iron", "Xbox-dna"]
PANEL = ["DTT", "H2O2", "heat", "NaCl", "MG132", "glucose_starvation"]
DOSES = (0.25, 0.5, 1.0, 2.0)

VARIANTS = [name for name in family_names() if name != "baseline"]


def clean(family, **overrides) -> np.ndarray:
    """Noise-free readings, so a comparison between families is about the families."""
    settings = dict(reporters=READERS, stressors=PANEL, doses=DOSES, replicates=1,
                    noise_cv=0.0)
    return family_dataset(family, **(settings | overrides)).readings


def scalar_residual(observed: np.ndarray, reference: np.ndarray) -> float:
    """Relative error left after the best single rescaling of `reference`.

    Zero means the two configurations differ by a change of units, which no subspace method
    can see. Anything else is a change of direction, which is what a family is for.
    """
    a, b = observed.ravel(), reference.ravel()
    scale = float(a @ b) / float(b @ b)
    return float(np.linalg.norm(a - scale * b) / np.linalg.norm(a))


class TestBaselineIsTheSameGenerator:
    """The origin the other families are measured from."""

    def test_it_reproduces_panel_dataset_reading_for_reading(self):
        """Not approximately: the same arithmetic in the same order, noise draws included.

        A drift here would leave every per-family delta measured against a generator that is
        not the one the reported numbers came from.
        """
        settings = dict(reporters=READERS, stressors=PANEL, doses=DOSES, replicates=3,
                        noise_cv=MEASURED_ACTIVITY_CV,
                        growth_rate_se=MEASURED_GROWTH_RATE_SE, seed=0,
                        combinations=RECOMMENDED_DESIGN)
        reference = panel_dataset(**settings)
        family = family_dataset(baseline_family(), **settings)

        np.testing.assert_allclose(family.readings, reference.readings, rtol=1e-12)
        np.testing.assert_allclose(family.modules, reference.modules, rtol=0, atol=1e-15)
        assert np.array_equal(family.labels, reference.labels)
        np.testing.assert_allclose(family.doses, reference.doses, rtol=1e-12)

    def test_the_linear_solve_agrees_with_the_sequential_cascade(self):
        """The literature cascade is acyclic, so one pass in topological order equals a solve.

        This is what licenses replacing `_propagate` with a solve at all: the solve is only a
        generalisation if it agrees with the special case.
        """
        from ystwin.generator.stress_panel import combination_response

        doses = {"DTT": 1.0, "heat": 4.0, "MG132": 20.0}
        solved = module_activity(baseline_family(), doses)
        sequential = combination_response(doses)

        for module in MODULES:
            assert solved[module] == pytest.approx(sequential[module], abs=1e-14)

    def test_an_undosed_well_reads_the_promoter_floor(self):
        readings = clean(baseline_family(), doses=(0.0,))
        floor = np.array([REPORTERS[r].basal for r in READERS])

        np.testing.assert_allclose(readings, np.tile(floor, (len(readings), 1)), rtol=1e-12)


class TestEveryFamilyIsStructural:
    @pytest.mark.parametrize("name", VARIANTS)
    def test_it_is_not_the_baseline_rescaled(self, name):
        """The load-bearing test: a change of direction, not a change of units."""
        residual = scalar_residual(clean(build_family(name, seed=1)),
                                   clean(baseline_family()))

        assert residual > 0.01, (
            f"{name} differs from the baseline by a factor of "
            f"{1.0 - residual:.4f} and almost nothing else, so any method invariant to "
            "scale would score identically on it")

    @pytest.mark.parametrize("name", VARIANTS)
    def test_it_moves_the_readings_at_all(self, name):
        assert not np.allclose(clean(build_family(name, seed=1)), clean(baseline_family()))


class TestWhatEachFamilyPerturbs:
    def test_dropping_an_edge_removes_that_edge_and_no_other(self):
        baseline, dropped = baseline_family(), build_family("edge_dropped")
        removed = {(m, d) for m, edges in baseline.drivers.items() for d in edges} - {
            (m, d) for m, edges in dropped.drivers.items() for d in edges}

        assert removed == {("proteasome", "heat")}

    def test_dropping_the_heat_edge_lowers_the_proteasome_under_heat(self):
        """The edge exists to carry Hsf1 into Rpn4, so removing it must show up there."""
        doses = {"heat": 4.0}
        baseline = module_activity(baseline_family(), doses)
        dropped = module_activity(build_family("edge_dropped"), doses)

        assert dropped["proteasome"] < baseline["proteasome"]
        assert dropped["heat"] == pytest.approx(baseline["heat"])

    def test_adding_an_edge_creates_activity_where_the_panel_has_none(self):
        """Glucose withdrawal touches no oxidative target, so a non-zero reading is the edge.

        Cleaner than a magnitude comparison: under the literature topology the number is
        exactly zero, so anything above it can only have come through the added edge.
        """
        doses = {"glucose_starvation": 0.6}
        baseline = module_activity(baseline_family(), doses)
        added = module_activity(build_family("edge_added"), doses)

        assert baseline["oxidative"] == 0.0
        assert added["oxidative"] > 0.0
        assert added["proteasome"] > baseline["proteasome"]

    def test_a_module_specific_limb_is_not_a_shared_viability_factor(self):
        """A shared falling limb scales every module alike; a per-module one cannot.

        Measured above the turnover, where the limb dominates: if the ratios to baseline were
        equal across modules the family would be a rescaling and the geometry unchanged.
        """
        doses = {"DTT": 2.0}
        baseline = module_activity(baseline_family(), doses)
        biphasic = module_activity(build_family("biphasic"), doses)
        ratios = [biphasic[m] / baseline[m] for m in MODULES if abs(baseline[m]) > 1e-9]

        assert max(ratios) - min(ratios) > 0.05
        assert biphasic["ESR"] > baseline["ESR"]

    def test_an_adapting_response_decays_only_the_general_regulon(self):
        """The ESR is the transient one; a stressor-specific regulon stays where it was."""
        doses = {"DTT": 1.0}
        baseline = module_activity(baseline_family(), doses)
        adapting = module_activity(build_family("adapting"), doses)

        assert adapting["ESR"] < baseline["ESR"]
        assert adapting["UPR"] == pytest.approx(baseline["UPR"])

    def test_an_adapting_response_lowers_every_reporter_that_reads_the_esr(self):
        columns = [i for i, name in enumerate(READERS)
                   if REPORTERS[name].module == "ESR" or "ESR" in REPORTERS[name].also_reads]
        baseline = clean(baseline_family(), doses=(1.0,))
        adapting = clean(build_family("adapting"), doses=(1.0,))

        assert columns
        assert np.all(adapting[:, columns] <= baseline[:, columns] + 1e-12)
        assert np.any(adapting[:, columns] < baseline[:, columns] - 1e-9)

    def test_feedback_closes_a_loop_the_sequential_cascade_cannot_represent(self):
        """A cycle has a non-zero spectral radius; an acyclic cascade is nilpotent."""
        baseline = np.max(np.abs(np.linalg.eigvals(baseline_family().cascade_matrix())))
        looped = np.max(np.abs(np.linalg.eigvals(build_family("feedback").cascade_matrix())))

        assert baseline == pytest.approx(0.0, abs=1e-12)
        assert looped > 1e-6

    def test_negative_feedback_damps_the_arm_it_returns_to(self):
        doses = {"heat": 4.0}
        baseline = module_activity(baseline_family(), doses)
        looped = module_activity(build_family("feedback"), doses)

        assert looped["heat"] < baseline["heat"]
        assert looped["proteasome"] < baseline["proteasome"]

    def test_loading_noise_moves_the_promoters_and_leaves_the_network_alone(self):
        """It is the control: the readings must move without any module activity moving."""
        doses = {"DTT": 1.0, "heat": 4.0}
        family = build_family("loading_noise", seed=3)
        baseline = module_activity(baseline_family(), doses)
        drawn = module_activity(family, doses)

        for module in MODULES:
            assert drawn[module] == pytest.approx(baseline[module])
        assert not np.allclose(family.loadings(READERS), baseline_family().loadings(READERS))

    def test_loading_noise_stays_inside_the_bracket_it_claims(self):
        """A two-fold bracket at 95% means a handful of draws outside it and none absurd."""
        baseline = baseline_family().also_reads
        folds = []
        for seed in range(40):
            drawn = build_family("loading_noise", seed=seed).also_reads
            for reporter, edges in baseline.items():
                for module, weight in edges.items():
                    ratio = drawn[reporter][module] / weight
                    folds.append(max(ratio, 1.0 / ratio))

        assert np.mean(np.array(folds) <= 2.0) > 0.9
        assert max(folds) < 8.0

    def test_loading_noise_leaves_a_reporter_own_element_at_one(self):
        """That entry is a normalisation of the channel's units, not a measurement."""
        loadings = build_family("loading_noise", seed=5).loadings(READERS)
        order = list(MODULES)

        for row, name in enumerate(READERS):
            assert loadings[row, order.index(REPORTERS[name].module)] == 1.0

    def test_a_redraw_is_smaller_than_a_rewiring(self):
        """Why the control belongs in the set at all.

        If redrawing unmeasured constants moved the readings as far as removing an edge, the
        structural families would be measuring parameter sensitivity and the cross-family
        comparison would say nothing about topology.
        """
        baseline = clean(baseline_family())
        redraw = np.mean([scalar_residual(clean(build_family("loading_noise", seed=s)),
                                          baseline) for s in range(5)])
        rewire = scalar_residual(clean(build_family("edge_dropped")), baseline)

        assert redraw < rewire


class TestTheRegistry:
    def test_it_declares_what_it_perturbs_and_why(self):
        for name in family_names():
            family = build_family(name)
            assert len(family.perturbs) > 20
            assert len(family.plausible_because) > 80, (
                f"{name} does not make a case for itself, so surviving it means nothing")

    def test_no_two_families_claim_the_same_perturbation(self):
        claims = [build_family(name).perturbs for name in family_names()]

        assert len(set(claims)) == len(claims)

    def test_the_baseline_comes_first_so_a_report_has_an_origin(self):
        assert family_names()[0] == "baseline"

    def test_a_family_is_built_by_name_rather_than_hardcoded(self):
        assert set(family_names()) == set(FAMILY_BUILDERS)
        for name in family_names():
            assert build_family(name).name == name

    def test_building_families_does_not_mutate_the_panel_they_derive_from(self):
        """The footgun this file exists next to: `edge_dropped` builds its cascade by
        popping an edge, and `loading_noise` by rescaling weights. Both operate on copies.
        If either ever stopped copying, it would corrupt `MODULES` or `REPORTERS` for the
        whole process -- every later caller in the same session, tests included, would read
        a panel nobody wrote, and the symptom would surface as an unrelated test failing
        only when this module was imported first.
        """
        before = (copy.deepcopy({n: m.driven_by for n, m in MODULES.items()}),
                  copy.deepcopy({n: r.also_reads for n, r in REPORTERS.items()}))

        for name in family_names():
            for seed in range(3):
                build_family(name, seed=seed)
        rng = np.random.default_rng(0)
        for _ in range(20):
            sample_family(rng)

        assert ({n: m.driven_by for n, m in MODULES.items()},
                {n: r.also_reads for n, r in REPORTERS.items()}) == before

    def test_sampling_reaches_every_family(self):
        """The distribution the bound is stated over has to actually cover the registry."""
        rng = np.random.default_rng(0)
        drawn = {sample_family(rng).name for _ in range(200)}

        assert drawn == set(family_names())

    def test_sampling_is_reproducible_from_its_seed(self):
        first = [sample_family(np.random.default_rng(7)) for _ in range(1)][0]
        again = sample_family(np.random.default_rng(7))

        assert (first.name, first.draw) == (again.name, again.draw)

    def test_a_stochastic_family_differs_between_draws(self):
        assert (build_family("loading_noise", seed=1).also_reads
                != build_family("loading_noise", seed=2).also_reads)

    def test_a_deterministic_family_ignores_the_seed(self):
        """So a caller can seed uniformly without knowing which families are stochastic."""
        assert (build_family("edge_dropped", seed=1).drivers
                == build_family("edge_dropped", seed=99).drivers)


class TestItRefusesRatherThanGuesses:
    def test_an_unknown_family_names_the_ones_there_are(self):
        with pytest.raises(KeyError, match="baseline"):
            build_family("no_such_family")

    def test_a_family_with_no_case_for_itself_is_refused(self):
        with pytest.raises(ValueError, match="what it perturbs"):
            PanelFamily(name="bare", perturbs="", plausible_because="")

    def test_an_edge_to_a_module_that_does_not_exist_is_refused(self):
        with pytest.raises(KeyError, match="do not exist"):
            PanelFamily(name="bad", perturbs="x" * 30, plausible_because="y" * 100,
                        drivers={"ESR": {"telepathy": 0.5}})

    def test_a_reporter_that_does_not_exist_is_refused(self):
        with pytest.raises(KeyError, match="do not exist"):
            PanelFamily(name="bad", perturbs="x" * 30, plausible_because="y" * 100,
                        also_reads={"NoSuchReporter": {"ESR": 0.5}})

    def test_an_adaptation_above_one_is_refused(self):
        """Above one is an induction, and belongs in a weight where it can be seen."""
        with pytest.raises(ValueError, match="adaptation outside"):
            PanelFamily(name="bad", perturbs="x" * 30, plausible_because="y" * 100,
                        adaptation={"ESR": 1.4})

    def test_a_non_positive_falling_limb_is_refused(self):
        with pytest.raises(ValueError, match="falling-limb"):
            PanelFamily(name="bad", perturbs="x" * 30, plausible_because="y" * 100,
                        lethal_scale={"ESR": 0.0})

    def test_a_self_sustaining_loop_is_refused(self):
        """At unit loop gain the dose stops deciding the answer and the network decides it."""
        with pytest.raises(ValueError, match="loop gain"):
            PanelFamily(name="bad", perturbs="x" * 30, plausible_because="y" * 100,
                        drivers={"heat": {"proteasome": 1.0}, "proteasome": {"heat": 1.0}})

    def test_an_unknown_stressor_is_refused_by_the_response(self):
        with pytest.raises(KeyError, match="no stressor"):
            module_activity(baseline_family(), {"absinthe": 1.0})

    def test_an_unknown_stressor_is_refused_by_the_dataset(self):
        with pytest.raises(KeyError, match="no such stressor"):
            family_dataset(baseline_family(), stressors=["absinthe"])

    def test_an_unknown_reporter_is_refused_by_the_dataset(self):
        with pytest.raises(KeyError, match="no reporter"):
            family_dataset(baseline_family(), reporters=["mCherry-vibes"])


class TestTheDataset:
    def test_it_returns_a_reading_per_well_and_channel(self):
        data = family_dataset(build_family("adapting"), reporters=READERS, stressors=PANEL,
                              doses=DOSES, replicates=2, noise_cv=0.05)

        assert data.readings.shape == (len(PANEL) * len(DOSES) * 2, len(READERS))
        assert set(data.labels) == set(PANEL)

    def test_co_dosed_treatments_appear_as_their_own_label(self):
        data = family_dataset(baseline_family(), reporters=READERS, stressors=PANEL,
                              doses=DOSES, replicates=1, noise_cv=0.0,
                              combinations=[("DTT", "H2O2")])

        assert "DTT+H2O2" in set(data.labels)

    def test_it_is_reproducible_from_its_seed_and_moves_with_it(self):
        """Both halves, in the form `scripts/audit_determinism.py` asks for them.

        The second is the one that catches a seed threaded halfway: a generator that accepts
        a seed and never draws from it reproduces perfectly and makes every replicate of a
        spread the same number.
        """
        settings = dict(reporters=READERS, stressors=["DTT"], doses=DOSES, replicates=2,
                        noise_cv=MEASURED_ACTIVITY_CV)
        first = family_dataset(baseline_family(), seed=11, **settings).readings

        assert np.array_equal(first, family_dataset(baseline_family(), seed=11,
                                                    **settings).readings)
        assert not np.allclose(first, family_dataset(baseline_family(), seed=12,
                                                     **settings).readings)

    def test_replicates_are_identical_without_noise_and_differ_with_it(self):
        settings = dict(reporters=READERS, stressors=["DTT"], doses=(1.0,), replicates=4)
        quiet = family_dataset(baseline_family(), noise_cv=0.0, **settings).readings
        noisy = family_dataset(baseline_family(), noise_cv=0.1, seed=0, **settings).readings

        assert np.allclose(quiet, quiet[0])
        assert not np.allclose(noisy, noisy[0])
