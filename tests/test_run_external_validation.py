"""The external-validation script: its condition parser, its labels, and its leak guard.

``scripts/run_external_validation.py`` is the only place in the repository where the latent
state meets a measurement it was not fitted to, and it writes the seven ``external_*``
tables ``docs/EXTERNAL_VALIDATION.md`` quotes. Almost everything load-bearing in it is a
pure function over free text or over JSON, and none of that needs the 96 arrays:

* ``assign_conditions`` decides which arrays are evidence, from GEO sample titles, by
  regex, in order. Loosening one pattern silently changes the denominator of every number
  in the file;
* ``load_regulons`` builds the labels. It has to refuse the measurement paper's own
  records, keep only records where the factor is the regulator, and strip every channel
  gene from every regulon -- the last is the leak guard that stops a channel predicting a
  module whose truth contains that channel;
* ``peer_count_measured`` is a thresholded count, and the threshold is quoted in prose.

The real Gasch data is touched only by the final class, which skips without it.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
import json
import importlib.util
import os
import pathlib
from queue import SimpleQueue
from threading import Barrier, get_ident
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "run_external_validation.py"
_FETCH_SCRIPT = _REPO / "scripts" / "fetch_external_validation_data.py"


@pytest.fixture(scope="module")
def script(tmp_path_factory):
    """The script as a module, with its output directory redirected away from outputs/."""
    redirect = tmp_path_factory.mktemp("external_outputs")
    previous = os.environ.get("YSTWIN_OUTPUTS")
    os.environ["YSTWIN_OUTPUTS"] = str(redirect)
    try:
        spec = importlib.util.spec_from_file_location("run_external_validation", _SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if previous is None:
            os.environ.pop("YSTWIN_OUTPUTS", None)
        else:
            os.environ["YSTWIN_OUTPUTS"] = previous
    return module


def _meta(*titles) -> pd.DataFrame:
    return pd.DataFrame({"title": list(titles)},
                        index=[f"GSM{i}" for i in range(len(titles))])


class TestTheConditionParserIsTheDenominator:
    """Every array the parser accepts becomes evidence and every array it drops leaves the
    study. The titles below are the real GSE18 strings, so a pattern that stopped matching
    one of them would shrink a block without failing anything."""

    @pytest.mark.parametrize("title,expected", [
        ("1 mM Menadione (50 min)redo", "menadione"),
        ("1.5 mM diamide (5 min)", "diamide"),
        ("2.5mM DTT 045 min", "DTT"),
        ("dtt 015 min", "DTT"),
        ("constant 0.32 mM H2O2 (30 min) redo", "H2O2"),
        ("1M sorbitol - 60 min", "sorbitol"),
        ("heat shock 15 min", "heat"),
        ("Heat Shock 005 minutes", "heat"),
        ("heat shock 33 to 37, 20 minutes", "heat"),
        ("Nitrogen Depletion 4 h", "rapamycin"),
        ("YPD 2 h (25C)", "glucose_starvation"),
        ("YPD_4_h_30C", "glucose_starvation"),
    ])
    def test_a_real_title_lands_on_the_stressor_the_rules_claim(self, script, title,
                                                                expected):
        assert script.assign_conditions(_meta(title)).stressor.iloc[0] == expected

    def test_an_adapted_sorbitol_steady_state_is_dropped_rather_than_dosed(self, script):
        """The sorbitol pattern is anchored at both ends for this reason. "1M sorbitol vs.
        YPD" is a culture grown in sorbitol, not one shocked with it, and a pattern
        loosened to ``^1M sorbitol`` would fold five adapted arrays into the shock block."""
        assigned = script.assign_conditions(_meta("1M sorbitol vs. YPD 30 minutes"))

        assert assigned.stressor.iloc[0] is None
        assert "steady state" in assigned.note.iloc[0]

    def test_a_genotype_perturbation_is_dropped_with_its_reason(self, script):
        """GSE18 also holds deletion strains. They are not stressor arrays and the reason
        has to survive into ``external_conditions.csv``, where a reader checks it."""
        assigned = script.assign_conditions(_meta("DBY7286 + 0.3 mM H2O2 (20 min)"))

        assert assigned.stressor.iloc[0] is None
        assert "genotype" in assigned.note.iloc[0]

    def test_a_title_matching_nothing_says_unmatched_rather_than_guessing(self, script):
        assigned = script.assign_conditions(_meta("something nobody has a rule for"))

        assert assigned.stressor.iloc[0] is None
        assert assigned.note.iloc[0] == "unmatched"

    def test_a_matched_array_carries_the_rules_own_note(self, script):
        """The note is what ``external_conditions.csv`` shows for an array that was used,
        so it has to be the rule's stated reading and not a drop reason."""
        assigned = script.assign_conditions(_meta("1M sorbitol - 5 min"))

        assert assigned.note.iloc[0] == "1 M sorbitol hyper-osmotic shock"

    def test_the_first_matching_rule_wins(self, script):
        """The rules are an ordered tuple and several patterns are loose enough to overlap.
        A title naming two agents takes the earlier rule, deterministically."""
        both = "1 mM menadione plus 1.5 mM diamide"
        order = [name for name, _, _ in script.CONDITION_RULES]

        got = script.assign_conditions(_meta(both)).stressor.iloc[0]

        assert got == min({"menadione", "diamide"}, key=order.index)

    def test_the_input_rows_are_neither_added_to_nor_dropped(self, script):
        """Assignment annotates; the exclusion happens later and is recorded per row. A
        parser that dropped here would lose the audit trail for the dropped arrays."""
        meta = _meta("heat shock 15 min", "nonsense", "dtt 015 min")

        assigned = script.assign_conditions(meta)

        assert list(assigned.index) == list(meta.index)

    def test_every_rule_names_a_real_stressor(self, script):
        """The assigned label is used to index ``STRESSORS`` when own-target modules are
        looked up, so a name that is not a stressor is a KeyError far from here."""
        from ystwin.generator.stress_panel import STRESSORS

        assert {name for name, _, _ in script.CONDITION_RULES} <= set(STRESSORS)

    def test_no_stressor_is_claimed_by_two_rules(self, script):
        """The note lookup takes the first rule with a matching name, so two rules for one
        stressor would attach the wrong reading to half its arrays."""
        names = [name for name, _, _ in script.CONDITION_RULES]

        assert len(set(names)) == len(names)


class TestTheLabelsAreIndependentOfTheMeasurements:
    """The regulons are the truth the transfer scores are computed against. Two things make
    them independent of Gasch 2000: a binding-only evidence source, and a real raise if the
    measurement paper turns out to contribute records after all."""

    FACTOR = "MSN2"

    @staticmethod
    def _record(script, factor, target, pmid, regulation_of="transcription",
                regulator=None):
        return {"regulation_of": regulation_of,
                "reference": {"pubmed_id": pmid},
                "locus1": {"display_name": factor if regulator is None else regulator},
                "locus2": {"format_name": target}}

    @pytest.fixture
    def sgd(self, script, tmp_path, monkeypatch):
        """An SGD directory holding one file per factor, all empty until a test fills one."""
        directory = tmp_path / "regulation"
        directory.mkdir()
        for factors in script.MODULE_FACTORS.values():
            for factor in factors:
                (directory / f"{factor}.json").write_text("[]")
        monkeypatch.setattr(script, "SGD", directory)
        return directory

    @staticmethod
    def _expression(*genes):
        return pd.DataFrame(index=list(genes))

    def test_a_record_from_the_measurement_paper_stops_the_run(self, script, sgd):
        """If the labels came partly from Gasch's own clustering the validation would be
        scoring the panel against the paper that parameterised it. A raise and not an
        assert, because ``python -O`` drops asserts and this check is load-bearing."""
        (sgd / f"{self.FACTOR}.json").write_text(json.dumps(
            [self._record(script, self.FACTOR, "YFL014W", script.GASCH_PMID)]))

        with pytest.raises(RuntimeError) as refusal:
            script.load_regulons(self._expression("YFL014W"), excluded=set())

        assert str(script.GASCH_PMID) in str(refusal.value)

    def test_only_the_binding_only_source_contributes_genes(self, script, sgd):
        """The wide set runs to 500+ promiscuous ChIP hits per factor, whose regulon mean is
        then the array mean. Restricting to MacIsaac's conserved motifs is what keeps the
        module means module-specific."""
        (sgd / f"{self.FACTOR}.json").write_text(json.dumps([
            self._record(script, self.FACTOR, "YFL014W", script.REGULON_SOURCE_PMID),
            self._record(script, self.FACTOR, "YGR088W", 999999),
        ]))

        regulons = script.load_regulons(self._expression("YFL014W", "YGR088W"), set())

        assert regulons["ESR"] == ["YFL014W"]

    def test_a_record_where_the_factor_is_the_target_is_not_its_own_regulon(self, script,
                                                                           sgd):
        """SGD's per-gene files hold records where the gene is somebody else's target.
        Taking ``locus2`` unconditionally put CAT8 into its own regulon and produced a
        one-gene regulon consisting entirely of CAT8."""
        (sgd / f"{self.FACTOR}.json").write_text(json.dumps([
            self._record(script, self.FACTOR, self.FACTOR, script.REGULON_SOURCE_PMID,
                         regulator="SOMEONE_ELSE"),
        ]))

        regulons = script.load_regulons(self._expression(self.FACTOR), set())

        assert regulons["ESR"] == []

    def test_regulation_of_something_other_than_transcription_is_ignored(self, script,
                                                                        sgd):
        (sgd / f"{self.FACTOR}.json").write_text(json.dumps([
            self._record(script, self.FACTOR, "YFL014W", script.REGULON_SOURCE_PMID,
                         regulation_of="protein activity"),
        ]))

        regulons = script.load_regulons(self._expression("YFL014W"), set())

        assert regulons["ESR"] == []

    def test_a_gene_with_no_probe_on_the_arrays_is_not_in_the_regulon(self, script, sgd):
        """The regulon mean is taken over the expression block, so a gene absent from it
        would be a KeyError rather than a missing row."""
        (sgd / f"{self.FACTOR}.json").write_text(json.dumps([
            self._record(script, self.FACTOR, "YFL014W", script.REGULON_SOURCE_PMID),
            self._record(script, self.FACTOR, "NOT_ON_THE_ARRAY",
                         script.REGULON_SOURCE_PMID),
        ]))

        regulons = script.load_regulons(self._expression("YFL014W"), set())

        assert regulons["ESR"] == ["YFL014W"]

    def test_a_channel_gene_is_stripped_from_every_regulon_and_not_only_its_own(
            self, script, sgd):
        """The leak guard. A channel that is inside a module's truth predicts that module
        by being part of it, and the module need not be the one the channel reads: KAR2 is
        the UPR channel, and leaving it in the ESR regulon would let it predict ESR."""
        (sgd / f"{self.FACTOR}.json").write_text(json.dumps([
            self._record(script, self.FACTOR, "YJL034W", script.REGULON_SOURCE_PMID),
            self._record(script, self.FACTOR, "YFL014W", script.REGULON_SOURCE_PMID),
        ]))

        regulons = script.load_regulons(self._expression("YJL034W", "YFL014W"),
                                        excluded={"YJL034W"})

        assert regulons["ESR"] == ["YFL014W"]

    def test_a_module_with_no_records_is_reported_empty_rather_than_dropped(self, script,
                                                                           sgd):
        """An absent module would make ``regulons.get(m, [])`` silently unmeasurable; an
        empty one is visible in ``external_module_map.csv`` with a zero gene count."""
        regulons = script.load_regulons(self._expression("YFL014W"), set())

        assert set(regulons) == set(script.MODULE_FACTORS)
        assert all(genes == [] for genes in regulons.values())

    def test_the_factor_table_is_the_one_the_bridge_defines(self, script):
        """It was duplicated here once, and a divergence between the two would be
        invisible: both would run and quietly disagree about which genes a module owns."""
        from ystwin.bridge import regulation

        assert script.MODULE_FACTORS is regulation.MODULE_FACTORS


class TestThePeerCountsAreCountsOfOtherStressors:
    """The redundancy law regresses recovery on peer count, so the count is half the
    statistic. Two definitions are computed, and neither may look at the held-out block."""

    def test_a_stressor_is_not_its_own_peer(self, script):
        """DTT and tunicamycin are the only two agents that drive UPR. Held out of a panel
        holding neither of the others, UPR has no peers -- and counting DTT itself would
        empty the zero-peer bin, which is the bin simulation says recovers nothing."""
        assert script.peer_count("UPR", "DTT", ["DTT", "H2O2", "heat"]) == 0

    def test_another_agent_driving_the_same_module_is_a_peer(self, script):
        assert script.peer_count("UPR", "DTT", ["DTT", "tunicamycin", "heat"]) == 1

    def test_own_targets_are_restricted_to_what_the_data_can_see(self, script):
        """A module with fewer than ``MIN_REGULON`` genes has no truth column worth
        scoring, and a pair scored against a column of zeros is not a transfer failure."""
        targets = script.own_targets("DTT", measurable=["UPR"])

        assert targets == ["UPR"]

    def test_own_targets_of_an_agent_are_the_panels_own_list(self, script):
        from ystwin.generator.stress_panel import MODULES, STRESSORS

        targets = script.own_targets("DTT", measurable=list(MODULES))

        assert targets == [m for m in STRESSORS["DTT"].targets if m in MODULES]

    @staticmethod
    def _dataset(script, values):
        """One array per stressor, with a chosen value in the ``UPR`` column."""
        from ystwin.generator.panel_experiment import PanelDataset
        from ystwin.generator.stress_panel import MODULES

        index = list(MODULES).index("UPR")
        modules = np.zeros((len(values), len(MODULES)))
        for row, value in enumerate(values.values()):
            modules[row, index] = value
        return PanelDataset(readings=np.zeros((len(values), 1)),
                            labels=np.array(list(values)), doses=np.zeros(len(values)),
                            modules=modules, reporters=["one"])

    def test_a_module_exactly_on_the_threshold_is_not_driven(self, script):
        """``DRIVEN_LOG2`` is quoted in prose as "about 19%", and the comparison is strict.
        A move of exactly the threshold is the boundary of the claim, not inside it."""
        dataset = self._dataset(script, {"DTT": 0.0, "H2O2": script.DRIVEN_LOG2})

        assert script.peer_count_measured(dataset, "UPR", "DTT", ["DTT", "H2O2"]) == 0

    def test_a_module_just_past_the_threshold_is_driven(self, script):
        dataset = self._dataset(script, {"DTT": 0.0, "H2O2": script.DRIVEN_LOG2 + 1e-6})

        assert script.peer_count_measured(dataset, "UPR", "DTT", ["DTT", "H2O2"]) == 1

    def test_a_repressed_module_counts_as_driven(self, script):
        """The count is over ``abs``: a regulon pushed down by an agent is as much evidence
        that the agent reaches the module as one pushed up."""
        dataset = self._dataset(script, {"DTT": 0.0, "H2O2": -1.0})

        assert script.peer_count_measured(dataset, "UPR", "DTT", ["DTT", "H2O2"]) == 1

    def test_the_held_out_block_is_never_inspected(self, script):
        """Every other stressor is in training when this one is held out, which is what
        makes the measured peer count legitimate. Counting the held-out agent's own block
        would put the test set into the predictor."""
        dataset = self._dataset(script, {"DTT": 5.0, "H2O2": 0.0})

        assert script.peer_count_measured(dataset, "UPR", "DTT", ["DTT", "H2O2"]) == 0

    def test_the_threshold_is_the_fold_change_the_comment_claims(self, script):
        """0.25 in log2 is a 19% move. The prose is what a reader checks the threshold
        against, so the arithmetic behind it is checked here."""
        assert 2 ** script.DRIVEN_LOG2 - 1 == pytest.approx(0.19, abs=0.005)


class TestAnUnscoreableSweepIsNotAZero:
    def test_a_table_with_no_pairs_scores_nan_rather_than_zero(self, script,
                                                              monkeypatch):
        """``--independent`` can leave a stressor with no measurable own-target module. A
        median of nothing is not a transfer score of nothing, and reporting zero would put
        a fabricated failure into the null comparison."""
        monkeypatch.setattr(script, "transfer_pairs",
                            lambda *args, **kwargs: pd.DataFrame())

        assert np.isnan(script.median_pair_score(None, [], [], 3))


class TestTheKnownLoadingsPathInvertsWhatItIsGiven:
    """``x = y (L+)'`` with no fitting anywhere. It is the one path whose failure would be
    arithmetic rather than statistical, so the arithmetic is pinned."""

    def test_an_invertible_loading_matrix_is_inverted_exactly(self, script):
        loadings = np.array([[2.0, 0.0], [0.0, 4.0]])
        modules = np.array([[1.0, 3.0], [-2.0, 0.5]])

        recovered = script.loadings_readout(modules @ loadings.T, loadings)

        assert recovered == pytest.approx(modules)

    def test_more_modules_than_channels_leaves_modules_unreachable(self, script):
        """17 channels against 24 modules, so seven directions are unreachable by
        construction -- the Ledermann-adjacent point, met on real readings."""
        from ystwin.generator.stress_panel import MODULES, transcriptional_reporters
        from ystwin.generator.stress_panel import reporter_loadings

        loadings = reporter_loadings(transcriptional_reporters())

        assert np.linalg.matrix_rank(loadings) < len(MODULES)

    def test_agreement_is_nan_when_no_module_is_measurable(self, script):
        """Rather than a median over an empty list, which is a warning and a NaN by
        accident instead of a NaN on purpose."""
        readings = np.ones((4, 2))
        truth = np.ones((4, 2))

        assert np.isnan(script.loadings_agreement(readings, truth, np.eye(2), []))


class TestTheChannelsAndTheDroppedOnes:
    """``CHANNEL_GENES`` is the map from a panel promoter to the marker gene the arrays
    carry, and ``SPARSE_CHANNELS`` is what could not be read. Both are literals."""

    def test_every_channel_is_a_reporter_the_panel_defines(self, script):
        from ystwin.generator.stress_panel import REPORTERS

        assert set(script.CHANNEL_GENES) <= set(REPORTERS)

    def test_every_dropped_channel_is_one_of_the_channels(self, script):
        """A name that is not a channel would drop nothing and the sparse probe would
        still poison its column."""
        assert set(script.SPARSE_CHANNELS) <= set(script.CHANNEL_GENES)

    def test_the_circularity_restricted_stressors_are_ones_the_parser_can_assign(self,
                                                                                 script):
        """``--independent`` drops them by name from the assigned labels. A name no rule
        produces would drop nothing while claiming to."""
        assigned = {name for name, _, _ in script.CONDITION_RULES}

        assert set(script.GASCH_DERIVED_STRESSORS) <= assigned

    def test_the_platforms_read_here_are_the_ones_the_fetcher_downloads(self, script):
        """Two scripts, two literal tuples, one dataset. GPL63 is listed by the fetcher and
        deliberately not read here, because GEO publishes no annotation table for it."""
        spec = importlib.util.spec_from_file_location("fetch_external", _FETCH_SCRIPT)
        fetch = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fetch)

        read = set(script.PLATFORMS) | set(script.PLATFORMS_EXTRA)
        assert read == set(fetch.PLATFORMS) - set(fetch.NO_ANNOTATION)


@pytest.mark.integration
class TestTheDatasetTheRealArraysProduce:
    """The 81 arrays the published tables are computed on. Magnitudes are not pinned; the
    shape of the evidence and the leak guard are."""

    @pytest.fixture(scope="class")
    def built(self, script):
        if not (script.GASCH / "GPL51.annot.gz").exists():
            pytest.skip("Gasch GSE18 not fetched; run "
                        "scripts/fetch_external_validation_data.py")
        if not (script.SGD / "MSN2.json").exists():
            pytest.skip("SGD regulation records not fetched; run "
                        "scripts/fetch_external_validation_data.py")
        return script.build(independent=False)

    def test_the_evidence_is_the_eighty_one_arrays_the_docstring_claims(self, built):
        dataset, kept, *_ = built

        assert len(kept) == 81
        assert len(set(dataset.labels)) == 8

    def test_no_reading_is_missing(self, built):
        """``train_stress_model`` centres with a plain mean, so one missing entry poisons a
        whole channel. Arrays with a gap are excluded by name rather than imputed."""
        dataset = built[0]

        assert np.isfinite(dataset.readings).all()

    def test_no_channel_gene_is_inside_any_module_regulon(self, built, script):
        """The leak guard, on the real records. A channel inside a module's truth predicts
        that module by being part of it. Every channel is stripped from every regulon --
        including the two dropped for sparse probes, whose genes would otherwise return to
        the truth after leaving the readings."""
        _, _, _, _, _, regulons, _ = built
        symbols = {gene for genes in script.CHANNEL_GENES.values() for gene in genes}
        # The regulons are ORFs and CHANNEL_GENES are symbols, so compare through the same
        # map the script builds ``excluded`` with.
        _, _, symbol_to_orf = script.load_expression()
        channels = {symbol_to_orf[s] for s in symbols if s in symbol_to_orf}
        sparse = {symbol_to_orf[gene] for name in script.SPARSE_CHANNELS
                  for gene in script.CHANNEL_GENES[name] if gene in symbol_to_orf}

        assert sparse <= channels
        for module, genes in regulons.items():
            assert not channels & set(genes), f"{module} contains a channel gene"

    def test_no_dose_is_fitted_anywhere(self, built):
        """Gasch's time courses are a time axis, not a concentration ladder. Doses are
        supplied as zeros so the panel's dose-response machinery is not exercised, and a
        nonzero one here would mean it silently was."""
        dataset = built[0]

        assert not dataset.doses.any()

    def test_every_excluded_array_says_why(self, built):
        """``external_conditions.csv`` is the audit trail for the 15 arrays that did not
        make it. An exclusion with an empty reason is indistinguishable from an inclusion."""
        _, _, meta, *_ = built
        excluded = meta[meta["excluded"] != ""]

        assert len(excluded) == len(meta) - 81
        assert excluded["excluded"].str.len().gt(0).all()

    def test_the_zero_dose_arrays_are_inside_their_blocks(self, built):
        """A documented decision, not an oversight: they are the bottom rung of their own
        ladders, and keeping them pulls the DTT and heat block means down. Recorded in
        docs/EXTERNAL_VALIDATION.md sections 5 and 8. Pinned because a regex tightened to
        exclude them would move every number in the file with nothing to point at."""
        _, kept, *_ = built
        zero_dose = kept[kept.title.str.contains("000", case=False)]

        assert set(zero_dose.stressor) == {"DTT", "heat"}

    def test_dropping_the_gasch_derived_agents_drops_their_arrays_entirely(self, script,
                                                                          built):
        """``--independent`` is the circularity control: the panel it leaves carries nothing
        traceable to the dataset except the ESR arm."""
        restricted = script.build(independent=True)
        labels = set(restricted[0].labels)

        assert labels == set(built[0].labels) - set(script.GASCH_DERIVED_STRESSORS)
        assert len(restricted[1]) < len(built[1])


def test_invalid_null_budget_is_refused_before_loading_or_writing(script, tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(script, "build", lambda *args: calls.append(1))
    output = tmp_path / "not_created"
    with pytest.raises(ValueError, match="at least"):
        script.main(False, script._MIN_DRAWS - 1, None, output)
    assert not calls
    assert not output.exists()


def test_transfer_pairs_uses_nested_training_only_widths_by_default(script):
    from ystwin.generator.panel_experiment import panel_dataset

    data = panel_dataset(reporters=["STRE-general", "UPRE-ER", "TRX2-oxidative"],
                         stressors=["DTT", "H2O2", "heat"],
                         doses=(0.1, 0.5, 1.0, 2.0), replicates=3, seed=4)
    table = script.transfer_pairs(data, ["DTT", "H2O2", "heat"], ["UPR", "ESR", "oxidative", "heat"], None)
    assert set(table.held_out) == {"DTT", "H2O2", "heat"}
    assert table.selection.str.contains("outer training only").all()
    assert table.n_states.between(1, 3).all()
    assert (table.n_train == 24).all()
    assert (table.n_test == 12).all()


def test_constant_loadings_estimates_are_unscorable_without_a_runtime_warning(script):
    value = script.loadings_agreement(np.ones((4, 2)), np.ones((4, 2)), np.eye(2), ["ESR", "UPR"])
    assert np.isnan(value)


@pytest.fixture(scope="module")
def synthetic_folds(script):
    rng = np.random.default_rng(1901)
    labels = np.repeat(["heat", "DTT", "H2O2", "NaCl", "DTT+H2O2"], 4)
    readings = rng.normal(size=(len(labels), 3))
    modules = (readings @ rng.normal(size=(3, len(script.MODULES)))
               + rng.normal(scale=0.1, size=(len(labels), len(script.MODULES))))
    modules[:, list(script.MODULES).index("UPR")] = 0.0
    data = script.PanelDataset(readings, labels, np.zeros(len(labels)), modules,
                               ["STRE-general", "UPRE-ER", "TRX2-oxidative"])
    for values in (data.readings, data.labels, data.doses, data.modules):
        values.setflags(write=False)
    return data


@pytest.mark.parametrize("width", [None, 2])
@pytest.mark.parametrize("seed", [0, 19])
def test_parallel_pairs_preserve_every_value_and_key(script, synthetic_folds, width, seed):
    panel = ["heat", "DTT", "H2O2", "NaCl", "heat"]
    measurable = list(script.MODULES)
    serial = script.transfer_pairs(synthetic_folds, panel, measurable, width, seed)
    parallel = script.transfer_pairs(synthetic_folds, panel, measurable, width, seed, workers=3)

    pd.testing.assert_frame_equal(serial, parallel, check_exact=True)
    assert serial.to_csv(index=False).encode() == parallel.to_csv(index=False).encode()
    assert list(parallel.columns) == [
        "held_out", "module", "score", "n_states", "n_train", "n_test", "selection",
        "baseline", "score_method", "peers", "peers_measured", "truth_sd", "truth_span",
    ]
    assert list(parallel[["held_out", "module"]].itertuples(index=False, name=None)) == [
        (held, module) for held in panel for module in script.own_targets(held, measurable)
    ]
    assert parallel.loc[parallel.module == "UPR", "score"].isna().all()
    assert parallel.score.notna().any()
    assert parallel.loc[parallel.held_out.isin(["DTT", "H2O2"]), "n_train"].eq(12).all()
    assert parallel.loc[parallel.held_out.isin(["heat", "NaCl"]), "n_train"].eq(16).all()
    assert parallel.n_test.eq(4).all()
    expected_selection = ("nested leave-stressor-out on outer training only" if width is None
                          else "fixed width supplied before outer scoring")
    assert parallel.selection.eq(expected_selection).all()
    assert script.median_pair_score(synthetic_folds, panel, measurable, width, seed,
                                    workers=3) == float(np.nanmedian(serial.score))


def test_parallel_pairs_preserve_full_validation_results_and_bound_threads(
        script, synthetic_folds, monkeypatch):
    panel = ["heat", "DTT", "H2O2"]
    validate = script.validate_modules
    expected = {held: validate(synthetic_folds, held, seed=19) for held in panel}
    completed = SimpleQueue()
    barrier = Barrier(len(panel), timeout=10)
    bounds = []

    def executor(*, max_workers):
        bounds.append(max_workers)
        return ThreadPoolExecutor(max_workers=max_workers)

    def inspect(data, held, **kwargs):
        barrier.wait()
        result = validate(data, held, **kwargs)
        completed.put((held, os.getpid(), get_ident(), result))
        return result

    monkeypatch.setattr(script, "ThreadPoolExecutor", executor)
    monkeypatch.setattr(script, "validate_modules", inspect)
    script.transfer_pairs(synthetic_folds, panel, list(script.MODULES), None, 19, workers=99)
    results = [completed.get_nowait() for _ in panel]
    assert bounds == [len(panel)]
    assert {pid for _, pid, _, _ in results} == {os.getpid()}
    assert len({thread for _, _, thread, _ in results}) == len(panel)
    assert get_ident() not in {thread for _, _, thread, _ in results}
    for held, _, _, result in results:
        reference = expected[held]
        assert result.n_states == reference.n_states
        assert result.selection == reference.selection
        assert result.inner_scores == reference.inner_scores
        assert list(result.scores) == list(reference.scores)
        np.testing.assert_array_equal(list(result.scores.values()), list(reference.scores.values()))
        np.testing.assert_array_equal(result.prediction, reference.prediction)
        np.testing.assert_array_equal(result.baseline, reference.baseline)
        for field in ("train", "test", "excluded"):
            np.testing.assert_array_equal(getattr(result.fold, field), getattr(reference.fold, field))


@pytest.mark.parametrize("failure,width,message", [
    ("nonfinite_targets", 2, "training targets must be finite"),
    ("constant_targets", None, "no width has a defined score"),
    ("missing_reading", None, "no width has a defined score"),
    ("invalid_width", 0, "width must be a positive integer"),
])
def test_parallel_pairs_preserve_refusals(script, synthetic_folds, failure, width, message):
    data = synthetic_folds
    if failure == "nonfinite_targets":
        modules = data.modules.copy()
        modules[:, 0] = np.nan
        data = replace(data, modules=modules)
    elif failure == "constant_targets":
        data = replace(data, modules=np.zeros_like(data.modules))
    elif failure == "missing_reading":
        readings = data.readings.copy()
        readings[0, 0] = np.nan
        data = replace(data, readings=readings)
    messages = []
    for workers in (1, 3):
        with pytest.raises(ValueError, match=message) as refused:
            script.transfer_pairs(data, ["heat", "DTT", "H2O2"], list(script.MODULES),
                                  width, workers=workers)
        messages.append(str(refused.value))
    assert messages[0] == messages[1]


def test_parallel_pairs_preserve_missing_reading_scores(script, synthetic_folds):
    readings = synthetic_folds.readings.copy()
    readings[0, 0] = np.nan
    data = replace(synthetic_folds, readings=readings)
    panel = ["heat", "DTT", "H2O2"]
    serial = script.transfer_pairs(data, panel, list(script.MODULES), 2)
    parallel = script.transfer_pairs(data, panel, list(script.MODULES), 2, workers=3)
    pd.testing.assert_frame_equal(serial, parallel, check_exact=True)
    assert serial.to_csv(index=False).encode() == parallel.to_csv(index=False).encode()
    assert parallel.loc[parallel.held_out == "heat", "score"].isna().all()


@pytest.mark.parametrize("workers", [1, 3])
def test_parallel_pairs_raise_the_first_requested_fold_error(script, synthetic_folds, workers):
    with pytest.raises(KeyError, match="absent"):
        script.transfer_pairs(synthetic_folds, ["absent", "DTT"], ["UPR"], 0, workers=workers)


@pytest.mark.parametrize("panel,measurable", [([], ["UPR"]), (["heat", "DTT"], [])])
def test_parallel_pairs_preserve_empty_tables(script, synthetic_folds, panel, measurable):
    serial = script.transfer_pairs(synthetic_folds, panel, measurable, 2)
    parallel = script.transfer_pairs(synthetic_folds, panel, measurable, 2, workers=3)
    pd.testing.assert_frame_equal(serial, parallel, check_exact=True)
    assert parallel.empty
    assert np.isnan(script.median_pair_score(synthetic_folds, panel, measurable, 2, workers=3))


@pytest.mark.parametrize("workers", [0, -1, 1.5, True, np.bool_(True), "3"])
def test_invalid_workers_are_refused_before_loading_or_writing(script, tmp_path, monkeypatch, workers):
    calls = []
    monkeypatch.setattr(script, "build", lambda *args: calls.append(1))
    monkeypatch.setattr(script, "validate_modules", lambda *args, **kwargs: calls.append(1))
    output = tmp_path / "not_created"
    with pytest.raises(ValueError, match="workers must be a positive integer"):
        script.main(False, 200, None, output, workers=workers)
    with pytest.raises(ValueError, match="workers must be a positive integer"):
        script.transfer_pairs(None, ["DTT"], ["UPR"], None, workers=workers)
    assert not calls
    assert not output.exists()


def test_default_pairs_do_not_start_a_thread_pool(script, synthetic_folds, monkeypatch):
    def forbidden(**kwargs):
        pytest.fail("the default execution path must remain serial")

    monkeypatch.setattr(script, "ThreadPoolExecutor", forbidden)
    script.transfer_pairs(synthetic_folds, ["heat", "DTT"], ["UPR", "ESR"], 2)


def test_main_routes_workers_only_to_fold_statistics(script, synthetic_folds, tmp_path, monkeypatch):
    data = synthetic_folds
    kept = pd.DataFrame({"stressor": data.labels})
    meta = kept.assign(excluded="")
    readings = pd.DataFrame(data.readings, columns=data.reporters)
    truth = pd.DataFrame(data.modules, columns=list(script.MODULES))
    measurable = list(script.MODULES)
    panel = ["heat", "DTT", "H2O2", "NaCl"]
    table = script.transfer_pairs(data, panel, measurable, 2)
    monkeypatch.setattr(script, "build", lambda independent: (
        data, kept, meta, readings, truth, {}, measurable))
    monkeypatch.setattr(script, "select_dimension", lambda centred, max_states: 2)
    calls, budgets = [], []

    def pairs(dataset, panel, measurable, width, seed=0, *, workers):
        calls.append((width, seed, workers))
        return table.copy()

    def compare(statistic, data, surrogate, *, n_draws, greater_is_better, label, seed):
        budgets.append((n_draws, seed))
        observed = statistic(data)
        return SimpleNamespace(observed=observed, null_median=0.0, p_value=1.0,
                               n_draws=n_draws, beats_null=False, summary=lambda: "test null")

    monkeypatch.setattr(script, "transfer_pairs", pairs)
    monkeypatch.setattr(script, "compare_to_null", compare)
    script.main(False, 200, None, tmp_path, workers=3)
    assert calls == [(width, 0, 3) for width in [1, 2, 3, None, None, None, None]]
    assert budgets == [(200, 0)] * 6
    assert len(list(tmp_path.glob("external_*.csv"))) == 7
