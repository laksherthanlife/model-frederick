"""The declared table and what it does to yeast-GEM, checked by running it.

Every number here was measured against yeast-GEM v9.0.2 under the model's SHIPPED DEFAULT
bounds, where growth is 0.085844 /h, except the one test that names the reference batch regime.
The metal-transport tests are the ones that matter: they are the predictions the repairs move.
"""

from __future__ import annotations

import pytest

from ystwin import paths
from ystwin.fba import physiology, solver
from ystwin.fba.stress_genes import (
    DECLARATION_PATH,
    NEW_REACTION_IDS,
    SGD_GAF_PATH,
    VERDICTS,
    additions,
    gene_label_map,
    load_declarations,
    new_reaction_targets,
    refusals,
    restore_stress_genes,
    sgd_annotations,
    unsupported_claims,
    verify_absence,
)

SHIPPED_DEFAULT_GROWTH = 0.085844
REFERENCE_BATCH_GROWTH = 0.345984


@pytest.fixture(scope="module")
def declarations():
    return load_declarations()


@pytest.fixture(scope="module")
def _base_model_template(yeast_gem_factory):
    """A private copy, so pinning the solver cannot leak into the session-scoped fixture."""
    model = yeast_gem_factory()
    solver.configure(model)
    return model


@pytest.fixture
def base_model(_base_model_template, model_copy):
    return model_copy(_base_model_template)


@pytest.fixture(scope="module")
def _fixed_model_template(_base_model_template, model_copy):
    model = restore_stress_genes(model_copy(_base_model_template))
    solver.configure(model)
    return model


@pytest.fixture
def fixed_model(_fixed_model_template, model_copy):
    return model_copy(_fixed_model_template)


def test_table_parses_and_every_verdict_is_known(declarations):
    assert len(declarations) == 41
    assert {r.verdict for r in declarations} <= set(VERDICTS)
    assert len({r.orf for r in declarations}) == len(declarations)


def test_every_claim_carries_a_go_term_and_an_sgd_evidence_line(declarations):
    for record in declarations:
        if record.verdict == "ALREADY_PRESENT":
            continue
        assert record.go_term.startswith("GO:"), record.orf
        assert ":" in record.evidence and record.evidence != "-", record.orf
        assert len(record.basis) > 40, record.orf


def test_the_refusal_list_is_the_larger_half(declarations):
    refused = refusals(declarations)
    assert len(refused) == 20
    assert len(additions(declarations)) == 19
    # The genes the audit singled out as correctly absent must all be refused.
    for orf in ("YER044C", "YMR038C", "YGR105W", "YHR060W", "YKL119C"):
        assert orf in {r.orf for r in refused}


def test_the_reaction_id_constant_is_the_table_and_not_a_second_opinion(declarations):
    assert new_reaction_targets(declarations) == NEW_REACTION_IDS


def test_new_reaction_rows_repeat_their_definition_and_are_checked(tmp_path, declarations):
    lines = DECLARATION_PATH.read_text().splitlines()
    broken = [
        line.replace("0,1000", "0,999", 1) if line.startswith("YML115C") else line
        for line in lines
    ]
    corrupted = tmp_path / "corrupted.tsv"
    corrupted.write_text("\n".join(broken) + "\n")
    with pytest.raises(ValueError, match="disagree"):
        load_declarations(corrupted)


def test_stoichiometry_and_bounds_parse(declarations):
    for record in additions(declarations):
        if record.verdict != "NEW_REACTION":
            continue
        assert record.parsed_bounds == (0.0, 1000.0)
        assert len(record.parsed_stoichiometry) == 4


def test_every_go_term_and_evidence_line_is_a_real_sgd_annotation(declarations):
    """The gate that caught the one invented citation: SMF2's manganese term is IBA, not IDA."""
    annotations = sgd_annotations()
    if annotations is None:
        pytest.skip(f"SGD GAF not vendored here ({SGD_GAF_PATH})")
    assert unsupported_claims(declarations, annotations) == ()


@pytest.mark.integration
def test_absence_was_verified_by_orf_not_by_symbol():
    path = paths.yeast_gem()
    if path is None:
        pytest.skip("yeast-GEM not present; set YSTWIN_YEAST_GEM")
    labels = gene_label_map(path)
    assert len(labels) == 1161
    # The alias trap: "CTR1" as a bare symbol is SGD's synonym for HNM1, which IS in the model.
    assert labels["YGL077C"] == "HNM1"
    assert "YPR124W" not in labels
    assert labels["YDL103C"] == "QRI1" and labels["YPL053C"] == "KTR6"


@pytest.mark.integration
def test_no_declared_absence_is_actually_present(base_model, declarations):
    assert verify_absence(base_model, declarations) == ()
    present = [r.orf for r in declarations if r.verdict == "ALREADY_PRESENT"]
    assert present and all(base_model.genes.has_id(orf) for orf in present)


@pytest.mark.integration
def test_every_added_gene_reaches_the_model_it_was_missing_from(fixed_model, declarations):
    """The false negative this exists to fix: each added ORF is now a gene carrying reactions."""
    for record in additions(declarations):
        assert fixed_model.genes.has_id(record.orf), record.orf
        assert fixed_model.genes.get_by_id(record.orf).reactions, record.orf


@pytest.mark.integration
def test_new_reactions_are_mass_and_charge_balanced(fixed_model):
    for rid in NEW_REACTION_IDS:
        reaction = fixed_model.reactions.get_by_id(rid)
        assert reaction.check_mass_balance() == {}, rid
        assert sum(m.charge * c for m, c in reaction.metabolites.items()) == 0, rid


@pytest.mark.integration
def test_growth_is_unchanged_in_both_regimes(base_model, fixed_model):
    """Loose by 1e-6 on purpose: the repo pins GLPK at tolerance 1e-7, so anything
    tighter would be testing the simplex rather than the additions."""
    before, after = solver.growth_or_none(base_model), solver.growth_or_none(fixed_model)
    assert before == pytest.approx(SHIPPED_DEFAULT_GROWTH, abs=1e-6)
    assert after == pytest.approx(before, abs=1e-6)
    for model in (base_model, fixed_model):
        with model:
            physiology.aerobic_batch_constraints(model, physiology.REFERENCE_AEROBIC_BATCH)
            assert solver.growth_or_none(model) == pytest.approx(
                REFERENCE_BATCH_GROWTH, abs=1e-6
            )


@pytest.mark.integration
def test_manganese_and_copper_lethality_calls_are_corrected(base_model, fixed_model):
    """SMF1 and ALR1 were the only genes on the copper, calcium and manganese transports."""
    for orf in ("YOL122C", "YOL130W"):
        with base_model:
            base_model.genes.get_by_id(orf).knock_out()
            assert solver.growth_or_none(base_model) == pytest.approx(0.0, abs=1e-9), orf
        with fixed_model:
            fixed_model.genes.get_by_id(orf).knock_out()
            assert solver.growth_or_none(fixed_model) == pytest.approx(
                SHIPPED_DEFAULT_GROWTH, abs=1e-6
            ), orf


@pytest.mark.integration
def test_the_dpm1_side_effect_is_declared_rather_than_silent(base_model, fixed_model, declarations):
    """MANPOL_G makes DPM1 dispensable here, which is wrong biology and is why it is written down."""
    with base_model:
        base_model.genes.get_by_id("YPR183W").knock_out()
        assert solver.growth_or_none(base_model) == pytest.approx(0.0, abs=1e-9)
    with fixed_model:
        fixed_model.genes.get_by_id("YPR183W").knock_out()
        assert solver.growth_or_none(fixed_model) == pytest.approx(
            SHIPPED_DEFAULT_GROWTH, abs=1e-6
        )
    basis = next(r.basis for r in declarations if r.orf == "YPL050C")
    assert "DPM1" in basis and "YPR183W" in basis


@pytest.mark.integration
def test_a_regulon_naming_tsa1_now_reaches_r_0550(base_model, fixed_model):
    assert not base_model.genes.has_id("YML028W")
    assert [r.id for r in fixed_model.genes.get_by_id("YML028W").reactions] == ["r_0550"]
    for orf, expected in (("YKR050W", "r_1249"), ("YER145C", "r_1178"),
                          ("YPR124W", "r_4589"), ("YLR411W", "r_4589"),
                          ("YHR050W", "r_4590"), ("YNL291C", "r_4587")):
        assert [r.id for r in fixed_model.genes.get_by_id(orf).reactions] == [expected]


@pytest.mark.integration
def test_tsa2_knockout_no_longer_silences_the_peroxiredoxin(base_model, fixed_model):
    from cobra.flux_analysis import flux_variability_analysis

    ranges = {}
    for label, model in (("before", base_model), ("after", fixed_model)):
        with model:
            model.genes.get_by_id("YDR453C").knock_out()
            fva = flux_variability_analysis(
                model, reaction_list=["r_0550"], fraction_of_optimum=0.0,
                processes=solver.FVA_PROCESSES,
            )
            ranges[label] = float(fva.maximum.iloc[0])
    assert ranges["before"] == pytest.approx(0.0, abs=1e-9)
    assert ranges["after"] > 1.0


@pytest.mark.integration
def test_manpol_carries_flux_and_alg3_is_declared_blocked(fixed_model, declarations):
    """ALG3_ER is exact chemistry into a branch that dead-ends at s_3887; the table says so."""
    from cobra.flux_analysis import flux_variability_analysis

    fva = flux_variability_analysis(
        fixed_model, reaction_list=["MANPOL_G", "ALG3_ER"], fraction_of_optimum=0.0,
        processes=solver.FVA_PROCESSES,
    )
    assert fva.maximum["MANPOL_G"] > 1e-3
    assert fva.maximum["ALG3_ER"] == pytest.approx(0.0, abs=1e-9)
    basis = next(r.basis for r in declarations if r.orf == "YBL082C")
    assert "s_3887" in basis and "zero flux" in basis


@pytest.mark.integration
def test_the_only_flux_ranges_the_additions_move_are_the_new_route_s_own_transports(
    base_model, fixed_model
):
    """A whole-model FVA moves four ranges and only these two by more than 1e-7.

    The full sweep costs five minutes, so the suite pins the two that carry the change: GDP and
    GDP-mannose transport could only run backwards before, and now reach the new route's own
    maximum. Everything else, including all 454 reactions already at |1000|, is untouched."""
    from cobra.flux_analysis import flux_variability_analysis

    ranges = {}
    for label, model in (("before", base_model), ("after", fixed_model)):
        ranges[label] = flux_variability_analysis(
            model, reaction_list=["r_1801", "r_1803"], fraction_of_optimum=0.0,
            processes=solver.FVA_PROCESSES,
        )
    for rid in ("r_1801", "r_1803"):
        assert ranges["before"].maximum[rid] == pytest.approx(0.0, abs=1e-9), rid
        assert ranges["after"].maximum[rid] == pytest.approx(0.060266, abs=1e-6), rid
        assert ranges["after"].minimum[rid] == pytest.approx(
            ranges["before"].minimum[rid], abs=1e-6), rid


@pytest.mark.integration
def test_the_refusals_cost_four_regulon_genes_and_the_docstring_says_which(base_model, regulons):
    """A correct refusal still leaves a module short, so the four it costs are pinned here."""
    from ystwin.fba import stress_genes as module

    refused = {r.orf for r in refusals(load_declarations())}
    stranded = {(mod, orf) for mod, regulon in regulons.items()
                for orf in regulon.genes if orf in refused}
    assert stranded == {("UPR", "YOL030W"), ("cell_wall", "YGR189C"),
                        ("hypoxia", "YER044C"), ("hypoxia", "YNL111C"),
                        ("oxidative", "YMR038C")}
    assert base_model.genes.has_id("YJR104C")
    assert [r.id for r in base_model.genes.get_by_id("YJR104C").reactions] == ["r_4270"]
    assert "GAS5 in UPR" in module.__doc__ and "costs it nothing" in module.__doc__


@pytest.mark.integration
def test_the_tok1_refusal_is_the_measurement_it_claims(base_model, fixed_model, declarations):
    """The refusal that was argued and is now measured: adding TOK1's K+ efflux moves nothing.

    The row used to say the reaction would be a free proton pump. It is not: r_2020 and r_1832
    are both (-1000, 1000), so r_1249 can already carry protons inward against nothing."""
    import cobra

    for model in (base_model, fixed_model):
        with model:
            efflux = cobra.Reaction("TOK1_K", lower_bound=0.0, upper_bound=1000.0)
            model.add_reactions([efflux])
            efflux.add_metabolites({model.metabolites.s_1373: -1, model.metabolites.s_1374: 1})
            assert solver.growth_or_none(model) == pytest.approx(
                SHIPPED_DEFAULT_GROWTH, abs=1e-6)
            with model:
                for reaction in model.exchanges:
                    reaction.bounds = (-1000.0, 1000.0)
                assert solver.growth_or_none(model) == pytest.approx(77.556090, abs=1e-6)
            with model:
                model.reactions.get_by_id("r_1714").bounds = (-1.0, -1.0)
                model.reactions.get_by_id("r_1992").lower_bound = -1000.0
                model.reactions.get_by_id("r_4046").bounds = (0.0, 1000.0)
                model.objective = model.reactions.get_by_id("r_4046")
                assert solver.growth_or_none(model) == pytest.approx(19.192, abs=1e-3)
    for reaction in ("r_2020", "r_1832"):
        assert base_model.reactions.get_by_id(reaction).bounds == (-1000.0, 1000.0)
    basis = next(r.basis for r in declarations if r.orf == "YJL093C")
    assert "r_2020" in basis and "moves no number" in basis


@pytest.mark.integration
def test_exactly_three_deletion_calls_move_and_they_are_the_declared_three(
    base_model, fixed_model
):
    """The ledger, over every gene a new reaction or a rewritten rule could reach.

    A full 1161/1180-gene sweep was run once and agrees -- three moved, lethal count 188 to 185
    -- but it costs three minutes, so the suite pins the neighbourhood the change can act in."""
    candidates = set()
    for rid in NEW_REACTION_IDS:
        for metabolite in fixed_model.reactions.get_by_id(rid).metabolites:
            for reaction in metabolite.reactions:
                candidates.update(gene.id for gene in reaction.genes)
    for record in additions(load_declarations()):
        if record.verdict == "GPR_EXTEND":
            candidates.update(g.id for g in base_model.reactions.get_by_id(record.target).genes)
            candidates.update(g.id for g in fixed_model.reactions.get_by_id(record.target).genes)
    moved = set()
    for orf in sorted(candidates):
        if not base_model.genes.has_id(orf):
            continue
        with base_model:
            base_model.genes.get_by_id(orf).knock_out()
            before = solver.growth_or_none(base_model)
        with fixed_model:
            fixed_model.genes.get_by_id(orf).knock_out()
            after = solver.growth_or_none(fixed_model)
        if abs((before or 0.0) - (after or 0.0)) > 1e-6:
            moved.add(orf)
    assert len(candidates) > 50
    assert moved == {"YOL122C", "YOL130W", "YPR183W"}


@pytest.mark.integration
def test_which_added_genes_are_load_bearing_is_measured_not_asserted(fixed_model, regulons):
    """Criterion (a) in one test: twelve of the nineteen move no single-gene observable.

    SMF2, MID1 and CCH1 look inert alone and are not -- they carry the rescues, which is why
    the double knockouts are here beside the single ones."""
    named = {orf for regulon in regulons.values() for orf in regulon.genes}
    added = [r.orf for r in additions(load_declarations())]
    assert {orf for orf in added if orf in named} == {
        "YML028W", "YPR124W", "YER145C", "YER001W"}
    inert = []
    for orf in added:
        with fixed_model:
            fixed_model.genes.get_by_id(orf).knock_out()
            if solver.growth_or_none(fixed_model) == pytest.approx(
                    SHIPPED_DEFAULT_GROWTH, abs=1e-6) and orf not in named:
                inert.append(orf)
    assert len(inert) == 15
    for partners in (("YOL122C", "YNL291C"), ("YOL122C", "YGR217W"),
                     ("YOL130W", "YHR050W"), ("YOL122C", "YPR124W", "YLR411W")):
        with fixed_model:
            for orf in partners:
                fixed_model.genes.get_by_id(orf).knock_out()
            assert solver.growth_or_none(fixed_model) == pytest.approx(0.0, abs=1e-9), partners
    assert sorted(set(inert) - {"YHR050W", "YNL291C", "YGR217W"}) == [
        "YBL082C", "YBR015C", "YDR245W", "YEL036C", "YGL257C", "YIL014W", "YJL183W",
        "YJL186W", "YKR050W", "YLR411W", "YML115C", "YPL050C"]


@pytest.mark.integration
def test_installing_twice_is_refused(fixed_model):
    with pytest.raises(ValueError, match="already installed"):
        restore_stress_genes(fixed_model)


@pytest.mark.integration
def test_a_rule_that_has_drifted_is_refused_rather_than_overwritten(base_model):
    drifted = base_model.copy()
    drifted.reactions.get_by_id("r_0550").gene_reaction_rule = "YDR453C"
    with pytest.raises(ValueError, match="Refusing to overwrite"):
        restore_stress_genes(drifted)


@pytest.fixture(scope="module")
def regulons():
    """SGD's curated regulation records, or a skip: they are gitignored, like the GAF."""
    from ystwin.bridge import regulation

    try:
        return regulation.load_regulons()
    except FileNotFoundError as exc:
        pytest.skip(f"SGD regulation records not vendored here ({exc})")


def _eflux_scales(model, activity, regulons):
    from ystwin.bridge import regulation

    layer = regulation.eflux_layer(model, activity, regulons)
    return layer, {row.reaction: row.scale for row in layer.scales}


@pytest.mark.integration
def test_the_oxidative_regulon_could_not_reach_the_peroxiredoxin_before(
    base_model, fixed_model, regulons
):
    """The false negative in one number, through the repo's own gene-level layer.

    YAP1/SKN7 name TSA1 and TRX2 but not TSA2, so before this module every clause of
    r_0550's rule contained a gene the regulon says nothing about."""
    before_layer, before = _eflux_scales(base_model, {"oxidative": 0.5}, regulons)
    after_layer, after = _eflux_scales(fixed_model, {"oxidative": 0.5}, regulons)
    assert before["r_0550"] == pytest.approx(1.0)
    assert after["r_0550"] == pytest.approx(1.5)
    assert "YML028W" in before_layer.genes_absent
    assert "YML028W" in after_layer.genes_reached
    assert len(after_layer.genes_absent) == len(before_layer.genes_absent) - 2


@pytest.mark.integration
def test_the_copper_regulon_reaches_the_copper_transport_only_after(
    base_model, fixed_model, regulons
):
    _, before = _eflux_scales(base_model, {"copper": 0.5}, regulons)
    _, after = _eflux_scales(fixed_model, {"copper": 0.5}, regulons)
    assert "r_4589" not in before
    assert after["r_4589"] == pytest.approx(1.5)


@pytest.mark.integration
def test_the_ftr1_complex_costs_the_hypoxia_regulon_its_iron_bound(
    base_model, fixed_model, regulons, declarations
):
    """Writing FET3 and FTR1 as one complex loses reach where a regulon names only FET3.
    Declared in the table rather than discovered later."""
    _, before = _eflux_scales(base_model, {"hypoxia": 0.5}, regulons)
    _, after = _eflux_scales(fixed_model, {"hypoxia": 0.5}, regulons)
    assert before["r_1178"] == pytest.approx(1.5)
    assert after["r_1178"] == pytest.approx(1.0)
    basis = next(r.basis for r in declarations if r.orf == "YER145C")
    assert "hypoxia" in basis and "1.500 to 1.000" in basis


@pytest.mark.integration
def test_the_additions_are_not_a_free_lunch(base_model, fixed_model):
    """The probe that has caught every serious bug on this workflow, run on this module.

    A new reaction is the classic way to make a model grow on nothing, so the two regimes
    that would show it are pinned: unlimited exchanges, and a closed model with the forced
    NGAM released. 77.556090 /h is what an all-open yeast-GEM v9.0.2 reports on GLPK."""
    import cobra

    opened = {}
    for label, model in (("before", base_model), ("after", fixed_model)):
        with model:
            for reaction in model.exchanges:
                reaction.bounds = (-1000.0, 1000.0)
            opened[label] = solver.growth_or_none(model)
    assert opened["before"] == pytest.approx(77.556090, abs=1e-6)
    assert opened["after"] == pytest.approx(opened["before"], abs=1e-6)

    drains = {
        "ATP": (("s_0434", -1), ("s_0803", -1), ("s_0394", 1), ("s_1322", 1), ("s_0794", 1)),
        "NADPH": (("s_1212", -1), ("s_1207", 1), ("s_0794", 1)),
        "NADH": (("s_1203", -1), ("s_1198", 1), ("s_0794", 1)),
    }
    for label, stoichiometry in drains.items():
        for model in (base_model, fixed_model):
            with model:
                _close(model)
                probe = cobra.Reaction("FREE_LUNCH_PROBE", lower_bound=0, upper_bound=1000)
                model.add_reactions([probe])
                probe.add_metabolites(
                    {model.metabolites.get_by_id(i): c for i, c in stoichiometry})
                model.objective = probe
                assert solver.growth_or_none(model) == pytest.approx(0.0, abs=1e-9), label

    from cobra.flux_analysis import flux_variability_analysis

    with fixed_model:
        _close(fixed_model)
        fva = flux_variability_analysis(
            fixed_model, reaction_list=list(NEW_REACTION_IDS), fraction_of_optimum=0.0,
            processes=solver.FVA_PROCESSES,
        )
    assert fva.abs().to_numpy().max() == pytest.approx(0.0, abs=1e-9)


def _close(model):
    """Shut every exchange. r_4046 is a maintenance flux FORCED at 0.7, so a closed model
    is infeasible until it is released, and the drain probe would read as a false zero."""
    for reaction in model.exchanges:
        reaction.bounds = (0.0, 0.0)
    model.reactions.get_by_id("r_4046").lower_bound = 0.0


@pytest.mark.integration
def test_manpol_enters_two_more_regulons_than_the_reaction_it_parallels(fixed_model, regulons):
    """MNN1 is in the cell_wall and oxidative regulons, so the new reaction is scaled by
    them. Inert at an FVA maximum of 0.06, but it is a bound this module moves."""
    for module in ("cell_wall", "oxidative"):
        _, after = _eflux_scales(fixed_model, {module: 0.5}, regulons)
        assert after["MANPOL_G"] == pytest.approx(1.5), module
    basis = next(r.basis for r in load_declarations() if r.orf == "YER001W")
    assert "cell_wall" in basis and "1000 -> 1500" in basis


@pytest.mark.integration
def test_only_the_declared_reactions_and_genes_are_added(base_model, fixed_model, declarations):
    added_reactions = {r.id for r in fixed_model.reactions} - {r.id for r in base_model.reactions}
    added_genes = {g.id for g in fixed_model.genes} - {g.id for g in base_model.genes}
    assert added_reactions == set(NEW_REACTION_IDS)
    assert added_genes == {r.orf for r in additions(declarations)}


@pytest.mark.integration
def test_no_refused_gene_ever_enters_the_model(fixed_model, declarations):
    """The refusal list is enforced, not decorative."""
    for record in refusals(declarations):
        assert not fixed_model.genes.has_id(record.orf), record.orf


@pytest.mark.integration
def test_a_refused_edit_leaves_the_input_model_untouched(base_model):
    """The drift is on the sixth row, so earlier rows had already been applied to the copy."""
    drifted = base_model.copy()
    drifted.reactions.get_by_id("r_4590").gene_reaction_rule = "YOL122C"
    with pytest.raises(ValueError, match="Refusing to overwrite"):
        restore_stress_genes(drifted)
    assert drifted.reactions.get_by_id("r_1249").gene_reaction_rule == "YDR456W or YJL129C"
    assert not drifted.reactions.has_id("MANPOL_G")
    assert not drifted.genes.has_id("YKR050W")
