"""The regulation layer: does module activity reach the right reactions, and only those?

`bridge/latent_bridge.py` collapses 24 modules into one ATP scalar, so there was nothing
to test about specificity -- a single number cannot be wrong about which reaction it hit.
A per-reaction layer can be, in both directions: it can miss reactions its regulon
controls, and it can constrain reactions it has no evidence about. Most of what follows
pins the second, because that is the failure that would make the layer look informative
while being arbitrary.

The structural property is `test_eflux_never_raises_a_lower_bound_above_zero`. D1's
relative width is `(max - min) / max` with `min = 0`, so anything that cannot lift a floor
cannot move that number, and asserting it here means the finding is pinned rather than
re-derived by whoever reads the table next.

Toy models rather than the GSMM for everything except the two integration tests: a
five-reaction network makes "exactly these reactions and no others" checkable by hand, and
a 8144-reaction one makes it checkable by nobody.
"""

from __future__ import annotations

import json

import cobra
import pytest

from ystwin.bridge.regulation import (
    COLIJN_PMID,
    EFLUX_LIFTS_LOWER_BOUND,
    MODULE_FACTORS,
    InfeasibleRegulation,
    MetabolicTask,
    MissingMeasurement,
    Regulon,
    apply_eflux,
    apply_task_bounds,
    eflux_layer,
    gene_scales,
    load_regulons,
    reference_capacity,
    task_priority_order,
)


def _reaction(model, rid, stoichiometry, gpr="", bounds=(0.0, 10.0)):
    reaction = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([reaction])
    reaction.add_metabolites(
        {model.metabolites.get_by_id(k): v for k, v in stoichiometry.items()})
    reaction.gene_reaction_rule = gpr
    return reaction


@pytest.fixture
def toy():
    """A-> B-> C with two routes to C, one reversible arm, and one gene-free step.

    Shaped so every fold in a gene-reaction rule is exercised by exactly one reaction:
    ``R1`` is a lone gene, ``R2`` is ``or`` with one regulated and one unregulated
    isozyme, ``R3`` is ``and``, ``R_side`` has no gene at all. ``R1`` is the only route
    out of A, which is what lets one repression make the model infeasible on purpose.
    """
    model = cobra.Model("toy")
    model.add_metabolites([cobra.Metabolite(m) for m in ("A_e", "A", "B", "C")])
    _reaction(model, "EX_A", {"A_e": -1}, bounds=(-10.0, 1000.0))
    _reaction(model, "R_in", {"A_e": -1, "A": 1})
    _reaction(model, "R1", {"A": -1, "B": 1}, gpr="g1")
    _reaction(model, "R2", {"B": -1, "C": 1}, gpr="g2 or g3")
    _reaction(model, "R3", {"B": -1, "C": 1}, gpr="g4 and g5")
    _reaction(model, "R_rev", {"B": -1, "C": 1}, gpr="g1", bounds=(-10.0, 10.0))
    _reaction(model, "R_side", {"C": -1})
    biomass = _reaction(model, "BIO", {"C": -1}, bounds=(0.0, 1000.0))
    model.objective = biomass
    return model


@pytest.fixture
def regulons():
    """Two modules, plus one that is empty and one whose gene is not in the toy model."""
    return {
        "m1": Regulon("m1", ("F1",), ("g1", "g2"), "test"),
        "m2": Regulon("m2", ("F2",), ("g4",), "test"),
        "m_empty": Regulon("m_empty", ("F3",), (), "test"),
        "m_offmodel": Regulon("m_offmodel", ("F4",), ("gZZZ",), "test"),
    }


def _layer(model, activity, regulons, **kwargs):
    return eflux_layer(model, activity, regulons, **kwargs)


def _bounds(model):
    return {r.id: r.bounds for r in model.reactions}


# --------------------------------------------------------------------------- #
# the structural property D1's number turns on
# --------------------------------------------------------------------------- #


def test_eflux_never_raises_a_lower_bound_above_zero(toy, regulons):
    """The whole reason E-Flux cannot rescue D1: a capacity bound never requires flux.

    D1's relative width is `(max - min) / max` and its `min` is zero. If no E-Flux bound
    can be positive, no E-Flux layer can move that number, whatever it does to the
    ceiling. Asserted over both signs and both magnitudes so it reads as structural.
    """
    for activity in ({"m1": -1.0}, {"m1": -0.5}, {"m1": 0.5}, {"m1": 4.0}):
        layer = _layer(toy, activity, regulons)
        assert all(row.lower_bound <= 0.0 for row in layer.scales)
    assert EFLUX_LIFTS_LOWER_BOUND is False


def test_a_reaction_with_a_zero_floor_still_has_a_zero_floor_after_regulation(toy, regulons):
    """The same property read off the model rather than off the layer."""
    before = toy.reactions.R_side.lower_bound
    layer = _layer(toy, {"m1": -0.9, "m2": -0.9}, regulons)
    with toy as scoped:
        apply_eflux(scoped, layer, biomass_reaction="BIO", acknowledge_infeasible=True)

        assert scoped.reactions.R_side.lower_bound == before == 0.0
        assert all(r.lower_bound <= 0.0 for r in scoped.reactions)


# --------------------------------------------------------------------------- #
# specificity
# --------------------------------------------------------------------------- #


def test_an_inactive_state_changes_no_bound(toy, regulons):
    """Zero activity is the unstressed reference, so it must be a no-op exactly.

    Not approximately: `module_response` returns exact zeros for untouched modules, and a
    layer that nudged every bound at zero activity would make every subsequent difference
    partly an artefact of being applied at all.
    """
    before = _bounds(toy)
    layer = _layer(toy, {"m1": 0.0, "m2": 0.0}, regulons)
    with toy as scoped:
        report = apply_eflux(scoped, layer, biomass_reaction="BIO")

        assert _bounds(scoped) == before
        assert report.bounds_changed == 0


def test_a_fully_repressive_module_reduces_exactly_the_reactions_its_regulon_controls(
        toy, regulons):
    """Specificity in both directions: these reactions and no others.

    ``m1`` binds g1 and g2. R1 and R_rev carry g1 alone and must fall to zero. R2 is
    ``g2 or g3`` and must NOT move, because g3 is unregulated and an isozyme at reference
    capacity keeps the reaction's capacity. Everything else must be untouched.
    """
    before = _bounds(toy)
    layer = _layer(toy, {"m1": -1.0}, regulons)
    with toy as scoped:
        apply_eflux(scoped, layer, biomass_reaction="BIO", acknowledge_infeasible=True)
        after = _bounds(scoped)

    assert after["R1"] == (0.0, 0.0)
    assert after["R_rev"] == (0.0, 0.0)
    for untouched in ("R2", "R3", "R_in", "R_side", "EX_A", "BIO"):
        assert after[untouched] == before[untouched], untouched


def test_a_gene_the_panel_says_nothing_about_is_not_treated_as_repressed(toy, regulons):
    """Silence is not evidence of repression.

    A 19-module panel names a few hundred genes. If an unnamed gene scored zero instead of
    one, the panel would constrain the whole proteome and every narrowing would be an
    artefact of the panel's coverage rather than of its content.
    """
    layer = _layer(toy, {"m1": -1.0}, regulons)
    by_id = {row.reaction: row for row in layer.scales}

    assert by_id["R2"].scale == 1.0  # g2 repressed, g3 unregulated, `or` keeps capacity
    assert "R3" not in by_id  # g4/g5 are in no active regulon, so not even a candidate


def test_a_complex_is_limited_by_its_scarcest_subunit(toy, regulons):
    """``and`` folds to the minimum: one missing subunit stops the whole complex."""
    layer = _layer(toy, {"m2": -1.0}, regulons)
    by_id = {row.reaction: row for row in layer.scales}

    assert by_id["R3"].scale == 0.0  # g4 repressed to zero, g5 unregulated at 1.0


# --------------------------------------------------------------------------- #
# direction
# --------------------------------------------------------------------------- #


def test_the_signed_direction_is_honoured_rather_than_absolute_valued(toy, regulons):
    """A negative weight must lower a bound. Taking |w| would silently invert it.

    ``STRESSORS[...].targets`` carries negative weights for modules a stressor represses,
    and `module_response` multiplies them through, so the sign arrives here. The two
    assertions are deliberately mirror images: equal magnitudes, opposite effects.
    """
    down = {row.reaction: row.scale for row in _layer(toy, {"m1": -0.5}, regulons).scales}
    up = {row.reaction: row.scale for row in _layer(toy, {"m1": 0.5}, regulons).scales}

    assert down["R1"] == pytest.approx(0.5)
    assert up["R1"] == pytest.approx(1.5)


def test_activities_on_one_gene_add_before_becoming_a_scale(toy):
    """Two modules on one promoter add, matching ``combination_response``'s convention."""
    shared = {
        "a": Regulon("a", ("F",), ("g1",), "test"),
        "b": Regulon("b", ("F",), ("g1",), "test"),
    }

    assert gene_scales({"a": -0.3, "b": -0.4}, shared)["g1"] == pytest.approx(0.3)


def test_repression_past_total_saturates_at_zero_rather_than_going_negative(toy):
    """A negative capacity is not a quantity, so the scale clips instead of inverting."""
    single = {"a": Regulon("a", ("F",), ("g1",), "test")}

    assert gene_scales({"a": -5.0}, single)["g1"] == 0.0


# --------------------------------------------------------------------------- #
# refusals
# --------------------------------------------------------------------------- #


def test_a_module_with_an_empty_regulon_is_refused_rather_than_doing_nothing(toy, regulons):
    """Silently dropping it would report a specificity the layer does not have.

    This is the failure `latent_bridge.py` has by construction -- every module's
    specificity discarded, nothing said about it -- so the replacement has to name it.
    """
    with pytest.raises(ValueError, match="cannot reach any reaction"):
        _layer(toy, {"m_empty": -0.8}, regulons)


def test_a_module_with_no_transcription_factor_at_all_is_refused_by_a_distinct_message(
        toy, regulons):
    """A metabolite pool has no regulon in principle, which is a different fact.

    The five pool modules cannot be transcribed; `calcium` has a factor and no targets at
    this evidence level. Both must be declared, and a reader needs to be able to tell
    "there is nothing to find" from "we did not find it".
    """
    with pytest.raises(ValueError, match="not a transcriptional module"):
        _layer(toy, {"redox": -0.9}, regulons)


def test_declaring_an_unmappable_module_lets_the_rest_of_the_state_through(toy, regulons):
    """The escape hatch is explicit and per-module, so the exclusion is on the record."""
    layer = _layer(toy, {"m_empty": -0.8, "m1": -0.5}, regulons,
                   allow_unmapped={"m_empty"})

    assert layer.modules_excluded == ("m_empty",)
    assert {row.reaction for row in layer.binding} == {"R1", "R_rev"}


def test_a_zero_activity_module_with_no_regulon_needs_no_declaration(toy, regulons):
    """A full ``module_response`` dict is mostly zeros; refusing on those would be noise."""
    layer = _layer(toy, {"m_empty": 0.0, "m1": -0.5}, regulons)

    assert layer.modules_excluded == ()


def test_a_gene_absent_from_the_gem_is_recorded_and_skipped_rather_than_crashing(
        toy, regulons):
    """Two thirds of these regulons' genes are not metabolic; that is normal, not an error."""
    layer = _layer(toy, {"m_offmodel": -1.0}, regulons)

    assert layer.genes_absent == ("gZZZ",)
    assert layer.genes_reached == ()
    assert layer.scales == ()


def test_regulation_that_makes_the_model_infeasible_says_so(toy, regulons):
    """An infeasible LP returns ``None``, which a caller reads as a flux of zero.

    R1 is the only route out of A, so repressing g1 to nothing stops growth. The raise
    names the tightest bounds, because "infeasible" on its own does not say which
    constraint to look at.
    """
    toy.reactions.R_side.lower_bound = 1.0
    layer = _layer(toy, {"m1": -1.0}, regulons)
    with toy as scoped:
        with pytest.raises(InfeasibleRegulation, match="unable to grow"):
            apply_eflux(scoped, layer, biomass_reaction="BIO")


def test_infeasibility_can_be_measured_instead_of_raised_when_asked_for_explicitly(
        toy, regulons):
    """Some analyses are about where the layer breaks the model, so the gate is a flag."""
    layer = _layer(toy, {"m1": -1.0}, regulons)
    with toy as scoped:
        report = apply_eflux(scoped, layer, biomass_reaction="BIO",
                            acknowledge_infeasible=True)

    assert report.max_growth == 0.0


# --------------------------------------------------------------------------- #
# applying it
# --------------------------------------------------------------------------- #


def test_applying_the_same_layer_twice_changes_nothing_the_second_time(toy, regulons):
    """Idempotent because bounds are assigned from the layer, never multiplied into.

    A layer that multiplied the current bound would compound on re-application, and the
    second call is easy to make by accident -- a retry, a loop over growth fractions, a
    nested context.
    """
    layer = _layer(toy, {"m1": -0.5, "m2": -0.5}, regulons)
    with toy as scoped:
        first = apply_eflux(scoped, layer, biomass_reaction="BIO")
        once = _bounds(scoped)
        second = apply_eflux(scoped, layer, biomass_reaction="BIO")

        assert _bounds(scoped) == once
        assert second.bounds_changed == 0
        assert first.bounds_changed > 0


def test_computing_a_layer_does_not_touch_the_model(toy, regulons):
    """The layer is a value, so it can be inspected and diffed before anything solves."""
    before = _bounds(toy)
    _layer(toy, {"m1": -1.0, "m2": -1.0}, regulons)

    assert _bounds(toy) == before


def test_a_reversible_reaction_has_both_directions_scaled(toy, regulons):
    """E-Flux bounds ``|v|``, so the reverse capacity scales by the same factor.

    Scaling only the upper bound would leave a repressed reversible reaction free to run
    backwards at full rate, which is not a capacity constraint at all.
    """
    layer = _layer(toy, {"m1": -0.5}, regulons)
    row = next(r for r in layer.scales if r.reaction == "R_rev")

    assert (row.lower_bound, row.upper_bound) == pytest.approx((-5.0, 5.0))


def test_the_layer_records_which_reference_capacity_it_used(toy, regulons):
    """The reference is the load-bearing choice, so a result must carry which one."""
    implied = _layer(toy, {"m1": -0.5}, regulons)
    supplied = _layer(toy, {"m1": -0.5}, regulons, capacities={"R1": 2.0})

    assert implied.reference == "model-implied"
    assert supplied.reference == "caller-supplied capacities"
    assert next(r for r in supplied.scales if r.reaction == "R1").upper_bound == 1.0


def test_the_citation_travels_with_the_layer(toy, regulons):
    """A bound whose method is not named beside it is a number with no provenance."""
    assert str(COLIJN_PMID) in _layer(toy, {"m1": -0.5}, regulons).citation


# --------------------------------------------------------------------------- #
# reference capacity
# --------------------------------------------------------------------------- #


def test_an_enzyme_draw_is_capped_by_the_protein_budget_it_draws_from(toy):
    """GECKO's own arithmetic: the most of one enzyme the pool can buy is cap / MW.

    Read off the model rather than asserted, which is what keeps this layer free of a
    second ``_MAINTENANCE_PER_ACTIVITY``.
    """
    toy.add_metabolites([cobra.Metabolite("prot_pool[c]"), cobra.Metabolite("prot_E1[c]")])
    _reaction(toy, "prot_pool_exchange", {"prot_pool[c]": 1}, bounds=(0.0, 0.1))
    draw = _reaction(toy, "draw_prot_E1", {"prot_pool[c]": -25.0, "prot_E1[c]": 1},
                     gpr="g1", bounds=(0.0, float("inf")))

    assert reference_capacity(draw, (0.1, "prot_pool[c]")) == pytest.approx((0.004, 0.0))


def test_a_reaction_with_no_finite_ceiling_has_no_reference_and_is_reported_as_skipped(toy):
    """Scaling infinity yields infinity, so pretending to bound it would overstate reach."""
    unbounded = _reaction(toy, "R_open", {"C": -1}, gpr="g1",
                          bounds=(0.0, float("inf")))

    assert reference_capacity(unbounded, None) is None

    layer = eflux_layer(toy, {"m": -1.0},
                        {"m": Regulon("m", ("F",), ("g1",), "test")})

    assert "R_open" not in {row.reaction for row in layer.scales}
    assert layer.unbounded_skipped >= 1


# --------------------------------------------------------------------------- #
# module -> gene, from SGD
# --------------------------------------------------------------------------- #


def _record(regulator, target, pmid=16522208, kind="transcription"):
    return {
        "regulation_of": kind,
        "reference": {"pubmed_id": pmid},
        "locus1": {"display_name": regulator, "format_name": f"ORF_{regulator}"},
        "locus2": {"display_name": target, "format_name": f"ORF_{target}"},
    }


def test_only_records_where_the_factor_is_the_regulator_enter_its_regulon(tmp_path):
    """SGD's files are per-gene, not per-role, so direction has to be read off the record.

    ``CAT8.json`` really does carry two conserved-motif records in which CAT8 is somebody
    else's target and none in which it is the regulator. Taking ``locus2`` unconditionally
    turns that into a one-gene regulon consisting of CAT8 itself.
    """
    (tmp_path / "F1.json").write_text(json.dumps([
        _record("F1", "TARGET"),          # kept
        _record("OTHER", "F1"),           # F1 is the target here, not the regulator
        _record("F1", "WRONG_PMID", pmid=11102521),
        _record("F1", "WRONG_KIND", kind="protein activity"),
    ]))

    regulons = load_regulons(tmp_path, factors={"m": ("F1",)})

    assert regulons["m"].genes == ("ORF_TARGET",)


def test_a_factor_with_no_records_at_this_evidence_level_gets_an_empty_regulon(tmp_path):
    """Empty and visible, not dropped: the gap is a result about the evidence."""
    (tmp_path / "F1.json").write_text(json.dumps([_record("F1", "T", pmid=999)]))

    regulons = load_regulons(tmp_path, factors={"m": ("F1",)})

    assert regulons["m"].genes == ()
    assert "m" in regulons


def test_a_missing_factor_file_is_an_error_and_not_an_empty_regulon(tmp_path):
    """An unknown regulon read as empty would silently leave that module unregulated."""
    with pytest.raises(FileNotFoundError, match="ABSENT"):
        load_regulons(tmp_path, factors={"m": ("ABSENT",)})


def test_every_named_module_factor_is_a_module_the_panel_defines():
    """Guards the one hand-written table here against drifting from ``MODULES``."""
    from ystwin.generator.stress_panel import MODULES

    assert set(MODULE_FACTORS) <= set(MODULES)
    pools = {"redox", "peroxide", "atp", "ph", "nadh"}
    assert set(MODULES) - set(MODULE_FACTORS) == pools


# --------------------------------------------------------------------------- #
# lexicographic task efficiency -- the mechanism that can lift a floor
# --------------------------------------------------------------------------- #


def test_a_task_with_no_measured_flux_is_refused_and_the_message_names_the_units(toy):
    """The refusal is the deliverable: no strain carries the pathway, so there is no q_p.

    Substituting a simulated flux would put the FBA maximum in both numerator and
    denominator, make every efficiency 1.0, and produce a narrow interval that is an
    artefact of the substitution rather than a result.
    """
    with pytest.raises(MissingMeasurement, match="mmol/gDW/h"):
        task_priority_order(toy, [MetabolicTask("product", "R_side")])


def test_a_measured_flux_without_a_source_is_refused_too(toy):
    """An unsourced number is indistinguishable from an assumption, which is the whole risk."""
    with pytest.raises(MissingMeasurement, match="no source"):
        task_priority_order(toy, [MetabolicTask("product", "R_side", 1.0)])


def test_tasks_are_ranked_by_measured_over_maximum_and_fixed_in_that_order(toy):
    """The ordering is lexicographic: fixing the leader is what re-ranks the rest."""
    order = task_priority_order(toy, [
        MetabolicTask("growth", "BIO", 8.0, "test"),
        MetabolicTask("side", "R_side", 1.0, "test"),
    ])

    assert [p.name for p in order] == ["growth", "side"]
    assert order[0].efficiency == pytest.approx(0.8)
    assert order[0].rank == 1


def test_a_measurement_the_model_cannot_reach_is_a_contradiction_not_a_priority(toy):
    """Efficiency above 1 means the model and the measurement disagree; say that.

    This fires on the real models: van Hoek's 0.40 /h aerobic batch growth exceeds
    ecYeastGEM's own maximum of 0.3768 /h, so biomass cannot be ranked at all.
    """
    with pytest.raises(ValueError, match="contradiction"):
        task_priority_order(toy, [MetabolicTask("growth", "BIO", 50.0, "test")])


def test_ranking_tasks_leaves_the_model_as_it_was_found(toy):
    """The ordering fixes bounds while it works and must not leak them to the caller."""
    before = _bounds(toy)
    task_priority_order(toy, [MetabolicTask("growth", "BIO", 8.0, "test")])

    assert _bounds(toy) == before


def test_an_empty_task_list_is_refused(toy):
    """A priority order over nothing is not a result."""
    with pytest.raises(ValueError, match="no tasks"):
        task_priority_order(toy, [])


def test_task_bounds_do_lift_lower_bounds_which_is_the_contrast_with_eflux(toy):
    """The one mechanism here that can move a floor, and it moves it to a measurement."""
    order = task_priority_order(toy, [MetabolicTask("growth", "BIO", 8.0, "test")])
    with toy as scoped:
        report = apply_task_bounds(scoped, order)

        assert scoped.reactions.BIO.lower_bound == pytest.approx(8.0)
        assert report.lifted_lower_bounds == 1


def test_only_the_higher_priority_tasks_can_be_fixed_when_asked(toy):
    """Asks what the ranked tasks imply about one whose own rate is not measured."""
    order = task_priority_order(toy, [
        MetabolicTask("growth", "BIO", 8.0, "test"),
        MetabolicTask("side", "R_side", 1.0, "test"),
    ])
    with toy as scoped:
        apply_task_bounds(scoped, order, up_to_rank=1)

        assert scoped.reactions.BIO.lower_bound == pytest.approx(8.0)
        assert scoped.reactions.R_side.lower_bound == 0.0


# --------------------------------------------------------------------------- #
# against the real models
# --------------------------------------------------------------------------- #


@pytest.mark.integration
def test_the_real_mapping_reaches_metabolism_through_both_gems(ec_yeast_gem, yeast_gem):
    """A narrowing computed on a mapping that reaches nothing is a number about nothing.

    Pinned loosely -- these are counts over vendored data, and a tight number here would
    fail on the next model version for no reason worth a test failure.
    """
    regulons = load_regulons()
    activity = {module: 0.5 for module in MODULE_FACTORS}
    for model, floor in ((ec_yeast_gem, 200), (yeast_gem, 200)):
        layer = eflux_layer(model, activity, regulons,
                           allow_unmapped={"calcium", "zinc", "xenobiotic", "UPR"})

        assert len(layer.genes_reached) > 100
        assert len(layer.scales) + layer.unbounded_skipped > floor


@pytest.mark.integration
def test_the_upr_regulon_reaches_no_metabolic_reaction_in_the_gem(ec_yeast_gem):
    """The module the wet lab's own constructs report on cannot constrain metabolism.

    UPRE1 and UPRE2 are two of the four constructs on the plates, and Hac1's conserved-motif
    regulon is seven genes, none of which the GSMM carries. Pinned because it is the single
    most consequential gap in the mapping and it should fail loudly if it ever closes.

    It is *named*, not raised on: a full ``module_response`` dict is mostly non-metabolic
    regulons and refusing on those would make the layer unusable. Being unnamed is the
    failure; being empty is a fact about the mapping.
    """
    regulons = load_regulons()
    layer = eflux_layer(ec_yeast_gem, {"UPR": 0.5}, regulons)

    assert len(regulons["UPR"].genes) == 7
    assert not set(regulons["UPR"].genes) & {g.id for g in ec_yeast_gem.genes}
    assert layer.scales == ()
    assert layer.modules_unreached == ("UPR",)
    assert "reached no reaction: ['UPR']" in layer.summary()


def test_reverse_stoichiometric_protein_pool_supply_sets_reference_capacity(toy):
    from ystwin.bridge.regulation import _protein_pool_cap

    toy.add_metabolites([cobra.Metabolite("prot_pool[c]")])
    _reaction(toy, "pool_supply", {"prot_pool[c]": -2.0}, bounds=(-0.05, 0.0))
    assert _protein_pool_cap(toy) == pytest.approx((0.1, "prot_pool[c]"))


def test_measured_task_is_a_floor_not_an_equality(toy):
    order = task_priority_order(toy, [MetabolicTask("growth", "BIO", 8.0, "test")])
    with toy as scoped:
        report = apply_task_bounds(scoped, order, biomass_reaction="BIO")
        assert scoped.reactions.BIO.bounds == (8.0, 1000.0)
        assert report.max_growth == pytest.approx(10.0)
        assert report.status == "optimal"


def test_coupled_task_can_exceed_its_measurement_without_false_infeasibility(toy):
    toy.add_cons_vars([toy.problem.Constraint(toy.reactions.BIO.flux_expression -
                                              toy.reactions.R_side.flux_expression, lb=0.0, ub=0.0)])
    order = task_priority_order(toy, [MetabolicTask("growth", "BIO", 4.0, "test"),
                                      MetabolicTask("side", "R_side", 1.0, "test")])
    with toy as scoped:
        report = apply_task_bounds(scoped, order, biomass_reaction="BIO")
        solution = scoped.optimize()
        assert report.max_growth == pytest.approx(5.0)
        assert solution.fluxes["R_side"] == pytest.approx(5.0)
        assert scoped.reactions.R_side.lower_bound == 1.0
        assert scoped.reactions.R_side.upper_bound == 10.0


def test_task_order_sets_maximization_explicitly(toy):
    toy.objective_direction = "min"
    order = task_priority_order(toy, [MetabolicTask("growth", "BIO", 8.0, "test")])
    assert order[0].fba_maximum == pytest.approx(10.0)
    assert toy.objective_direction == "min"


@pytest.mark.parametrize("value", [-1.0, float("nan"), float("inf"), True])
def test_invalid_measured_task_rate_is_refused_before_solving(toy, value):
    with pytest.raises(ValueError, match="finite nonnegative"):
        task_priority_order(toy, [MetabolicTask("growth", "BIO", value, "test")])


def test_duplicate_tasks_are_not_silently_dropped(toy):
    with pytest.raises(ValueError, match="duplicate"):
        task_priority_order(toy, [
            MetabolicTask("growth", "BIO", 1.0, "test"),
            MetabolicTask("growth again", "BIO", 2.0, "test"),
        ])


def test_feasible_zero_task_and_infeasible_model_are_distinct(toy, regulons):
    toy.reactions.BIO.upper_bound = 0.0
    order = task_priority_order(toy, [MetabolicTask("growth", "BIO", 0.0, "test")])
    assert order[0].efficiency is None
    report = apply_task_bounds(toy, order, biomass_reaction="BIO")
    assert report.status == "optimal"
    assert report.max_growth == 0.0
    toy.reactions.R_side.lower_bound = 1.0
    layer = _layer(toy, {"m1": -1.0}, regulons)
    report = apply_eflux(toy, layer, biomass_reaction="BIO", acknowledge_infeasible=True)
    assert report.status == "infeasible"
    assert report.max_growth is None


@pytest.mark.parametrize("coefficient, bounds", [(-2.0, (-10.0, 1000.0)), (2.0, (0.0, 10.0))])
def test_uptake_uses_stoichiometric_magnitude_not_reaction_sign_heuristics(coefficient, bounds):
    from ystwin.bridge.regulation import constrain_uptake

    model = cobra.Model("scaled_exchange")
    model.add_metabolites([cobra.Metabolite("carbon", compartment="e"),
                           cobra.Metabolite("carbon_c", compartment="c")])
    _reaction(model, "feed", {"carbon": coefficient}, bounds=bounds)
    _reaction(model, "transport", {"carbon": -1.0, "carbon_c": 1.0}, bounds=(0.0, 1000.0))
    _reaction(model, "BIO", {"carbon_c": -1.0}, bounds=(0.0, 1000.0))
    model.objective = "BIO"
    constrain_uptake(model, "feed", 3.0, source="synthetic uptake")
    assert model.slim_optimize() == pytest.approx(3.0)


def test_split_uptake_discovers_all_boundary_routes_from_stoichiometry():
    from ystwin.bridge.regulation import constrain_uptake, exchange_coefficients

    model = cobra.Model("split")
    model.add_metabolites([cobra.Metabolite("carbon", compartment="e"),
                           cobra.Metabolite("carbon_c", compartment="c")])
    _reaction(model, "EX", {"carbon": -1.0}, bounds=(0.0, 0.0))
    _reaction(model, "unrelated_import_name", {"carbon": 2.0}, bounds=(0.0, 1000.0))
    _reaction(model, "another_import", {"carbon": 1.0}, bounds=(0.0, 1000.0))
    _reaction(model, "transport", {"carbon": -1.0, "carbon_c": 1.0}, bounds=(0.0, 1000.0))
    _reaction(model, "BIO", {"carbon_c": -1.0}, bounds=(0.0, 1000.0))
    model.objective = "BIO"
    coefficients = exchange_coefficients(model, "EX")
    assert coefficients["unrelated_import_name"] == 2.0
    constrain_uptake(model, "EX", 3000.0, unit="umol/gDW/h", source="synthetic uptake")
    assert model.reactions.EX.bounds == (0.0, 0.0)
    assert model.slim_optimize() == pytest.approx(3.0)


def test_uptake_constraint_has_explicit_cap_and_measured_modes(toy):
    from ystwin.bridge.regulation import constrain_uptake

    with toy as scoped:
        constrain_uptake(scoped, "EX_A", 4.0, mode="fixed", source="synthetic uptake")
        scoped.objective = "EX_A"
        scoped.objective_direction = "max"
        assert scoped.slim_optimize() == pytest.approx(-4.0)
    with toy as scoped:
        constrain_uptake(scoped, "EX_A", 4.0, source="synthetic uptake")
        scoped.objective = "EX_A"
        scoped.objective_direction = "max"
        assert scoped.slim_optimize() == pytest.approx(0.0)


def test_eflux_does_not_erase_an_existing_required_task_floor(toy, regulons):
    toy.reactions.R1.lower_bound = 8.0
    layer = eflux_layer(toy, {"m1": -0.5}, regulons)
    report = apply_eflux(toy, layer, biomass_reaction="BIO", acknowledge_infeasible=True)
    assert report.status == "infeasible"
    assert report.max_growth is None
    assert toy.reactions.R1.lower_bound == 8.0


def test_task_measurement_condition_and_host_must_match(toy):
    task = MetabolicTask("growth", "BIO", 1.0, "test", unit="1/h",
                         condition_id="culture-a", host="yeast")
    with pytest.raises(ValueError, match="condition"):
        task_priority_order(toy, [task], condition_id="culture-b", host="yeast")
    with pytest.raises(ValueError, match="host"):
        task_priority_order(toy, [task], condition_id="culture-a", host="CHO")
