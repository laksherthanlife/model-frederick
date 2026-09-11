"""The GEM's honest job: audit a flux that arrived from elsewhere, and never supply one.

Capping every reaction of the installed pathway at 1e-3, then at 1e-4 -- the measured flux
magnitude -- gives FVA ranges of [0, 1e-3] and [0, 1e-4]. The ceiling becomes exactly the
cap and the floor never leaves zero, because a capacity constraint is an upper bound and an
upper bound cannot make a flux mandatory. The first class here pins that floor, since it is
the entire reason ``pathway/flux.py`` exists at all: if the GEM could ever be persuaded to
return a positive minimum, the expression law would be redundant.

The second thing this file guards is the growth cost, which is regime-dependent, and was
wrong for that reason. The cost was being computed against the model's own default uptake
while the envelope it was subtracted from came from ``glucose_uptake``, so the subtraction
compared two different cultures -- and reported a 90% growth cost for a flux of 9.5e-4
mmol/gDCW/h, four orders of magnitude below a ceiling it had just cleared by 106x. Both
sides now use the same uptake, ``glucose_uptake`` is required rather than defaulted, and it
is carried out in the result. A burden number without its uptake regime is not a number.

Every solve here is shared through module-scoped fixtures: the SBML parse alone is seconds,
and the point of these tests is the arithmetic around the LP, not the LP.
"""

from __future__ import annotations

import pytest

from ystwin.fba.audit import audit_predicted_flux, precursor_floor
from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID, add_beta_carotene_pathway
from ystwin.pathway.spec import load_pathway

pytestmark = pytest.mark.integration

# The flux the expression law predicts for a strain at relative CrtE expression 1.0, which
# is what the audit layer is actually asked about.
PREDICTED = 1.0085e-3
# A bound nobody measured, and the measured uptake of Elizondo's slow chemostat. The whole
# point of the pair is that the cost is not the same number under the two.
LIBRARY_BOUND, MEASURED_UPTAKE = 10.0, 1.22


@pytest.fixture(scope="module")
def gsmm():
    """Yeast9, or a skip naming the variable that would find it.

    ``paths._resolve`` returns whatever an override points at as long as it exists, so an
    override aimed at a directory resolves to the directory and the load raises "the SBML
    model is not valid" instead of skipping. Running this suite without its assets is a
    supported thing to do, so the guard checks for a file and not merely for a path.
    """
    from ystwin import paths

    path = paths.yeast_gem()
    if path is None or not path.is_file():
        pytest.skip("yeast-GEM not present; set YSTWIN_YEAST_GEM to yeast-GEM.xml")
    import cobra

    return cobra.io.read_sbml_model(str(path))


@pytest.fixture(scope="module")
def carotenoid_gsmm(gsmm):
    """Yeast9 carrying crtE/crtYB/crtI. Module-scoped: the audits do not modify it."""
    return add_beta_carotene_pathway(gsmm)


@pytest.fixture(scope="module")
def audit_at_library_bound(carotenoid_gsmm):
    return audit_predicted_flux(carotenoid_gsmm, PRODUCT_DEMAND_ID, PREDICTED,
                                glucose_uptake=LIBRARY_BOUND)


@pytest.fixture(scope="module")
def audit_at_measured_uptake(carotenoid_gsmm):
    return audit_predicted_flux(carotenoid_gsmm, PRODUCT_DEMAND_ID, PREDICTED,
                                glucose_uptake=MEASURED_UPTAKE)


@pytest.fixture(scope="module")
def audit_above_the_ceiling(carotenoid_gsmm):
    return audit_predicted_flux(carotenoid_gsmm, PRODUCT_DEMAND_ID, 1.0,
                                glucose_uptake=LIBRARY_BOUND)


class TestTheFeasibleMinimumIsZero:
    def test_the_network_never_requires_any_beta_carotene(self, audit_at_library_bound):
        """The result that sends product forecasting off FBA and onto kinetics.

        Growth is held at 90% of maximum and the demand reaction is still free to sit at
        zero, because nothing in the stoichiometry makes a heterologous product mandatory.
        Every enzyme-constrained variant inherits it: GECKO's ``v <= kcat*[E]`` under a
        shared protein pool has the zero vector feasible too. A positive minimum here would
        mean this model is not the model the claim was established on.
        """
        assert audit_at_library_bound.feasible_min == pytest.approx(0.0, abs=1e-9)

    def test_the_ceiling_is_far_above_anything_the_strains_make(self,
                                                               audit_at_library_bound):
        """Stoichiometry alone permits about a hundred times the measured flux.

        Which is the other half of the same point: an envelope this wide bounds the
        prediction without informing it.
        """
        assert audit_at_library_bound.headroom > 50

    def test_a_predicted_flux_inside_the_envelope_is_reported_as_inside(
            self, audit_at_library_bound):
        assert audit_at_library_bound.within_envelope is True

    def test_a_flux_inside_the_envelope_draws_no_complaint(self, audit_at_library_bound):
        assert audit_at_library_bound.notes == ()


class TestTheGrowthCostIsAPropertyOfTheUptake:
    def test_the_same_flux_costs_more_growth_under_the_measured_uptake(
            self, audit_at_library_bound, audit_at_measured_uptake):
        """The reason ``glucose_uptake`` is a required argument and not a defaulted one.

        One flux, one model, one pathway -- and the cost moves by nearly an order of
        magnitude with an input that a caller could easily leave unstated. Reporting a
        burden without naming the regime it was computed under is reporting the default.
        """
        assert (audit_at_measured_uptake.growth_cost_fraction
                > 3 * audit_at_library_bound.growth_cost_fraction)

    def test_a_flux_far_below_the_ceiling_costs_almost_no_growth(
            self, audit_at_library_bound):
        """The defect this file was written after, stated as a number.

        The cost used to be computed against the model's own default uptake while the
        maximum it was subtracted from came from ``glucose_uptake``, so the subtraction
        compared two different cultures. It reported a 90% growth cost for a flux of
        9.5e-4 mmol/gDCW/h -- a flux the very same audit had just found 106x below the
        ceiling. A number that large for a flux this small is arithmetically impossible,
        and that is what this bound says.
        """
        assert audit_at_library_bound.growth_cost_fraction < 0.01

    def test_the_cost_is_still_a_real_cost_and_not_a_rounding(self,
                                                             audit_at_measured_uptake):
        """Guards the other direction: a cost pinned at zero would hide the burden entirely.

        Carotenoid carbon leaves the mevalonate pathway and does not come back, so at a
        measured uptake of 1.22 the drain is visible -- under a percent, but not nothing.
        """
        assert audit_at_measured_uptake.growth_cost_fraction > 1e-4

    def test_the_uptake_regime_is_carried_out_in_the_result(self,
                                                            audit_at_measured_uptake):
        assert audit_at_measured_uptake.glucose_uptake == MEASURED_UPTAKE

    def test_the_summary_names_the_uptake_it_was_computed_under(self,
                                                               audit_at_measured_uptake):
        """A burden quoted in prose is the form most likely to travel without its regime."""
        assert "glucose <= 1.22" in audit_at_measured_uptake.summary()

    def test_a_tighter_uptake_lowers_the_ceiling_too(self, audit_at_library_bound,
                                                     audit_at_measured_uptake):
        assert audit_at_measured_uptake.feasible_max < audit_at_library_bound.feasible_max


class TestAFluxAboveTheCeilingIsRefuted:
    def test_it_is_reported_as_outside_the_envelope(self, audit_above_the_ceiling):
        """Stoichiometry alone refutes it, with no kinetics and no fitted constant.

        This is the one verdict the GEM can deliver on its own, and it has to be a flag a
        caller can branch on rather than a note somebody has to read.
        """
        assert audit_above_the_ceiling.within_envelope is False

    def test_it_arrives_with_a_note_naming_the_envelope_it_left(self,
                                                               audit_above_the_ceiling):
        assert any("cannot carry" in note for note in audit_above_the_ceiling.notes)

    def test_the_note_quotes_the_flux_that_was_asked_for(self, audit_above_the_ceiling):
        assert any("predicted 1" in note for note in audit_above_the_ceiling.notes)


class TestWhatTheAuditRefuses:
    def test_the_uptake_regime_cannot_be_omitted(self, carotenoid_gsmm):
        """A defaulted uptake is how the growth cost went wrong in the first place.

        The signature refuses rather than choosing 10.0, so a caller with a measured uptake
        cannot accidentally audit against a library default, and a caller without one has
        to decide what regime the burden is being quoted for.
        """
        with pytest.raises(TypeError):
            audit_predicted_flux(carotenoid_gsmm, PRODUCT_DEMAND_ID, PREDICTED)

    def test_a_non_positive_uptake_is_refused(self, carotenoid_gsmm):
        with pytest.raises(ValueError, match="glucose_uptake must be positive"):
            audit_predicted_flux(carotenoid_gsmm, PRODUCT_DEMAND_ID, PREDICTED,
                                 glucose_uptake=0.0)

    def test_a_negative_flux_is_refused(self, carotenoid_gsmm):
        """A negative product flux is the pathway running backwards, which it does not.

        Refused before any LP runs, so the caller sees the input error rather than an
        infeasibility that reads like a claim about the network.
        """
        with pytest.raises(ValueError, match="predicted_flux must be non-negative"):
            audit_predicted_flux(carotenoid_gsmm, PRODUCT_DEMAND_ID, -1e-6,
                                 glucose_uptake=LIBRARY_BOUND)


class TestThePrecursorBudgetFailsByName:
    def test_a_metabolite_the_model_does_not_carry_raises_by_name(self,
                                                                  carotenoid_gsmm):
        """Returning zero would read as "no native demand on this node", which is a claim.

        The precursor id comes from the pathway spec, so a typo there or a metabolite
        renamed by a Yeast9 release both land here, and both need the name to be diagnosed.
        """
        with pytest.raises(KeyError, match="s_9999"):
            precursor_floor(carotenoid_gsmm, "s_9999")

    def test_the_precursor_named_by_the_shipped_spec_is_in_this_model(self,
                                                                      carotenoid_gsmm):
        """``precursor_metabolite`` is a Yeast9 id declared in TOML and resolved here.

        Nothing between the two files checks that the id is real, so a Yeast9 release that
        renumbers GGPP would surface as the KeyError above during a run rather than as a
        spec that no longer matches its host. Asserted without an LP: reading the floor
        itself is an FVA over every consumer of the node.
        """
        spec = load_pathway("beta_carotene")

        assert carotenoid_gsmm.metabolites.has_id(spec.precursor_metabolite)


class TestAnInfeasibleSolveIsNotSilentlyNan:
    """`slim_optimize()` signals failure with **nan**, not None, and both obvious guards
    are silently wrong about nan:

        value is None   ->  False, so the "infeasible" branch is dead code
        value <= 0      ->  False, so a "cannot grow" check lets it straight through

    Three call sites in `ystwin.fba` had one or the other. The audit's "refuted by
    stoichiometry rather than merely expensive" note could never fire; a nan propagated into
    a reported growth cost instead, which reads as a missing measurement rather than a
    refuted prediction. `error_value=None` is not the fix either -- it makes cobra raise
    `Infeasible` rather than return anything -- so the check has to be an explicit isnan.
    """

    @staticmethod
    def _infeasible():
        """A model whose only producer overproduces what its only consumer can take."""
        import cobra

        model = cobra.Model("infeasible")
        species = cobra.Metabolite("a")
        source = cobra.Reaction("source")
        source.add_metabolites({species: 1})
        source.bounds = (5.0, 5.0)
        sink = cobra.Reaction("sink")
        sink.add_metabolites({species: -1})
        sink.bounds = (0.0, 1.0)
        model.add_reactions([source, sink])
        model.objective = "sink"
        return model

    def test_cobra_really_does_return_nan_rather_than_none(self):
        """The premise. If a cobra upgrade changed this, the helper would be dead weight
        and the tests below would pass for the wrong reason."""
        import math

        value = self._infeasible().slim_optimize()

        assert isinstance(value, float) and math.isnan(value)

    def test_and_neither_obvious_guard_catches_it(self):
        value = self._infeasible().slim_optimize()

        assert (value is None) is False
        assert (value <= 0) is False

    def test_the_helper_turns_it_into_none(self):
        from ystwin.fba.solver import growth_or_none

        assert growth_or_none(self._infeasible()) is None

    def test_a_feasible_model_still_returns_its_value(self):
        """The other half: a guard that returned None for everything would also pass the
        test above."""
        import cobra

        from ystwin.fba.solver import growth_or_none

        model = cobra.Model("feasible")
        reaction = cobra.Reaction("r")
        reaction.bounds = (0.0, 3.0)
        model.add_reactions([reaction])
        model.objective = "r"

        assert growth_or_none(model) == pytest.approx(3.0)


class TestRegulationIsWiredAndProvablyInert:
    """E-Flux now reaches the audit, and the honest result is that it cannot do anything.

    Wiring it here rather than into the prediction was the right *place*: scaling a ceiling
    cannot set a flux, but it can narrow the envelope a prediction is checked against, and
    that is a real question. Having asked it, the answer is that for every stressor at every
    dose, E-Flux tightens **zero** reaction bounds in Yeast9. Two independent reasons, and
    both are structural rather than a missing constant:

    1. **Every transcriptional module a stressor touches is INDUCED.** `gene_scales` computes
       ``scale = max(0, 1 + activity)``, so positive activity gives a scale at or above 1,
       and E-Flux applies UPPER bounds only. Raising a ceiling of 1000 to 1500 changes
       nothing about a flux of 10.
    2. **Every repressing signal sits on a metabolite pool** -- redox under DTT, nadh under
       menadione, atp under glucose starvation and antimycin A, ph under acetic acid. Those
       are not transcriptional, they carry no regulon, and `gene_scales` refuses them rather
       than pretending. So the one direction that could tighten a bound is exactly the
       direction with no route to a gene.

    This is worth a test rather than a paragraph, because "regulation is connected" reads as
    "regulation does something" and here it does not. If a repressed transcriptional module
    is ever added, or a measured reference capacity replaces the model-implied one, these
    tests fail and that is the intended signal to go and look.
    """

    def test_no_stressor_at_any_dose_represses_a_transcriptional_module(self):
        from ystwin.generator.stress_panel import (
            METABOLITE_POOLS,
            STRESSORS,
            module_response,
        )

        repressed = set()
        for stressor in STRESSORS:
            for dose in (0.1, 0.5, 1.0, 2.0, 5.0):
                try:
                    activity = module_response(stressor, dose)
                except Exception:
                    continue
                repressed |= {m for m, v in activity.items() if v < -1e-9}

        assert repressed <= METABOLITE_POOLS

    def test_and_the_pools_are_the_ones_with_no_regulon(self, carotenoid_gsmm):
        """The other half. If a pool ever gained a regulon this would stop being a reason."""
        from ystwin.bridge.regulation import load_regulons
        from ystwin.generator.stress_panel import METABOLITE_POOLS

        regulons = load_regulons()

        assert not (set(regulons) & METABOLITE_POOLS)

    def test_eflux_tightens_nothing_for_a_real_stressor(self, carotenoid_gsmm):
        from ystwin.bridge.regulation import eflux_layer
        from ystwin.generator.stress_panel import METABOLITE_POOLS, module_response

        layer = eflux_layer(carotenoid_gsmm, module_response("DTT", 1.0),
                            allow_unmapped=METABOLITE_POOLS)

        assert layer.binding == ()

    def test_it_computes_scales_all_the_same_so_this_is_not_a_wiring_fault(self, carotenoid_gsmm):
        """The distinction that makes the result meaningful. Zero bindings because the
        machinery never ran would be a bug; zero bindings out of sixty computed scales is a
        property of the model and the modules."""
        from ystwin.bridge.regulation import eflux_layer
        from ystwin.generator.stress_panel import METABOLITE_POOLS, module_response

        layer = eflux_layer(carotenoid_gsmm, module_response("DTT", 1.0),
                            allow_unmapped=METABOLITE_POOLS)

        assert len(layer.scales) > 10

    def test_the_audit_says_so_rather_than_returning_the_same_number_quietly(self, carotenoid_gsmm):
        """A stressed audit that silently equals the unstressed one is the same trap as a
        stressed prediction that silently equals the unstressed one."""
        from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID
        from ystwin.fba.audit import audit_predicted_flux
        from ystwin.generator.stress_panel import module_response

        stressed = audit_predicted_flux(
            carotenoid_gsmm, PRODUCT_DEMAND_ID, 9.46e-4, glucose_uptake=10.0,
            module_activity=module_response("DTT", 1.0))

        assert stressed.regulation == "no reaction reached"

    def test_an_unstressed_audit_records_that_it_asked_the_unstressed_question(self, carotenoid_gsmm):
        from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID
        from ystwin.fba.audit import audit_predicted_flux

        plain = audit_predicted_flux(carotenoid_gsmm, PRODUCT_DEMAND_ID, 9.46e-4, glucose_uptake=10.0)

        assert plain.regulation == "unstressed"
