"""A pathway is a file, so the file is where a wrong pathway has to be caught.

``pathway/spec.py`` moved beta-carotene out of Python and into TOML precisely so that
adding a product is an edit to ``data/pathways/`` rather than a diff against the solver.
That trade only pays if the loader refuses a spec the rest of the stack cannot honour,
because a malformed spec does not crash: it walks into :func:`solve_pathway` and comes back
out as a number.

The refusal that matters most is the secreted one. Every balance downstream closes on
``d[X]/dt = v_in - v_out - mu*[X] = 0``, and a species leaving through a transporter has no
``mu*[X]`` term at all -- the arithmetic still runs, still returns a content, and the
content is meaningless. So the message has to name the missing term rather than say "not
supported", and these tests pin that.

The rest guard the smaller assumptions the solver makes without checking: that the product
is the last node it walks, that the last node weighs something (content is reported in
mg/gDCW), that nothing consumes the product, and that a spec nobody can calibrate says so
AT LOAD rather than waiting to be asked -- `phb` sat in `data/pathways/` looking fittable
for as long as its non-calibratability lived in a property nobody knew to read.
"""

from __future__ import annotations

import dataclasses
import logging

import pytest

from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS
from ystwin.pathway.solve import NodeKinetics, solve_pathway
from ystwin.pathway.spec import (
    Fate,
    Node,
    PathwaySpec,
    RateLaw,
    available_pathways,
    load_pathway,
)

# One valid spec, so every negative case below differs from it by exactly one field.
VALID = """
[pathway]
product = "squalene"
organism = "Saccharomyces cerevisiae"
entry_enzyme = "ERG9"
precursor_metabolite = "s_0190"

[[node]]
name = "presqualene"
enzyme = "ERG9"
rate_law = "passthrough"
molar_mass_g_per_mol = 426.72

[[node]]
name = "squalene"
enzyme = ""
rate_law = "passthrough"
molar_mass_g_per_mol = 410.72
"""


def _write(directory, name, body):
    """Put one spec on disk and hand back the directory ``load_pathway`` should read."""
    (directory / f"{name}.toml").write_text(body)
    return directory


def _load(directory, name, body):
    return load_pathway(name, _write(directory, name, body))


# A chain whose middle node holds a pool nobody measures: alpha and gamma are measurable,
# beta is not. Two measurable nodes, which is what the old `calibratable` counted.
HIDDEN_POOL = PathwaySpec(
    product="gamma",
    organism="Saccharomyces cerevisiae",
    nodes=(Node("alpha", RateLaw.SATURATING, enzyme="E1", measurable=True),
           Node("beta", RateLaw.SATURATING, enzyme="E2", measurable=False),
           Node("gamma", RateLaw.PASSTHROUGH, molar_mass_g_per_mol=100.0)),
    entry_enzyme="E1",
    precursor_metabolite="s_0001")

MU = 0.1
ALPHA_CONTENT = 2.0
GAMMA_CONTENT = 1.0

# Two ``(vmax_per_growth, km)`` pairs for the hidden node, both able to carry the flux it
# must carry. The second sits just above that flux, where ``km*v/(vmax - v)`` blows up --
# which is the shape of the ridge and not a pathological choice: nothing measurable
# distinguishes them, so nothing rules it out.
TIGHT_HIDDEN_NODE = (5.0, 1e-4)
SLACK_HIDDEN_NODE = (1.001, 1.0)


def _chain_reaching(vmax_hidden, km_hidden):
    """Entry flux and kinetics putting the SAME contents on both measurable nodes.

    The hidden node absorbs the difference. For any ``vmax`` above the flux beta must
    carry, the pool ``km*v/(vmax - v)`` delivers that flux exactly, so the free pair
    ``(vmax_hidden, km_hidden)`` buys an entry flux without moving anything observable.
    """
    out_of_beta = MU * GAMMA_CONTENT
    beta_content = km_hidden * GAMMA_CONTENT / (vmax_hidden - GAMMA_CONTENT)
    out_of_alpha = out_of_beta + MU * beta_content
    km_alpha = 1e-3
    vmax_alpha = out_of_alpha * (km_alpha + ALPHA_CONTENT) / (MU * ALPHA_CONTENT)
    entry_flux = out_of_alpha + MU * ALPHA_CONTENT
    return entry_flux, {"alpha": NodeKinetics(vmax_per_growth=vmax_alpha, km=km_alpha),
                        "beta": NodeKinetics(vmax_per_growth=vmax_hidden, km=km_hidden)}


class TestASecretedProductIsRefused:
    """Split on 2026-09-01, and the split is the point rather than a rearrangement.

    A secreted TERMINAL node is one measured capacity short of a closable balance, so its
    refusal moved to `solve_pathway`, where the missing number can be named -- see
    `tests/test_secretion_outlet.py`. A secreted INTERMEDIATE is refused here still, and no
    measurement lifts it: the walk carries one flux from each node to the next, and an
    intermediate leaking sideways breaks the walk rather than a coefficient in it.

    So the two tests that used `squalene` (the product) now assert that it LOADS, and the
    one that used `presqualene` (an intermediate) is unchanged.
    """

    def test_a_secreted_terminal_node_now_loads(self, tmp_path):
        """The relocation, asserted at the place it moved from."""
        body = VALID.replace('name = "squalene"', 'name = "squalene"\nfate = "secreted"')

        spec = _load(tmp_path, "squalene", body)

        assert spec.nodes[-1].fate == Fate.SECRETED

    def test_a_secreted_intermediate_is_still_refused_at_load(self, tmp_path):
        """The half that did NOT move, and the reason is the chain's shape.

        The whole steady state here is ``v_in - v_out - mu*[X] = 0``, and the walk closes it
        one node at a time by handing ``v_out`` forward. An intermediate that also leaves
        the cell loses the ``mu*[X]`` term AND breaks that handoff: the flux arriving at the
        next node is not the flux this one passed on. Naming the term is what tells a reader
        the refusal is structural rather than a policy that could be waived with a constant.
        """
        body = VALID.replace('name = "presqualene"',
                             'name = "presqualene"\nfate = "secreted"')

        with pytest.raises(ValueError) as raised:
            _load(tmp_path, "squalene", body)

        assert "mu*[X]" in str(raised.value)

    def test_a_degraded_intermediate_is_refused_at_load_too(self, tmp_path):
        """The guard nobody had, and it was live until 2026-09-01.

        When the DEGRADED refusal moved from load time to solve time it landed inside the
        TERMINAL arm of `solve_pathway`, which narrowed it from every node to the last one.
        A chain with a degraded intermediate and no rate constant then loaded, solved, and
        returned the diluted answer with its largest outlet silently dropped -- refused at
        neither place. The fix is here rather than in the solver because an intermediate
        outlet is not a missing measurement; it is a chain this walk cannot describe.
        """
        body = VALID.replace('name = "presqualene"',
                             'name = "presqualene"\nfate = "degraded"')

        with pytest.raises(ValueError, match="presqualene"):
            _load(tmp_path, "squalene", body)

    def test_the_refusal_names_the_offending_species(self, tmp_path):
        body = VALID.replace('name = "presqualene"',
                             'name = "presqualene"\nfate = "secreted"')

        with pytest.raises(ValueError, match="presqualene"):
            _load(tmp_path, "squalene", body)

    def test_an_unknown_fate_is_refused_rather_than_treated_as_dilution(self, tmp_path):
        """``fate`` is the field that decides whether this module may model the species.

        A typo silently falling back to ``diluted`` would put a secreted product through
        the intracellular balance, which is the one outcome the module exists to prevent.
        """
        body = VALID.replace('name = "squalene"', 'name = "squalene"\nfate = "exported"')

        with pytest.raises(ValueError, match="fate must be one of"):
            _load(tmp_path, "squalene", body)


class TestTheTerminalNodeCarriesTheProduct:
    def test_a_product_that_is_not_the_last_node_is_refused(self, tmp_path):
        """The solver walks the chain in order and reports whatever comes off the end.

        If the product is declared somewhere in the middle, that walk still succeeds and
        the reported content belongs to a different species.
        """
        body = VALID.replace('product = "squalene"', 'product = "presqualene"')

        with pytest.raises(ValueError, match="must be the terminal node"):
            _load(tmp_path, "presqualene", body)

    def test_a_terminal_node_with_no_molar_mass_is_refused(self, tmp_path):
        """Content is reported in mg/gDCW, and the alternative is a constant in reporting.

        The one number needed to get there is the product's molar mass. Refusing here keeps
        it beside the chemistry instead of hard-coded wherever the unit conversion lands.
        """
        body = VALID.replace("\nmolar_mass_g_per_mol = 410.72", "")

        with pytest.raises(ValueError, match="molar_mass_g_per_mol"):
            _load(tmp_path, "squalene", body)

    def test_a_terminal_node_naming_an_enzyme_is_refused(self, tmp_path):
        """An enzyme on the last node is a chain that was truncated by accident.

        Nothing consumes the product, so the solver ignores the enzyme -- and the reader is
        left believing a step is modelled that never runs.
        """
        body = VALID.replace('name = "squalene"\nenzyme = ""',
                             'name = "squalene"\nenzyme = "ERG1"')

        with pytest.raises(ValueError, match="nothing consumes the product"):
            _load(tmp_path, "squalene", body)

    def test_a_non_positive_molar_mass_is_refused(self, tmp_path):
        body = VALID.replace("molar_mass_g_per_mol = 410.72",
                             "molar_mass_g_per_mol = 0.0")

        with pytest.raises(ValueError, match="molar mass must be positive"):
            _load(tmp_path, "squalene", body)


class TestAChainMustBeAChain:
    def test_a_repeated_node_name_is_refused(self, tmp_path):
        """``PathwaySolution.node`` returns the first match, so a duplicate hides a pool.

        Two nodes sharing a name also make ``kinetics`` ambiguous: one entry would silently
        parameterise both steps.
        """
        body = VALID.replace('name = "presqualene"', 'name = "squalene"')

        with pytest.raises(ValueError, match="repeats a node name"):
            _load(tmp_path, "squalene", body)

    def test_a_spec_with_no_nodes_is_refused_by_the_loader(self, tmp_path):
        body = VALID.split("[[node]]")[0]

        with pytest.raises(ValueError, match=r"declares no \[\[node\]\] entries"):
            _load(tmp_path, "squalene", body)

    def test_an_empty_node_tuple_is_refused_by_the_spec_itself(self):
        with pytest.raises(ValueError, match="declares no nodes"):
            PathwaySpec(product="squalene", organism="Saccharomyces cerevisiae",
                        nodes=(), entry_enzyme="ERG9", precursor_metabolite="s_0190")

    def test_an_unknown_rate_law_is_refused(self, tmp_path):
        """The solver has one branch per rate law and no default.

        A law it does not know would either fall through to the unreachable arm or, worse,
        be read as the saturating default and quietly acquire two fitted constants.
        """
        body = VALID.replace('rate_law = "passthrough"', 'rate_law = "hill"', 1)

        with pytest.raises(ValueError, match="rate_law must be one of"):
            _load(tmp_path, "squalene", body)

    def test_a_non_positive_precursor_stoichiometry_is_refused(self, tmp_path):
        body = VALID.replace('precursor_metabolite = "s_0190"',
                             'precursor_metabolite = "s_0190"\nprecursor_stoichiometry = 0.0')

        with pytest.raises(ValueError, match="precursor_stoichiometry must be positive"):
            _load(tmp_path, "squalene", body)


class TestCountingMeasurableNodesIsTheWrongTest:
    """The predicate `calibratable` used to be ``>= 2 measurable nodes``, and a count is
    wrong in both directions.

    The entry flux is ``mu`` times the sum of every pool in the chain, so what decides
    whether it can be recovered from measurements is WHICH nodes hold a pool -- a property
    of ``rate_law`` -- and not how many nodes happen to be quantified. A count says yes to
    a chain hiding a pool between two measured ones, and no to a chain that holds no pool
    at all and whose entry flux is therefore exact.
    """

    def test_two_measurable_nodes_with_a_pool_hidden_between_them_are_not_enough(self):
        assert sum(1 for n in HIDDEN_POOL.nodes if n.measurable) == 2  # what was counted

        assert HIDDEN_POOL.calibratable is False

    def test_the_hidden_pool_buys_an_entry_flux_three_hundred_fold_wide(self):
        """The reason the count is not merely imprecise. This is the ridge, run.

        Both parameter sets put 2.0 and 1.0 on the two MEASURED nodes -- identical to nine
        decimals -- while implying entry fluxes 334x apart. Nothing measurable moves,
        so no fit can choose between them, and a `calibratable` that returned True here
        would licence reporting whichever one the optimiser started nearest.
        ``docs/research/KINETIC_FIT.md`` section 3 reports the same ridge for the
        desaturase step, over three orders of magnitude in vmax.
        """
        low_flux, _ = _chain_reaching(*TIGHT_HIDDEN_NODE)
        high_flux, _ = _chain_reaching(*SLACK_HIDDEN_NODE)

        assert high_flux / low_flux > 300.0

    def test_and_neither_of_those_entry_fluxes_moves_a_measured_content(self):
        low_flux, low_kinetics = _chain_reaching(*TIGHT_HIDDEN_NODE)
        high_flux, high_kinetics = _chain_reaching(*SLACK_HIDDEN_NODE)
        low = solve_pathway(HIDDEN_POOL, low_flux, MU, low_kinetics)
        high = solve_pathway(HIDDEN_POOL, high_flux, MU, high_kinetics)

        for measured in ("alpha", "gamma"):
            assert low.node(measured).content_mmol_per_gdcw == pytest.approx(
                high.node(measured).content_mmol_per_gdcw, rel=1e-9)

    def test_the_verdict_names_the_node_that_would_lift_it(self):
        """A refusal to fit is useless without the one measurement that would end it."""
        verdict = HIDDEN_POOL.calibratability

        assert verdict.unmeasured_pools == ("beta",)

    def test_a_hidden_pool_is_an_identifiability_failure_and_says_which_kind(self):
        verdict = HIDDEN_POOL.calibratability

        assert verdict.flux_is_recoverable is False

    def test_beta_carotene_stops_being_calibratable_if_phytoene_holds_a_pool(self):
        """The declaration that actually carries the identification, made load-bearing.

        `beta_carotene.toml` declares phytoene ``passthrough`` and its own notes call that
        known-wrong -- Chen 2016 measures phytoene at 3.99% of total carotenoid. Flip the
        rate law to what that implies and the entry flux stops being ``mu*([lyc] + [bcar])``
        and acquires the desaturase ridge, while the count of measurable nodes does not
        move at all. A predicate that cannot see this cannot protect the fit that rests on
        it.
        """
        shipped = load_pathway("beta_carotene")
        accumulating = dataclasses.replace(
            shipped,
            nodes=(dataclasses.replace(shipped.nodes[0], rate_law=RateLaw.SATURATING),
                   *shipped.nodes[1:]))
        assert sum(1 for n in accumulating.nodes if n.measurable) == 2  # unchanged

        assert accumulating.calibratable is False


class TestASpecWithNoPoolButTheProductIsNotCalibratableEither:
    """The other direction the count got wrong, and `phb` is the case.

    `phb` is passthrough end to end, so its only pool is the product and the entry flux is
    ``mu*[PHB]`` with no fitted parameter in it. The old predicate called it
    non-calibratable for a stated reason -- "the entry flux and the branch capacity trade
    off along a ridge" -- that this chain structurally cannot have, because it has no
    branch capacity. Right answer, wrong mechanism, and the mechanism is what a reader
    acts on: it sends somebody looking for a second node when what is missing is an
    entry-enzyme expression series.
    """

    def test_a_chain_holding_no_pool_but_the_product_is_not_calibratable(self, tmp_path):
        spec = _load(tmp_path, "squalene", VALID)

        assert spec.calibratable is False

    def test_but_its_entry_flux_is_recoverable_exactly(self, tmp_path):
        """Which is the opposite of a ridge, and is why the old reason was wrong."""
        spec = _load(tmp_path, "squalene", VALID)

        assert spec.calibratability.flux_is_recoverable is True

    def test_the_solver_returns_its_own_input_on_such_a_chain(self, tmp_path):
        """What "nothing can contradict it" means, run rather than argued.

        Every node above the product is passthrough, so ``solve_pathway`` needs no kinetics
        and returns ``content = v_in/mu``. No parameter exists to be wrong, so no measured
        content can disagree with any fitted scalar -- which is why the flux law cannot be
        SCORED here even though it can be evaluated.
        """
        spec = _load(tmp_path, "squalene", VALID)
        solved = solve_pathway(spec, entry_flux=7e-4, growth_rate=0.18, kinetics={})
        content = solved.terminal.content_mmol_per_gdcw

        assert content == pytest.approx(7e-4 / 0.18, rel=1e-12)

    def test_the_shipped_phb_spec_is_the_real_case(self):
        assert load_pathway("phb").calibratable is False

    def test_and_its_reason_does_not_blame_a_ridge_it_cannot_have(self):
        """The finding this class exists for.

        `phb`'s chain contains no saturating node, so there is no branch capacity for an
        entry flux to trade off against. A message claiming otherwise sends the reader to
        measure a CoA thioester, which would change nothing.
        """
        reason = load_pathway("phb").calibratability.reason

        assert "no branch capacity" in reason

    def test_and_it_names_the_entry_enzyme_whose_expression_series_is_missing(self):
        """What would actually lift it. Kocharin 2012's four strains carry one plasmid, so
        relative PhaA expression is 1.0 by construction and there is no axis to regress on
        -- a dataset property, which is why no spec can carry it."""
        reason = load_pathway("phb").calibratability.reason

        assert "PhaA" in reason


class TestASpecSaysWhetherItCanBeCalibrated:
    def test_the_shipped_beta_carotene_spec_is_calibratable(self):
        """Lycopene holds a saturating pool and beta-carotene is the product; both come off
        one HPLC run, so the entry flux is ``mu*([lycopene] + [beta_carotene])`` with no
        fitted constant in it. Losing the lycopene measurement would take the fit below
        identifiability without changing a line of Python."""
        spec = load_pathway("beta_carotene")

        assert spec.calibratable is True

    def test_a_passthrough_node_is_not_counted_as_a_pool(self):
        """It is the one rate law whose steady state is exactly zero, so it contributes
        nothing to ``mu * sum(pools)`` and needs no measurement to be pinned."""
        spec = load_pathway("beta_carotene")

        assert [n.name for n in spec.pool_bearing_nodes] == ["lycopene", "beta_carotene"]

    def test_the_terminal_node_always_holds_a_pool_whatever_it_declares(self):
        """`beta_carotene` declares the product ``passthrough``, and the solver ignores that
        and gives it ``content = flux_in/mu`` regardless. If the two disagreed, a spec could
        declare its product away and the verdict would call an unanchored chain fitted."""
        spec = load_pathway("beta_carotene")
        assert spec.nodes[-1].rate_law == RateLaw.PASSTHROUGH  # the declaration ignored

        assert spec.pool_bearing_nodes[-1] is spec.nodes[-1]

    def test_a_declared_zero_nobody_can_check_is_recorded_rather_than_forgotten(self):
        """Phytoene is ``passthrough`` AND unmeasurable, so the whole identification rests
        on a pool nothing can ever weigh. That is not a refusal -- the alternative,
        inventing CrtI kinetics, is worse -- but it must not vanish into a True."""
        verdict = load_pathway("beta_carotene").calibratability

        assert verdict.unfalsifiable_zeros == ("phytoene",)

    def test_and_the_calibratable_verdict_carries_that_caveat_in_its_reason(self):
        verdict = load_pathway("beta_carotene").calibratability
        assert verdict.calibratable is True  # the caveat rides on a clean verdict

        assert "phytoene" in verdict.reason


class TestTheVerdictIsAnnouncedAtLoadAndIsNotARefusal:
    """`calibratable` was discoverable only by asking, and nobody asks a property they do
    not know exists.

    Loading `phb` printed nothing, said nothing and returned a spec indistinguishable from
    a fittable one. The record is a log line rather than a `warnings.warn` for a reason
    that belongs to this repository: `pyproject.toml` sets ``filterwarnings = ["error"]``,
    so a warning here would raise in every test that loads `phb` -- which is exactly the
    refusal these tests pin that it must not be.
    """

    def test_loading_a_non_calibratable_spec_records_the_reason(self, caplog):
        with caplog.at_level(logging.WARNING, logger="ystwin.pathway.spec"):
            load_pathway("phb")

        assert "CANNOT BE CALIBRATED" in caplog.text

    def test_the_record_names_the_file_so_a_reader_knows_which_spec(self, caplog):
        with caplog.at_level(logging.WARNING, logger="ystwin.pathway.spec"):
            load_pathway("phb")

        assert "phb.toml" in caplog.text

    def test_the_record_says_what_would_lift_it(self, caplog):
        """A note that says only "cannot be calibrated" gets the flag deleted rather than
        the measurement taken."""
        with caplog.at_level(logging.WARNING, logger="ystwin.pathway.spec"):
            load_pathway("phb")

        assert "What lifts it" in caplog.text

    def test_a_calibratable_spec_loads_silently(self, caplog):
        """A note on every load is a note nobody reads. `beta_carotene` is loaded by the
        prediction chain, four test modules and two scripts."""
        with caplog.at_level(logging.WARNING, logger="ystwin.pathway.spec"):
            load_pathway("beta_carotene")

        assert caplog.text == ""

    def test_a_hidden_pool_is_recorded_at_load_too(self, tmp_path, caplog):
        body = VALID.replace('name = "presqualene"\nenzyme = "ERG9"\n'
                             'rate_law = "passthrough"',
                             'name = "presqualene"\nenzyme = "ERG9"\n'
                             'rate_law = "saturating"\nmeasurable = false')

        with caplog.at_level(logging.WARNING, logger="ystwin.pathway.spec"):
            _load(tmp_path, "squalene", body)

        assert "presqualene" in caplog.text

    def test_a_non_calibratable_spec_still_loads(self):
        """The capability the refusal would have deleted. `phb` is the only second product
        in the tree and the only one with an environment series, so it exercises the
        solver, the growth-rate guard and the environment layer on a second paper."""
        spec = load_pathway("phb")

        assert spec.product == "phb"

    def test_and_still_solves(self):
        """A record that stopped the spec working would be a refusal wearing a note."""
        solved = solve_pathway(load_pathway("phb"), entry_flux=1e-3, growth_rate=0.1)

        assert solved.terminal.content_mmol_per_gdcw > 0

    def test_require_calibratable_is_the_gate_a_fitting_path_calls(self):
        """Separate from loading, because loading must not refuse and fitting must.

        Nothing calls it yet: `scripts/fit_pathway_flux.py` reads `data/carotenoid/*.tsv`
        straight into `fit_flux_law` and never loads a spec at all, which is why the
        load-time record is the only thing between a non-calibratable spec and a fit.
        """
        with pytest.raises(ValueError, match="CANNOT BE CALIBRATED"):
            load_pathway("phb").require_calibratable()

    def test_and_it_passes_on_a_pathway_that_can_be_fitted(self):
        assert load_pathway("beta_carotene").require_calibratable() is None


class TestTheLoaderFailsByName:
    def test_a_missing_spec_names_what_is_available(self, tmp_path):
        """A typo in a product name must not read as an unsupported product.

        The two failures need different responses -- write the TOML, or fix the string --
        and only the message distinguishes them.
        """
        _write(tmp_path, "squalene", VALID)

        with pytest.raises(FileNotFoundError) as raised:
            load_pathway("farnesene", tmp_path)

        assert "squalene" in str(raised.value)

    def test_a_missing_spec_raises_file_not_found_and_not_value_error(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_pathway("farnesene", tmp_path)

    def test_an_empty_directory_still_reports_the_available_set(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="none"):
            load_pathway("farnesene", tmp_path)

    def test_a_file_with_no_pathway_table_is_refused(self, tmp_path):
        body = VALID.replace("[pathway]", "[header]")

        with pytest.raises(ValueError, match=r"no \[pathway\] table"):
            _load(tmp_path, "squalene", body)

    def test_an_unknown_node_field_is_refused_rather_than_ignored(self, tmp_path):
        """A field the loader drops is a claim the reader believes and the solver ignores.

        ``Node(**entry)`` is what enforces this, so a renamed field fails at load rather
        than reverting to the dataclass default several layers down.
        """
        body = VALID + "\nkm = 1e-4\n"

        with pytest.raises(TypeError):
            _load(tmp_path, "squalene", body)


class TestTheShippedSpecMatchesWhatTheCalibrationsExpect:
    def test_beta_carotene_ships_and_is_discoverable(self):
        assert "beta_carotene" in available_pathways()

    def test_the_entry_enzyme_is_the_gene_the_flux_law_was_fitted_on(self):
        """A calibration applied to a different gene's numbers is silently wrong.

        ``FluxCalibration`` carries its gene for exactly this reason; the spec names the
        same gene, and nothing else checks that the two agree.
        """
        spec = load_pathway("beta_carotene")

        assert spec.entry_enzyme == BETA_CAROTENE_FLUX.entry_enzyme

    def test_every_calibrated_node_exists_in_the_spec(self):
        """``solve_pathway`` looks kinetics up by node name and defaults to empty.

        A rename on either side would not raise: the saturating node would ask for
        parameters that are no longer keyed to it, and the failure would surface as a
        refusal about missing vmax rather than about a renamed node.
        """
        spec = load_pathway("beta_carotene")

        assert all(spec.node(name) for name in BETA_CAROTENE_KINETICS)

    def test_the_calibrated_node_is_the_one_declared_saturating(self):
        spec = load_pathway("beta_carotene")

        assert all(spec.node(name).rate_law == RateLaw.SATURATING
                   for name in BETA_CAROTENE_KINETICS)

    def test_the_product_is_the_terminal_node(self):
        spec = load_pathway("beta_carotene")

        assert spec.nodes[-1].name == spec.product

    def test_the_product_carries_its_molar_mass(self):
        spec = load_pathway("beta_carotene")

        assert spec.nodes[-1].molar_mass_g_per_mol == pytest.approx(536.87)

    def test_no_node_of_the_shipped_spec_is_secreted(self):
        spec = load_pathway("beta_carotene")

        assert all(node.fate == Fate.DILUTED for node in spec.nodes)

    def test_the_precursor_metabolite_is_declared_for_the_audit_layer(self):
        """``audit_predicted_flux`` reads the precursor budget off this id.

        It is a Yeast9 metabolite id, not a name, and an empty string here would send the
        audit layer looking up ``''`` in the model.
        """
        spec = load_pathway("beta_carotene")

        assert spec.precursor_metabolite == "s_0189"

    def test_two_ggpp_are_declared_per_phytoene(self):
        spec = load_pathway("beta_carotene")

        assert spec.precursor_stoichiometry == 2.0

    def test_phytoene_is_declared_unmeasurable_so_the_assumption_stays_visible(self):
        """Verwaal 2007 reports phytoene accumulating in this pathway in S. cerevisiae.

        The spec declares it passthrough anyway, which is an assumption rather than a
        measurement, and ``measurable = false`` is the record that no phytoene number
        exists for these strains. If that flag flips without a measurement arriving, the
        branch would look identifiable when it is not.
        """
        spec = load_pathway("beta_carotene")

        assert spec.node("phytoene").measurable is False

    def test_an_unknown_node_name_fails_by_name(self):
        with pytest.raises(KeyError, match="phytofluene"):
            load_pathway("beta_carotene").node("phytofluene")


class TestADegradedProductIsUnderParameterisedNotOutOfScope:
    """The second fate the balance cannot represent *without one more constant*.

    The solver used to give a terminal node exactly two outlets -- the next step and growth
    dilution -- so every native storage compound was refused at load: glycogen has GPH1 and
    SGA1, trehalose has NTH1/NTH2/ATH1, triacylglycerol has TGL3/4/5, and declaring one of
    those `diluted` returns ``content = flux/mu`` while omitting the term doing most of the
    work. How much work: Ferreira 2018 raises triacylglycerol from 129 to 218 mg/gDCW by
    DELETING the lipases alone.

    **The refusal moved on 2026-09-01; it did not go away.** `solve_pathway` carries a third
    outlet now, so the missing thing is a rate constant rather than a whole class of
    chemistry -- and a rate constant is a fitted number, which belongs in `kinetics` beside
    a citation and not in a spec. These tests assert the same guarantee at its new location:
    no degradation rate, no number.
    """

    STORAGE = ('[pathway]\nproduct = "storage"\norganism = "S. cerevisiae"\n'
               'entry_enzyme = "E1"\nprecursor_metabolite = "s_0001"\n\n'
               '[[node]]\nname = "storage"\nenzyme = ""\nrate_law = "passthrough"\n'
               'fate = "degraded"\nmolar_mass_g_per_mol = 100.0\n')

    def _spec(self, tmp_path):
        (tmp_path / "storage.toml").write_text(self.STORAGE)
        return load_pathway("storage", tmp_path)

    def test_a_degraded_fate_now_loads(self, tmp_path):
        """Being one constant short is not the same as being outside the arithmetic."""
        spec = self._spec(tmp_path)

        assert spec.product == "storage"
        assert spec.nodes[-1].fate == "degraded"

    def test_solving_without_a_degradation_rate_is_refused(self, tmp_path):
        """The guarantee the load-time refusal protected: no silent `flux/mu`."""
        from ystwin.pathway.solve import solve_pathway

        with pytest.raises(ValueError, match="degradation_rate_per_h"):
            solve_pathway(self._spec(tmp_path), 0.5, 0.1, {})

    def test_the_refusal_still_says_how_large_the_omitted_term_can_be(self, tmp_path):
        """A refusal without a magnitude reads as pedantry. 129 to 218 mg/gDCW is a
        measurement."""
        from ystwin.pathway.solve import solve_pathway

        with pytest.raises(ValueError, match="129 to 218"):
            solve_pathway(self._spec(tmp_path), 0.5, 0.1, {})

    def test_the_refusal_names_what_would_lift_it(self, tmp_path):
        from ystwin.pathway.solve import solve_pathway

        with pytest.raises(ValueError, match="measured on"):
            solve_pathway(self._spec(tmp_path), 0.5, 0.1, {})

    def test_a_supplied_rate_changes_the_answer_and_by_how_much(self, tmp_path):
        """The whole reason the outlet exists. At mu = 0.1 and k_deg = 0.4 the pool is a
        fifth of what growth dilution alone would predict."""
        from ystwin.pathway.solve import NodeKinetics, solve_pathway

        solved = solve_pathway(self._spec(tmp_path), 0.5, 0.1,
                               {"storage": NodeKinetics(degradation_rate_per_h=0.4)})

        assert solved.terminal.content_mmol_per_gdcw == pytest.approx(0.5 / 0.5)
        assert solved.terminal.content_mmol_per_gdcw < 0.5 / 0.1

    def test_a_diluted_node_is_untouched_by_the_new_outlet(self, tmp_path):
        """Every existing result must be bit-for-bit what it was: the degradation term is
        exactly zero unless a node declares the fate."""
        from ystwin.pathway.solve import solve_pathway

        (tmp_path / "plain.toml").write_text(self.STORAGE
                                             .replace("storage", "plain")
                                             .replace('fate = "degraded"', 'fate = "diluted"'))
        solved = solve_pathway(load_pathway("plain", tmp_path), 0.5, 0.1, {})

        assert solved.terminal.content_mmol_per_gdcw == 0.5 / 0.1

    def test_the_shipped_glycogen_spec_loads_and_refuses_at_solve(self):
        """`data/pathways/glycogen.toml` documents its own GPH1 and SGA1 reaction ids. It is
        the case that prompted the fate, and it is still not answerable -- nobody has
        measured that turnover in a producing strain."""
        from ystwin.pathway.solve import solve_pathway

        spec = load_pathway("glycogen")

        assert spec.nodes[-1].fate == "degraded"
        with pytest.raises(ValueError, match="degradation_rate_per_h"):
            solve_pathway(spec, 0.5, 0.1, {})


class TestTerminalLossMustBeMeasuredToRecoverEntryFlux:
    @staticmethod
    def _spec(fate, terminal_law=RateLaw.PASSTHROUGH):
        return PathwaySpec(
            product="product_pool", organism="S. cerevisiae",
            nodes=(Node("intermediate", RateLaw.PROPORTIONAL, enzyme="E1"),
                   Node("product_pool", terminal_law, fate=fate,
                        molar_mass_g_per_mol=100.0)),
            entry_enzyme="E0", precursor_metabolite="precursor")

    def test_glycogen_gross_synthesis_includes_the_degradation_loss(self):
        spec = load_pathway("glycogen")
        mu, k_deg, entry_flux = 0.025, 0.1, 1e-3
        solved = solve_pathway(
            spec, entry_flux, mu, {"glycogen": NodeKinetics(degradation_rate_per_h=k_deg)})
        content = solved.terminal.content_mmol_per_gdcw
        loss = solved.terminal.loss_flux_mmol_per_gdcw_h

        assert spec.nodes[-1].fate == Fate.DEGRADED
        assert content == pytest.approx(entry_flux / (mu + k_deg), rel=1e-12)
        assert loss == pytest.approx(k_deg * content, rel=1e-12)
        assert mu * content + loss == pytest.approx(entry_flux, rel=1e-12)
        assert entry_flux / (mu * content) == pytest.approx(5.0, rel=1e-12)
        verdict = spec.calibratability
        assert not verdict.flux_is_recoverable
        assert verdict.required_loss_parameters == ("degradation_rate_per_h",)
        assert "mu*sum(pools) + degradation_rate_per_h*[glycogen]" in verdict.reason
        assert "mu*[glycogen] and there is no branch capacity" not in verdict.reason

    def test_the_same_glycogen_pool_does_not_identify_its_gross_synthesis(self):
        spec = load_pathway("glycogen")
        mu, measured_pool = 0.025, 0.01
        synthesis = []
        for k_deg in (0.025, 0.225):
            entry_flux = (mu + k_deg) * measured_pool
            solved = solve_pathway(
                spec, entry_flux, mu,
                {"glycogen": NodeKinetics(degradation_rate_per_h=k_deg)})
            assert solved.terminal.content_mmol_per_gdcw == pytest.approx(
                measured_pool, rel=1e-12)
            synthesis.append(entry_flux)

        assert synthesis[1] / synthesis[0] == pytest.approx(5.0, rel=1e-12)
        with pytest.raises(ValueError, match="degradation_rate_per_h"):
            spec.require_calibratable()

    def test_secretion_is_also_a_loss_beyond_the_measured_pool_sum(self):
        spec = self._spec(Fate.SECRETED)
        mu, entry_flux = 0.1, 1e-3
        vmax, km = 2e-3, 5e-3
        solved = solve_pathway(spec, entry_flux, mu, {
            "intermediate": NodeKinetics(rate_constant=0.4),
            "product_pool": NodeKinetics(
                secretion_vmax_mmol_per_gdcw_h=vmax, secretion_km_mmol_per_gdcw=km,
                growth_rate_range=(mu, mu)),
        })
        content = solved.terminal.content_mmol_per_gdcw
        diluted = mu * sum(node.content_mmol_per_gdcw for node in solved.nodes)
        loss = solved.terminal.loss_flux_mmol_per_gdcw_h

        assert loss == pytest.approx(vmax * content / (km + content), rel=1e-12)
        assert diluted + loss == pytest.approx(entry_flux, rel=1e-12)
        assert diluted < entry_flux
        assert not spec.calibratability.flux_is_recoverable
        assert "secretion_vmax_mmol_per_gdcw_h" in spec.calibratability.reason
        assert "secretion_km_mmol_per_gdcw" in spec.calibratability.reason
        assert "growth_rate_range" in spec.calibratability.reason

    @pytest.mark.parametrize("terminal_law", RateLaw.ALL)
    @pytest.mark.parametrize("fate,required_parameters", [
        (Fate.DILUTED, ()),
        (Fate.DEGRADED, ("degradation_rate_per_h",)),
        (Fate.SECRETED, ("secretion_vmax_mmol_per_gdcw_h", "secretion_km_mmol_per_gdcw",
                         "growth_rate_range")),
    ])
    def test_recovery_follows_the_fate_not_the_terminal_rate_law_or_product_name(
            self, fate, required_parameters, terminal_law):
        spec = self._spec(fate, terminal_law)
        verdict = spec.calibratability

        assert verdict.has_a_fitted_pool
        assert verdict.flux_is_recoverable is (fate == Fate.DILUTED)
        assert verdict.calibratable is (fate == Fate.DILUTED)
        assert verdict.required_loss_parameters == required_parameters
        if fate != Fate.DILUTED:
            with pytest.raises(ValueError, match=fate.upper()):
                spec.require_calibratable()

    @pytest.mark.parametrize("fate", [Fate.DEGRADED, Fate.SECRETED])
    def test_hidden_pools_do_not_hide_the_additional_loss_requirement(self, fate):
        spec = self._spec(fate)
        spec = dataclasses.replace(
            spec, nodes=tuple(dataclasses.replace(node, measurable=False)
                              for node in spec.nodes))
        verdict = spec.calibratability

        assert verdict.unmeasured_pools == ("intermediate", "product_pool")
        assert not verdict.flux_is_recoverable
        assert verdict.required_loss_parameters
        assert all(parameter in verdict.reason for parameter in verdict.required_loss_parameters)
        assert fate.upper() in verdict.reason
        assert "no kinetic parameter" not in verdict.reason
        assert "content = v_in/mu whatever" not in verdict.reason
