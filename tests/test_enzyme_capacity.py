"""`vmax = kcat * [E]` audits the fitted ceiling and may never replace it.

The fitted `vmax_per_growth` is scored leave-one-strain-out against six measured states. A
derived capacity is scored against nothing, so the whole point of this layer is that it
CONTRADICTS without being able to overrule. These tests pin both halves: that it can
disagree loudly, and that disagreeing changes no number.

The interesting case is the shipped one. At the slowest turnover in the host model's own
isoprenoid branch the derived capacity sits far above the fitted ceiling, in the same
direction and the same order as Arhar 2024's measurement -- which is the point of building
it, because that refutation had to come from a paper and this reaches it from the enzyme.
"""

from __future__ import annotations

import functools
import gzip
import math
from xml.etree import ElementTree

import pytest

from ystwin import paths
from ystwin.pathway.enzyme_capacity import (
    ISOPRENOID_BRANCH_GENES,
    ISOPRENOID_KCAT_RANGE_PER_S,
    EnzymeCapacityUnmeasured,
    derived_capacity,
)
from ystwin.pathway.proteome import enzyme_content_from_mass_fraction

CRTYB_MW = 74736.0
FITTED_CAPACITY = 2.3252e-3        # mmol/gDCW, the shipped ceiling
MU = 0.101
ERG1_KCAT = 0.076                  # 1/s, squalene epoxidase, the neighbour nearest CrtI
#: Arhar 2024's measured 79 mg/gDCW as a capacity, at beta-carotene's 536.87 g/mol. The
#: measurement that refuted the fitted ceiling by 63x, in the units this module works in.
ARHAR_CAPACITY = 79.0 / 536.87


def _content(fraction):
    """A mass fraction of total protein as a content, mmol/gDCW.

    Named rather than inlined because this module no longer accepts a fraction: the
    convention now lives in the function name, in `pathway/proteome.py`.
    """
    return enzyme_content_from_mass_fraction(fraction, CRTYB_MW)


#: The second axis of the module's claim, declared here because the claim is quantified over
#: it. A heterologous cassette enzyme somewhere between a weak single integration and a strong
#: high-copy one, as a MASS fraction of total protein. Unmeasured for every strain in this
#: repository -- that is the whole reason the layer refuses.
MASS_FRACTION_WINDOW = (0.001, 0.05)

#: Every corner of the two-parameter box, plus the interior point the old single-corner tests
#: used, so a claim can never again be written from one of them. `_derived` HAS NO DEFAULT
#: FRACTION on purpose: until 2026-09-04 it defaulted to 0.005 and the module docstring
#: quantified over kcat alone while the dominant axis sat pinned at that value.
BOX = [(kcat, fraction)
       for kcat in ISOPRENOID_KCAT_RANGE_PER_S
       for fraction in (*MASS_FRACTION_WINDOW, 0.005)]


def _derived(kcat, fraction):
    return derived_capacity("crtYB", kcat, _content(fraction))


_SBML = "{http://www.sbml.org/sbml/level3/version1/core}"
_FBC = "{http://www.sbml.org/sbml/level3/version1/fbc/version2}"

_needs_ecmodel = pytest.mark.skipif(
    paths.ec_yeast_gem() is None or not paths.ec_yeast_gem().is_file(),
    reason="ecYeastGEM not vendored")


@functools.lru_cache(maxsize=1)
def _ecmodel():
    """The vendored GECKO model's `<model>` element.

    Parsed with the standard library rather than through cobrapy: what is read here is the
    RAW protein stoichiometry, which is where GECKO stores 1/kcat in hours, and a solver
    wrapper would hide it.
    """
    with gzip.open(paths.ec_yeast_gem(), "rb") as handle:
        return ElementTree.parse(handle).getroot().find(_SBML + "model")


def _gene_labels(model):
    return {gp.get(_FBC + "id"): gp.get(_FBC + "label")
            for gp in model.iter(_FBC + "geneProduct")}


def _protein_coefficients(reaction):
    """The `prot_*` reactant stoichiometries of one reaction, excluding the pool."""
    reactants = reaction.find(_SBML + "listOfReactants")
    if reactants is None:
        return {}
    return {ref.get("species"): float(ref.get("stoichiometry"))
            for ref in reactants
            if ref.get("species").startswith("prot_")
            and not ref.get("species").startswith("prot_pool")}


def _is_forward_arm(reaction_id):
    """The direction rule declared beside `ISOPRENOID_BRANCH_GENES`, as one predicate."""
    return not (reaction_id.startswith("draw_prot_")
                or reaction_id.startswith("arm_")
                or "_REV" in reaction_id)


@functools.lru_cache(maxsize=1)
def _branch_kcats_per_s():
    """`{gene: {reaction id: kcat /s}}` over `ISOPRENOID_BRANCH_GENES`, read from the file.

    GECKO encodes a turnover as the protein species' stoichiometric coefficient, which is
    `1/kcat` in HOURS -- so the conversion is `1/(coeff * 3600)`.
    """
    model = _ecmodel()
    labels = _gene_labels(model)
    declared = set(ISOPRENOID_BRANCH_GENES)
    found = {}
    for reaction in model.find(_SBML + "listOfReactions"):
        reaction_id = reaction.get("id")
        if not _is_forward_arm(reaction_id):
            continue
        association = reaction.find(_FBC + "geneProductAssociation")
        if association is None:
            continue
        genes = {labels.get(ref.get(_FBC + "geneProduct"))
                 for ref in association.iter(_FBC + "geneProductRef")} & declared
        if not genes:
            continue
        for coefficient in _protein_coefficients(reaction).values():
            for gene in genes:
                found.setdefault(gene, {})[reaction_id] = 1.0 / (coefficient * 3600.0)
    return found


@functools.lru_cache(maxsize=1)
def _gpr_kcat_coverage():
    """`(reactions carrying a turnover, GPR reactions)` over the whole vendored model.

    Both classes of GECKO pseudo-reaction are out of the denominator: `draw_prot_*` draws
    protein from the pool and `arm_*` pools isoenzymes, and neither carries a turnover by
    construction. Including `arm_*` gives 80.7% and is the error to avoid here.
    """
    model = _ecmodel()
    total = with_kcat = 0
    for reaction in model.find(_SBML + "listOfReactions"):
        reaction_id = reaction.get("id")
        if reaction_id.startswith("draw_prot_") or reaction_id.startswith("arm_"):
            continue
        if reaction.find(_FBC + "geneProductAssociation") is None:
            continue
        total += 1
        with_kcat += bool(_protein_coefficients(reaction))
    return with_kcat, total


def _sig3(value):
    """Three significant figures, which is the precision the shipped constant carries."""
    return float(f"{value:.3g}")


class TestItRefusesRatherThanDefaulting:
    """A default here would BE the result, so there are none."""

    @pytest.mark.parametrize("field", ["kcat_per_s", "enzyme_mmol_per_gdcw"])
    def test_every_required_input_refuses_by_name(self, field):
        kwargs = {"enzyme": "crtYB", "kcat_per_s": ERG1_KCAT,
                  "enzyme_mmol_per_gdcw": _content(0.005)}
        kwargs[field] = None

        with pytest.raises(EnzymeCapacityUnmeasured, match=field):
            derived_capacity(**kwargs)

    def test_the_refusal_names_the_experiment_and_not_just_the_field(self):
        """A refusal naming a missing number without naming how to get it sends a reader
        nowhere. Proteome fraction comes from one heavy peptide per enzyme."""
        with pytest.raises(EnzymeCapacityUnmeasured) as raised:
            derived_capacity("crtYB", None, None)

        assert "targeted proteomics" in str(raised.value)

    @pytest.mark.parametrize("bad", [0.0, -1.0, math.inf, math.nan])
    def test_a_non_finite_or_non_positive_input_is_refused(self, bad):
        with pytest.raises(ValueError, match="finite and positive"):
            _derived(kcat=bad, fraction=0.005)

    def test_a_percentage_mistaken_for_a_fraction_is_refused(self):
        """`0.5` meaning half the proteome is almost certainly `0.5%` mistyped, and the
        error is a hundredfold in the dominant term.

        The guard moved to `pathway/proteome.py` on 2026-09-04 along with the conversion
        itself -- it belongs where a FRACTION is a fraction, not where a content is."""
        with pytest.raises(ValueError, match="FRACTION"):
            _derived(kcat=ERG1_KCAT, fraction=5.0)


class TestTheArithmetic:
    def test_vmax_is_kcat_times_enzyme_content(self):
        d = _derived(ERG1_KCAT, 0.005)
        expected_mmol = 0.49974 * 0.005 / CRTYB_MW * 1000.0

        assert d.enzyme_mmol_per_gdcw == pytest.approx(expected_mmol, rel=1e-12)
        assert d.vmax_mmol_per_gdcw_h == pytest.approx(
            ERG1_KCAT * 3600.0 * expected_mmol, rel=1e-12)

    def test_capacity_divides_by_growth_because_the_fitted_one_is_per_unit_growth(self):
        d = _derived(ERG1_KCAT, 0.005)

        assert d.capacity_mmol_per_gdcw(MU) == pytest.approx(
            d.vmax_mmol_per_gdcw_h / MU, rel=1e-12)

    def test_a_non_positive_growth_rate_is_refused(self):
        with pytest.raises(ValueError, match="growth_rate"):
            _derived(ERG1_KCAT, 0.005).capacity_mmol_per_gdcw(0.0)

    def test_the_content_scales_linearly_in_both_inputs(self):
        """Neither input is privileged: doubling kcat and doubling the enzyme do the same
        thing, which is why the proteome fraction being unmeasured is as bad as kcat being
        unmeasured."""
        base = _derived(ERG1_KCAT, 0.005).vmax_mmol_per_gdcw_h

        assert (_derived(kcat=2 * ERG1_KCAT, fraction=0.005).vmax_mmol_per_gdcw_h
                == pytest.approx(2 * base))
        assert (_derived(kcat=ERG1_KCAT, fraction=0.010).vmax_mmol_per_gdcw_h
                == pytest.approx(2 * base))


class TestItDisagreesWithTheFitInTheDirectionArharMeasured:
    """The reason to build it. This is a claim about the SHIPPED calibration, so if the
    calibration is ever refitted these numbers move and the test should be re-derived rather
    than relaxed."""

    def test_erg1_the_nearest_neighbour_overshoots_the_fitted_ceiling_by_more_than_10x(self):
        """ERG1 is a NAMED NEIGHBOUR, not the branch minimum, and the name said otherwise
        until 2026-09-04.

        It read `test_the_slowest_host_isoprenoid_enzyme_...`, which stopped being true when
        `ISOPRENOID_KCAT_RANGE_PER_S` was recomputed over a declared gene set: the slowest is
        ERG6 at 0.011 /s, and ERG3, ERG5 and ERG24 are also below ERG1's 0.076. The 10x
        margin asserted here is reachable only BECAUSE ERG1 is not the minimum -- at the true
        minimum and this same 0.5% of protein the fold is 5.64x, so this assertion would
        fail at the enzyme the old name claimed it was testing.

        What holds at the real minimum is the weaker, directional claim, and it is made by
        `test_it_is_above_the_fitted_ceiling_at_every_corner_of_the_box`: over the whole box
        the fold stays above 1.0, and at the conservative corner (ERG6 at 0.1% of protein)
        it is 1.127x. Two different claims, each pinned where it is true.
        """
        fold = _derived(ERG1_KCAT, 0.005).disagreement_with(FITTED_CAPACITY, MU)

        assert fold > 10.0, f"derived is only {fold:.3g}x the fit"

    @pytest.mark.parametrize("kcat, fraction", BOX)
    def test_it_is_above_the_fitted_ceiling_at_every_corner_of_the_box(self, kcat, fraction):
        """The claim the module docstring is now allowed to make, and the only one.

        DIRECTIONAL over the whole (kcat, mass-fraction) box rather than a value at one
        point in it. The corner that matters is the conservative one -- the slowest branch
        enzyme at the weakest expression -- and even there the derived capacity is above the
        fit. A test written at a single interior point cannot say that, which is what the
        docstring's earlier "at ANY turnover" claim rested on.
        """
        fold = _derived(kcat, fraction).disagreement_with(FITTED_CAPACITY, MU)

        assert fold > 1.0, f"kcat {kcat:g} /s at {fraction:.1%} of protein gives {fold:.4g}x"

    def test_and_the_box_brackets_arhars_measurement_rather_than_landing_on_it(self):
        """What the module may NOT claim. The same box that is uniformly above the FIT spans
        Arhar's measurement by six orders of magnitude in each direction, so "lands on the
        measurement" is a statement about one corner and not about the method."""
        folds = [_derived(kcat, fraction).capacity_mmol_per_gdcw(MU) / ARHAR_CAPACITY
                 for kcat, fraction in BOX]

        assert min(folds) < 1.0 < max(folds)
        assert max(folds) / min(folds) > 1e6

    def test_and_it_overshoots_the_more_the_faster_the_enzyme(self):
        slow = _derived(kcat=ISOPRENOID_KCAT_RANGE_PER_S[0], fraction=0.005)
        fast = _derived(kcat=ISOPRENOID_KCAT_RANGE_PER_S[1], fraction=0.005)

        assert (fast.disagreement_with(FITTED_CAPACITY, MU)
                > slow.disagreement_with(FITTED_CAPACITY, MU))

    def test_reproducing_the_fitted_ceiling_needs_a_turnover_below_the_host_range(self):
        """The sharp version, and the one that says the FIT is the odd number rather than
        the method. Solving `kcat` for the fitted capacity puts it under the slowest enzyme
        in the host model's own isoprenoid branch."""
        enzyme_mmol = _derived(ERG1_KCAT, 0.005).enzyme_mmol_per_gdcw
        needed_kcat = FITTED_CAPACITY * MU / (3600.0 * enzyme_mmol)

        assert needed_kcat < ISOPRENOID_KCAT_RANGE_PER_S[0]

    def test_the_summary_says_which_way_the_disagreement_runs(self):
        text = _derived(ERG1_KCAT, 0.005).summary(FITTED_CAPACITY, MU)

        assert "above the fit" in text


class TestTheCheckCanActuallyFail:
    """An audit that cannot disagree is decoration."""

    def test_a_capacity_matching_the_fit_reports_agreement_not_disagreement(self):
        enzyme_mmol = _derived(ERG1_KCAT, 0.005).enzyme_mmol_per_gdcw
        matching_kcat = FITTED_CAPACITY * MU / (3600.0 * enzyme_mmol)
        d = derived_capacity("crtYB", matching_kcat, _content(0.005))

        assert d.disagreement_with(FITTED_CAPACITY, MU) == pytest.approx(1.0, rel=1e-9)
        assert "below the fit" not in d.summary(FITTED_CAPACITY, MU)


@_needs_ecmodel
class TestTheHostRangeIsRecomputedFromTheVendoredModel:
    """The constant is a min/max over a declared gene set, so this recomputes it.

    REPLACES `test_the_host_range_is_a_real_interval_and_not_a_placeholder`, which asserted
    only `0 < low < high` and `high/low > 100`. That is satisfied by almost any pair, and it
    is what let `(0.076, 46.0)` ship: the pair was two NAMED enzymes -- ERG1 and HMG1 -- with
    the gene set and direction rule that would have made it a min/max recorded nowhere. Five
    branch enzymes are slower than ERG1 and five faster than 46, two of them core MVA steps,
    so the upper bound was ~650x low.

    A derived number frozen as a literal cannot be repaired by retyping it, because the
    vendored model can be refreshed at any time. These tests reopen the file.
    """

    def test_the_shipped_interval_is_the_min_and_max_of_the_declared_branch(self):
        kcats = _branch_kcats_per_s()
        low, high = ISOPRENOID_KCAT_RANGE_PER_S
        values = [k for pairs in kcats.values() for k in pairs.values()]

        assert _sig3(min(values)) == _sig3(low), (
            f"branch minimum is {min(values):.6g} /s, constant says {low}")
        assert _sig3(max(values)) == _sig3(high), (
            f"branch maximum is {max(values):.6g} /s, constant says {high}")

    def test_the_interval_contains_every_enzyme_of_the_declared_branch(self):
        """The endpoints are rounded OUTWARD, and this is why. Rounding them to nearest put
        ERG6 and IDI1 -- the two enzymes that DEFINE the interval -- outside it, and
        `predict.py`'s layer then reported IDI1 as unlike anything the ecModel carries."""
        low, high = ISOPRENOID_KCAT_RANGE_PER_S
        for gene, reactions in _branch_kcats_per_s().items():
            for reaction_id, kcat in reactions.items():
                assert low <= kcat <= high, f"{gene} {reaction_id} at {kcat:.6g} /s"

    def test_every_declared_gene_is_actually_in_the_vendored_model(self):
        """A gene set is only data if a typo in it fails. A silently absent gene would
        narrow the interval without narrowing the claim."""
        kcats = _branch_kcats_per_s()

        assert set(kcats) == set(ISOPRENOID_BRANCH_GENES), (
            f"declared but not found: {sorted(set(ISOPRENOID_BRANCH_GENES) - set(kcats))}")

    def test_the_old_pair_was_two_named_enzymes_and_not_the_extremes(self):
        """The retraction, kept executable. ERG1 and HMG1 are real branch members and their
        turnovers are real -- what was false is that they bounded the branch."""
        kcats = _branch_kcats_per_s()
        erg1 = min(kcats["ERG1"].values())
        hmg1 = max(kcats["HMG1"].values())

        assert _sig3(erg1) == _sig3(0.076)
        assert round(hmg1, 1) == 46.2
        slower = [g for g, r in kcats.items() if min(r.values()) < erg1]
        faster = [g for g, r in kcats.items() if max(r.values()) > hmg1]

        assert len(slower) >= 3, slower
        assert len(faster) >= 3, faster

    def test_the_direction_rule_excludes_the_reverse_arms_and_the_pseudo_reactions(self):
        """The rule is declared beside the gene set, so it has to be checked too. IDI1 is
        the case that matters: its forward arm sets the branch maximum at 29900 /s and its
        reverse arm sits at 0.061 /s, so including reverse arms would move both ends."""
        kcats = _branch_kcats_per_s()

        assert not any(rid.startswith("draw_prot_") or rid.startswith("arm_")
                       or "_REV" in rid
                       for reactions in kcats.values() for rid in reactions)
        assert _sig3(max(kcats["IDI1"].values())) == _sig3(29900.0)

    def test_the_range_is_far_too_wide_to_be_a_prior(self):
        """Recorded as an assertion because it is what the width COSTS. A branch spanning
        six orders of magnitude rejects a typo and nothing finer, which is why the module's
        disagreement with the fitted ceiling is read at a NAMED enzyme's turnover."""
        low, high = ISOPRENOID_KCAT_RANGE_PER_S

        assert high / low > 1e6


@_needs_ecmodel
class TestTheGprCoverageClaimIsAlsoRecomputed:
    """The neighbouring statistic in the same docstring, pinned the same way.

    "a turnover number for 90.6% of its GPR reactions" is TRUE -- it was checked and
    survived on 2026-09-04 -- but it is derived prose frozen beside a derived interval that
    was not, and it was un-recomputable for exactly the same reason. The denominator is the
    subtlety: GECKO's `arm_*` isoenzyme-pooling reactions carry no protein stoichiometry by
    construction, so counting them gives 80.7% and is the wrong denominator.
    """

    def test_the_coverage_is_the_shipped_percentage(self):
        with_kcat, total = _gpr_kcat_coverage()

        assert total == 4258
        assert round(100.0 * with_kcat / total, 1) == 90.6


class TestThePredictLayerReportsTheRangeItActuallyHolds:
    """`predict.py`'s enzyme-capacity layer flags a kcat outside the branch's interval.

    Pinned against `ISOPRENOID_KCAT_RANGE_PER_S` and never against a literal. The NOTE fired
    falsely for every kcat in (46, 29900] while the constant was wrong -- ERG10 at 900 /s
    and IDI1 at 29900 /s, both core MVA enzymes IN the branch, were reported as "an enzyme
    unlike the ones the vendored ecModel carries". A test written against `46.0` would have
    survived that, which is why this one reads the constant.
    """

    @staticmethod
    def _detail(kcat):
        from ystwin.predict import _enzyme_capacity_layer

        return _enzyme_capacity_layer(
            {"enzyme": "crtYB", "kcat_per_s": kcat,
             "enzyme_mmol_per_gdcw": _content(0.005)}, FITTED_CAPACITY, MU).detail

    def test_it_stays_quiet_just_inside_each_end(self):
        low, high = ISOPRENOID_KCAT_RANGE_PER_S

        assert "NOTE" not in self._detail(low * 1.001)
        assert "NOTE" not in self._detail(high * 0.999)

    def test_it_fires_just_outside_each_end(self):
        low, high = ISOPRENOID_KCAT_RANGE_PER_S

        assert "NOTE" in self._detail(low * 0.999)
        assert "NOTE" in self._detail(high * 1.001)

    def test_the_note_does_not_fire_for_a_core_mva_enzyme_of_the_branch(self):
        """The specific false positive the old constant produced. ERG10 and IDI1 are both
        in `ISOPRENOID_BRANCH_GENES`, so nothing may call them unlike it."""
        kcats = _branch_kcats_per_s() if paths.ec_yeast_gem() else None
        if kcats is None:
            pytest.skip("ecYeastGEM not vendored")
        for gene in ("ERG10", "IDI1"):
            for kcat in kcats[gene].values():
                assert "NOTE" not in self._detail(kcat), gene
