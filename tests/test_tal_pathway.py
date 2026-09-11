"""Triacetic acid lactone was the top-ranked second product, and it does not ship. This is why.

`docs/SECOND_DATASET_HUNT.md` searched eight product classes for a dataset that could
calibrate the flux law a second time, and put **TAL in *Kluyveromyces marxianus*** first of
four: eighteen states, six genotypes by three temperatures, a biomass-normalised product,
and the only candidate with any growth rate attached. Lang X, Besada-Lombana PB, Li M,
Da Silva NA, Wheeldon I. 2020, *Metab Eng Commun* 11:e00145, PMID 32995271. The PMC
identifier, the doi and the article's erratum are in `data/tal/SOURCE.md`, which is the
provenance document for this dataset.

The paper was retrieved and read in full -- JATS XML from Europe PMC plus the `mmc1.docx`
supplement -- and **no `data/pathways/tal.toml` can be written from it honestly.** Three
blockers, recorded at length in `data/tal/SOURCE.md` and pinned here so that the negative is
executable rather than a paragraph somebody has to find:

1. **The product axis is a bar chart.** Fig. 4B plots `TAL (mg/L)/OD600` for all eighteen
   states and prints no values; Fig. S7 repeats it at stationary phase. The article states
   exactly one number about that panel -- a 17.8-fold spread between its extremes -- which
   is one constraint on eighteen unknowns. There is no data-availability statement, and the
   supplement holds primers and promoter sequences only. Digitising the bars is what this
   repository refuses to do, so the fit was never attempted.

2. **TAL is assayed in the supernatant**, by both the HPLC and the plate-reader method. It
   is extracellular, so its terminal node's fate is `Fate.SECRETED`, which `spec.py` refuses
   at load. That refusal is the right answer and not an obstacle -- the tests below show
   that the *only* thing between it and a meaningless mg/gDCW content is which word the
   author types into the file.

3. **Expression is numeric at 30 degC only**, and the authors report that the rank
   correspondence between promoter strength and TAL production is *lost* at 37 and 41 degC.
   So the promoter-identity pairing that `SECOND_DATASET_HUNT.md` flagged as this
   dataset's catch is not merely untested; it is reported broken on twelve of the eighteen
   states by the people who took the measurements.

**What is shipped** is the two axes that *are* numeric, so nobody repeats the retrieval:
`data/tal/lang2020_promoter_strength_xylose_30c.tsv` and
`data/tal/lang2020_growth_rates_xylose.tsv`. Both are quoted from article prose, not read
off a figure, and the tests here re-derive the article's own percentages from the first one
to prove which carbon source it belongs to.

**And one thing worth stating positively.** Six genotypes is a better-resolved design than
the three the flux law is actually calibrated on: the smallest attainable permutation p goes
from 1/3! = 0.167, which cannot reach 0.05 at any alpha, to 1/6! = 0.00139, which can. The
missing eighteen numbers are therefore worth asking the authors for, and
`test_six_genotypes_would_clear_the_floor_that_three_cannot` is the arithmetic that says so.
"""

from __future__ import annotations

import math

import pytest

from ystwin import paths
from ystwin.analysis.multiplicity import sign_test_floor
from ystwin.pathway.flux import fit_flux_law
from ystwin.pathway.solve import solve_pathway
from ystwin.pathway.spec import Fate, Node, PathwaySpec, RateLaw, available_pathways

TAL_DIR = paths.data_dir() / "tal"

# The six promoters that drive 2-PS in Fig. 4B, and the only expression numbers the paper
# prints. Fig. 3B caption, verbatim: "The absolute value of RFU/OD for each promoter is
# 25 +/- 18 (ADH1), 47 +/- 7 (HHF1), 95 +/- 14 (NC1), 10 +/- 1 (PGK), 15 +/- 11 (SSA3) and
# 86 +/- 3 (TEF3) at 30 degC." Background-subtracted, biological triplicates.
#
# THESE ARE EGFP REPORTER STRAINS, NOT TAL PRODUCERS. Separate transformants on the same
# low-copy backbone. Expression is joined to the product panel by promoter identity alone
# and was never measured on a cell that makes TAL.
PROMOTER_RFU_PER_OD = {"PGK": 10.0, "SSA3": 15.0, "ADH1": 25.0,
                       "HHF1": 47.0, "TEF3": 86.0, "NC1": 95.0}
PROMOTER_RFU_SD = {"PGK": 1.0, "SSA3": 11.0, "ADH1": 18.0,
                   "HHF1": 7.0, "TEF3": 3.0, "NC1": 14.0}

# Results, verbatim: "(0.28 h-1 at 30 degC, 0.34 h-1 at 37 degC, and 0.35 h-1 at 41 degC)".
# Measured on the BLANK-vector strain, one value per temperature, none per genotype.
MU_PER_H_BY_TEMPERATURE = {30: 0.28, 37: 0.34, 41: 0.35}

# The TAL chain as the paper draws it (Fig. 4A): acetyl-CoA + malonyl-CoA, condensed twice
# by the Gerbera hybrida type III polyketide synthase 2-PS with two decarboxylations, to
# 4-hydroxy-6-methyl-2-pyrone. Written out here ONLY so the refusal below has something
# real to refuse -- it is deliberately not a TOML file in data/pathways/, because there is
# no product measurement to go with it.
TAL_MOLAR_MASS_G_PER_MOL = 126.11  # PubChem CID 54675757, C6H6O3, MolecularWeight 126.11, IUPACName
# "4-hydroxy-6-methylpyran-2-one". Fetched from the PUG REST property endpoint.


def _tal_chain(product_fate: str) -> PathwaySpec:
    """The TAL chain with the terminal node's fate as the caller declares it."""
    return PathwaySpec(
        product="tal",
        organism="Kluyveromyces marxianus",
        entry_enzyme="2PS",
        precursor_metabolite="",
        nodes=(
            Node(name="tal", enzyme="", rate_law=RateLaw.PASSTHROUGH,
                 fate=product_fate, molar_mass_g_per_mol=TAL_MOLAR_MASS_G_PER_MOL),
        ),
    )


class TestTheSpecIsAbsentOnPurpose:
    """Not an oversight, and not a to-do. There is no product column to calibrate against."""

    def test_no_tal_pathway_spec_ships(self) -> None:
        """`data/pathways/` names no TAL, because Fig. 4B is a bitmap.

        If this ever fails, the person adding `tal.toml` owes the file a `SOURCE.md` line
        saying where the eighteen per-strain numbers came from -- and they cannot have come
        from this paper, which does not contain them.
        """
        assert "tal" not in available_pathways()

    def test_the_shipped_tables_carry_no_product_column(self) -> None:
        """Nothing in `data/tal/` claims a TAL measurement, at any precision.

        The guard is against the one failure mode this directory invites: someone reads
        `SOURCE.md`, decides the bars are readable enough, and adds a `tal_mg_per_l_od`
        column. A digitised value and a quoted value look identical in a TSV, so the
        invariant is enforced on the schema rather than on anyone's memory.
        """
        for table in sorted(TAL_DIR.glob("*.tsv")):
            header = table.read_text(encoding="utf-8").splitlines()[0].split("\t")
            product_columns = [c for c in header
                               if "tal" in c.lower() or "titer" in c.lower()
                               or "titre" in c.lower() or "product" in c.lower()]
            assert product_columns == [], (
                f"{table.name} declares {product_columns}. The per-strain TAL values exist "
                "only as bar heights in Fig. 4B / Fig. S7 of PMID 32995271; if a real "
                "source for them has been found, say so in data/tal/SOURCE.md and delete "
                "this test rather than working around it")


class TestASecretedProductIsWhyThisChainCannotBeModelledHere:
    """Blocker 2, executable. The assay says supernatant, and the loader says no."""

    def test_the_secreted_declaration_is_refused_and_names_the_missing_term(self) -> None:
        """`Fate.SECRETED` on the terminal node is refused, and the message says why.

        The message has to name `mu*[X]` rather than say "unsupported", because the reason
        is arithmetic and not policy: TAL leaves through the medium at a rate the cell sets,
        so the balance the whole solver closes on does not describe it.

        The refusal MOVED on 2026-09-01, from the loader to the solver, and the blocker is
        unchanged in force: the chain now loads, and asking it for a number still fails,
        because an export capacity for TAL has never been measured. What moved is what the
        message can say -- it now names the constants and the experiment that would produce
        them, instead of declaring the pathway out of scope.
        """
        chain = _tal_chain(Fate.SECRETED)

        with pytest.raises(ValueError, match=r"mu\*\[X\]"):
            solve_pathway(chain, entry_flux=1e-3, growth_rate=0.1)

    def test_the_identical_chain_loads_the_moment_the_fate_is_typed_diluted(self) -> None:
        """One word is the whole difference, which is why the word must come from the assay.

        This is the failure the refusal prevents, demonstrated rather than described:
        nothing about the chemistry, the molar mass or the chain changes. Declare `diluted`
        and the spec is accepted, the solver will return `content = v_in/mu` in mg/gDCW, and
        every one of those numbers is about an intracellular pool that this paper never
        measured. Lang 2020's Methods are explicit -- both TAL assays read the supernatant --
        so `secreted` is the reading the data supports and `diluted` is a typo with no
        symptom.
        """
        spec = _tal_chain(Fate.DILUTED)
        assert spec.nodes[-1].fate == Fate.DILUTED
        assert spec.nodes[-1].molar_mass_g_per_mol == TAL_MOLAR_MASS_G_PER_MOL


class TestTheExpressionAxisIsRealAndSmallerThanItLooks:
    """Blocker 3, and the provenance check that pins the six numbers to a carbon source."""

    def test_the_absolutes_are_the_xylose_panel_not_the_glucose_one(self) -> None:
        """The Fig. 3B caption does not say which carbon source, and the panel compares two.

        The article's prose decides it: "In glucose, P_PGK was found to be a medium level
        promoter, reaching 28% of P_TEF3, but in xylose expression was reduced to less than
        12% of P_TEF3. Growth in xylose had the opposite effect on P_ADH1, increasing
        expression to 28% of P_TEF3."

        PGK/TEF3 from the quoted absolutes is 11.6%, which is "less than 12%" and is not
        28%; ADH1/TEF3 is 29.1%, which is the ~28% quoted for xylose. So the absolutes are
        the xylose set -- the same carbon source as the TAL panel -- and this is a derivation
        from the paper rather than an assumption made here.
        """
        tef3 = PROMOTER_RFU_PER_OD["TEF3"]
        pgk_fraction = PROMOTER_RFU_PER_OD["PGK"] / tef3
        adh1_fraction = PROMOTER_RFU_PER_OD["ADH1"] / tef3
        assert pgk_fraction < 0.12, "xylose: 'reduced to less than 12% of P_TEF3'"
        assert abs(adh1_fraction - 0.28) < 0.02, "xylose: 'increasing expression to 28%'"
        # And the glucose reading is excluded, not merely less likely.
        assert abs(pgk_fraction - 0.28) > 0.15, "glucose would put PGK at 28% of TEF3"

    def test_the_shipped_table_matches_those_absolutes(self) -> None:
        """`lang2020_promoter_strength_xylose_30c.tsv` is the caption, not a re-derivation."""
        rows = TAL_DIR.joinpath("lang2020_promoter_strength_xylose_30c.tsv").read_text(
            encoding="utf-8").splitlines()
        header = rows[0].split("\t")
        table = {r.split("\t")[header.index("promoter")]: r.split("\t") for r in rows[1:]}
        assert set(table) == set(PROMOTER_RFU_PER_OD)
        for promoter, row in table.items():
            assert float(row[header.index("rfu_per_od600")]) == PROMOTER_RFU_PER_OD[promoter]
            assert float(row[header.index("rfu_per_od600_sd")]) == PROMOTER_RFU_SD[promoter]
            derived = PROMOTER_RFU_PER_OD[promoter] / PROMOTER_RFU_PER_OD["TEF3"]
            assert abs(float(row[header.index("relative_to_tef3")]) - derived) < 5e-5

    def test_two_of_the_six_expression_values_are_barely_resolved(self) -> None:
        """SSA3 and ADH1 carry relative standard deviations above 70% on triplicates.

        Worth pinning because it survives even the optimistic reading of this dataset. A law
        fitted in log space is most sensitive at the low end, and the lower one-sigma edge of
        SSA3 is 4 RFU/OD against a mean of 15 -- so the weakest two of the six promoters,
        which are exactly the ones that set the span of the expression axis, are the two
        whose values are least certain.
        """
        noisy = {p for p in PROMOTER_RFU_PER_OD
                 if PROMOTER_RFU_SD[p] / PROMOTER_RFU_PER_OD[p] > 0.70}
        assert noisy == {"SSA3", "ADH1"}
        assert PROMOTER_RFU_PER_OD["SSA3"] - PROMOTER_RFU_SD["SSA3"] == 4.0


class TestGrowthRateCarriesNoGenotypeInformationHere:
    """The denominator this chain needs exists, and it does not vary the way it must."""

    def test_three_growth_rates_for_six_genotypes(self) -> None:
        """mu is per TEMPERATURE, so within a temperature it is a shared constant.

        `docs/SECOND_DATASET_HUNT.md` ranks TAL first because it is "the only [candidate]
        with any mu at all", which is true and is the right reason to rank it. But the pools
        layer is `X = v_in / mu`, and dividing six genotypes at one temperature by one number
        rescales all six identically -- it cannot separate them. Contrast the beta-carotene
        chemostat, where the dilution rate IS the state: two per strain, set by the operator,
        and the thing the growth-exponent question is asked of.
        """
        assert len(MU_PER_H_BY_TEMPERATURE) == 3
        assert len(PROMOTER_RFU_PER_OD) == 6
        assert len(set(MU_PER_H_BY_TEMPERATURE.values())) < len(PROMOTER_RFU_PER_OD)

    def test_the_temperature_span_of_mu_is_small_and_would_rescale_a_derived_q(self) -> None:
        """1.25-fold from 30 to 41 degC, which matters only because `q ~ mu * (P/X)`.

        Fig. 4B's axis is `TAL (mg/L)/OD600` -- an extracellular concentration over an
        optical density, with no time dimension. It is a specific TITRE, not the specific
        productivity `docs/SECOND_DATASET_HUNT.md` credits it as. Converting one to the
        other needs `q ~ mu * (P/X)`, which holds only under balanced exponential growth with
        q constant from inoculation. That assumption is not neutral across this panel: it
        multiplies the 41 degC states by 1.25 relative to the 30 degC ones, so it moves the
        temperature axis it would be invoked to interpret.
        """
        span = MU_PER_H_BY_TEMPERATURE[41] / MU_PER_H_BY_TEMPERATURE[30]
        assert abs(span - 1.25) < 0.005

    def test_the_shipped_growth_table_matches_and_marks_45c_out_of_scope(self) -> None:
        rows = TAL_DIR.joinpath("lang2020_growth_rates_xylose.tsv").read_text(
            encoding="utf-8").splitlines()
        header = rows[0].split("\t")
        table = {int(r.split("\t")[header.index("temperature_c")]): r.split("\t")
                 for r in rows[1:]}
        for temperature, mu in MU_PER_H_BY_TEMPERATURE.items():
            assert float(table[temperature][header.index("mu_per_h")]) == mu
            assert table[temperature][header.index("in_scope_for_the_2ps_panel")] == "yes"
        # 45 degC exists in the paper's growth data and in no TAL panel.
        assert float(table[45][header.index("mu_per_h")]) == 0.14
        assert table[45][header.index("in_scope_for_the_2ps_panel")] == "no"


class TestTheStatisticsThisDatasetWouldSupport:
    """The honest reason to keep asking for the eighteen numbers, stated as arithmetic."""

    def test_six_genotypes_would_clear_the_floor_that_three_cannot(self) -> None:
        """1/6! = 0.00139 against 1/3! = 0.167, and only one of those can reach 0.05.

        The fixed CrtE comparison now has skill +0.745 against training-only baselines.
        Permuting expression profiles across three strains still has only 3! = 6
        assignments, so even rank 1 cannot reach 0.05. Renaming strain labels alone would
        not change the partition and is not this null. Six genotypes permit a finer
        permutation rank; neither sample size nor a good score guarantees significance.
        """
        beta_carotene_floor = 1.0 / math.factorial(3)
        tal_floor = 1.0 / math.factorial(6)
        assert abs(beta_carotene_floor - 0.1667) < 1e-3
        assert abs(tal_floor - 0.0013889) < 1e-6
        assert beta_carotene_floor > 0.05 > tal_floor

    def test_the_temperature_axis_has_its_own_floor_and_it_is_a_quarter(self) -> None:
        """Three temperatures is three clusters, and an exact sign test floors at 0.25.

        Recorded because the temperature contrast is the comparison this paper's own
        Discussion leans on ("specific TAL titers at 37 and 41 degC were significantly
        greater than those achieved at 30 degC"). Treated as a distribution-free test over
        the three temperature levels it could never return anything below 0.25 -- the same
        shape of limit `multiplicity.py` found at three plates. Their claim is not made that
        way, and is not being disputed here; what is pinned is that anyone re-testing it
        across the temperature axis alone inherits this floor.
        """
        assert sign_test_floor(len(MU_PER_H_BY_TEMPERATURE)) == 0.25
        assert sign_test_floor(3) > 0.05

    def test_a_fit_of_the_expression_axis_against_itself_scores_a_perfect_one(self) -> None:
        """What a circular fit looks like, on the numbers this dataset actually hands you.

        This is the specific temptation `data/tal/` creates. Six real expression values are
        sitting here and the product column is not, so the cheapest way to make
        `fit_flux_law` run is to give it something else in place of a flux. Do that with the
        expression itself and it returns alpha = 1, zero leave-one-strain-out error and skill
        +1.000 -- a better score than any real calibration in this repository, from a fit
        that contains no information at all.

        `flux.py` exists because exactly this happened once already, with the product
        supplying its own predictor. The signature is a skill indistinguishable from 1.
        """
        promoters = sorted(PROMOTER_RFU_PER_OD)
        values = [PROMOTER_RFU_PER_OD[p] for p in promoters]
        circular = fit_flux_law(values, values, promoters, "2PS", "not a fit; see docstring")
        assert circular.loso_skill == pytest.approx(1.0)
        assert circular.loso_rmse_log == pytest.approx(0.0, abs=1e-12)
        assert circular.alpha == pytest.approx(1.0)

    def test_the_real_calibration_scores_far_below_that(self) -> None:
        """The real fixed-model comparison is imperfect and was selected on these data."""
        from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX

        assert BETA_CAROTENE_FLUX.loso_skill == pytest.approx(0.7447, abs=5e-4)
        assert BETA_CAROTENE_FLUX.loso_rmse_log == pytest.approx(0.2002, abs=5e-4)
        assert BETA_CAROTENE_FLUX.loso_skill < 1.0
