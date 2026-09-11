"""The pathway audit, re-measured rather than trusted.

Every number asserted here was produced by running this file against yeast-GEM v9.0.2 through
GLPK with cobra 0.32.1: growth 0.085844 /h under the model's shipped default bounds and
0.345984 /h under ``physiology.REFERENCE_AEROBIC_BATCH``. The floor throughout is
``generator.panel_experiment.MEASURED_GROWTH_RATE_SE`` = 0.0117 /h, which is the standard error
of a growth rate read off a real optical-density trace; the plate reader's activity CV is
deliberately NOT borrowed for it, because it measures a different quantity.

Three tests are tripwires rather than checks, and they are the point of the file. The
zero-reach test pins the audit's headline finding -- four declared panel modules whose regulon
reaches no reaction at all -- so that a later repair to ``bridge/regulation.py`` turns this file
red and forces the table to be re-audited instead of silently going stale. The symbol-trap test
pins the fact that makes the whole audit possible, that cobrapy leaves every ``Gene.name`` empty.
And the neutrality test pins that installing the coverage groups changes no optimum, which is
the claim that this module declares and does not model.
"""

from __future__ import annotations

import csv
import gzip
import re

import cobra
import pytest

from ystwin import paths
from ystwin.fba import physiology, stress_genes, stress_ros
from ystwin.fba.solver import configure
from ystwin.fba.stress_pathways import (
    COVERAGE_PATH,
    LAYERS,
    REFUSED_REPRESSION_COUPLING,
    UNREACHED_MODULES,
    VERDICTS,
    PathwayItem,
    apply_pathway_coverage,
    coverage_by_layer,
    gaps,
    group_name,
    load_coverage,
    panel_module_reactions,
    pathway_essentiality,
    refusals,
    represented,
    unrepresented_pathways,
    verify_coverage,
)
from ystwin.generator.panel_experiment import MEASURED_GROWTH_RATE_SE
from ystwin.generator.stress_panel import MODULES, STRESSORS
from ystwin.mech.params import FREE_TAGS, RefusedValue

SHIPPED_DEFAULT_GROWTH = 0.085844
REFERENCE_BATCH_GROWTH = 0.345984

#: The owner's diagram, transcribed from the assignment. Every name here must appear in the
#: table under some verdict; the table may be a superset but never a subset.
DIAGRAM_SENSORS = (
    "WSC1", "MID2", "MTL1", "SLN1", "SHO1", "IRE1", "YAP1", "SKN7", "HSF1", "SNF1", "MEC1",
    "TEL1", "GPR1", "TOR", "PMA1", "VATPASE", "RIM101", "SPS",
)
DIAGRAM_CASCADES = (
    "HOG", "CWI", "UPR", "cAMP", "TORC1", "SNF1", "RIM101_proteolysis", "calcineurin",
)
DIAGRAM_TFS = (
    "MSN2", "MSN4", "HAC1", "HOT1", "SKO1", "RIM101", "HSF1", "YAP1", "SKN7", "CRZ1", "RLM1",
    "AFT1", "AFT2", "ZAP1", "MET4", "PDR1", "PDR3", "RTG1", "RTG3", "ADR1", "CAT8", "GLN3",
    "GAT1", "ROX1", "UPC2", "HAP1", "ACE1", "MAC1", "RFX1", "DOT6", "TOD6",
)
DIAGRAM_ADAPTATIONS = (
    "chaperones", "er_expansion", "glycerol", "trehalose", "glycogen", "antioxidant_sod",
    "antioxidant_catalase", "antioxidant_gpx", "antioxidant_trx", "antioxidant_tsa1",
    "antioxidant_glutathione", "cell_wall_chitin", "cell_wall_b13_glucan",
    "cell_wall_b16_glucan", "cell_wall_mannan", "ergosterol", "ion_ena1", "ion_nha1",
    "ion_trk", "autophagy_glycophagy", "autophagy_bulk_protein", "dna_repair_rad",
    "dna_repair_rnr", "ribosome_reduction", "flux_reprogramming",
    "manganese_homeostasis", "n_glycosylation", "cell_wall_crosslinking",
)


#: Every PMID the table cites, each resolved against PubMed E-utilities on 2026-09-05 and
#: checked first-author-and-year against the claim it sits beside. A first pass of this table
#: carried 61 PMIDs that resolved to unrelated papers -- a dopamine-receptor cloning paper for
#: Kre6, an otitis-media study for Sod1 -- with every author and year correct and only the
#: digits drifted, which is exactly why reading the prose never catches it. This set is the
#: tripwire: a new or edited citation trips it, and the fix is to re-resolve, never to widen it.
VERIFIED_PMIDS = frozenset({
    "1456112", "1657642", "1715094", "1837148", "2072919", "3005867", "3037314", "3043194",
    "8643535",
    "3044613", "3526554", "7507493", "7528927", "7624781", "7629054", "7681220", "7891685",
    "7926756", "7961686", "8078477", "8262047", "8293472", "8313910", "8422683", "8599111",
    "8602515", "8622686", "8641288", "8654575", "8662189", "8670839", "8808622", "8898193",
    "8947468", "9017390", "9032238", "9118942", "9171333", "9199298", "9271382", "9371781",
    "8842708", "9407035", "9409150", "9510529", "9741624", "9858577", "9920864", "10330137",
    "10361302", "10409737", "10476026", "10594829", "10604478", "10844700", "10856230",
    "11095737", "11102521", "11113201", "11152943", "11154269", "11212295", "11309404",
    "11321575", "11445588", "11448968", "11486014", "11533229", "11698381", "11918814",
    "10809732", "11027260", "12408816", "12437921", "12456786", "12676948", "13034889",
    "15944456", "16027116",
    "16524922", "17397106", "17981722", "18694928", "19901341", "19948500", "21730963",
    "24762745", "28857745",
})

@pytest.fixture(scope="module")
def items():
    return load_coverage()


@pytest.fixture(scope="module")
def _base_model_template(yeast_gem_factory):
    """A private copy, so pinning the solver cannot leak into the session-scoped fixture."""
    model = yeast_gem_factory()
    configure(model)
    return model


@pytest.fixture
def base_model(_base_model_template, model_copy):
    return model_copy(_base_model_template)


@pytest.fixture(scope="module")
def _restored_model_template(_base_model_template, model_copy):
    model = stress_genes.restore_stress_genes(model_copy(_base_model_template))
    configure(model)
    return model


@pytest.fixture
def restored_model(_restored_model_template, model_copy):
    return model_copy(_restored_model_template)


def _row(items: tuple[PathwayItem, ...], key: str) -> PathwayItem:
    return next(r for r in items if r.item == key)


# --------------------------------------------------------------------------- #
# the table itself
# --------------------------------------------------------------------------- #


def test_table_parses_and_covers_four_layers(items):
    assert len(items) == 97
    assert {r.layer for r in items} == set(LAYERS)
    assert {r.verdict for r in items} <= set(VERDICTS)
    assert len({r.item for r in items}) == len(items)


def test_every_item_of_the_owners_diagram_has_a_row(items):
    """The roll-call. A diagram item missing from the table is the failure this audit exists
    to prevent, and it cannot be caught by any check that only reads the table."""
    keys = {r.item for r in items}
    for layer, names in (
        ("sensor", DIAGRAM_SENSORS),
        ("cascade", DIAGRAM_CASCADES),
        ("tf", DIAGRAM_TFS),
        ("adaptation", DIAGRAM_ADAPTATIONS),
    ):
        for name in names:
            matches = [k for k in keys if k.startswith(f"{layer}.{name}")]
            assert matches, f"{layer}.{name} is on the owner's diagram and has no row"


def test_the_gap_list_is_six_and_not_fiftyseven(items):
    """NOT_METABOLIC is the verdict that carries the audit: 57 rows are absent from the GEM
    and CORRECTLY so, and counting them as gaps would give a 57-item list where the real
    one is six."""
    census = {v: sum(1 for r in items if r.verdict == v) for v in VERDICTS}
    assert census["NOT_METABOLIC"] == 57
    assert census["IN_GEM"] == 24
    assert census["IN_GEM_RESTORED"] == 8
    assert census["IN_MODULE"] == 2
    assert census["REFUSED"] == 6
    assert census["GAP"] == 0, "a closable GAP would need a route in its basis, not a refusal"
    assert {r.item for r in gaps()} == {
        "tf.DOT6",
        "tf.TOD6",
        "adaptation.xenobiotic_drug_efflux",
        "adaptation.autophagy_bulk_protein",
        "adaptation.transcriptional_repression",
        "adaptation.cell_wall_crosslinking",
    }
    assert len(represented(items)) == 34
    assert {r.verdict for r in refusals(items)} == {"REFUSED"}


def test_the_panel_layer_holes_are_three_and_are_not_the_gap_list(items):
    """The second kind of hole, and the reason it needs its own accessor. A NOT_METABOLIC
    row is correctly absent from the MODEL; its pathway can still be absent from the PANEL,
    which is the only layer that could ever carry a GPCR or a kinase. Two pathways, three
    rows, and disjoint from gaps() by construction."""
    holes = unrepresented_pathways(items)
    assert {r.item for r in holes} == {"sensor.GPR1", "cascade.PKA_kinases", "sensor.SPS"}
    assert {r.verdict for r in holes} == {"NOT_METABOLIC"}
    assert not {r.item for r in holes} & {r.item for r in gaps(items)}
    for record in holes:
        assert record.panel_module in ("", "-") and record.code_module in ("", "-")


def test_the_two_inert_rows_no_module_reads_are_not_called_holes(items):
    """The audit's own padding check. Both rows have their chemistry in the GEM and this
    file measures both inert at the optimum, so neither is a gap at either layer -- but a
    definition written on the module columns alone would sweep both in."""
    unread = {
        r.item for r in items
        if r.panel_module in ("", "-") and r.code_module in ("", "-")
        and r.verdict in ("IN_GEM", "IN_GEM_RESTORED")
    }
    assert unread == {"sensor.VATPASE", "cascade.cAMP_second_messenger"}
    assert not unread & {r.item for r in unrepresented_pathways(items)}


def test_not_metabolic_rows_never_claim_a_reaction(items):
    for record in items:
        if record.verdict in ("NOT_METABOLIC", "GAP", "REFUSED", "IN_MODULE"):
            assert not record.reaction_ids, f"{record.item} claims reactions under {record.verdict}"


def test_every_row_carries_a_basis_and_a_source(items):
    for record in items:
        assert len(record.basis) > 40, f"{record.item}: basis is too thin to be a ground"
        assert record.source, f"{record.item}: no source"


def test_no_citation_enters_the_table_without_being_resolved(items):
    """The guard on the repository's worst failure mode. Adding a PMID here must be a
    deliberate act that includes resolving it, not a paste that nobody can check offline."""
    cited = {m for r in items for m in re.findall(r"PMID[: ]*(\d{5,8})", r.source + " " + r.basis)}
    assert cited == VERIFIED_PMIDS, (
        f"unverified: {sorted(cited - VERIFIED_PMIDS)}; gone: {sorted(VERIFIED_PMIDS - cited)}. "
        "Resolve each against PubMed and check the title against the claim beside it")


def test_every_row_names_an_author_and_a_year_beside_its_identifier(items):
    """A bare PMID cannot be checked by eye. The author and year are what make the drift
    visible in a diff."""
    for record in items:
        if "PMID" in record.source:
            assert re.search(r"[A-Z][\w'-]+.*\d{4}.*PMID", record.source), record.item


def test_every_orf_is_sgds_primary_systematic_name(items):
    """The ORF-side twin of the PMID guard, and it caught one. A merged secondary systematic
    name resolves in SGD and is absent from yeast-GEM exactly like the primary, so nothing
    else in this file would notice: ``cascade.RIM101_proteolysis`` carried Rim8 as YGL046W,
    the name merged into YGL045W, and both are absent from the GEM so every verdict held.
    SGD's GAF puts the primary systematic name FIRST in the synonym field, which is the
    authoritative signal; the aliases behind it are history.
    """
    gaf = paths.data_dir() / "external" / "sgd" / "gene_association.sgd.gaf.gz"
    if not gaf.exists():
        pytest.skip(f"{gaf} is gitignored and absent")
    primary, secondary = set(), {}
    with gzip.open(gaf, "rt") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            orfs = [p for p in fields[10].split("|") if re.fullmatch(r"Y[A-Z]{2}\d{3}[WC](-[A-Z])?", p)]
            if not orfs:
                continue
            primary.add(orfs[0])
            for alias in orfs[1:]:
                secondary.setdefault(alias, orfs[0])
    named = {orf for record in items for orf in record.orf_ids}
    assert named, "the table names no ORF at all"
    merged = {o: secondary[o] for o in named if o in secondary and o not in primary}
    assert not merged, f"merged secondary names; use the primary: {merged}"
    unknown = {o for o in named if o not in primary and o not in secondary}
    assert not unknown, f"not an SGD systematic name at all: {sorted(unknown)}"


def test_every_gene_the_restore_adds_lands_in_some_pathway_row(items):
    """The seam between the two halves of this assignment, and the invariant that keeps them
    from drifting apart. ``fba/stress_genes.py`` puts 19 genes into the model; if one of them
    belongs to no pathway item then a gene was restored for a function this table does not
    represent, which is exactly the state the audit was commissioned to end. It held only
    after this pass added manganese, N-glycosylation and the four late mannan transferases.

    The 20 gene-layer REFUSALS are treated differently and deliberately: seven of them are the
    cell-wall cross-linking cluster, a distinct adaptation that now has its own REFUSED row,
    and the other thirteen are accessory subunits or paralogues of items already represented.
    They are pinned by name rather than given rows, because a row apiece would be the padding
    the NOT_METABOLIC verdict exists to prevent.
    """
    genes = paths.data_dir() / "gem" / "absent_stress_genes.tsv"
    with genes.open(newline="") as handle:
        rows = list(csv.DictReader(
            (ln for ln in handle if not ln.startswith("#")), delimiter="\t"))
    named = {orf for record in items for orf in record.orf_ids}
    added = {r["orf"]: r["symbol"] for r in rows
             if r["verdict"] in ("GPR_EXTEND", "NEW_REACTION")}
    assert len(added) == 19
    assert not {s for o, s in added.items() if o not in named}, "restored for no represented item"
    refused = {r["orf"]: r["symbol"] for r in rows if r["verdict"] == "REFUSED"}
    assert len(refused) == 20
    assert {s for o, s in refused.items() if o in named} == {
        "GAS1", "GAS2", "GAS3", "GAS4", "GAS5", "CRH1", "UTR2"}
    assert {s for o, s in refused.items() if o not in named} == {
        "CCS1", "CYB5", "DAP1", "ERG28", "GLC3", "KHA1", "MNN4", "TOK1", "VCX1", "VMA21",
        "VMA22", "VNX1", "VPH2"}


def test_every_named_panel_module_and_code_module_exists(items):
    """A row pointing at a module that is not there would be coverage on paper only."""
    for record in items:
        if record.panel_module and record.panel_module != "-":
            assert record.panel_module in MODULES, record.item
        if record.code_module and record.code_module != "-":
            __import__(f"ystwin.{record.code_module}")


@pytest.mark.parametrize(
    "row,message",
    [
        ({"layer": "NOWHERE"}, "unknown layer"),
        ({"verdict": "MAYBE"}, "unknown verdict"),
        ({"basis": ""}, "basis is required"),
        ({"source": ""}, "source is required"),
        ({"verdict": "IN_GEM", "gem_reactions": "-"}, "names no reaction"),
        ({"verdict": "NOT_METABOLIC", "gem_reactions": "r_0005"}, "names reactions"),
        ({"verdict": "IN_MODULE", "gem_reactions": "-", "panel_module": "-",
          "code_module": "-"}, "names no module"),
    ],
)
def test_loader_refuses_a_self_inconsistent_row(tmp_path, items, row, message):
    """Half-applying an inconsistent declaration is worse than refusing it."""
    good = _row(items, "cascade.CWI_RHO1")
    fields = {c: getattr(good, c) for c in
              ("layer", "item", "label", "orfs", "verdict", "gem_reactions",
               "panel_module", "code_module", "basis", "source")}
    fields.update(row)
    path = tmp_path / "bad.tsv"
    path.write_text("\t".join(fields) + "\n" + "\t".join(fields[c] for c in fields) + "\n")
    with pytest.raises(ValueError, match=message):
        load_coverage(path)


def test_loader_refuses_a_duplicated_item(tmp_path, items):
    header = COVERAGE_PATH.read_text().splitlines()
    body = [ln for ln in header if not ln.startswith("#")]
    path = tmp_path / "dup.tsv"
    path.write_text("\n".join([body[0], body[1], body[1]]) + "\n")
    with pytest.raises(ValueError, match="appears twice"):
        load_coverage(path)


def test_coverage_by_layer_partitions_the_table(items):
    split = coverage_by_layer(items)
    assert sum(len(v) for v in split.values()) == len(items)
    assert [len(split[layer]) for layer in LAYERS] == [18, 10, 34, 35]


# --------------------------------------------------------------------------- #
# the refusal
# --------------------------------------------------------------------------- #


def test_the_repression_coupling_is_refused_and_raises_on_float():
    """Criterion (b). The repository refuses rather than inventing, and the sentinel has to
    raise, not return a placeholder a caller can read past."""
    assert str(REFUSED_REPRESSION_COUPLING.tag) == "REFUSED"
    assert REFUSED_REPRESSION_COUPLING.value is None
    assert REFUSED_REPRESSION_COUPLING.tag in FREE_TAGS
    with pytest.raises(RefusedValue):
        float(REFUSED_REPRESSION_COUPLING)


def test_the_refusal_names_the_measurement_that_would_close_it():
    assert "PAC-motif" in REFUSED_REPRESSION_COUPLING.missing
    assert "19901341" in REFUSED_REPRESSION_COUPLING.source  # Lippman & Broach 2009
    assert "21730963" in REFUSED_REPRESSION_COUPLING.source  # Huber 2011


def test_repression_is_absent_from_the_whole_panel_which_is_why_it_is_refused():
    """The measured ground for adaptation.transcriptional_repression, re-run here: across all
    25 stressors there are exactly five negative target weights and every one is on a
    metabolite pool, which has no transcription factor and so no regulon to repress."""
    negative = [
        (name, target) for name, spec in STRESSORS.items()
        for target, weight in spec.targets.items() if weight < 0
    ]
    assert len(STRESSORS) == 25
    assert len(negative) == 5
    assert {t for _, t in negative} == {"redox", "nadh", "atp", "ph"}
    assert "DOT6" not in MODULES and "TOD6" not in MODULES


# --------------------------------------------------------------------------- #
# against the real model
# --------------------------------------------------------------------------- #


@pytest.mark.integration
def test_every_claim_in_the_table_survives_the_real_model(base_model, restored_model, items):
    """The whole table, re-checked. Empty is the pass, and the message names the row."""
    assert verify_coverage(base_model, items, restored=restored_model) == ()


@pytest.mark.integration
def test_cobrapy_leaves_every_gene_name_empty(base_model):
    """The trap the audit is built to avoid. Not one of the 1161 genes carries a symbol, so a
    naive symbol sweep returns nothing and the entire diagram reads as absent."""
    assert len(base_model.genes) == 1161
    assert [g.id for g in base_model.genes if g.name] == []
    assert not base_model.genes.has_id("GPX3")
    assert not base_model.genes.has_id("RHO1")


@pytest.mark.integration
def test_the_two_genes_a_symbol_sweep_wrongly_calls_absent(base_model):
    """GPX3 is in the model as HYR1, and RHO1 is a required conjunct of every clause of the
    glucan synthase. Both are found only by ORF."""
    assert base_model.genes.has_id("YIR037W")
    assert {r.id for r in base_model.genes.get_by_id("YIR037W").reactions} == {"r_0483"}
    assert base_model.genes.has_id("YPR165W")
    rule = base_model.reactions.get_by_id("r_0005").gene_reaction_rule
    assert rule.count("YPR165W") == len(rule.split(" or "))


@pytest.mark.integration
def test_four_declared_panel_modules_reach_no_reaction_at_all(base_model):
    """TRIPWIRE, and the audit's headline finding. If a later repair to bridge/regulation.py
    gives any of these a reaction, this test SHOULD fail: the table's basis columns for UPR,
    calcium, zinc and xenobiotic are then stale and the audit needs re-running."""
    regulation = pytest.importorskip("ystwin.bridge.regulation")
    try:
        regulons = regulation.load_regulons()
    except Exception as exc:  # noqa: BLE001 -- absent SGD data is a skip, not a failure
        pytest.skip(f"SGD regulation records unavailable: {exc}")
    reach = {}
    for name, regulon in regulons.items():
        present = [g for g in regulon.genes if base_model.genes.has_id(g)]
        reach[name] = len({r.id for g in present for r in base_model.genes.get_by_id(g).reactions})
    assert set(UNREACHED_MODULES) == {"UPR", "calcium", "zinc", "xenobiotic"}
    for name in UNREACHED_MODULES:
        assert reach[name] == 0, f"{name} now reaches {reach[name]} reactions; re-run the audit"
    assert reach["cell_wall"] == 1 and reach["alkaline_ph"] == 1
    assert len(regulons["calcium"].genes) == 0, "CRZ1 has no target at this evidence level"


@pytest.mark.integration
def test_the_curated_map_gives_the_unreached_modules_reactions(base_model, restored_model):
    """What the table is FOR. The regulon route reaches nothing for these four; the declared
    map reaches the effectors that were there all along."""
    mapped = panel_module_reactions(base_model)
    assert {m: len(mapped[m]) for m in UNREACHED_MODULES} == {
        "UPR": 3, "calcium": 4, "zinc": 3, "xenobiotic": 6,
    }
    assert mapped["cell_wall"] == frozenset({"r_0005", "r_0006", "r_0272", "r_0362"})
    assert len(panel_module_reactions(restored_model)["cell_wall"]) == 5
    for reactions in mapped.values():
        assert all(base_model.reactions.has_id(r) for r in reactions)


@pytest.mark.integration
def test_installing_the_groups_changes_no_optimum(base_model, items):
    """The declaration is made addressable and nothing else. No bound, no rule, no
    stoichiometry, and the same optimum to the solver's tolerance."""
    reference = base_model.copy()
    grouped = apply_pathway_coverage(base_model, items)
    configure(reference)
    configure(grouped)
    before = reference.slim_optimize()
    assert grouped.slim_optimize() == pytest.approx(before, abs=1e-9)
    assert before == pytest.approx(SHIPPED_DEFAULT_GROWTH, abs=1e-6)
    assert len(grouped.reactions) == len(base_model.reactions)
    assert len(grouped.genes) == len(base_model.genes)
    assert len(grouped.metabolites) == len(base_model.metabolites)
    installed = {g.id for g in grouped.groups} - {g.id for g in base_model.groups}
    assert len(installed) == 31
    assert installed == {group_name(r) for r in items
                         if any(base_model.reactions.has_id(x) for x in r.reaction_ids)}
    with pytest.raises(ValueError, match="already installed"):
        apply_pathway_coverage(grouped, items)


@pytest.mark.integration
def test_a_group_answers_which_reactions_implement_an_adaptation(base_model, items):
    grouped = apply_pathway_coverage(base_model, items)
    group = next(g for g in grouped.groups if g.id == group_name(_row(items, "adaptation.trehalose")))
    assert {m.id for m in group.members} == {"r_0195", "r_1051", "r_0194", "r_0193"}


# --------------------------------------------------------------------------- #
# criterion (a): the ablation
# --------------------------------------------------------------------------- #


@pytest.mark.integration
def test_undosed_every_row_is_either_lethal_or_below_the_floor(base_model, items):
    """Criterion (a), the null result, and it is the honest reading of this table. Thirty-one
    rows carry reactions in the SHIPPED model; nineteen are lethal, and the other twelve move
    growth by at most 0.000086 /h = 0.0074x the 0.0117 /h growth SE. Nothing in between."""
    scored = pathway_essentiality(base_model, items)
    assert len(scored) == 31
    lethal = {k for k, v in scored.items() if v > SHIPPED_DEFAULT_GROWTH - 1e-6}
    assert len(lethal) == 19
    assert "adaptation.manganese_homeostasis" in lethal, (
        "r_4590 carries SMF1 in the shipped model, so manganese is scored before the restore")
    assert len([k for k in lethal if k.startswith("adaptation.")]) == 18
    assert "cascade.CWI_RHO1" in lethal
    survivors = {k: v for k, v in scored.items() if k not in lethal}
    assert max(survivors.values()) / MEASURED_GROWTH_RATE_SE < 0.01
    assert sum(1 for v in survivors.values() if abs(v) < 1e-7) == 9


@pytest.mark.integration
def test_lethality_is_regime_independent_which_is_what_a_biomass_requirement_looks_like(
    base_model, items
):
    """The reason none of the eighteen counts as a stress-response effect: a stress response
    binds under stress, and these bind identically under both regimes because they are fixed
    biomass components. Chitin and the two glucans are literal coefficients of r_4048."""
    shipped = pathway_essentiality(base_model, items)
    with base_model as scoped:
        physiology.aerobic_batch_constraints(scoped)
        assert scoped.slim_optimize() == pytest.approx(REFERENCE_BATCH_GROWTH, abs=1e-6)
        batch = pathway_essentiality(scoped, items)
    lethal_shipped = {k for k, v in shipped.items() if v > SHIPPED_DEFAULT_GROWTH - 1e-6}
    lethal_batch = {k for k, v in batch.items() if v > REFERENCE_BATCH_GROWTH - 1e-6}
    assert lethal_shipped == lethal_batch
    carbohydrate = base_model.reactions.get_by_id("r_4048")
    coefficients = {m.name: -c for m, c in carbohydrate.metabolites.items() if c < 0}
    assert coefficients["chitin"] == pytest.approx(0.02361, abs=1e-5)
    assert coefficients["(1->3)-beta-D-glucan"] == pytest.approx(0.73914, abs=1e-5)
    assert coefficients["(1->6)-beta-D-glucan"] == pytest.approx(0.24696, abs=1e-5)
    assert coefficients["trehalose"] == pytest.approx(0.13655, abs=1e-5)
    assert coefficients["glycogen"] == pytest.approx(0.35689, abs=1e-5)


@pytest.mark.integration
def test_catalase_clears_the_growth_floor_only_once_peroxide_is_dosed(base_model, items):
    """CRITERION (a), and the one row in this table that passes it. Undosed, blocking catalase
    moves growth by 1.1e-6 /h. Under stress_ros.dose_peroxide at 1.0 mmol/gDCW/h it moves
    0.017238 /h = 1.47x the 0.0117 /h growth SE, and at 5.0 it moves 0.086187 /h = 7.37x.

    The conditionality is the finding: a stress pathway is invisible until the stressor it
    answers is applied, which is how every other row in this table should be tested.
    """
    model = stress_ros.add_ros_module(base_model)
    configure(model)
    with model:
        physiology.aerobic_batch_constraints(model)
        deltas = {}
        for dose in (0.0, 1.0, 5.0):
            with model as scoped:
                stress_ros.dose_peroxide(scoped, dose)
                deltas[dose] = pathway_essentiality(scoped, items)[
                    "adaptation.antioxidant_catalase"]
    assert deltas[0.0] < MEASURED_GROWTH_RATE_SE
    assert deltas[1.0] == pytest.approx(0.017238, abs=1e-5)
    assert deltas[1.0] / MEASURED_GROWTH_RATE_SE == pytest.approx(1.47, abs=0.01)
    assert deltas[5.0] == pytest.approx(0.086187, abs=1e-5)
    assert deltas[5.0] / MEASURED_GROWTH_RATE_SE == pytest.approx(7.37, abs=0.01)


# --------------------------------------------------------------------------- #
# the audit's own false positive, kept as a test so it cannot come back
# --------------------------------------------------------------------------- #


@pytest.mark.integration
def test_the_v_atpase_is_not_an_unpriced_proton_escape(base_model):
    """A first pass called r_1085/r_1086 a second escape that stress_ph fails to close.
    Measured, that is wrong: with the budget closed the price is 1.0000 ATP either way, both
    reactions carry zero flux at the optimum, and blocking them costs 0.000000 /h. The
    stoichiometric difference is real -- they consume a cytosolic proton and r_0227 does not --
    but it costs nothing, so stress_ph's omission is correct.
    """
    stress_ph = pytest.importorskip("ystwin.fba.stress_ph")
    proton_c = base_model.metabolites.get_by_id("s_0794")
    assert proton_c not in base_model.reactions.get_by_id("r_0227").metabolites
    for rid in ("r_1085", "r_1086"):
        assert base_model.reactions.get_by_id(rid).metabolites[proton_c] == -1.0
    with base_model as model:
        physiology.aerobic_batch_constraints(model)
        stress_ph.constrain_proton_budget(model)
        open_price = stress_ph.proton_price(model)
        open_growth = model.slim_optimize()
        fluxes = model.optimize().fluxes
        assert fluxes["r_1085"] == pytest.approx(0.0, abs=1e-9)
        assert fluxes["r_1086"] == pytest.approx(0.0, abs=1e-9)
        with model as blocked:
            for rid in ("r_1085", "r_1086"):
                blocked.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            assert stress_ph.proton_price(blocked) == pytest.approx(open_price, abs=1e-6)
            assert blocked.slim_optimize() == pytest.approx(open_growth, abs=1e-9)
    assert open_price == pytest.approx(1.0, abs=1e-4)


@pytest.mark.integration
def test_the_camp_second_messenger_is_six_reactions_and_no_module_reads_it(base_model, items):
    """Finding 3. The chemistry is fully present and inert; the gap is at the panel layer."""
    record = _row(items, "cascade.cAMP_second_messenger")
    assert record.reaction_ids == ("r_0147", "r_0033", "r_0034", "r_0035", "r_0036", "r_0037")
    assert base_model.reactions.get_by_id("r_0147").name == "adenylate cyclase"
    for rid in ("r_0033", "r_0034", "r_0035", "r_0036", "r_0037"):
        assert "phosphodiesterase" in base_model.reactions.get_by_id(rid).name
    carried = " ".join(m.transcription_factor for m in MODULES.values())
    assert "cAMP" not in carried and "PKA" not in carried and "Tpk" not in carried
    assert pathway_essentiality(base_model, (record,))[record.item] < MEASURED_GROWTH_RATE_SE


@pytest.mark.integration
def test_glycophagy_is_the_only_autophagy_reaction_and_has_an_empty_rule(base_model):
    """Why adaptation.autophagy_bulk_protein is REFUSED rather than GAP: there is no reaction
    to attach an ATG gene to, and the one autophagy reaction that exists has no GPR at all."""
    assert base_model.reactions.get_by_id("r_1174").gene_reaction_rule == ""
    autophagy = [r.id for r in base_model.reactions if "autophagy" in r.name.lower()]
    assert autophagy == ["r_1174"]
    for orf in ("YGL180W", "YBL078C", "YPR185W"):
        assert not base_model.genes.has_id(orf)


@pytest.mark.integration
def test_pathway_essentiality_omits_rows_the_model_lacks_rather_than_scoring_them_zero(
    base_model, items
):
    """A missing row reported as 0.0 would read as a measured null, which is the failure mode
    this audit is trying to remove from the repository, not add to it."""
    scored = pathway_essentiality(base_model, items)
    addressable = {r.item for r in items
                   if any(base_model.reactions.has_id(x) for x in r.reaction_ids)}
    assert set(scored) == addressable
    assert all(r.item not in scored for r in items if not r.reaction_ids)
    absent_only = cobra.Model("absent_only")
    assert pathway_essentiality(absent_only, items) == {}


@pytest.mark.integration
def test_pathway_essentiality_does_not_modify_the_model(base_model, items):
    before = {r.id: r.bounds for r in base_model.reactions}
    pathway_essentiality(base_model, items)
    assert {r.id: r.bounds for r in base_model.reactions} == before


def test_apply_pathway_coverage_skips_reactions_the_model_lacks():
    """Safe on any model, which is what lets the same table serve the shipped and the
    restored one."""
    empty = cobra.Model("empty")
    assert len(apply_pathway_coverage(empty).groups) == 0
