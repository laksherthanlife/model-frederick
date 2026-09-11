"""The published-strain survey has a reader now, and these are the claims it can make.

The file had none. Five module docstrings made categorical claims about
`data/carotenoid_batch/published_batch_titres.tsv` -- "every cassette is qualitative", "every
carotenoid row is glucose or galactose" -- and `grep -rn published_batch_titres` returned two
prose mentions and nothing else. Every one had gone false, and the file itself is what
falsified them: it gained `copies_*` and `copy_evidence` columns on 2026-09-03, three days
after the oldest of those claims was written.

So the tests here are not about the numbers. They are about the two properties that stop the
claims coming back: the table has a DECLARED DOMAIN a violation of which fails at load, and
every claim about it is DERIVED at the point of use rather than copied into a docstring.
"""

from __future__ import annotations

import csv
import pathlib

import pytest

from ystwin import paths
from ystwin.pathway.published_cassettes import (
    CARBON_SPECIES,
    EPISOMAL_2U,
    HIGH_COPY,
    UNRECORDED,
    CassetteDomainError,
    crtyb_dosage_unfittable_reason,
    numeric_crtyb_dosages,
    published_cassettes,
)

_TABLE = paths.data_dir() / "carotenoid_batch" / "published_batch_titres.tsv"

pytestmark = pytest.mark.skipif(not _TABLE.is_file(), reason="the survey is not vendored")


class TestEveryRowLoadsUnderTheDeclaredDomain:
    """The nag. A row nobody has looked at closely fails here rather than downstream."""

    def test_the_whole_table_loads(self):
        assert len(published_cassettes()) == 29

    def test_no_column_holds_a_free_text_non_value(self):
        """The specific defect this domain was written for: `promoter` held the string
        "not recorded here" in a column otherwise holding promoter names, contradicted by its
        own row's `copy_evidence`. A reader cannot tell that from a promoter."""
        for cassette in published_cassettes():
            for value in (cassette.promoter, *cassette.copies.values()):
                if isinstance(value, str):
                    assert "not " not in value.lower(), (
                        f"{cassette.reference}: {value!r} is a sentence, not a value")

    def test_a_dosage_is_an_integer_or_one_of_the_declared_literals(self):
        for cassette in published_cassettes():
            for gene, value in cassette.copies.items():
                assert isinstance(value, int) or value in (
                    EPISOMAL_2U, HIGH_COPY, UNRECORDED), f"{cassette.reference} {gene}"

    def test_an_episomal_dosage_is_not_silently_a_number(self):
        """`2u` is a distribution over a population, not a copy number. Converting it would
        manufacture the quantitative dosage the survey is being asked whether it has."""
        episomal = [c for c in published_cassettes()
                    if EPISOMAL_2U in c.copies.values()]

        assert episomal
        for cassette in episomal:
            assert not cassette.states_a_numeric_crtyb_dosage or \
                cassette.copies["crtYB"] != EPISOMAL_2U


class TestThePopulationRuleIsEnforcedRatherThanObserved:
    """`promoter` is filled exactly on `dfba_ready = yes` rows. A real convention that was
    written down nowhere until `data/carotenoid_batch/SOURCE.md`."""

    def test_promoter_is_populated_exactly_on_the_dfba_ready_rows(self):
        for cassette in published_cassettes():
            assert (cassette.promoter is not None) == cassette.dfba_ready, cassette.reference

    def test_the_rule_holds_in_both_directions_and_is_not_vacuous(self):
        ready = [c for c in published_cassettes() if c.dfba_ready]

        assert 0 < len(ready) < len(published_cassettes())

    def test_a_promoter_on_a_non_ready_row_is_refused(self, tmp_path):
        """The rule has teeth, or it is a description rather than a convention."""
        rows = _survey_rows()
        broken = next(r for r in rows if r["dfba_ready"] == "no")
        broken["promoter"] = "TDH3p"
        assert _reload(tmp_path, rows) is not None


class TestTheCarbonClaimIsDerivedAndSharper:
    """What replaced "every carotenoid row is glucose or galactose", which was false."""

    def test_no_row_has_ethanol_in_the_feed(self):
        """The reason `pathway/environment_flux.py` refuses beta-carotene. The fitted
        descriptor is the feed's ETHANOL carbon fraction, so z = 0 on every row here and the
        coefficient is unidentified by construction -- a reason that survives a xylose row
        being added, which the old enumeration did not."""
        with_ethanol = [c.reference for c in published_cassettes()
                        if "ethanol" in c.carbon_species]

        assert with_ethanol == []

    def test_ethanol_is_in_the_vocabulary_so_the_absence_is_a_measurement(self):
        """Otherwise the test above would pass because nothing could ever match."""
        assert "ethanol" in CARBON_SPECIES

    def test_the_old_enumeration_was_already_false(self):
        """Recorded as an assertion so the retracted sentence cannot be written again."""
        beyond = {species for cassette in published_cassettes()
                  for species in cassette.carbon_species} - {
                      "glucose", "galactose", "undeclared_rich_medium"}

        assert beyond, "the old 'glucose or galactose' claim would be true again"
        assert {"sucrose", "xylose"} <= beyond

    def test_a_medium_naming_no_known_species_is_refused(self, tmp_path):
        rows = _survey_rows()
        rows[0]["carbon_source"] = "grown on something nobody declared"

        assert _reload(tmp_path, rows) is not None


class TestTheRefusalReasonIsComputedNotCopied:
    def test_it_names_the_papers_that_state_a_number(self):
        reason = crtyb_dosage_unfittable_reason()

        for paper in ("Verwaal 2007", "Xie 2014"):
            assert paper in reason

    def test_it_says_confounding_rather_than_absent_numbers(self):
        assert "CONFOUNDING" in crtyb_dosage_unfittable_reason()

    def test_it_moves_when_the_table_does(self, tmp_path):
        """The property the old prose lacked. Adding a third dosage must change what the
        reason says, or it is a literal with extra steps."""
        before = crtyb_dosage_unfittable_reason()
        rows = _survey_rows()
        third = next(r for r in rows if r["copies_crtYB"] == "1")
        third["copies_crtYB"] = "5"

        assert _reload(tmp_path, rows, expect_error=False) != before

    def test_the_dosage_map_carries_the_papers_and_not_just_the_counts(self):
        for dosage, papers in numeric_crtyb_dosages().items():
            assert papers, dosage


def _survey_rows() -> list[dict]:
    """The survey as raw dicts, with the handle closed. Kept in one place so no test leaks a
    file object -- an unclosed handle surfaces as `PytestUnraisableExceptionWarning`, which
    this suite turns into an error."""
    with _TABLE.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _reload(tmp_path: pathlib.Path, rows: list[dict], expect_error: bool = True):
    """Rewrite the survey into ``tmp_path``, re-read it, and return the error or the reason.

    The loader is cached on the vendored path, so a mutation test has to point the module at
    a different file. Done by monkeying the module-level path and clearing the cache, which
    is restored by the fixture's teardown of ``tmp_path`` and by the explicit reset below.
    """
    import ystwin.pathway.published_cassettes as module

    path = tmp_path / "published_batch_titres.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    original = module._TABLE
    module._TABLE = path
    module.published_cassettes.cache_clear()
    try:
        if expect_error:
            with pytest.raises(CassetteDomainError) as raised:
                module.published_cassettes()
            return str(raised.value)
        return module.crtyb_dosage_unfittable_reason()
    finally:
        module._TABLE = original
        module.published_cassettes.cache_clear()
