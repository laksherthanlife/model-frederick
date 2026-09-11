"""The lowercased Bio-Rad export, and whether the repair leaves a workbook a reader opens.

``scripts/repair_qpcr_export.py`` is the only reason replicate 1 of the qPCR anchor exists
in readable form: the instrument wrote its archive entries in lowercase and OOXML readers
match those case-sensitively. So the test that matters is not that the rename table looks
right, it is that openpyxl opens the repaired file and refuses the original -- a table that
renames five of six entries produces a workbook that still cannot be opened, and the names
alone do not show it.

Every unrecognised entry has to survive untouched. The rename table names six paths and a
real export carries thirty, so a rule that dropped or renamed what it did not recognise
would corrupt the archive while every test on the six passed.
"""

from __future__ import annotations

import importlib.util
import io
import pathlib
import zipfile

import openpyxl
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "repair_qpcr_export.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("repair_qpcr_export", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def lowercased(tmp_path):
    """A workbook rewritten the way the instrument wrote the 2026-07-24 export.

    Every archive name lowercased and nothing else touched. The reference rewriting inside
    the XML is exercised separately, because openpyxl writes inline strings and so emits no
    shared-strings part to refer to.
    """
    original = tmp_path / "original.xlsx"
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.title = "Quantification Cq Results"
    sheet.append(["Well", "Target", "Cq"])
    sheet.append(["A1", "UBC", 21.5])
    book.save(original)

    broken = tmp_path / "lowercased.xlsx"
    with zipfile.ZipFile(original) as source, \
            zipfile.ZipFile(broken, "w", zipfile.ZIP_DEFLATED) as destination:
        for item in source.infolist():
            destination.writestr(item.filename.lower(), source.read(item.filename))
    return original, broken


class TestTheRenameTableIsAppliedWhereItApplies:
    @pytest.mark.parametrize("lowered,expected", [
        ("[content_types].xml", "[Content_Types].xml"),
        ("xl/sharedstrings.xml", "xl/sharedStrings.xml"),
        ("xl/_rels/workbook.xml.rels", "xl/_rels/workbook.xml.rels"),
        ("xl/worksheets/sheet1.xml", "xl/worksheets/sheet1.xml"),
    ])
    def test_a_lowercased_name_comes_back_in_its_canonical_case(self, script, lowered,
                                                               expected):
        assert script.canonical_name(lowered) == expected

    def test_a_name_already_in_canonical_case_is_left_alone(self, script):
        """The repair has to be safe to run twice: an export that is already correct must
        come out byte-identical rather than mangled by a second pass."""
        for name in script.CANONICAL.values():
            assert script.canonical_name(name) == name

    def test_a_multi_digit_worksheet_keeps_its_own_number(self, script):
        """The regex captures the index. Rebuilding the name from the match group rather
        than the captured digits would collapse sheet12 onto sheet1 and silently point two
        relationships at one sheet."""
        assert script.canonical_name("xl/worksheets/sheet12.xml") == \
            "xl/worksheets/sheet12.xml"

    def test_windows_separators_become_archive_separators(self, script):
        assert script.canonical_name("xl\\worksheets\\sheet1.xml") == \
            "xl/worksheets/sheet1.xml"

    @pytest.mark.parametrize("name", [
        "docProps/core.xml", "xl/theme/theme1.xml",
        "xl/worksheets/_rels/sheet1.xml.rels", "xl/calcChain.xml",
    ])
    def test_an_entry_the_table_does_not_name_passes_through_unchanged(self, script, name):
        """Six paths are named and a real export carries thirty. Anything unrecognised is
        somebody else's part of the archive, and renaming or dropping it would corrupt a
        file whose data is not wrong."""
        assert script.canonical_name(name) == name


class TestTheRepairedFileOpens:
    def test_the_lowercased_original_cannot_be_opened_at_all(self, script, lowercased):
        """The premise. If openpyxl were to open the broken file this script would have no
        reason to exist and every test below would pass without testing anything."""
        _, broken = lowercased

        with pytest.raises(Exception):
            openpyxl.load_workbook(io.BytesIO(broken.read_bytes()))

    def test_and_the_repaired_copy_reads_back_the_values_that_were_in_it(self, script,
                                                                        lowercased,
                                                                        tmp_path):
        """Nothing is wrong with the data, so every cell has to survive the rebuild."""
        _, broken = lowercased

        repaired = script.repair(broken, tmp_path / "out" / "repaired.xlsx")

        book = openpyxl.load_workbook(io.BytesIO(repaired.read_bytes()))
        try:
            sheet = book.active
            title = sheet.title
            rows = [[c.value for c in row] for row in sheet.iter_rows(max_row=2)]
        finally:
            book.close()
        assert title == "Quantification Cq Results"
        assert rows == [["Well", "Target", "Cq"], ["A1", "UBC", 21.5]]

    def test_no_archive_entry_is_lost_on_the_way(self, script, lowercased, tmp_path):
        """A part silently dropped can leave a workbook that still opens but has lost a
        sheet, so the count is checked rather than inferred from the file opening."""
        _, broken = lowercased

        repaired = script.repair(broken, tmp_path / "repaired.xlsx")

        with zipfile.ZipFile(broken) as before, zipfile.ZipFile(repaired) as after:
            assert len(after.namelist()) == len(before.namelist())

    def test_the_destination_directory_is_created_rather_than_required(self, script,
                                                                      lowercased,
                                                                      tmp_path):
        _, broken = lowercased

        repaired = script.repair(broken, tmp_path / "made" / "up" / "path.xlsx")

        assert repaired.exists()


class TestTheReferenceInsideTheXmlIsRewrittenOnlyWhereItIsAReference:
    """The half a rename table cannot do, and the half it must not do.

    The relationship XML names its target in lowercase too, so renaming the entry alone
    leaves an archive whose parts are all correctly named and whose references point at
    nothing. But the same string appearing in a worksheet is cell text, and rewriting that
    would edit the measurements -- which is the one thing this script exists not to do.
    """

    @pytest.fixture
    def rebuilt(self, script, tmp_path):
        source = tmp_path / "referring.xlsx"
        with zipfile.ZipFile(source, "w") as archive:
            archive.writestr("xl/_rels/workbook.xml.rels",
                             b'<Relationship Target="sharedstrings.xml"/>')
            archive.writestr("xl/workbook.xml", b'<workbook ref="sharedstrings.xml"/>')
            archive.writestr("xl/worksheets/sheet1.xml", b"<v>sharedstrings.xml</v>")
        destination = script.repair(source, tmp_path / "rebuilt.xlsx")
        with zipfile.ZipFile(destination) as archive:
            return {name: archive.read(name) for name in archive.namelist()}

    @pytest.mark.parametrize("part", ["xl/_rels/workbook.xml.rels", "xl/workbook.xml"])
    def test_a_part_that_points_at_the_shared_strings_is_rewritten(self, rebuilt, part):
        assert rebuilt[part] == rebuilt[part].replace(b"sharedstrings.xml",
                                                      b"sharedStrings.xml")
        assert b"sharedStrings.xml" in rebuilt[part]

    def test_and_the_same_text_inside_a_worksheet_is_left_exactly_as_it_was(self, rebuilt):
        assert rebuilt["xl/worksheets/sheet1.xml"] == b"<v>sharedstrings.xml</v>"
