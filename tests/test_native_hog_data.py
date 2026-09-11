from __future__ import annotations

import shutil

import pytest

from ystwin import paths
from ystwin.analysis.hog_data import load_native_hog_data


@pytest.fixture
def native_source(tmp_path):
    source = paths.data_dir() / "hog2013"
    for filename in ("sources.json", "observations.xls", "methods.pdf"):
        shutil.copyfile(source / filename, tmp_path / filename)
    return tmp_path


def test_native_reader_never_requires_or_opens_preliminary_data(native_source):
    data = load_native_hog_data(native_source)
    assert data.observations.supplement_id.eq("s002").all()
    assert len(data.observations) == 317
    assert set(data.metadata["sources"]) == {"observations", "methods"}
    assert not (native_source / "western_blots.xls").exists()
    wild_type = data.observations.loc[data.observations.genotype.eq("wild_type")]
    assert len(wild_type.loc[wild_type.nacl_molar.eq(0.4)]) == 54
    assert len(wild_type.loc[wild_type.nacl_molar.eq(0.0)]) == 28


def test_native_reader_preserves_original_roles_and_provenance(native_source):
    data = load_native_hog_data(native_source)
    assert data.observations["split"].eq("fit").sum() == 13
    assert data.observations.source_cell.notna().all()
    assert data.observations.unit.eq("mol/L").all()
    assert data.metadata["native_only"] is True
    assert data.observations.loc[data.observations.genotype.eq("gpd1_del")
                                 & data.observations.observable_id.eq("glycerol_measured")].empty


def test_native_reader_rejects_changed_training_bytes(native_source):
    path = native_source / "observations.xls"
    path.write_bytes(path.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="SHA-256"):
        load_native_hog_data(native_source)
