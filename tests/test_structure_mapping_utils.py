from pathlib import Path


def test_primary_covariance_model_mapping_imports_its_helpers():
    from src.map_rfam_cm_to_pdb import FAMILIES, ROOT

    assert set(FAMILIES) == {"RF00001", "RF00005"}
    assert (ROOT / "data" / "raw" / "rf00001_seed.sto").is_file()
    assert (ROOT / "data" / "raw" / "rf00005_seed.sto").is_file()
