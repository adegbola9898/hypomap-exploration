import json

import pytest

from hypomap_exploration.run_extraction import run_extraction


@pytest.fixture
def project(tmp_path):
    source = tmp_path / "data/raw/atlas.h5ad"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"Not an HDF5 file: extraction must not reach this.")

    config = tmp_path / "config/analysis.json"
    config.parent.mkdir(parents=True)
    config.write_text(json.dumps({
        "source_file": "data/raw/atlas.h5ad",
        "source_sha256": "0" * 64,
    }))

    return tmp_path, config


def test_checksum_mismatch_creates_no_final_outputs(project):
    root, config = project

    with pytest.raises(ValueError, match="checksum"):
        run_extraction(config, root, max_nnz_per_chunk=2)

    assert not (root / "data/derived/all_nuclei_gene_counts.tsv.gz").exists()
    assert not (root / "provenance/extraction_manifest.json").exists()


@pytest.mark.parametrize("existing_path", [
    "data/derived/all_nuclei_gene_counts.tsv.gz",
    "provenance/extraction_manifest.json",
])
def test_existing_artifact_is_preserved(project, existing_path):
    root, config = project
    existing = root / existing_path
    existing.parent.mkdir(parents=True, exist_ok=True)
    original = b"Existing artifact must remain unchanged."
    existing.write_bytes(original)

    with pytest.raises(FileExistsError):
        run_extraction(config, root, max_nnz_per_chunk=2)

    assert existing.read_bytes() == original


@pytest.fixture
def valid_project(tmp_path):
    import hashlib
    import h5py
    import numpy as np

    source = tmp_path / "data/raw/atlas.h5ad"
    source.parent.mkdir(parents=True)

    with h5py.File(source, "w") as f:
        matrix = f.create_group("raw/X")
        matrix.attrs["encoding-type"] = "csr_matrix"
        matrix.attrs["shape"] = [2, 2]
        matrix.create_dataset("data", data=[1., 2., 3.])
        matrix.create_dataset("indices", data=[0, 1, 1])
        matrix.create_dataset("indptr", data=[0, 2, 3])

        strings = h5py.string_dtype("utf-8")
        var = f.create_group("raw/var")
        var.create_dataset(
            "_index", data=["DENND1B", "MC4R"], dtype=strings
        )
        obs = f.create_group("obs")
        obs.create_dataset("_index", data=["cell-b", "cell-a"], dtype=strings)
        obs.create_dataset("region", data=["NA", "Vent"], dtype=strings)
        obs.create_dataset("nCount_RNA", data=[3., np.nan])

    config = {
        "source_file": "data/raw/atlas.h5ad",
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "matrix_path": "raw/X",
        "feature_index_path": "raw/var/_index",
        "cell_index_path": "obs/_index",
        "genes": ["DENND1B", "MC4R"],
        "extraction": {
            "population": "all_nuclei",
            "expected_nuclei": 2,
            "preserve_obs_fields": ["region", "nCount_RNA"],
        },
    }
    config_path = tmp_path / "config/analysis.json"
    config_path.parent.mkdir()
    config_path.write_text(json.dumps(config))
    return tmp_path, config_path


def test_success_writes_table_and_matching_manifest(valid_project):
    import hashlib
    import pandas as pd

    root, config_path = valid_project
    source = root / "data/raw/atlas.h5ad"
    source_before = source.read_bytes()

    returned = run_extraction(config_path, root, max_nnz_per_chunk=1)

    output = root / "data/derived/all_nuclei_gene_counts.tsv.gz"
    manifest_path = root / "provenance/extraction_manifest.json"
    saved = json.loads(manifest_path.read_text())

    assert returned == saved
    assert saved["status"] == "complete"
    assert saved["source_sha256"] == hashlib.sha256(source_before).hexdigest()
    assert saved["config_sha256"] == hashlib.sha256(
        config_path.read_bytes()
    ).hexdigest()
    assert saved["output_sha256"] == hashlib.sha256(
        output.read_bytes()
    ).hexdigest()
    assert saved["n_nuclei"] == 2
    assert saved["gene_columns"] == {"DENND1B": 0, "MC4R": 1}
    assert saved["missing_metadata_counts"] == {
        "region": 0, "nCount_RNA": 1
    }
    assert saved["table_read_options"] == {
        "sep": "\t",
        "keep_default_na": False,
        "na_values": ["__MISSING__"],
    }
    assert source.read_bytes() == source_before

    table = pd.read_csv(output, **saved["table_read_options"])
    assert table["cell_id"].tolist() == ["cell-b", "cell-a"]
    assert table["region"].tolist() == ["NA", "Vent"]
    assert pd.isna(table.loc[1, "nCount_RNA"])
    assert table[["DENND1B_count", "MC4R_count"]].values.tolist() == [
        [1, 2], [0, 3]
    ]

    for directory, expected in [
        (root / "data/derived", {output.name}),
        (root / "provenance", {manifest_path.name}),
    ]:
        assert {p.name for p in directory.iterdir()} == expected


def test_partial_table_write_is_cleaned_up(valid_project, monkeypatch):
    import pandas as pd

    root, config_path = valid_project

    def broken_write(self, destination, *args, **kwargs):
        from pathlib import Path

        Path(destination).write_bytes(b"incomplete output")
        raise OSError("simulated disk failure")

    monkeypatch.setattr(pd.DataFrame, "to_csv", broken_write)

    with pytest.raises(OSError, match="simulated disk failure"):
        run_extraction(config_path, root, max_nnz_per_chunk=1)

    for directory in [root / "data/derived", root / "provenance"]:
        assert not directory.exists() or not list(directory.iterdir())


def test_manifest_publication_failure_removes_new_table(
    valid_project, monkeypatch
):
    import importlib

    module = importlib.import_module("hypomap_exploration.run_extraction")
    root, config_path = valid_project
    original_publish = module._publish_without_overwrite
    published = []

    def fail_on_manifest(temporary, destination):
        if destination.name == "extraction_manifest.json":
            assert (
                root / "data/derived/all_nuclei_gene_counts.tsv.gz"
            ).exists()
            raise OSError("simulated manifest publication failure")
        original_publish(temporary, destination)
        published.append(destination.name)

    monkeypatch.setattr(
        module, "_publish_without_overwrite", fail_on_manifest
    )

    with pytest.raises(OSError, match="manifest publication failure"):
        run_extraction(config_path, root, max_nnz_per_chunk=1)

    assert published == ["all_nuclei_gene_counts.tsv.gz"]

    for directory in [root / "data/derived", root / "provenance"]:
        assert not list(directory.iterdir())
