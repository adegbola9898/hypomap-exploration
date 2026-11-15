from copy import deepcopy

import h5py
import numpy as np
import pytest

from hypomap_exploration.extraction import extract_atlas_counts


@pytest.fixture
def atlas(tmp_path):
    path = tmp_path / "atlas.h5ad"
    text_dtype = h5py.string_dtype("utf-8")

    with h5py.File(path, "w") as f:
        matrix = f.create_group("raw/X")
        matrix.attrs["encoding-type"] = "csr_matrix"
        matrix.attrs["shape"] = [3, 3]

        # Columns: MC4R, OTHER, DENND1B
        # Rows:
        # cell-z: [2, 9, 1]
        # cell-a: [0, 0, 3]
        # cell-m: [4, 0, 0]
        matrix.create_dataset("data", data=[2., 9., 1., 3., 4.])
        matrix.create_dataset("indices", data=[0, 1, 2, 2, 0])
        matrix.create_dataset("indptr", data=[0, 3, 4, 5])

        var = f.create_group("raw/var")
        var.create_dataset(
            "_index",
            data=["MC4R", "OTHER", "DENND1B"],
            dtype=text_dtype,
        )

        obs = f.create_group("obs")
        obs.create_dataset(
            "_index", data=["cell-z", "cell-a", "cell-m"], dtype=text_dtype
        )
        obs.create_dataset(
            "C0_named", data=["Neuron", "Glia", "Neuron"], dtype=text_dtype
        )
        obs.create_dataset(
            "Donor_ID", data=["D2", "D1", "D2"], dtype=text_dtype
        )
        obs.create_dataset("nCount_RNA", data=[12., np.nan, 4.])

    config = {
        "matrix_path": "raw/X",
        "feature_index_path": "raw/var/_index",
        "cell_index_path": "obs/_index",
        "genes": ["DENND1B", "MC4R"],
        "selection_field": "C0_named",
        "selection_value": "Neuron",
        "extraction": {
            "population": "all_nuclei",
            "expected_nuclei": 3,
            "preserve_obs_fields": ["C0_named", "Donor_ID", "nCount_RNA"],
        },
    }
    return path, config


def test_preserves_all_rows_metadata_and_gene_order(atlas):
    path, config = atlas
    original_config = deepcopy(config)

    with h5py.File(path, "r") as f:
        table, mapping = extract_atlas_counts(
            f, config, max_nnz_per_chunk=2
        )

    assert table.columns.tolist() == [
        "cell_id", "source_row", "C0_named", "Donor_ID", "nCount_RNA",
        "DENND1B_count", "MC4R_count",
    ]
    assert table["cell_id"].tolist() == ["cell-z", "cell-a", "cell-m"]
    assert table["source_row"].tolist() == [0, 1, 2]
    assert table["C0_named"].tolist() == ["Neuron", "Glia", "Neuron"]
    assert table["Donor_ID"].tolist() == ["D2", "D1", "D2"]
    assert np.isnan(table.loc[1, "nCount_RNA"])
    np.testing.assert_array_equal(
        table[["DENND1B_count", "MC4R_count"]].to_numpy(),
        [[1, 2], [3, 0], [0, 4]],
    )
    assert mapping == {"DENND1B": 2, "MC4R": 0}
    assert config == original_config


@pytest.mark.parametrize("problem", ["missing", "ambiguous"])
def test_rejects_missing_or_ambiguous_gene(atlas, problem):
    path, config = atlas
    with h5py.File(path, "r+") as f:
        if problem == "missing":
            f["raw/var/_index"][2] = "OTHER2"
        else:
            f["raw/var/_index"][1] = "DENND1B"

    with h5py.File(path, "r") as f:
        with pytest.raises(ValueError, match="gene"):
            extract_atlas_counts(f, config, max_nnz_per_chunk=2)


def test_rejects_duplicate_cell_ids(atlas):
    path, config = atlas
    with h5py.File(path, "r+") as f:
        f["obs/_index"][1] = "cell-z"

    with h5py.File(path, "r") as f:
        with pytest.raises(ValueError, match="cell"):
            extract_atlas_counts(f, config, max_nnz_per_chunk=2)


def test_rejects_misaligned_metadata(atlas):
    path, config = atlas
    with h5py.File(path, "r+") as f:
        del f["obs/nCount_RNA"]
        f["obs"].create_dataset("nCount_RNA", data=[12., 3.])

    with h5py.File(path, "r") as f:
        with pytest.raises(ValueError, match="metadata"):
            extract_atlas_counts(f, config, max_nnz_per_chunk=2)


@pytest.mark.parametrize("genes", [[], ["MC4R", "MC4R"], [""]])
def test_rejects_invalid_gene_requests(atlas, genes):
    path, config = atlas
    config["genes"] = genes
    with h5py.File(path, "r") as f:
        with pytest.raises(ValueError, match="gene"):
            extract_atlas_counts(f, config, max_nnz_per_chunk=2)


def test_rejects_unexpected_nucleus_count(atlas):
    path, config = atlas
    config["extraction"]["expected_nuclei"] = 4
    with h5py.File(path, "r") as f:
        with pytest.raises(ValueError, match="expected_nuclei"):
            extract_atlas_counts(f, config, max_nnz_per_chunk=2)


def test_rejects_missing_metadata(atlas):
    path, config = atlas
    config["extraction"]["preserve_obs_fields"].append("Missing_Field")
    with h5py.File(path, "r") as f:
        with pytest.raises(ValueError, match="metadata"):
            extract_atlas_counts(f, config, max_nnz_per_chunk=2)


def test_rejects_metadata_column_collision(atlas):
    path, config = atlas
    config["extraction"]["preserve_obs_fields"].append("MC4R_count")
    with h5py.File(path, "r") as f:
        with pytest.raises(ValueError, match="collide"):
            extract_atlas_counts(f, config, max_nnz_per_chunk=2)


@pytest.mark.parametrize("layout", ["legacy", "group"])
def test_rejects_undecoded_categorical_metadata(atlas, layout):
    path, config = atlas
    with h5py.File(path, "r+") as f:
        del f["obs/Donor_ID"]

        if layout == "legacy":
            f["obs"].create_dataset("Donor_ID", data=[1, 0, 1])
            categories = f["obs"].create_group("__categories")
        else:
            categories = f["obs"].create_group("Donor_ID")
            categories.attrs["encoding-type"] = "categorical"
            categories.create_dataset("codes", data=[1, 0, 1])

        categories.create_dataset(
            "Donor_ID" if layout == "legacy" else "categories",
            data=["D1", "D2"],
            dtype=h5py.string_dtype("utf-8"),
        )

    with h5py.File(path, "r") as f:
        with pytest.raises(ValueError, match="categorical"):
            extract_atlas_counts(f, config, max_nnz_per_chunk=2)


@pytest.mark.parametrize("axis", ["features", "cells"])
def test_rejects_identifier_dimension_mismatch(atlas, axis):
    path, config = atlas
    target = (
        config["feature_index_path"]
        if axis == "features"
        else config["cell_index_path"]
    )
    with h5py.File(path, "r+") as f:
        del f[target]
        f.create_dataset(
            target,
            data=["MC4R", "DENND1B"] if axis == "features" else ["a", "b"],
            dtype=h5py.string_dtype("utf-8"),
        )

    with h5py.File(path, "r") as f:
        with pytest.raises(ValueError, match="gene|cell"):
            extract_atlas_counts(f, config, max_nnz_per_chunk=2)
