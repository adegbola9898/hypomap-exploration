import h5py
import numpy as np
import pytest

from hypomap_exploration.extraction import read_csr_columns


@pytest.fixture
def matrix(tmp_path):
    """Five rows, four columns, including empty rows and duplicates."""
    path = tmp_path / "fixture.h5ad"

    with h5py.File(path, "w") as handle:
        group = handle.create_group("raw/X")
        group.attrs["encoding-type"] = "csr_matrix"
        group.attrs["shape"] = [5, 4]

        # Dense interpretation:
        # [0, 2, 0, 5]
        # [0, 0, 0, 0]
        # [7, 0, 0, 4]  <- duplicate column 3 entries: 1 + 3
        # [0, 0, 0, 0]  <- explicitly stored zero
        # [0, 6, 0, 0]
        group.create_dataset(
            "data", data=np.array([2, 5, 7, 1, 3, 0, 6], dtype=float)
        )
        group.create_dataset(
            "indices", data=np.array([1, 3, 0, 3, 3, 2, 1], dtype=np.int64)
        )
        group.create_dataset(
            "indptr", data=np.array([0, 2, 2, 5, 6, 7], dtype=np.int64)
        )
        yield group


@pytest.mark.parametrize("chunk_size", [1, 2, 3, 100])
def test_requested_order_empty_rows_and_duplicates(matrix, chunk_size):
    actual = read_csr_columns(
        matrix, [3, 1], max_nnz_per_chunk=chunk_size
    )
    expected = np.array([
        [5, 2],
        [0, 0],
        [4, 0],
        [0, 0],
        [0, 6],
    ])
    np.testing.assert_array_equal(actual, expected)


def test_single_column_preserves_two_dimensions(matrix):
    actual = read_csr_columns(matrix, [0], max_nnz_per_chunk=2)
    np.testing.assert_array_equal(actual, [[0], [0], [7], [0], [0]])


def test_zero_column_is_retained(matrix):
    actual = read_csr_columns(matrix, [2], max_nnz_per_chunk=2)
    np.testing.assert_array_equal(actual, np.zeros((5, 1)))


@pytest.mark.parametrize(
    "columns",
    [[], [1, 1], [-1], [4], [1.5], [True]],
)
def test_rejects_invalid_columns(matrix, columns):
    with pytest.raises(ValueError, match="column"):
        read_csr_columns(matrix, columns, max_nnz_per_chunk=2)


@pytest.mark.parametrize("chunk_size", [0, -1, 1.5, True])
def test_rejects_invalid_chunk_size(matrix, chunk_size):
    with pytest.raises(ValueError, match="chunk"):
        read_csr_columns(matrix, [1], max_nnz_per_chunk=chunk_size)


@pytest.mark.parametrize("bad_value", [np.nan, np.inf, -1.0, 1.5])
def test_rejects_invalid_counts_even_in_unselected_columns(matrix, bad_value):
    # Position 2 belongs to column 0, which is not requested.
    matrix["data"][2] = bad_value
    with pytest.raises(ValueError, match="count"):
        read_csr_columns(matrix, [3, 1], max_nnz_per_chunk=2)


@pytest.mark.parametrize(
    "pointers",
    [
        [1, 2, 2, 5, 6, 7],  # Must start at zero.
        [0, 2, 1, 5, 6, 7],  # Must be nondecreasing.
        [0, 2, 2, 5, 6, 6],  # Must end at number of stored entries.
        [0, 2, 2, 8, 8, 7],  # Out-of-range intermediate pointers.
    ],
)
def test_rejects_invalid_pointer_values(matrix, pointers):
    matrix["indptr"][:] = pointers
    with pytest.raises(ValueError, match="indptr"):
        read_csr_columns(matrix, [1], max_nnz_per_chunk=2)


def test_rejects_wrong_pointer_length(matrix):
    del matrix["indptr"]
    matrix.create_dataset("indptr", data=[0, 2, 2, 5, 7])
    with pytest.raises(ValueError, match="indptr"):
        read_csr_columns(matrix, [1], max_nnz_per_chunk=2)


@pytest.mark.parametrize("bad_index", [-1, 4])
def test_rejects_out_of_range_stored_indices(matrix, bad_index):
    matrix["indices"][2] = bad_index
    with pytest.raises(ValueError, match="indices"):
        read_csr_columns(matrix, [1], max_nnz_per_chunk=2)


def test_rejects_mismatched_array_lengths(matrix):
    del matrix["indices"]
    matrix.create_dataset("indices", data=[1, 3])
    with pytest.raises(ValueError, match="length"):
        read_csr_columns(matrix, [1], max_nnz_per_chunk=2)
