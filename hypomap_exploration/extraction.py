"""Bounded-memory extraction of gene columns from HDF5 CSR matrices."""

import h5py
import numpy as np


def _is_integer(value):
    return isinstance(value, (int, np.integer)) and not isinstance(
        value, (bool, np.bool_)
    )


def read_csr_columns(matrix_group, column_indices, *, max_nnz_per_chunk):
    """Validate a count-valued CSR group and extract requested columns."""
    if not _is_integer(max_nnz_per_chunk) or max_nnz_per_chunk <= 0:
        raise ValueError("chunk size must be a positive integer")

    shape = np.asarray(matrix_group.attrs.get("shape", []))
    if (
        shape.shape != (2,)
        or shape.dtype.kind not in "iu"
        or np.any(shape < 0)
    ):
        raise ValueError("matrix shape must contain two nonnegative integers")
    n_rows, n_columns = map(int, shape)

    try:
        columns = list(column_indices)
    except TypeError as exc:
        raise ValueError("column indices must be an iterable") from exc

    if not columns or any(not _is_integer(c) for c in columns):
        raise ValueError("column indices must be nonempty integers")
    if len(set(columns)) != len(columns):
        raise ValueError("column indices must be unique")
    if any(c < 0 or c >= n_columns for c in columns):
        raise ValueError("column indices are outside the matrix bounds")

    for name in ("data", "indices", "indptr"):
        if name not in matrix_group:
            raise ValueError(f"missing CSR array: {name}")
        obj = matrix_group[name]
        if not isinstance(obj, h5py.Dataset) or obj.ndim != 1:
            raise ValueError(f"{name} must be a one-dimensional dataset")

    data = matrix_group["data"]
    indices = matrix_group["indices"]
    pointer_data = matrix_group["indptr"]

    if len(data) != len(indices):
        raise ValueError("data and indices lengths must match")
    if indices.dtype.kind not in "iu":
        raise ValueError("indices must have integer dtype")
    if pointer_data.dtype.kind not in "iu":
        raise ValueError("indptr must have integer dtype")
    if data.dtype.kind not in "iuf":
        raise ValueError("count values must have real numeric dtype")
    if len(pointer_data) != n_rows + 1:
        raise ValueError("indptr length must equal number of rows plus one")

    pointers = pointer_data[:]
    if (
        pointers[0] != 0
        or pointers[-1] != len(data)
        or np.any(pointers < 0)
        or np.any(pointers > len(data))
        or np.any(pointers[1:] < pointers[:-1])
    ):
        raise ValueError("indptr contains invalid offsets")

    output = np.zeros((n_rows, len(columns)), dtype=np.float64)

    for start in range(0, len(data), int(max_nnz_per_chunk)):
        stop = min(start + int(max_nnz_per_chunk), len(data))
        chunk_columns = indices[start:stop]
        chunk_values = data[start:stop]

        if np.any(chunk_columns < 0) or np.any(chunk_columns >= n_columns):
            raise ValueError(f"indices outside matrix bounds in chunk {start}")

        if (
            not np.all(np.isfinite(chunk_values))
            or np.any(chunk_values < 0)
            or np.any(chunk_values != np.floor(chunk_values))
        ):
            raise ValueError(
                f"count values must be finite nonnegative integers "
                f"in chunk {start}"
            )

        for output_column, source_column in enumerate(columns):
            offsets = np.flatnonzero(chunk_columns == source_column)
            if offsets.size == 0:
                continue

            positions = start + offsets
            rows = np.searchsorted(pointers, positions, side="right") - 1
            np.add.at(
                output[:, output_column],
                rows,
                chunk_values[offsets],
            )

    if not np.all(np.isfinite(output)):
        raise ValueError("count accumulation produced non-finite values")

    return output
