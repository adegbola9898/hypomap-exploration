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


def _read_vector(handle, path, label):
    """Read a plain HDF5 vector without silently interpreting categories."""
    if path not in handle:
        raise ValueError(f"{label}: missing field {path}")

    obj = handle[path]
    if not isinstance(obj, h5py.Dataset) or obj.ndim != 1:
        raise ValueError(
            f"{label}: expected a plain one-dimensional dataset at {path}; "
            "categorical groups require explicit decoding"
        )

    legacy_categories = f"{obj.parent.name}/__categories/{obj.name.rsplit('/', 1)[-1]}"
    if (
        legacy_categories in handle
        or "categories" in obj.attrs
        or obj.attrs.get("encoding-type") in ("categorical", b"categorical")
    ):
        raise ValueError(f"{label}: categorical encoding requires decoding")

    if h5py.check_string_dtype(obj.dtype) is not None:
        return obj.asstr()[:]
    return obj[:]


def extract_atlas_counts(handle, config, *, max_nnz_per_chunk):
    """Return an aligned count/metadata table and gene-column mapping."""
    import pandas as pd

    genes = config["genes"]
    if (
        not isinstance(genes, list)
        or not genes
        or any(not isinstance(gene, str) or not gene for gene in genes)
        or len(set(genes)) != len(genes)
    ):
        raise ValueError("gene names must be a nonempty list of unique strings")

    matrix_path = config["matrix_path"]
    if matrix_path not in handle:
        raise ValueError(f"missing matrix: {matrix_path}")
    matrix = handle[matrix_path]

    shape = np.asarray(matrix.attrs.get("shape", []))
    if (
        shape.shape != (2,)
        or shape.dtype.kind not in "iu"
        or np.any(shape < 0)
    ):
        raise ValueError("matrix shape must contain two nonnegative integers")
    n_rows, n_columns = map(int, shape)

    extraction = config["extraction"]
    if extraction["population"] != "all_nuclei":
        raise ValueError("extraction population must be all_nuclei")
    expected = extraction["expected_nuclei"]
    if not _is_integer(expected) or expected < 0 or n_rows != expected:
        raise ValueError("source cell count differs from expected_nuclei")

    features = _read_vector(handle, config["feature_index_path"], "gene")
    if len(features) != n_columns:
        raise ValueError("gene feature count differs from matrix columns")

    mapping = {}
    for gene in genes:
        matches = np.flatnonzero(features == gene)
        if len(matches) != 1:
            raise ValueError(
                f"gene {gene!r} must match exactly one feature; "
                f"found {len(matches)}"
            )
        mapping[gene] = int(matches[0])

    cells = _read_vector(handle, config["cell_index_path"], "cell")
    if len(cells) != n_rows:
        raise ValueError("cell identifier count differs from matrix rows")
    if any(not isinstance(cell, str) or not cell for cell in cells):
        raise ValueError("cell identifiers must be nonempty strings")
    if not pd.Index(cells).is_unique:
        raise ValueError("cell identifiers must be unique")

    fields = extraction["preserve_obs_fields"]
    if (
        not isinstance(fields, list)
        or any(not isinstance(field, str) or not field for field in fields)
        or len(set(fields)) != len(fields)
    ):
        raise ValueError("metadata fields must be unique nonempty strings")

    generated = {"cell_id", "source_row", *(f"{gene}_count" for gene in genes)}
    if generated.intersection(fields):
        raise ValueError("metadata fields collide with generated columns")

    values = {
        "cell_id": cells,
        "source_row": np.arange(n_rows, dtype=np.int64),
    }
    for field in fields:
        vector = _read_vector(handle, f"obs/{field}", "metadata")
        if len(vector) != n_rows:
            raise ValueError(f"metadata length mismatch for {field}")
        values[field] = vector

    counts = read_csr_columns(
        matrix,
        list(mapping.values()),
        max_nnz_per_chunk=max_nnz_per_chunk,
    )
    for position, gene in enumerate(genes):
        values[f"{gene}_count"] = counts[:, position]

    return pd.DataFrame(values), mapping
