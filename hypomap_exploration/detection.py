"""Descriptive gene detection summaries."""

import numpy as np
import pandas as pd


def summarise_detection(table, *, min_count):
    """Summarise detection within an already-selected population."""
    if (
        isinstance(min_count, (bool, np.bool_))
        or not isinstance(min_count, (int, np.integer))
        or min_count <= 0
    ):
        raise ValueError("min_count must be a positive integer")

    required = ["DENND1B_count", "MC4R_count"]
    if not isinstance(table, pd.DataFrame):
        raise ValueError("table must be a pandas DataFrame")
    if not table.columns.is_unique:
        raise ValueError("table column names must be unique")
    if any(column not in table.columns for column in required):
        raise ValueError("table must contain both required count columns")

    for column in required:
        values = table[column].to_numpy()
        if values.dtype.kind not in "iuf":
            raise ValueError(f"{column}: counts must have real numeric dtype")
        if (
            not np.isfinite(values).all()
            or (values < 0).any()
            or (values != np.floor(values)).any()
        ):
            raise ValueError(
                f"{column}: counts must be finite nonnegative integers"
            )

    dennd1b = table["DENND1B_count"].to_numpy() >= min_count
    mc4r = table["MC4R_count"].to_numpy() >= min_count

    n_nuclei = len(table)
    n_dennd1b = int(dennd1b.sum())
    n_mc4r = int(mc4r.sum())
    n_both = int((dennd1b & mc4r).sum())

    def percentage(numerator, denominator):
        return 100.0 * numerator / denominator if denominator else None

    return {
        "n_nuclei": n_nuclei,
        "n_dennd1b": n_dennd1b,
        "n_mc4r": n_mc4r,
        "n_both": n_both,
        "pct_dennd1b": percentage(n_dennd1b, n_nuclei),
        "pct_mc4r": percentage(n_mc4r, n_nuclei),
        "pct_both": percentage(n_both, n_nuclei),
        "pct_dennd1b_given_mc4r": percentage(n_both, n_mc4r),
    }
