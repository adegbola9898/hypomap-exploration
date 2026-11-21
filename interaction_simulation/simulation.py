"""Synthetic data for the genotype-by-group interaction exercise."""

from numbers import Integral, Real

import numpy as np
import pandas as pd


def simulate_dataset(settings, rng):
    """Return genotype G, group N and outcome Y from the specified model."""
    for field in ("n_intact", "n_neutered"):
        value = settings[field]
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, Integral)
            or value <= 0
        ):
            raise ValueError(f"{field} must be a positive integer")

    numeric_fields = (
        "allele_frequency",
        "intercept",
        "beta_g",
        "beta_n",
        "beta_gn",
        "residual_sd",
    )
    for field in numeric_fields:
        value = settings[field]
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, Real)
            or not np.isfinite(value)
        ):
            raise ValueError(f"{field} must be a finite number")

    if not 0 < settings["allele_frequency"] < 1:
        raise ValueError("allele_frequency must be strictly between 0 and 1")

    if settings["residual_sd"] < 0:
        raise ValueError("residual_sd must be nonnegative")

    n_intact = int(settings["n_intact"])
    n_neutered = int(settings["n_neutered"])
    n_total = n_intact + n_neutered

    group = np.concatenate([
        np.zeros(n_intact, dtype=np.int64),
        np.ones(n_neutered, dtype=np.int64),
    ])
    genotype = rng.binomial(
        n=2,
        p=settings["allele_frequency"],
        size=n_total,
    )
    residual = rng.normal(
        loc=0.0,
        scale=settings["residual_sd"],
        size=n_total,
    )

    outcome = (
        settings["intercept"]
        + settings["beta_g"] * genotype
        + settings["beta_n"] * group
        + settings["beta_gn"] * genotype * group
        + residual
    )

    return pd.DataFrame({
        "G": genotype,
        "N": group,
        "Y": outcome,
    })
