from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from interaction_simulation.simulation import simulate_dataset


@pytest.fixture
def settings():
    return {
        "n_intact": 40,
        "n_neutered": 60,
        "allele_frequency": 0.30,
        "intercept": 1.0,
        "beta_g": 0.5,
        "beta_n": 2.0,
        "beta_gn": 1.5,
        "residual_sd": 1.0,
    }


def test_structure_and_settings_unchanged(settings):
    original = deepcopy(settings)
    table = simulate_dataset(settings, np.random.default_rng(123))

    assert table.columns.tolist() == ["G", "N", "Y"]
    assert len(table) == 100
    assert table["N"].tolist() == [0] * 40 + [1] * 60
    assert pd.api.types.is_integer_dtype(table["G"])
    assert set(table["G"].unique()).issubset({0, 1, 2})
    assert pd.api.types.is_float_dtype(table["Y"])
    assert np.isfinite(table["Y"]).all()
    assert settings == original


def test_same_seed_reproduces_dataset(settings):
    first = simulate_dataset(settings, np.random.default_rng(456))
    second = simulate_dataset(settings, np.random.default_rng(456))

    pd.testing.assert_frame_equal(first, second)


def test_zero_noise_matches_known_equation(settings):
    settings["residual_sd"] = 0.0
    table = simulate_dataset(settings, np.random.default_rng(789))

    intact = table.loc[table["N"] == 0]
    neutered = table.loc[table["N"] == 1]

    # Intact: intercept 1, genotype effect 0.5.
    np.testing.assert_allclose(
        intact["Y"], 1.0 + 0.5 * intact["G"]
    )

    # Neutered: intercept 1 + 2, genotype effect 0.5 + 1.5.
    np.testing.assert_allclose(
        neutered["Y"], 3.0 + 2.0 * neutered["G"]
    )


@pytest.mark.parametrize(
    "field,value",
    [
        ("n_intact", 0),
        ("n_neutered", -1),
        ("n_intact", 2.5),
        ("n_neutered", True),
        ("allele_frequency", 0.0),
        ("allele_frequency", 1.0),
        ("allele_frequency", np.nan),
        ("intercept", np.inf),
        ("beta_g", np.nan),
        ("beta_n", -np.inf),
        ("beta_gn", np.inf),
        ("residual_sd", -1.0),
        ("residual_sd", np.nan),
    ],
)
def test_invalid_settings_rejected(settings, field, value):
    settings[field] = value

    with pytest.raises(ValueError, match=field):
        simulate_dataset(settings, np.random.default_rng(123))
