import pandas as pd
import pytest

from hypomap_exploration.detection import summarise_detection


@pytest.fixture
def counts():
    # At threshold 1: D=6, M=4, both=3, total=10.
    return pd.DataFrame({
        "DENND1B_count": [2, 1, 1, 1, 1, 1, 0, 0, 0, 0],
        "MC4R_count":    [2, 1, 1, 0, 0, 0, 1, 0, 0, 0],
    })


def test_known_detection_percentages(counts):
    original = counts.copy(deep=True)
    result = summarise_detection(counts, min_count=1)

    assert result == {
        "n_nuclei": 10,
        "n_dennd1b": 6,
        "n_mc4r": 4,
        "n_both": 3,
        "pct_dennd1b": 60.0,
        "pct_mc4r": 40.0,
        "pct_both": 30.0,
        "pct_dennd1b_given_mc4r": 75.0,
    }
    pd.testing.assert_frame_equal(counts, original)


def test_threshold_is_inclusive(counts):
    result = summarise_detection(counts, min_count=2)

    assert result["n_dennd1b"] == 1
    assert result["n_mc4r"] == 1
    assert result["n_both"] == 1
    assert result["pct_both"] == 10.0
    assert result["pct_dennd1b_given_mc4r"] == 100.0


def test_no_mc4r_detection_has_undefined_conditional_percentage():
    table = pd.DataFrame({
        "DENND1B_count": [1, 0],
        "MC4R_count": [0, 0],
    })
    result = summarise_detection(table, min_count=1)

    assert result["pct_dennd1b"] == 50.0
    assert result["pct_mc4r"] == 0.0
    assert result["pct_both"] == 0.0
    assert result["pct_dennd1b_given_mc4r"] is None


def test_empty_population_has_undefined_percentages():
    table = pd.DataFrame({
        "DENND1B_count": pd.Series(dtype=float),
        "MC4R_count": pd.Series(dtype=float),
    })
    result = summarise_detection(table, min_count=1)

    for key in ["n_nuclei", "n_dennd1b", "n_mc4r", "n_both"]:
        assert result[key] == 0
    for key in [
        "pct_dennd1b", "pct_mc4r", "pct_both",
        "pct_dennd1b_given_mc4r",
    ]:
        assert result[key] is None


@pytest.mark.parametrize("threshold", [0, -1, 1.5, True])
def test_rejects_invalid_threshold(counts, threshold):
    with pytest.raises(ValueError, match="min_count"):
        summarise_detection(counts, min_count=threshold)


@pytest.mark.parametrize("column", ["DENND1B_count", "MC4R_count"])
@pytest.mark.parametrize("bad_value", [float("nan"), float("inf"), -1., 1.5])
def test_rejects_invalid_counts(counts, column, bad_value):
    table = counts.astype(float)
    table.loc[0, column] = bad_value

    with pytest.raises(ValueError, match="counts"):
        summarise_detection(table, min_count=1)


def test_rejects_missing_count_column(counts):
    with pytest.raises(ValueError, match="count columns"):
        summarise_detection(counts.drop(columns="MC4R_count"), min_count=1)


def test_rejects_string_counts(counts):
    table = counts.astype(str)
    with pytest.raises(ValueError, match="numeric dtype"):
        summarise_detection(table, min_count=1)


def test_percentages_are_not_rounded():
    table = pd.DataFrame({
        "DENND1B_count": [1, 0, 0],
        "MC4R_count": [1, 1, 1],
    })
    result = summarise_detection(table, min_count=1)

    assert result["pct_dennd1b"] == pytest.approx(100 / 3)
    assert result["pct_both"] == pytest.approx(100 / 3)
    assert result["pct_dennd1b_given_mc4r"] == pytest.approx(100 / 3)
    assert result["pct_both"] != round(100 / 3, 2)
