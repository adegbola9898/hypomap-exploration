"""Generate primary detection summaries from the verified extract."""

import json
from pathlib import Path

import pandas as pd

from hypomap_exploration.detection import summarise_detection
from hypomap_exploration.run_extraction import sha256_file


def main():
    root = Path(__file__).resolve().parents[1]
    config_path = root / "config/analysis.json"
    config = json.loads(config_path.read_text())
    manifest = json.loads(
        (root / "provenance/extraction_manifest.json").read_text()
    )
    source = root / manifest["output_file"]

    if manifest["status"] != "complete":
        raise ValueError("Extraction manifest is not complete")
    if sha256_file(source) != manifest["output_sha256"]:
        raise ValueError("Extracted-table checksum mismatch")
    if sha256_file(config_path) != manifest["config_sha256"]:
        raise ValueError("Configuration differs from the extraction checkpoint")

    table = pd.read_csv(source, **manifest["table_read_options"])
    if len(table) != config["extraction"]["expected_nuclei"]:
        raise ValueError("Unexpected total nucleus count")
    if not table["cell_id"].is_unique:
        raise ValueError("Duplicate cell identifiers")

    neurons = table.loc[
        table[config["selection_field"]].eq(config["selection_value"])
    ]
    if len(neurons) != config["expected_selected_nuclei"]:
        raise ValueError("Unexpected neuronal nucleus count")

    group_field = config["group_field"]
    if set(neurons[group_field]) != set(config["group_order"]):
        raise ValueError("Unexpected neuronal population labels")

    threshold = config["primary_min_count"]

    def summarise(label, subset):
        return {
            "population": label,
            "min_count": threshold,
            **summarise_detection(subset, min_count=threshold),
        }

    overall = pd.DataFrame([
        summarise("All nuclei", table),
        summarise("C0 neurons", neurons),
    ])
    populations = pd.DataFrame([
        summarise(label, neurons.loc[neurons[group_field].eq(label)])
        for label in config["group_order"]
    ])

    # C1 populations must partition the selected neurons.
    for column in ["n_nuclei", "n_dennd1b", "n_mc4r", "n_both"]:
        if populations[column].sum() != overall.iloc[1][column]:
            raise ValueError(f"C1 totals do not reconcile for {column}")

    destination = root / "results/tables"
    destination.mkdir(parents=True, exist_ok=True)
    overall.to_csv(
        destination / "primary_overall.tsv",
        sep="\t", index=False, na_rep="__MISSING__",
    )
    populations.to_csv(
        destination / "primary_c1.tsv",
        sep="\t", index=False, na_rep="__MISSING__",
    )

    # Round only the terminal display; saved percentages retain precision.
    print("Overall summaries")
    print(overall.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print("\nC1 neuronal populations")
    print(populations.to_string(
        index=False, float_format=lambda x: f"{x:.4f}"
    ))
    print("\nC1 count reconciliation: PASS")


if __name__ == "__main__":
    main()
