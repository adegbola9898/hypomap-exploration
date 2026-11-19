"""Compare two-gene detection summaries against published Table S10."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from hypomap_exploration.detection import summarise_detection
from hypomap_exploration.run_extraction import sha256_file


def main():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (root / "provenance/extraction_manifest.json").read_text()
    )
    source = root / manifest["output_file"]

    if manifest["status"] != "complete":
        raise ValueError("Extraction is not marked complete")
    if sha256_file(source) != manifest["output_sha256"]:
        raise ValueError("Extracted-table checksum mismatch")

    checksums = {}
    for line in (root / "references/SHA256SUMS").read_text().splitlines():
        if line.strip():
            checksum, filename = line.split(maxsplit=1)
            checksums[filename] = checksum

    reference_path = root / "references/s10.xlsx"
    if sha256_file(reference_path) != checksums["references/s10.xlsx"]:
        raise ValueError("S10 checksum mismatch")

    table = pd.read_csv(source, **manifest["table_read_options"])
    if len(table) != manifest["n_nuclei"] or not table["cell_id"].is_unique:
        raise ValueError("Unexpected row count or duplicate cell identifiers")

    published = pd.read_excel(
        reference_path, sheet_name="S10", header=1, keep_default_na=False
    )
    published = published.loc[
        published["Cluster Number"].astype(str).str.strip().ne("")
    ].copy()

    if published["Cluster Number"].duplicated().any():
        raise ValueError("Duplicate S10 cluster identifiers")

    reference = published.set_index("Cluster Number")
    if len(reference) != 452 or set(table["C4"]) != set(reference.index):
        raise ValueError("Atlas and S10 cluster identifiers differ")

    # Tolerances are in percentage points.
    # Single-gene percentages are published to two decimal places.
    metrics = [
        ("pct_dennd1b", "% DENND1B+", 0.0050000001),
        ("pct_mc4r", "% MC4R+", 0.0050000001),
        ("pct_both", "% DENND1B+MC4R+", 1e-8),
    ]

    records = []
    for cluster, subset in table.groupby("C4", sort=True):
        published_row = reference.loc[cluster]

        if len(subset) != float(published_row["Total number of cells"]):
            raise ValueError(f"Cluster size mismatch: {cluster}")
        if set(subset["C4_named"]) != {published_row["Cluster Name"]}:
            raise ValueError(f"Cluster name mismatch: {cluster}")

        calculated = summarise_detection(subset, min_count=1)

        for metric, column, tolerance in metrics:
            reported = float(published_row[column])
            if not np.isfinite(reported) or not 0 <= reported <= 100:
                raise ValueError(f"Invalid reference value: {cluster}, {column}")

            difference = calculated[metric] - reported
            records.append({
                "cluster": cluster,
                "cluster_name": published_row["Cluster Name"],
                "n_nuclei": len(subset),
                "min_count": 1,
                "metric": metric,
                "calculated_percent": calculated[metric],
                "published_percent": reported,
                "difference_percentage_points": difference,
                "tolerance_percentage_points": tolerance,
                "matches": abs(difference) <= tolerance,
            })

    results = pd.DataFrame(records)
    output = root / "results/tables/s10_detection_comparison.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(output, sep="\t", index=False)

    print("Cluster identifiers, names and sizes: PASS")
    for metric, group in results.groupby("metric", sort=False):
        maximum = group["difference_percentage_points"].abs().max()
        print(
            f"{metric}: {int(group['matches'].sum())}/{len(group)} match; "
            f"maximum absolute difference = {maximum:.10g} percentage points"
        )

    mismatches = results.loc[~results["matches"]]
    print(f"\nTotal mismatches: {len(mismatches)}")
    if not mismatches.empty:
        print(mismatches.head(20).to_string(index=False))

    print(f"\nSaved: {output.relative_to(root)}")


if __name__ == "__main__":
    main()
