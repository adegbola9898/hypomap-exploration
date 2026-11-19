# S10 detection comparison findings

Using the public human HYPOMAP dataset, raw/X counts >= 1,
and the existing C4 cluster annotations, independently written code
recovered all selected S10 detection percentages within declared
reporting tolerances.

Cluster identifiers, names and nucleus counts matched for all
452 clusters, collectively containing 433,369 nuclei.

Results:
- DENND1B detection: 452/452 values matched.
- MC4R detection: 452/452 values matched.
- Joint DENND1B–MC4R detection: 452/452 values matched.
- Total: 1,356 comparisons; zero mismatches.

Comparison tolerances:
- Single-gene percentages: 0.0050000001 percentage points,
  accommodating two-decimal reporting and numerical precision.
- Joint-detection percentages: 1e-8 percentage points.

Maximum observed absolute differences:
- DENND1B: approximately 0.005 percentage points.
- MC4R: approximately 0.005 percentage points.
- Joint detection: approximately 4.796163466e-14 percentage points.

All percentages use the total nuclei within the cluster as denominator.

Interpretation:
This reproduces the selected published detection summaries using
the same underlying public atlas. It supports the count >= 1
operational definition for these results but does not establish
the authors' exact code.

Not reproduced:
- Average-expression columns.
- Other candidate genes or receptor pairs.
- Mouse results.
- Experimental findings or mechanistic conclusions.

This is computational reproduction, not independent biological
validation. The primary exploratory analysis retains its predefined
C0 neuronal selection and C1 grouping.

Reproduce with:
uv run python -m scripts.compare_s10

Evidence:
- references/s10.xlsx
- results/tables/s10_detection_comparison.tsv
- provenance/extraction_manifest.json
