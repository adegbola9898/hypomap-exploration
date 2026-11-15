# Detection summary contract

## Interface
summarise_detection(table, *, min_count)

Input:
- A DataFrame containing DENND1B_count and MC4R_count.
- Each row represents one nucleus in an already-selected population.
- min_count is a positive integer, excluding booleans.
- Counts must be finite, nonnegative integers.

Output dictionary:
- n_nuclei
- n_dennd1b
- n_mc4r
- n_both
- pct_dennd1b
- pct_mc4r
- pct_both
- pct_dennd1b_given_mc4r

Detection uses count >= min_count for each gene.

The first three percentages use all nuclei as their denominator.
The conditional percentage uses MC4R-positive nuclei.

Percentages are on a 0–100 scale and remain unrounded.
Undefined percentages return None:
- All percentages are undefined for an empty population.
- The conditional percentage is undefined when n_mc4r is zero.

Preserve the input table without modification.
Perform no population selection, grouping or significance testing.
