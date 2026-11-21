# Interaction simulation analysis plan

## Shared settings
- Genotype distribution: Binomial(2, 0.30).
- Intercept: 0.
- Group main effect: 0.30.
- Residual distribution: Normal(mean=0, standard deviation=1).
- Master random seed: 20261121.
- Repetitions: 1,000 per scenario.
- Technical pilot: 20 repetitions per scenario.
- Significance rule: two-sided p < 0.05.
- Confidence intervals: 95%.

These are illustrative settings, not estimates from canine research.
The significance threshold is for a single-variant demonstration,
not a genome-wide discovery analysis.

## Scenarios

| Scenario | Intact n | Neutered n | Intact genetic effect | Neutered genetic effect | Interaction |
|---|---:|---:|---:|---:|---:|
| A: Same effect, unequal samples | 200 | 2000 | 0.15 | 0.15 | 0 |
| B: Different effects, larger samples | 2000 | 2000 | 0.15 | 0.35 | 0.20 |
| C: Different effects, smaller samples | 200 | 200 | 0.15 | 0.35 | 0.20 |

## Models
Fit ordinary least squares:
1. Y ~ G within intact individuals.
2. Y ~ G within neutered individuals.
3. Y ~ G + N + G:N in the combined dataset.

Retain effect estimates, standard errors, confidence intervals
and p-values. Record group sizes and observed allele frequencies.

Check that the combined interaction estimate equals the difference
between the two group-specific slope estimates within numerical tolerance.

## Repeated-simulation summaries
For each scenario, report:
- Association rejection rate in each group.
- Frequency of significance in exactly one group.
- Interaction rejection rate.
- Mean interaction estimate and bias.
- Coverage of the true interaction by its 95% confidence interval.
- Monte Carlo standard errors for rejection proportions.

Interaction rejection measures false-positive rate in A and power
in B and C. Discordant group significance is not an interaction test.

## Figures
1. Effect estimates and confidence intervals from the first repetition
   of each scenario, with known true effects indicated.
2. Frequencies of discordant group significance and significant
   interaction across repeated simulations, with Monte Carlo uncertainty.

Do not search for a more dramatic illustrative repetition.

## Reproducibility
Use distinct reproducible random streams for scenarios and repetitions.
Record settings, seed, software versions, code version and fit failures.
Store per-repetition model summaries rather than all individual-level
synthetic datasets.

## Completion criteria
A tested implementation, checked results table, two figures and
a short interpretation explaining findings and limitations.

## Status
Plan documented; implementation has not started.

## Data-generator interface
simulate_dataset(settings, rng)

Inputs:
- settings: mapping containing n_intact, n_neutered, allele_frequency,
  intercept, beta_g, beta_n, beta_gn and residual_sd.
- rng: an explicitly supplied NumPy random Generator.

Output:
- A pandas DataFrame with columns G, N, Y.
- Exactly n_intact + n_neutered rows.
- Intact rows first (N=0), then neutered rows (N=1).
- Integer genotype dosage G in {0, 1, 2}.
- Floating-point outcome Y in arbitrary units.

Behaviour:
- Generate genotypes with the same allele frequency in both groups.
- Generate independent normally distributed residuals.
- Apply the documented data-generating equation.
- Do not mutate settings, write files or fit models.
- The supplied random generator advances as draws are consumed.

Validation:
- Group sizes must be positive integers, excluding booleans.
- Allele frequency must be finite and strictly between 0 and 1.
- Model coefficients must be finite.
- Residual standard deviation must be finite and nonnegative.
- Zero residual standard deviation is permitted for deterministic tests.
