# Does a stronger association signal mean a stronger genetic effect?

## Research question
How do differences in sample size and genetic effects produce different
association results between two groups, and what does a formal
interaction test add?

## Purpose
A reproducible simulation exercise to develop understanding of
group-specific genetic associations and genotype-by-group interaction.

The labels "intact" and "neutered" are illustrative.
No real canine data are analysed, and the simulation does not reproduce
the biology or population structure of dogs.

## Model
Y = beta_0 + beta_G*G + beta_N*N + beta_GN*G*N + error

- G: genotype dosage, 0, 1 or 2.
- N: group indicator, 0 = intact and 1 = neutered.
- Y: hypothetical continuous adiposity outcome in arbitrary units.
- Intact genetic effect: beta_G.
- Neutered genetic effect: beta_G + beta_GN.
- Interaction: the difference between those genetic effects.

## Initial scope
- One variant and two groups.
- Three predefined scenarios.
- Separate-group regressions and a combined interaction model.
- Repeated simulations with known underlying effects.
- One summary table, two figures and a short interpretation.

## Assumptions
Independent individuals, equal allele frequency across groups,
normally distributed errors and constant residual variance.
Genotype and group assignment are independent.

## Boundaries
No real GWAS, breed structure, relatedness, mixed models, binary
outcomes, hormone measurements or genome-wide testing.

## Working rules
Record scenarios and random seed before the full experiment.
Do not select examples because they produce striking results.
Do not interpret significance in only one group as an interaction test.
Do not silently discard failed model fits.

## Current status
Planning only. No simulation code implemented or results generated.
