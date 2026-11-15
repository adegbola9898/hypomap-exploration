# HYPOMAP exploratory analysis plan

## Purpose
Develop a reproducible exploratory analysis of public human HYPOMAP
data for discussion with a prospective supervisor.

## Primary question
How do DENND1B and MC4R detection, separately and together, vary
across predefined hypothalamic neuronal populations?

## Core scope
- One public human HYPOMAP dataset.
- Two genes: DENND1B and MC4R.
- Existing neuronal annotations at one documented clustering level.
- Descriptive detection and co-detection summaries.
- Sensitivity to the expression-detection threshold.
- Donor-level summaries if suitable metadata are available.

## Proposed detection rules
Primary: at least one measured raw count per gene.
Sensitivity: at least two measured raw counts per gene.
These rules require verification that raw counts are available.
They are exploratory choices, not presumed publication settings.

## Population selection
Choose populations using documented annotations before examining
co-detection results. Record the selection and its rationale.
Do not select populations because they show high co-detection.

## Exclusions
No genotype imputation, GWAS, atlas reconstruction, new clustering,
mouse comparison, hormone-receptor extension or significance tests
in the initial analysis.

## Interpretation
RNA co-detection does not establish protein interaction or causation.
Zero counts indicate non-detection, not proven biological absence.
Reanalysis of the same atlas is not independent biological replication.
Exact reproduction of published results remains unverified until
data versions, definitions and numerical outputs have been compared.

## Initial deliverables
- Documented data source and metadata audit.
- Tested detection-summary calculations.
- One primary results table and two figures.
- Threshold sensitivity assessment.
- A short supervisor discussion brief.

## Status
Environment preparation only. No biological data analysed.
