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
Environment setup and initial dataset audit completed. Primary population and count-matrix choices documented. Gene extraction and co-detection analysis have not started.

## Frozen primary population and matrix choices

Select C0_named == "C0-3 Neuron": 166,475 nuclei.
This includes both "Neurons" and "P2RX2_OTP" celltype annotations.
Group by the seven existing C1_named neuronal populations.
Retain every selected population, including zero-detection results.
Population selection was completed before gene co-detection analysis.

Use raw/X for count-based analysis. In 30 sampled rows, values
were finite nonnegative integers and totals matched nCount_RNA.
This supports the count interpretation but is not full-matrix validation.
Validate the extracted gene values before analysis.

Preserve Donor_ID, Dataset and existing QC metadata.
Donor coverage varies strongly by population; pooled percentages
describe sampled nuclei and are not donor-balanced estimates.
No new clustering or significance testing is planned.

The machine-readable choices are recorded in config/analysis.json.

## Amendment following inspection of supplementary tables S9 and S10

The original extraction proposal selected neurons immediately.
The supplied reference tables require broader coverage:
S9 contains all-cell summaries, and S10 includes non-neuronal clusters.

Extract the two target genes across all 433,369 nuclei and retain
C3/C4 identifiers and names alongside the existing metadata.
The primary exploratory analysis remains restricted to the previously
selected 166,475 C0 neurons and seven C1 populations.

S10 has 452 cluster rows with mixed C3/C4 identifiers. Mapping these
rows to atlas annotations requires verification before comparison.
Matching total nucleus counts alone does not establish identical
cell membership or analysis settings.

The first comparison concerns human DENND1B and MC4R detection
percentages only. Published detection thresholds and receptor-subgroup
definitions remain unverified. Average-expression comparisons are
deferred until their calculation is established.

This amendment follows inspection of reference tables, before
examining extracted DENND1B–MC4R co-detection results.
