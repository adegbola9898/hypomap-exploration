# Why was DENND1B prioritised?

## Purpose
I reviewed the published evidence at one canine obesity-associated locus
to understand how researchers move from a genetic association to a
candidate gene. This is a retrospective applicant exercise: I knew the
published conclusion before starting.

## What I established
S3 defines the region tagged by 7:5004016 as
CanFam3.1 chr7:3062016-5279201, approximately 2.2 Mb.
The main paper lists seven genes: NR5A2, PTPRC, ATP6V1G3, NEK7,
LHX9, DENND1B and CRB1.

I recorded this published candidate set and reviewed its coverage in
S5-S7. I have not independently reconstructed the set from genome
annotation. The methods describe supplementing canine annotation with
genes from corresponding human regions.

## Why DENND1B stands out in these sources
DENND1B has directly reported human BMI associations alongside
expression-based and regulatory evidence. Its reported SMR/HEIDI and
colocalisation results meet the authors' stated thresholds.
The ABC analysis links associated variants or their proxies to
DENND1B regulatory elements.

The rare-variant association and single homozygous patient observation
provide additional, qualified support. Cell perturbation experiments
support a role in MC4R regulation, but do not establish how the canine
lead variant acts.

These findings are not independent votes: S6 and S7 repeat some
association evidence, and G2G combines component analyses.

## What the alternatives show
S7 reports WHR and lipid evidence for NR5A2, and HDL evidence for
ATP6V1G3, NEK7 and PTPRC. The latter three share a reported HDL lead
signal. Those results should not be treated as equivalent to BMI
evidence or as independent association discoveries.

LHX9 and CRB1 have no matching rows in the inspected S5-S7 tables.
This is a reporting gap for this exercise, not evidence against them.
The tables do not support an exhaustive, equally informed ranking of
all seven candidates.

## Relationship to my HYPOMAP work
I previously recalculated selected DENND1B and MC4R detection metrics
and matched the published S10 values across 452 clusters within the
comparison tolerances. That checks those calculations using the same
atlas; it is not independent biological replication.

Co-detection supplies cellular context. It does not demonstrate protein
interaction, receptor trafficking or the causal gene at a GWAS locus.

## My conclusion
The inspected evidence gives a coherent reason to prioritise DENND1B
for functional investigation. My work reconstructs and critiques that
argument; it neither discovers a new obesity gene nor excludes the
other candidates.

## Questions for the supervisor
1. How should the seven-gene list be reconciled with the canine interval,
   human syntenic additions and shorter supplementary description?
2. Are complete gene-by-trait results available for comparing all seven
   candidates, including results omitted from S7?
3. What existing data could address the link between the canine
   associated variants and DENND1B regulation?

## Evidence trail
- locus_definition.md and candidates.tsv
- evidence_overview.tsv and dennd1b_evidence.md
- references/s3.xlsx through references/s7.xlsx
- Supplied supplementary methods and the main paper:
  https://doi.org/10.1126/science.ads2145
