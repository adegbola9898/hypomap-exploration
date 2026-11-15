# Initial S9 comparison findings

Five selected human S9 percentages match when calculations are rounded
to the publication's two decimal places.

All nuclei:
- DENND1B detection: 30.84%.
- MC4R detection: 1.16%.

Using celltype_annotation == "Neurons" (165,270 nuclei):
- DENND1B detection: 107,600 / 165,270 = 65.105585%, reported as 65.11%.
- MC4R detection: 4,796 / 165,270 = 2.901918%, reported as 2.90%.
- DENND1B among MC4R-positive neurons:
  3,834 / 4,796 = 79.941618%, reported as 79.94%.

Calculations use raw/X counts >= 1.

This supports the narrower neuronal annotation as a comparison
population for these S9 metrics. It does not establish the authors'
exact code or uniquely prove their detection threshold.

The primary exploratory population remains C0_named == "C0-3 Neuron",
containing 166,475 nuclei. Its additional 1,205 P2RX2_OTP nuclei explain
the difference between the two neuronal selections.

Comparison limits:
- Five selected human percentages checked at published precision.
- No mouse results, other genes, POMC/GHSR subgroups or average-expression
  results reproduced.
- S10 expression values have not yet been compared.
- Reanalysis of the same atlas is not independent biological validation.

Evidence:
- results/tables/primary_overall.tsv
- results/tables/s9_annotation_diagnostic.tsv
- references/s9.xlsx, sheet S9, rows 12, 13 and 19
