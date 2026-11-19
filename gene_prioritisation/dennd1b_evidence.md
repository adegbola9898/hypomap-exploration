# DENND1B: published evidence assessment

Status: retrospective assessment of the supplied tables.
This document does not reproduce the underlying association analyses.

| Evidence | Reported observation | Source | Contribution and limitation |
|---|---|---|---|
| Human common-variant BMI association | rs6702421: T allele associated with approximately +0.011 kg/m² BMI, p=9.42e-9; secondary signal rs1009188: approximately +0.012 kg/m² per T allele, p=7.15e-11 | S6 B3; S7 C6:L6 | Supports association in the human region. Association alone does not identify the causal gene. S6 and S7 overlap in their evidence. |
| Allele direction | S7 reports C as the effect allele at rs6702421, beta=-0.0114; S6 reports the opposite T allele with a positive effect | S6 B3; S7 G6:J6 | Directions agree after accounting for the effect allele and rounding. |
| Expression colocalisation | H4 posterior probability=0.979195408; source labelled eQTLGen | S7 U6:W6 | Supports a shared association signal under the colocalisation model. It is not a probability that DENND1B causes obesity and does not establish a hypothalamic mechanism. |
| SMR and HEIDI | SMR FDR-corrected p=1.72e-6; HEIDI p=0.03007395 | S7 X6:AA6 | Passes the authors' criteria: FDR-corrected SMR p < 0.05 and HEIDI p > 0.001 (supplementary methods, page 13). Supports their shared-signal interpretation under the method's assumptions; does not establish causality. |
| PoPS | DENND1B labelled the top PoPS gene; score=0.011524202 | S7 AB6:AC6 | Supports computational prioritisation. It is not independent experimental validation. |
| Regulatory links | ABC link marked yes; two elements; maximum score=0.182727; maximum-scoring element tissue=HAP1 | S7 AD6:AG6 | Supports a regulatory connection. The tissue listed describes the maximum-scoring element, not necessarily both elements. This row alone does not establish regulation in hypothalamic neurons. |
| Integrated G2G result | Genome-wide rank=238; percentile=96; causal predictor N=4 | S7 M6:P6 | Summarises component evidence. The percentile is not a causal probability; do not count the summary as an additional independent confirmation. |
| Rare-variant association | UKB ExWAS p=0.0087; beta=0.35 kg/m² | S6 C3 | Report as nominal evidence. This cell does not establish significance after correction for the full testing burden. |
| Patient observation | One patient homozygous for DENND1B p.R501C | S6 E3 | A potentially informative observation, not proof of pathogenicity. Assess segregation, phenotype, population frequency and functional evidence separately. |
| Prior biological findings | Asthma susceptibility and immune-system phenotypes | S5 C6:F6 | Provides biological context but does not by itself establish a role in appetite or BMI. |
| Unreported or unavailable entries | PheWAS=N/A; tissue enrichment blank; SCOOP-UK entry blank; S7 coding fields=NA | S5 B6/G6; S6 D3; S7 AH6:AL6 | Preserve these source states. They do not establish absence of a biological effect. |

## Provisional assessment

Among the seven candidates, DENND1B has the most directly documented
BMI-focused evidence in the inspected S5-S7 rows. This describes the
coverage and content of these tables, not an exhaustive comparison of
everything known about all seven genes.

The common-variant association, regulatory analyses and integrated G2G
result partly reuse evidence and must not be counted as independent votes.

The alternative candidates remain in the comparison. Their other-trait
results and missing BMI entries do not exclude a contribution to canine
body condition.

## Outstanding questions

- Can the specific ABC element coordinates and variant overlaps be independently verified?
- What additional evidence supports interpretation of p.R501C?
- Can the seven-gene published set be reconciled with interval annotation?
- Were all seven candidates assessed comparably for human BMI?

## Verified analysis criteria

Source: supplied supplementary methods PDF, page 13.

- SMR: FDR-corrected p < 0.05.
- HEIDI: p > 0.001.
- Coloc: H4 posterior probability > 0.75.

DENND1B meets all three stated criteria using the values in S7.
This is a check of published results against published thresholds,
not a rerun of SMR, HEIDI or colocalisation.

For ABC annotation, the methods considered independent signals and
closely linked variants (r² > 0.8) overlapping mapped regulatory elements.
The main paper identifies a DENND1B promoter in HAP1 cells and an enhancer active in bipolar iPSC neurons. These are published regulatory assignments, not independently verified element coordinates or evidence specific to hypothalamic MC4R-expressing neurons.

The methods also state that orthologous genes within each region were
queried for G2G evidence across the anthropometric and lipid GWAS.
The displayed S7 rows do not provide a complete gene-by-trait results
matrix, so an omitted BMI row cannot be interpreted as a negative result.

## Regulatory evidence source

Main paper, section "DENND1B is associated with canine and human obesity":
https://pmc.ncbi.nlm.nih.gov/articles/PMC7618706/

S7 reports only the tissue of the maximum-scoring ABC element.
The main text supplies the identities and cellular contexts of both
reported elements. The two descriptions are therefore complementary.

## Patient and functional evidence: distinct contributions

### Patient observation
S6 E3 reports one patient homozygous for DENND1B p.R501C.
The supplementary Fig. S4B caption adds that both normal-weight parents
were heterozygous and that exome sequencing found no alternative
genetic diagnosis.

This family pattern is compatible with a recessive hypothesis, but one
family does not establish penetrance or prove that the variant caused
the phenotype. Failure to identify another diagnosis does not exclude
all other genetic explanations.

### Gene-level perturbation
The main paper reports that DENND1B overexpression increased
ligand-induced MC4R internalisation and reduced cAMP signalling.
Knockdown produced the converse effects, with increased cAMP signalling
reported at maximal ligand concentrations.

These experiments support a role for DENND1B in MC4R regulation in the
tested cell system. They do not demonstrate that the canine lead SNP
changes DENND1B expression or acts through this mechanism.

### Variant-specific experiments
Supplementary Fig. S7D reports no significant difference in maximum
MC4R internalisation between wild-type and p.R501C DENND1B
overexpression under the tested conditions.

Fig. S7E reports differences in ligand-induced cAMP responses between
conditions; Fig. S7F reports no significant difference in log EC50.
Different assay endpoints must not be collapsed into a single claim
that every aspect of receptor function changed.

The p.R501C experiment must be distinguished from DENND1B knockdown.
The patient variant, human common BMI variants and canine lead variant
are different genetic observations.

### Relationship to this exercise
Genetic and regulatory evidence explains why DENND1B was a plausible
candidate to investigate. Perturbation experiments provide additional
functional support. This is a distinction between evidence types,
not a reconstruction of the chronological order of the research.

Our work reviews these published observations; it does not reproduce
the cell experiments or establish a clinical variant classification.

Sources:
- S6, cell E3.
- Supplied supplementary PDF, Fig. S4B and Fig. S7D-F captions.
- Main paper, DENND1B/MC4R functional results:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7618706/
