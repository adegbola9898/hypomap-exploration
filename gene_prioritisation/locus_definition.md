# Working locus definition

## Direct evidence from the supplied S3 workbook
Source: references/s3.xlsx, sheet S3, row 6.
The interpretation of region selection comes from the caption in row 1.

- Lead variant: 7:5004016.
- LD region (r² >= 0.8): 7:3062016-5279201.
- Haplotype entry: N/A.
- Published refined region: 7:3062016-5279201.
- Published approximate span: 2.2 Mb.

The refined interval equals the LD-defined interval for this row.
It is a published region of interest, not a demonstrated causal interval.

## Checks still required
- Confirm in the paper that this is the DENND1B-associated lead variant.
- Verify the reference genome assembly.
- Verify coordinate conventions before computational interval overlap.
- Establish the complete published candidate list and its selection rule.
- Check candidate coordinates using assembly-compatible annotation.

## Source provenance
Workbook supplied by the project owner.
Original download URL and download date have not been recorded.
The checksum identifies this supplied copy; it does not independently
authenticate the file against the publisher's version.
See source_checksums.sha256.

## Source review: supplementary methods and main paper

The following findings resolve some of the checks listed above.

- Supplementary methods, page 3: canine GWAS data used CanFam3.1,
  assembly accession GCA_000002285.2.
- Supplementary methods, page 8: initial regional gene lists used
  Ensembl BioMart archive release 104 (May 2021).
- The authors supplemented canine annotation with protein-coding genes
  from syntenic human GRCh38 regions and checked alternative canine builds.
- Supplementary text, page 19: rs24430444 lies within DENND1B;
  LHX9 is also mentioned in the mapped region.
- Main paper: the DENND1B intronic variant is 7:5004016:T>C,
  linking this association to S3 row 6.
- Main paper: seven genes are listed at the chromosome 7 locus:
  NR5A2, PTPRC, ATP6V1G3, NEK7, LHX9, DENND1B and CRB1.

Main paper:
https://pmc.ncbi.nlm.nih.gov/articles/PMC7618706/

## Working candidate set and limitations

Use the seven genes explicitly listed in the main paper as the
published comparison set.

This is a reconstruction of the published candidate list, not an
independent annotation-based verification of interval membership.

The supplementary paragraph mentions DENND1B and LHX9 without explicitly
stating that these are the only genes. The reason for the difference in
coverage between that paragraph and the main text remains unresolved.

S5/S6 table membership does not define the complete locus candidate set.
S7 results must retain their phenotype: BMI, WHR, HDL or triglycerides.
Absence from a table is not evidence that a gene has no relevant effect.

Coordinate conventions and per-gene annotation checks remain pending.
The supplementary PDF was reviewed separately and is not yet stored
in this repository.
