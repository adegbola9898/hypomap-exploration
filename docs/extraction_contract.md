# Gene count extraction contract

## Purpose
Extract DENND1B and MC4R counts for all source nuclei, preserving
metadata and source order. Apply neuronal selection downstream.

## Configuration
config/analysis.json defines source identity, matrix paths and genes.
Its extraction section defines extraction coverage and metadata.
Its selection and grouping fields define the exploratory analysis.
Publication comparison settings remain separate and provisional.

## Low-level interface
read_csr_columns(matrix_group, column_indices, *, max_nnz_per_chunk)

Inputs:
- An HDF5 CSR group containing data, indices and indptr.
- A nonempty list of unique valid integer column indices.
- A positive integer limit on stored entries read per chunk.

Output:
- A NumPy array with one row per source nucleus.
- Columns in exactly the requested order.
- Implicit sparse zeros preserved.
- Duplicate sparse coordinates summed.
- No normalisation or thresholding.

Memory behaviour:
- Read stored sparse entries in bounded chunks.
- Permit chunk boundaries within rows.
- Handle empty rows.
- Never create a dense version of the full expression matrix.
- Only requested columns may form a dense output.

Validation:
- Validate CSR dimensions, pointers, indices and array lengths.
- Reject invalid column requests and chunk sizes.
- Reject non-finite, negative or fractional stored count values.
- Raise informative errors instead of returning partial results.

## Atlas wrapper
- Verify source checksum against configuration.
- Require each requested gene to map to exactly one feature.
- Require unique cell identifiers and aligned metadata lengths.
- Require 433,369 source nuclei for the configured dataset.
- Extract all nuclei in source order.
- Preserve all configured metadata fields.
- Fail clearly if a required field is absent.
- Report and preserve missing metadata values without inventing
  replacements or silently dropping nuclei.
- Validate the extracted gene counts before writing output.

## Output
data/derived/all_nuclei_gene_counts.tsv.gz

Columns:
- cell_id
- source_row
- All configured preserve_obs_fields
- DENND1B_count
- MC4R_count

Also write provenance/extraction_manifest.json containing:
- Source and configuration checksums.
- Code commit and whether the working tree was dirty.
- Gene positions and extracted dimensions.
- Validation results and missing-metadata summary.
- Output checksum, software versions and completion time.

Write final artifacts only after successful validation.
Never modify the source atlas.

## Downstream exploratory analysis
Select C0_named == "C0-3 Neuron".
Require 166,475 selected nuclei.
Summarise the seven configured C1 populations.
Primary detection threshold: raw count >= 1.
Sensitivity threshold: raw count >= 2.
These are exploratory definitions, not verified publication settings.

## Publication comparison
S9 includes all-cell, neuronal and receptor-expressing subgroups.
Verify the published subgroup definitions before claiming agreement.

S10 contains 452 rows with both C3 and C4 identifiers.
Preserve both annotation levels. Establish and validate the mapping
before comparing results. Check cluster counts and cell membership
where possible; do not assume a single annotation column suffices.

Compare only the selected human DENND1B and MC4R detection metrics
initially. Average-expression columns and other genes are deferred.

## Required tests
- Known-answer CSR matrix with empty rows and explicit zeros.
- Requested columns supplied in reversed source order.
- Chunk boundaries within a row.
- Duplicate sparse coordinates.
- Invalid columns and chunk sizes.
- Malformed CSR pointers, indices and array lengths.
- Non-finite, negative and fractional values.
- Wrapper rejects absent or ambiguous genes and duplicate cell IDs.
- Wrapper preserves source order and metadata alignment.

## Exclusions
The extractor performs no neuronal filtering, detection thresholding,
percentage calculations, plotting, clustering or statistical inference.

## Wrapper implementation boundary

extract_atlas_counts(handle, config, *, max_nnz_per_chunk)

Inputs:
- An already-open, read-only HDF5 file.
- The parsed analysis configuration.
- The stored-entry chunk limit.

Output:
- A pandas DataFrame containing all source nuclei in source order.
- Columns ordered as cell_id, source_row, configured metadata,
  followed by configured gene count columns.
- A mapping from requested gene names to source column positions.

Validation before expression extraction:
- Requested genes are nonempty and unique.
- Feature identifiers match the matrix column dimension.
- Each requested gene has exactly one matching feature.
- Cell identifiers match the matrix row dimension and are unique.
- Source row count matches extraction.expected_nuclei.
- Every configured metadata field exists and has one value per nucleus.
- Metadata names do not collide with generated output column names.

This function does not apply neuronal selection or write files.
It must not mutate the supplied configuration.

The command-line orchestration layer will verify the source checksum
before calling this function, then write validated outputs and their
provenance. The wrapper alone does not authenticate the source file.

String metadata must be decoded to text. Missing values must remain
missing. Categorical encodings must be decoded explicitly or rejected
with a clear error; category codes must not silently become labels.
