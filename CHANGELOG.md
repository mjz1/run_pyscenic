# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

## [0.2.1] - 2026-03-20
### Added
- New regulon flattening and filtering step between ctx and AUCell.
- `flatten_regulons()` parses multi-header `regulons.csv` into a flat
  TF → target gene table with columns: TF, MotifID, AUC, NES,
  MotifSimilarityQvalue, OrthologousIdentity, Annotation, Context,
  RankAtMax, n_targets, target, weight.
- `filter_regulons()` applies configurable quality filters (NES ≥ 3.0,
  activating context, direct/orthologous-direct annotation, ≥ 10 targets)
  and can be re-invoked by users with custom parameters.
- `deduplicate_regulons()` collapses any flat table to one row per unique
  TF-target pair (keeping the best-NES motif) for clean gene-set scoring.
- Pipeline writes `regulons_flat.tsv` (unfiltered),
  `regulons_flat_filtered.tsv` (default-filtered),
  `regulons_flat_dedup.tsv` (deduplicated unfiltered), and
  `regulons_flat_filtered_dedup.tsv` (deduplicated filtered) to results
  directory.
- `--skip-flatten` CLI option to bypass the flattening step.
- 6 new tests covering unfiltered output, default and custom filtering,
  deduplication, file writing, and skip/overwrite logic.

## [0.2.0] - 2026-03-19
### Added
- Configurable cell subsampling for GRN inference and motif enrichment steps
  while AUCell scoring runs on all cells (two-matrix flow).
- New `subsample.py` module with three strategies: `random`, `stratified`,
  and `geosketch`.
- Stratified sampling preserves group proportions with configurable
  `--subsample-min-per-group` floor.
- Batch-aware geosketch embedding resolution via
  `_resolve_geosketch_embedding()` supporting Harmony and Scanorama
  integration methods.
- Embedding caching to `subsample_embedding.npy` for reproducibility.
- Subsampling metadata written to `subsample_metadata.json`.
- CLI options: `--subsample-method`, `--subsample-n`,
  `--subsample-annotation`, `--subsample-min-per-group`,
  `--subsample-embedding-key`, `--subsample-batch-correction-method`,
  `--subsample-batch-key`, `--subsample-n-pca-components`.
- Optional dependency groups for `geosketch`, `harmony`, `scanorama`, and
  a combined `subsampling` extra.
- Comprehensive test suite for all subsampling paths (37 new tests).

### Changed
- Refactored `write_expression_tsv()` into `_load_and_filter_adata()` and
  `_write_expression_tsv()` for reuse in the two-matrix flow.
- Pipeline now writes `expression_subsampled.tsv` (for GRN/ctx) and
  `expression_full.tsv` (for AUCell) when subsampling is enabled.

### Removed
- `--max-cells` CLI option (replaced by the subsampling system).

## [v0.1.2] - 2026-01-29
### Added
- Expose AUCell threshold options (`--rank-threshold`, `--auc-threshold`, `--nes-threshold`) in the wrapper.

### Tests
- Add coverage for AUCell threshold passthrough.

## [v0.1.1] - 2025-12-12
### Added
- Add binarization for AUCell matrices.

### Documentation
- Update README.

### Other
- Version bump to v0.1.1.

## [v0.1.0] - 2025-12-09
### Added
- Initial release.
- Add RegDiffusion support and Docker updates.
- Add uv lockfile.

### Documentation
- Add README and copy updates.
