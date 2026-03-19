# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

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
