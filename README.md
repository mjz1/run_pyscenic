# run_pyscenic

[![Tests](https://github.com/mjz1/run_pyscenic/actions/workflows/test.yaml/badge.svg)](https://github.com/mjz1/run_pyscenic/actions/workflows/test.yaml)
[![codecov](https://codecov.io/gh/mjz1/run_pyscenic/branch/main/graph/badge.svg)](https://codecov.io/gh/mjz1/run_pyscenic)

A containerized Python wrapper for the [pySCENIC](https://pyscenic.readthedocs.io/) workflow, enabling scalable gene regulatory network (GRN) inference from single-cell RNA-seq data.

## Overview

This project provides a complete, reproducible implementation of the pySCENIC workflow for inferring gene regulatory networks (GRNs) from single-cell expression matrices. It supports multiple GRN inference methods and includes automatic data quality checks.

**Key Features:**
- 🔄 **Dual GRN Methods**: Choose between [GRNBoost2](https://pyscenic.readthedocs.io/en/latest/tutorials.html) (via Singularity) or [RegDiffusion](https://tuftsbcb.github.io/RegDiffusion/)
- 🧬 **Smart Subsampling**: Run GRN inference on a representative cell subset (random, stratified, or geometry-preserving via [geosketch](https://github.com/brianhie/geosketch)) while scoring all cells with AUCell — with optional Harmony / Scanorama batch correction
- ✅ **Data Quality Checks**: Automatic validation of integer counts and fallback to alternative count matrices
- 🐳 **Fully Containerized**: Docker image with all dependencies and pre-downloaded resources
- 📊 **Complete Pipeline**: Expression matrix export → GRN inference → motif enrichment (ctx) → regulon flattening & filtering → AUCell scoring (per-cell regulon activity) → AUCell binarization for on/off regulon calls
- ⚙️ **Flexible Configuration**: Customizable resources, parameters, and workflow steps

## Installation

### Using Docker

```bash
docker run zatzmanm/run_pyscenic:latest
```

### Using Singularity (HPC)

Convert the Docker image to Singularity format:

```bash
singularity pull docker://zatzmanm/run_pyscenic:latest
```

### Local Installation (uv)

Recommended: use [uv](https://docs.astral.sh/uv/) for fast, reproducible installs (Python 3.10+):

```bash
git clone https://github.com/mjz1/run_pyscenic.git
cd run_pyscenic

# Install UV if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Ensure Python 3.10 is available (installs if missing)
uv python install 3.10
uv venv --python 3.10
source .venv/bin/activate  # Windows: .venv\Scripts\activate
uv pip install -e .

run-pyscenic --help
```

## Quick Start

### Basic Usage with Docker

See available parameters:

```bash
docker run zatzmanm/run_pyscenic:latest
```

Recommended (preferred): mount input data and an output directory so outputs persist

```bash
docker run \
  --gpus all \
  -v /path/to/data:/data \
  -v /path/to/results:/results \
  zatzmanm/run_pyscenic:latest \
    --anndata-path /data/your_data.h5ad \
    --results-dir /results \
    --n-cores 16
```

> Note the `--gpus all` flag is needed to enable GPU support for RegDiffusion in Docker (if available). Omit if not using RegDiffusion.

### Basic Usage with Singularity

```bash
# Pull the image (if not already done)
singularity pull docker://zatzmanm/run_pyscenic:latest
# Run the container
singularity run -e --nv \
  --bind /path/to/data:/data \
  run_pyscenic_latest.sif \
    --anndata-path /data/your_data.h5ad \
    --results-dir /data/results \
    --n-cores 16
```

> Note the `--nv` flag is needed to enable GPU support for RegDiffusion in Singularity (if available). Omit if not using RegDiffusion.

### Using RegDiffusion Instead of GRNBoost2

RegDiffusion offers a significantly faster alternative to GRNBoost2 for GRN inference, providing substantial speedups (especially with GPU acceleration) and making it ideal for large-scale single-cell datasets. For more details, see the [RegDiffusion method](https://tuftsbcb.github.io/RegDiffusion/). Note, these results may differ from GRNBoost2.

```bash
docker run -v /path/to/data:/data zatzmanm/run_pyscenic:latest \
    --anndata-path /data/your_data.h5ad \
    --results-dir /data/results \
    --grn-method regdiff
```

## Input Requirements

### AnnData File Format

Input must be an AnnData `.h5ad` file with:
- **X matrix**: Raw integer count data (required) or normalized data with integer counts in `adata.layers["counts"]`
- **var**: Gene annotations (must include `gene_names` or default index)
- **obs**: Cell metadata (optional)
- Genes are automatically filtered to match the ranking databases used

> ℹ️ **Note:** It is expected that your anndata object is already filtered for high-quality cells. Ensure this is done before running


### Output Files

The pipeline generates:
- `expression.tsv` - Filtered expression matrix (cells × genes)
- `adjacency.tsv` - GRN adjacency matrix (TF-target interactions)
- `regulons.csv` - Inferred regulons with motif evidence
- `regulons_flat.tsv` - Flattened TF → target gene table (one row per TF-motif-target triple, unfiltered)
- `regulons_flat_filtered.tsv` - Quality-filtered version of the flat table (NES ≥ 3, activating, direct annotation, ≥ 10 targets)
- `regulons_flat_dedup.tsv` - Deduplicated unfiltered table (one row per unique TF-target pair, best-NES motif retained)
- `regulons_flat_filtered_dedup.tsv` - Deduplicated filtered table (one row per unique TF-target pair, best-NES motif retained)
- `auc_mtx.csv` - Per-cell regulon activity scores (AUC matrix)
- `auc_mtx_binarized.csv` - Binarized per-cell regulon activity (0/1 calls)
- `auc_binarization_thresholds.csv` - Adaptive thresholds used for binarization (per regulon)

When [subsampling](#subsampling) is enabled, the single `expression.tsv` is replaced by two matrices and additional metadata:
- `expression_subsampled.tsv` - Subsampled expression matrix used for GRN + ctx
- `expression_full.tsv` - Full expression matrix used for AUCell scoring
- `subsample_metadata.json` - Run parameters, cell counts, per-group breakdowns (stratified), and embedding metadata (geosketch)
- `subsampled_barcodes.txt` - One selected cell barcode per line
- `subsample_embedding.npy` - Cached embedding array (geosketch only)

## Command-Line Options

```
--anndata-path PATH              Path to input .h5ad file [required]
--results-dir PATH               Output directory (default: ./pyscenic_results)
--resource-dir PATH              Directory with pySCENIC resources 
                                 (default: /opt/pyscenic_resources)
--n-cores N                      Number of CPU cores (default: 16)
--seed SEED                      Random seed for reproducibility (default: 42)
--log-level LEVEL                Logging verbosity: DEBUG|INFO|WARNING|ERROR|CRITICAL 
                                 (default: INFO)

GRN Inference:
--grn-method METHOD              grnboost2 (default) or regdiff
--skip-grn                       Skip GRN step (use existing adjacency.tsv)

RegDiffusion Options (with --grn-method regdiff):
--regdiff-percentile PERCENTILE  Weight threshold for edge retention (default: 50)

Motif Enrichment (ctx):
--skip-ctx                       Skip motif enrichment step
--ranking-files FILE             Custom ranking database files (can specify multiple)
--motif-f FILE                   Custom motif annotations file
--tf-file FILE                   Custom TF list file
--overwrite                      Regenerate outputs even if they exist

Subsampling:
--subsample-method METHOD         Subsampling strategy: none (default), random, stratified,
                                  geosketch
--subsample-n N                  Target cell count for subsampling (required if method != none)
--subsample-annotation COL       adata.obs column for stratified sampling
                                  (required if method = stratified)
--subsample-min-per-group N      Minimum cells per group for stratified sampling (default: 50)
--subsample-embedding-key KEY    adata.obsm key with pre-computed embedding for geosketch
                                  (overrides batch correction)
--subsample-batch-correction-method METHOD
                                  Batch correction for geosketch embedding: none (default),
                                  harmony, scanorama
--subsample-batch-key COL        adata.obs batch column (required if batch-correction != none)
--subsample-n-pca-components N   PCA components for geosketch embedding (default: 50)

Regulon Flattening:
--skip-flatten                   Skip regulon flattening/filtering step

AUCell (activity scoring):
--skip-aucell                    Skip AUCell scoring step (use existing auc_mtx.csv)

Binarization (AUCell on/off calls):
--skip-binarize                  Skip binarization step (use existing auc_mtx_binarized.csv)
```

## Workflow Steps

The pipeline runs these steps in sequence:

### 1. Expression Matrix Export (`write_expression_tsv`)
- Loads AnnData file
- Validates that count matrix contains integers (warns if normalized)
- Filters to genes present in all ranking databases
- Removes non-expressed genes
- When subsampling is enabled, writes two matrices:
  - `expression_subsampled.tsv` — subset of cells for GRN/ctx steps
  - `expression_full.tsv` — all cells for AUCell scoring
- Without subsampling, exports a single `expression.tsv`

### 2. GRN Inference
**Option A: GRNBoost2** (default)
- Runs via pySCENIC wrapper
- Requires Singularity with GRNBoost2 image
- Returns TF-target adjacency matrix

**Option B: RegDiffusion**
- Uses RegDiffusion for GRN inference (much faster with GPU acceleration)
- Filters edges by weight percentile
- Returns TF-target adjacency matrix

### 3. Motif Enrichment (ctx)
- Computes TF-motif associations
- Links TF-target interactions to known motifs
- Produces final regulon predictions with confidence scores

### 4. Regulon Flattening & Filtering (`run_flatten_regulons`)
- Parses the multi-header `regulons.csv` into a flat, analysis-ready table with one row per TF–motif–target triple
- Columns: TF, MotifID, AUC, NES, MotifSimilarityQvalue, OrthologousIdentity, Annotation, Context, RankAtMax, n_targets, target, weight
- Writes `regulons_flat.tsv` (unfiltered), `regulons_flat_filtered.tsv` (filtered with defaults: NES ≥ 3.0, activating context, direct/orthologous-direct annotation, ≥ 10 targets per TF-motif pair), `regulons_flat_dedup.tsv` (deduplicated unfiltered), and `regulons_flat_filtered_dedup.tsv` (deduplicated filtered) — dedup tables have one row per unique TF-target, keeping the best-NES motif
- `flatten_regulons()`, `filter_regulons()`, and `deduplicate_regulons()` are exposed as standalone functions for post-hoc re-filtering with custom parameters
- Can be skipped with `--skip-flatten`

### 5. AUCell (per-cell regulon activity)
- Computes per-cell regulon activity (AUC) from the `regulons.csv` produced by the `ctx` step
- Produces `auc_mtx.csv` (cells × regulons) containing activity scores usable for downstream visualization and clustering
- Can be skipped with `--skip-aucell` when you only need regulons or already have `auc_mtx.csv`

### 6. AUCell Binarization (per-cell regulon on/off)
- Converts continuous AUCell scores into binary on/off calls per regulon using adaptive thresholds
- Produces `auc_mtx_binarized.csv` plus `auc_binarization_thresholds.csv` (per-regulon thresholds used)
- Can be skipped with `--skip-binarize` when you already have binarized outputs or prefer continuous scores

## Subsampling

GRN inference (GRNBoost2/RegDiffusion) and motif enrichment (ctx) are the most
computationally expensive pipeline steps and their runtime scales with cell
count. For large datasets you can run these steps on a **representative
subset** of cells while AUCell scoring — which is comparatively fast —
still runs on **all** cells. This keeps regulon activity scores available for
every cell in your dataset.

Enable subsampling with `--subsample-method` and `--subsample-n`.

### How It Works (Two-Matrix Flow)

When subsampling is enabled the pipeline writes **two** expression matrices
instead of the usual single `expression.tsv`:

```
expression_subsampled.tsv   ← subset (used by GRN + ctx)
expression_full.tsv         ← all cells (used by AUCell)
```

1. The AnnData file is loaded and gene-filtered once.
2. The chosen subsampling strategy selects cell indices.
3. `expression_subsampled.tsv` (cells in the subset) feeds GRN inference and
   motif enrichment.
4. `expression_full.tsv` (all cells) feeds AUCell, so per-cell regulon
   activities are computed genome-wide.
5. Both matrices share the same gene columns, guaranteeing compatibility.

If `--subsample-n` is greater than or equal to the number of cells after
gene-filtering,  subsampling is skipped automatically and a single
`expression.tsv` is written (with a warning).

### Strategies

#### `random`

Uniform random sampling without replacement. Fast, simple, and reproducible
(controlled by `--seed`). Best when you have no strong cell-type imbalance.

#### `stratified`

Preserves group proportions from an `adata.obs` annotation column (e.g.
`cell_type`). The allocation algorithm:

1. Computes each group's ideal share proportional to its fraction of the total.
2. Guarantees a minimum of `--subsample-min-per-group` cells per group (default
   50), or the entire group if it is smaller than the minimum.
3. Redistributes leftover budget proportionally among larger groups.
4. Adjusts for rounding so the total matches `--subsample-n` exactly.

This avoids underrepresenting rare cell types. Per-group original and sampled
counts are logged and recorded in `subsample_metadata.json`.

#### `geosketch`

[Geometric sketching](https://github.com/brianhie/geosketch) selects a subset
that preserves the geometry of a low-dimensional embedding, so transcriptomic
diversity is maximally retained even when cell types form overlapping clusters.
Runs on an embedding matrix (cells × dimensions) — see
[Batch-Aware Geosketch Embedding](#batch-aware-geosketch-embedding) for how the
embedding is resolved.

### Usage Examples

**Random subsampling** — select 5 000 cells at random:

```bash
run-pyscenic \
    --anndata-path data.h5ad \
    --results-dir results/ \
    --subsample-method random \
    --subsample-n 5000
```

**Stratified subsampling** — 5 000 cells, preserving cell-type proportions
with at least 100 cells per type:

```bash
run-pyscenic \
    --anndata-path data.h5ad \
    --results-dir results/ \
    --subsample-method stratified \
    --subsample-n 5000 \
    --subsample-annotation cell_type \
    --subsample-min-per-group 100
```

**Geosketch** — geometry-preserving subset with a default naive PCA embedding:

```bash
run-pyscenic \
    --anndata-path data.h5ad \
    --results-dir results/ \
    --subsample-method geosketch \
    --subsample-n 5000
```

**Geosketch with Harmony batch correction:**

```bash
run-pyscenic \
    --anndata-path data.h5ad \
    --results-dir results/ \
    --subsample-method geosketch \
    --subsample-n 5000 \
    --subsample-batch-correction-method harmony \
    --subsample-batch-key sample_id
```

**Geosketch with Scanorama batch correction:**

```bash
run-pyscenic \
    --anndata-path data.h5ad \
    --results-dir results/ \
    --subsample-method geosketch \
    --subsample-n 5000 \
    --subsample-batch-correction-method scanorama \
    --subsample-batch-key sample_id
```

**Geosketch with a pre-computed embedding stored in `adata.obsm`:**

```bash
run-pyscenic \
    --anndata-path data.h5ad \
    --results-dir results/ \
    --subsample-method geosketch \
    --subsample-n 5000 \
    --subsample-embedding-key X_pca_harmony
```

### Batch-Aware Geosketch Embedding

When using `geosketch`, the embedding used for geometric sketching is resolved
with the following precedence:

1. **Pre-computed embedding** — if `--subsample-embedding-key` is set, the
   corresponding `adata.obsm` slot is used directly. This is the fastest
   option if you already have a batch-corrected or curated embedding in your
   AnnData object.
2. **In-pipeline batch-corrected embedding** — if
   `--subsample-batch-correction-method` is `harmony` or `scanorama`, the
   pipeline normalises a copy of the count matrix, runs PCA
   (`--subsample-n-pca-components`, default 50), then integrates across
   batches identified by `--subsample-batch-key`. The original `adata.X`
   (raw counts) is never modified.
3. **Naive PCA** — fallback when neither of the above is specified. The same
   normalise → log1p → PCA path runs, but without batch correction.

The resolved embedding is cached to `subsample_embedding.npy` in the results
directory. On subsequent runs the cache is reused automatically (pass
`--overwrite` to recompute).

> **Tip:** If `--subsample-embedding-key` is set and `--subsample-batch-correction-method` is also
> specified, the pre-computed key takes priority and a warning is logged.

### Subsampling Output Files

When subsampling is active the following additional files are written to
`--results-dir`:

| File | Description |
|---|---|
| `expression_subsampled.tsv` | Cells × genes matrix for the selected subset (input to GRN + ctx). |
| `expression_full.tsv` | Cells × genes matrix for all cells (input to AUCell). |
| `subsample_metadata.json` | JSON with method, target/actual cell counts, seed, total cells/genes, per-group breakdowns (stratified), and embedding source/dimensions (geosketch). |
| `subsampled_barcodes.txt` | One selected cell barcode per line, for external use or auditing. |
| `subsample_embedding.npy` | Cached embedding array — geosketch only. |

<details>
<summary>Example <code>subsample_metadata.json</code> (stratified)</summary>

```json
{
  "method": "stratified",
  "target_n": 5000,
  "actual_n": 5000,
  "seed": 42,
  "total_cells_original": 25000,
  "total_genes": 1800,
  "annotation_column": "cell_type",
  "per_group_counts": {
    "T_cell": 10000,
    "B_cell": 8000,
    "Monocyte": 5000,
    "NK": 2000
  },
  "per_group_subsampled": {
    "T_cell": 2000,
    "B_cell": 1600,
    "Monocyte": 1000,
    "NK": 400
  }
}
```

</details>

<details>
<summary>Example <code>subsample_metadata.json</code> (geosketch + harmony)</summary>

```json
{
  "method": "geosketch",
  "target_n": 5000,
  "actual_n": 5000,
  "seed": 42,
  "total_cells_original": 25000,
  "total_genes": 1800,
  "n_pca_components": 50,
  "embedding_source": "harmony",
  "batch_correction_method": "harmony",
  "batch_key": "sample_id",
  "embedding_dimensions": 50
}
```

</details>

### Validation and Warnings

The pipeline validates subsampling options before any computation:

| Check | Behaviour |
|---|---|
| `--subsample-n` missing when method ≠ `none` | Error |
| `--subsample-annotation` missing when method = `stratified` | Error |
| `--subsample-batch-key` missing when batch correction ≠ `none` | Error |
| `--subsample-batch-key` not in `adata.obs` | Error |
| `--subsample-embedding-key` not in `adata.obsm` | Error |
| `--subsample-n` < 10 000 | Warning (GRN quality may be impacted) |
| `--subsample-n` ≥ total cells after filtering | Warning; subsampling skipped |
| Batch correction options set with method ≠ `geosketch` | Warning (options ignored) |
| `--subsample-embedding-key` set alongside batch correction | Warning (embedding key takes priority) |
| Embedding has < 5 or > 100 dimensions | Warning |
| A stratified group falls below `--subsample-min-per-group` after sampling | Warning (per group) |

### Choosing a Strategy

- **`random`** is the simplest option and works well when cell types are
  roughly balanced.
- **`stratified`** is recommended when your dataset has significant cell-type
  imbalance and you want to ensure rare populations are represented in the GRN.
- **`geosketch`** provides the best transcriptomic coverage because it selects
  cells that span the full embedding space. Use this for heterogeneous datasets
  or multi-sample experiments, especially with batch correction enabled.

### Optional Dependencies

The base install supports `random` and `stratified` strategies (they use numpy
and scanpy, which are already required). Additional packages are needed only
for specific features:

| Package | Extra | Required for |
|---|---|---|
| `geosketch` | `pip install run-pyscenic[geosketch]` | `--subsample-method geosketch` |
| `harmonypy` | `pip install run-pyscenic[harmony]` | `--subsample-batch-correction-method harmony` |
| `scanorama` | `pip install run-pyscenic[scanorama]` | `--subsample-batch-correction-method scanorama` |

Install all subsampling extras at once:

```bash
pip install run-pyscenic[subsampling]
```

Or with uv:

```bash
uv pip install -e ".[subsampling]"
```

## Configuration

### Resource Management

The Docker image includes pre-downloaded resources for **human (hg38)**:
- Ranking databases: 2 feather format files (~1.3 GB each)
- Motif annotations: `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`
- TF list: `allTFs_hg38.txt`

**Custom Resources:**
To use different ranking databases or TF lists:

```bash
docker run -v /path/to/custom_resources:/custom_res \
  -v /path/to/data:/data \
  zatzmanm/run_pyscenic:latest \
    --anndata-path /data/data.h5ad \
    --resource-dir /custom_res \
    --ranking-files /custom_res/custom_rankings1.feather \
                    /custom_res/custom_rankings2.feather \
    --tf-file /custom_res/custom_tfs.txt \
    --motif-f /custom_res/custom_motifs.tbl
```


## References

- [pySCENIC Documentation](https://pyscenic.readthedocs.io/)
- [RegDiffusion Downstream Analysis](https://tuftsbcb.github.io/RegDiffusion/downstream_with_pyscenic.html)
- [CisTopic Resources](https://resources.aertslab.org/cistarget/)

## Citation

If you use this wrapper, please cite:

- **pySCENIC**: Aibar S, González-Blas CB, Moerman T, Huynh-Thu VA, Imrichova H, Hulselmans G, Rambow F, Marine JC, Geurts P, Aerts J, van den Oord J, Atak ZK, Wouters J, Aerts S. SCENIC: single-cell regulatory network inference and clustering. Nat Methods. 2017 Nov;14(11):1083-1086. doi: 10.1038/nmeth.4463. Epub 2017 Oct 9. PMID: 28991892; PMCID: PMC5937676.
- **GRNBoost2**: Moerman T, Aibar Santos S, Bravo González-Blas C, Simm J, Moreau Y, Aerts J, Aerts S. GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks. Bioinformatics. 2019 Jun 1;35(12):2159-2161. doi: 10.1093/bioinformatics/bty916. PMID: 30445495.
- **RegDiffusion**: Zhu H, Slonim D. From Noise to Knowledge: Diffusion Probabilistic Model-Based Neural Inference of Gene Regulatory Networks. J Comput Biol. 2024 Nov;31(11):1087-1103. doi: 10.1089/cmb.2024.0607. Epub 2024 Oct 10. PMID: 39387266; PMCID: PMC11698671.

