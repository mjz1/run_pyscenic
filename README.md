# run_pyscenic

A containerized Python wrapper for the [pySCENIC](https://pyscenic.readthedocs.io/) workflow, enabling scalable gene regulatory network (GRN) inference from single-cell RNA-seq data.

## Overview

This project provides a complete, reproducible implementation of the pySCENIC workflow for inferring gene regulatory networks (GRNs) from single-cell expression matrices. It supports multiple GRN inference methods and includes automatic data quality checks.

**Key Features:**
- 🔄 **Dual GRN Methods**: Choose between [GRNBoost2](https://pyscenic.readthedocs.io/en/latest/tutorials.html) (via Singularity) or [RegDiffusion](https://tuftsbcb.github.io/RegDiffusion/)
- ✅ **Data Quality Checks**: Automatic validation of integer counts and fallback to alternative count matrices
- 🐳 **Fully Containerized**: Docker image with all dependencies and pre-downloaded resources
- 📊 **Complete Pipeline**: Expression matrix export → GRN inference → motif enrichment (ctx) → AUCell scoring (per-cell regulon activity)
- ⚙️ **Flexible Configuration**: Customizable resources, parameters, and workflow steps

## Installation

### Using Docker

```bash
docker run zatzmanm/run_pyscenic:main
```

### Using Singularity (HPC)

Convert the Docker image to Singularity format:

```bash
singularity pull docker://zatzmanm/run_pyscenic:main
```

### Local Installation (uv)

Recommended: use [uv](https://docs.astral.sh/uv/) for fast, reproducible installs (Python 3.10+):

```bash
git clone https://github.com/mjz1/run_pyscenic.git
cd run_pyscenic

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

**Data Quality Notes:**
- The X matrix should contain **raw integer counts**, not normalized data
- If X contains normalized/float data, the script checks for an integer counts layer at `adata.layers["counts"]` and uses that instead
- Genes are automatically filtered to match the ranking databases used

### Output Files

The pipeline generates:
- `expression.tsv` - Filtered expression matrix (cells × genes)
- `adjacency.tsv` - GRN adjacency matrix (TF-target interactions)
- `regulons.csv` - Inferred regulons with motif evidence
- `auc_mtx.csv` - Per-cell regulon activity scores (AUC matrix)

## Command-Line Options

```
--anndata-path PATH              Path to input .h5ad file [required]
--results-dir PATH               Output directory (default: ./pyscenic_results)
--resource-dir PATH              Directory with pySCENIC resources 
                                 (default: /opt/pyscenic_resources)
--n-cores N                      Number of CPU cores (default: 16)
--seed SEED                      Random seed for reproducibility (default: 42)
--max-cells N                    Process only first N cells (optional)
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

AUCell (activity scoring):
--skip-aucell                    Skip AUCell scoring step (use existing auc_mtx.csv)
```

## Workflow Steps

The pipeline runs these steps in sequence:

### 1. Expression Matrix Export (`write_expression_tsv`)
- Loads AnnData file
- Validates that count matrix contains integers (warns if normalized)
- Filters to genes present in all ranking databases
- Optionally limits to N cells
- Removes non-expressed genes
- Exports as TSV format

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

