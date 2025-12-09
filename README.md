# run_pyscenic

A containerized Python wrapper for the [pySCENIC](https://pyscenic.readthedocs.io/) workflow, enabling scalable gene regulatory network (GRN) inference from single-cell RNA-seq data.

## Overview

This project provides a complete, reproducible implementation of the pySCENIC workflow for inferring gene regulatory networks (GRNs) from single-cell expression matrices. It supports multiple GRN inference methods and includes automatic data quality checks.

**Key Features:**
- 🔄 **Dual GRN Methods**: Choose between [GRNBoost2](https://pyscenic.readthedocs.io/en/latest/tutorials.html) (via Singularity) or [RegDiffusion](https://tuftsbcb.github.io/RegDiffusion/)
- ✅ **Data Quality Checks**: Automatic validation of integer counts and fallback to alternative count matrices
- 🐳 **Fully Containerized**: Docker image with all dependencies and pre-downloaded resources
- 📊 **Complete Pipeline**: Expression matrix export → GRN inference → motif enrichment (ctx)
- ⚙️ **Flexible Configuration**: Customizable resources, parameters, and workflow steps

## Installation

### Using Docker (Recommended)

Pull the pre-built image from Docker Hub:

```bash
docker pull zatzmanm/run_pyscenic:main
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

```bash
docker run -v /path/to/data:/data zatzmanm/run_pyscenic:main \
  run_pyscenic.py \
    --anndata-path /data/your_data.h5ad \
    --results-dir /data/results \
    --n-cores 16
```

### Basic Usage with Singularity

```bash
singularity run \
  --bind /path/to/data:/data \
  run_pyscenic.sif \
  run_pyscenic.py \
    --anndata-path /data/your_data.h5ad \
    --results-dir /data/results \
    --n-cores 16
```

### Using RegDiffusion Instead of GRNBoost2

```bash
docker run -v /path/to/data:/data zatzmanm/run_pyscenic:main \
  run_pyscenic.py \
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
  zatzmanm/run_pyscenic:main \
  run_pyscenic.py \
    --anndata-path /data/data.h5ad \
    --resource-dir /custom_res \
    --ranking-files /custom_res/custom_rankings1.feather \
                    /custom_res/custom_rankings2.feather \
    --tf-file /custom_res/custom_tfs.txt \
    --motif-f /custom_res/custom_motifs.tbl
```

### Performance Tuning

- **`--n-cores`**: Set to number of available CPU cores for parallelization
- **`--max-cells`**: Reduce dataset size for testing/quick runs
- **`--regdiff-percentile`**: Lower values (e.g., 25) = more edges; higher values (e.g., 75) = fewer edges

## Examples

### Run Complete Pipeline on Sample Data

```bash
docker run -v ~/data:/data zatzmanm/run_pyscenic:main \
  run_pyscenic.py \
    --anndata-path /data/pbmc_5k.h5ad \
    --results-dir /data/pyscenic_results \
    --n-cores 32 \
    --seed 123
```

### Process Only First 1000 Cells (Quick Test)

```bash
docker run -v ~/data:/data zatzmanm/run_pyscenic:main \
  run_pyscenic.py \
    --anndata-path /data/pbmc_5k.h5ad \
    --results-dir /data/test_results \
    --max-cells 1000 \
    --n-cores 8
```

### Use RegDiffusion with Custom Motif Enrichment Only (Skip GRN)

```bash
docker run -v ~/data:/data zatzmanm/run_pyscenic:main \
  run_pyscenic.py \
    --anndata-path /data/pbmc_5k.h5ad \
    --results-dir /data/results \
    --skip-grn \  # Use existing adjacency.tsv
    --ranking-files /data/my_rankings.feather \
    --motif-f /data/my_motifs.tbl
```

### Run on HPC with Singularity

```bash
#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=08:00:00

singularity run \
  --bind /scratch:/scratch \
  /path/to/run_pyscenic.sif \
  run_pyscenic.py \
    --anndata-path /scratch/data.h5ad \
    --results-dir /scratch/results \
    --n-cores $SLURM_CPUS_PER_TASK
```

## Troubleshooting

### Data Quality Issues

**Warning: "Main matrix does not contain integer counts"**
- The input data appears to be normalized (log-transformed, scaled, etc.)
- Solution: Check if `adata.layers["counts"]` exists; if yes, the script will use it automatically
- Alternative: Preprocess your data to obtain raw counts

**Error: "Missing required files"**
- Ranking databases or motif files not found
- Solution: Ensure resource files exist in the resource directory, or specify custom paths with `--ranking-files` and `--motif-f`

### Performance Issues

**Memory/CPU bottleneck?**
- Reduce `--n-cores` if hitting memory limits
- Use `--max-cells` to test on subset of data first
- For large datasets, increase Docker memory: `-m 256g`

**Singularity pull fails?**
- Ensure Singularity version >= 3.0
- Try specifying the tag explicitly: `singularity pull docker://zatzmanm/run_pyscenic:main`

## Building Locally

### Build Docker Image

```bash
git clone https://github.com/mjz1/run_pyscenic.git
cd run_pyscenic
docker build -t myrepo/run_pyscenic:latest .
```

The Dockerfile uses multi-stage builds and GitHub Actions for automated CI/CD.

## References

- [pySCENIC Documentation](https://pyscenic.readthedocs.io/)
- [RegDiffusion Downstream Analysis](https://tuftsbcb.github.io/RegDiffusion/downstream_with_pyscenic.html)
- [CisTopic Resources](https://resources.aertslab.org/cistarget/)

## Citation

If you use this wrapper, please cite:

- **pySCENIC**: Aibar et al. (2017). SCENIC: single-cell regulatory network inference and clustering. *Nature Methods*, 14(11), 1083-1086.
- **GRNBoost2**: Van de Sande et al. (2020). A robust gene regulatory network inference method...
- **RegDiffusion**: [See RegDiffusion paper](https://tuftsbcb.github.io/RegDiffusion/)

## License

This wrapper is open-source. Check the LICENSE file for details.

## Support

For issues or questions:
- Open a GitHub issue: https://github.com/mjz1/run_pyscenic/issues
- Check existing issues for common problems
- Include logs (use `--log-level DEBUG` for detailed output)

