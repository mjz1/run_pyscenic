#!/usr/bin/env python3
import logging
import os
import shutil
import subprocess
import sys
import tempfile

import click
import pandas as pd
import scanpy as sc

# If no arguments provided, show help
if len(sys.argv) == 1:
    sys.argv.append("--help")


def _atomic_move(src_path, dst_path):
    """Atomically move a file from src to dst using temporary file + rename.

    This ensures that the destination file is only updated when the source
    is complete, preventing partial/empty files from being treated as valid.

    Args:
        src_path: Source file path
        dst_path: Destination file path
    """
    # Ensure the destination directory exists
    os.makedirs(os.path.dirname(dst_path), exist_ok=True)
    # Atomic rename (on most filesystems, rename is atomic)
    shutil.move(src_path, dst_path)


def _ensure_files_exist(paths):
    missing = [p for p in paths if p and not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(f"Missing required files: {missing}")


def _real(path):
    return os.path.realpath(path) if path else path


def write_expression_tsv(
    anndata_path, results_dir, max_cells, overwrite, ranking_files
):
    """Load AnnData, filter to genes in ranking databases, export cells x genes TSV; return path. Skip if present unless overwrite.

    Args:
        anndata_path: Path to input .h5ad file
        results_dir: Directory to write output files
        max_cells: Maximum number of cells to process (None for all)
        overwrite: Whether to overwrite existing files
        ranking_files: Tuple/iterable of ranking database files to filter genes against
    """
    expr_path = os.path.join(results_dir, "expression.tsv")
    if os.path.exists(expr_path) and not overwrite:
        logging.info(
            "Found existing expression matrix at %s; skipping export (use --overwrite to regenerate)",
            expr_path,
        )
        return expr_path

    logging.info("Reading AnnData from %s", anndata_path)
    adata = sc.read_h5ad(anndata_path)
    logging.info("Initial gene count: %d", adata.n_vars)

    # Check if the main matrix contains integer counts
    import numpy as np

    test_data = adata.X[: min(100, adata.n_obs), : min(100, adata.n_vars)]
    if hasattr(test_data, "toarray"):
        test_data = test_data.toarray()

    is_integer = np.allclose(test_data, np.round(test_data))

    if not is_integer:
        logging.warning(
            "Main matrix (adata.X) does not appear to contain integer counts. "
            "This may be normalized data, which is not suitable for GRN inference."
        )
        # Check if counts layer exists
        if "counts" in adata.layers:
            logging.info(
                "Found 'counts' layer in adata.layers - checking if it contains integers"
            )
            test_counts = adata.layers["counts"][
                : min(100, adata.n_obs), : min(100, adata.n_vars)
            ]
            if hasattr(test_counts, "toarray"):
                test_counts = test_counts.toarray()

            if np.allclose(test_counts, np.round(test_counts)):
                logging.info(
                    "Using integer counts from adata.layers['counts'] instead of adata.X"
                )
                adata.X = adata.layers["counts"]
            else:
                logging.warning(
                    "adata.layers['counts'] also does not contain integers. "
                    "Proceeding with adata.X but results may be suboptimal."
                )
        else:
            logging.warning(
                "No 'counts' layer found in AnnData. "
                "Proceeding with adata.X but results may be suboptimal if data is normalized."
            )

    # Load ranking databases and filter to genes present in all of them
    logging.info("Loading ranking databases to identify common genes")
    db_genes_list = []
    for ranking_file in ranking_files:
        df = pd.read_feather(ranking_file)
        # Get all column names except the last one (which is typically a ranking/metadata column)
        db_genes = set(df.columns[:-1])
        db_genes_list.append(db_genes)
        logging.info(
            "Database %s has %d genes", os.path.basename(ranking_file), len(db_genes)
        )

    # Find intersection of all database genes
    if db_genes_list:
        common_genes = set.intersection(*db_genes_list)
        logging.info("Common genes across all databases: %d", len(common_genes))

        # Subset AnnData to genes present in all databases
        adata = adata[:, adata.var_names.isin(common_genes)]
        logging.info(
            "Filtered AnnData to genes in all databases: %d genes", adata.n_vars
        )

    # Slice to max_cells before removing non-expressed genes
    cell_slice = slice(None) if max_cells is None else slice(0, max_cells)
    adata = adata[cell_slice, :]
    logging.info("Sliced to %d cells", adata.n_obs)

    # Remove non-expressed genes (genes with zero counts across the sliced cells)
    expr = adata.X
    if hasattr(expr, "toarray"):
        expr_array = expr.toarray()
    else:
        expr_array = expr

    zero_genes = expr_array.sum(axis=0) == 0
    initial_n_vars = adata.n_vars
    adata = adata[:, ~zero_genes]
    logging.info(
        "Removed non-expressed genes: %d -> %d genes", initial_n_vars, adata.n_vars
    )

    expr = adata.X
    expr_dense = expr.toarray() if hasattr(expr, "toarray") else expr

    obs_names = adata.obs_names
    df_cells_x_genes = pd.DataFrame(
        expr_dense, index=obs_names, columns=adata.var_names
    )

    # Write to temporary file first, then atomically move to final location
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", dir=results_dir, delete=False
    ) as tmp_file:
        tmp_path = tmp_file.name

    try:
        df_cells_x_genes.to_csv(tmp_path, sep="\t", index=True)
        _atomic_move(tmp_path, expr_path)
        logging.info("Wrote expression matrix to %s", expr_path)
    except Exception as e:
        # Clean up temporary file on error
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise e

    return expr_path


def run_grnboost2(
    results_dir,
    resource_dir,
    tf_file,
    expression_path,
    n_cores,
    seed,
    overwrite,
):
    """Run pySCENIC GRNBoost2; skip if output exists unless overwrite."""
    tf_resolved = tf_file or os.path.join(resource_dir, "allTFs_hg38.txt")
    adjacency_path = os.path.join(results_dir, "adjacency.tsv")
    if os.path.exists(adjacency_path) and not overwrite:
        logging.info(
            "Found existing adjacency at %s; skipping grn (use --overwrite to rerun)",
            adjacency_path,
        )
        return adjacency_path

    # Use temporary file for output to avoid empty files being treated as valid
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", dir=results_dir, delete=False
    ) as tmp_file:
        tmp_adjacency_path = tmp_file.name

    try:
        cmd = [
            "pyscenic",
            "grn",
            "--num_workers",
            str(n_cores),
            "--seed",
            str(seed),
            "--method",
            "grnboost2",
            "--output",
            tmp_adjacency_path,
            expression_path,
            tf_resolved,
        ]

        logging.info("[GRN] Command: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error("pySCENIC grn failed. STDOUT: %s", result.stdout)
            logging.error("pySCENIC grn failed. STDERR: %s", result.stderr)
            raise subprocess.CalledProcessError(
                result.returncode, cmd, output=result.stdout, stderr=result.stderr
            )

        # Atomically move temporary file to final location
        _atomic_move(tmp_adjacency_path, adjacency_path)
        logging.info("[GRN] Completed; adjacency at %s", adjacency_path)
        return adjacency_path
    except Exception as e:
        # Clean up temporary file on error
        if os.path.exists(tmp_adjacency_path):
            os.remove(tmp_adjacency_path)
        raise e


def run_regdiff(
    results_dir,
    expression_path,
    n_cores,
    top_gene_percentile,
    overwrite,
    tf_file=None,
    resource_dir=None,
):
    """Run RegDiffusion GRN inference; skip if output exists unless overwrite.

    Args:
        results_dir: Directory to save output adjacency file.
        expression_path: Path to expression matrix TSV (cells x genes).
        n_cores: Number of cores/workers for edge extraction.
        top_gene_percentile: Percentile threshold for GRN edge weights (e.g., 50 for top 50%).
        overwrite: If True, regenerate adjacency even if it exists.
        tf_file: Path to TF list file (optional, defaults to resource_dir/allTFs_hg38.txt).
        resource_dir: Directory with pySCENIC resource files.

    Returns:
        Path to the generated adjacency.tsv file.
    """
    import numpy as np

    try:
        import regdiffusion as rd
    except ImportError:
        raise ImportError(
            "RegDiffusion is not installed. Install it with: pip install regdiffusion"
        )

    adjacency_path = os.path.join(results_dir, "adjacency.tsv")
    if os.path.exists(adjacency_path) and not overwrite:
        logging.info(
            "Found existing adjacency at %s; skipping regdiff (use --overwrite to rerun)",
            adjacency_path,
        )
        return adjacency_path

    tf_resolved = tf_file or os.path.join(resource_dir, "allTFs_hg38.txt")
    tf_df = pd.read_csv(tf_resolved, header=None, names=["TF"])

    logging.info("[RegDiff] Loading expression matrix from %s", expression_path)
    expr_df = pd.read_csv(expression_path, sep="\t", index_col=0)

    # Convert to numpy array and apply log transformation
    x = expr_df.values
    x = np.log(x + 1.0)

    logging.info(
        "[RegDiff] Training RegDiffusion model on %d cells x %d genes",
        x.shape[0],
        x.shape[1],
    )
    rd_trainer = rd.RegDiffusionTrainer(x)
    rd_trainer.train()

    logging.info(
        "[RegDiff] Extracting GRN with top %d percentile edges", top_gene_percentile
    )
    grn = rd_trainer.get_grn(
        gene_names=expr_df.columns.tolist(), top_gene_percentile=top_gene_percentile
    )

    # Extract all edges for each gene
    logging.info("[RegDiff] Extracting edgelist with %d workers", n_cores)
    adjacencies = grn.extract_edgelist(k=-1, workers=n_cores)
    adjacencies.columns = ["TF", "target", "importance"]

    # Restrict to TFs and sort by importance
    adjacencies = adjacencies[adjacencies["TF"].isin(tf_df["TF"])].sort_values(
        by="importance", ascending=False
    )

    # Write to temporary file first, then atomically move to final location
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", dir=results_dir, delete=False
    ) as tmp_file:
        tmp_path = tmp_file.name

    try:
        adjacencies.to_csv(tmp_path, sep="\t", index=False)
        _atomic_move(tmp_path, adjacency_path)
        logging.info("[RegDiff] Completed; adjacency at %s", adjacency_path)
        return adjacency_path
    except Exception as e:
        # Clean up temporary file on error
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise e


def run_ctx(
    results_dir,
    resource_dir,
    adjacency_path,
    expression_path,
    motif_f,
    ranking_files,
    n_cores,
    overwrite,
):
    """Run pySCENIC ctx (motif enrichment); skip if output exists unless overwrite.

    Args:
        ranking_files: Iterable of ranking database files (at least one required).
    """
    regulons_path = os.path.join(results_dir, "regulons.csv")
    if os.path.exists(regulons_path) and not overwrite:
        logging.info(
            "Found existing regulons at %s; skipping ctx (use --overwrite to rerun)",
            regulons_path,
        )
        return regulons_path

    # Use temporary file for output to avoid empty files being treated as valid
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".csv", dir=results_dir, delete=False
    ) as tmp_file:
        tmp_regulons_path = tmp_file.name

    try:
        cmd = [
            "pyscenic",
            "ctx",
            adjacency_path,
        ]
        cmd.extend(ranking_files)
        cmd.extend(
            [
                "--annotations_fname",
                motif_f,
                "--expression_mtx_fname",
                expression_path,
                "--mode",
                "custom_multiprocessing",
                "--output",
                tmp_regulons_path,
                "--num_workers",
                str(n_cores),
            ]
        )

        logging.info("[CTX] Command: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error("pySCENIC ctx failed. STDOUT: %s", result.stdout)
            logging.error("pySCENIC ctx failed. STDERR: %s", result.stderr)
            raise subprocess.CalledProcessError(
                result.returncode, cmd, output=result.stdout, stderr=result.stderr
            )

        # Atomically move temporary file to final location
        _atomic_move(tmp_regulons_path, regulons_path)
        logging.info("[CTX] Completed; regulons at %s", regulons_path)
        return regulons_path
    except Exception as e:
        # Clean up temporary file on error
        if os.path.exists(tmp_regulons_path):
            os.remove(tmp_regulons_path)
        raise e


def run_aucell(
    results_dir,
    expression_path,
    regulons_path,
    n_cores,
    overwrite,
):
    """Run pySCENIC aucell (AUCell activity scoring); skip if output exists unless overwrite.

    Args:
        results_dir: Directory to save output AUC matrix file.
        expression_path: Path to expression matrix TSV (cells x genes).
        regulons_path: Path to regulons CSV (from ctx step).
        n_cores: Number of workers for AUCell computation.
        overwrite: If True, regenerate AUC matrix even if it exists.

    Returns:
        Path to the generated AUC matrix CSV file.
    """
    auc_mtx_path = os.path.join(results_dir, "auc_mtx.csv")
    if os.path.exists(auc_mtx_path) and not overwrite:
        logging.info(
            "Found existing AUC matrix at %s; skipping aucell (use --overwrite to rerun)",
            auc_mtx_path,
        )
        return auc_mtx_path

    # Use temporary file for output to avoid empty files being treated as valid
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".csv", dir=results_dir, delete=False
    ) as tmp_file:
        tmp_auc_mtx_path = tmp_file.name

    try:
        cmd = [
            "pyscenic",
            "aucell",
            expression_path,
            regulons_path,
            "-o",
            tmp_auc_mtx_path,
            "--num_workers",
            str(n_cores),
        ]

        logging.info("[AUCell] Command: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error("pySCENIC aucell failed. STDOUT: %s", result.stdout)
            logging.error("pySCENIC aucell failed. STDERR: %s", result.stderr)
            raise subprocess.CalledProcessError(
                result.returncode, cmd, output=result.stdout, stderr=result.stderr
            )

        # Atomically move temporary file to final location
        _atomic_move(tmp_auc_mtx_path, auc_mtx_path)
        logging.info("[AUCell] Completed; AUC matrix at %s", auc_mtx_path)
        return auc_mtx_path
    except Exception as e:
        # Clean up temporary file on error
        if os.path.exists(tmp_auc_mtx_path):
            os.remove(tmp_auc_mtx_path)
        raise e


@click.command()
@click.option(
    "--anndata-path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to input .h5ad file",
)
@click.option(
    "--results-dir",
    default="./pyscenic_results",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Directory to write outputs",
)
@click.option(
    "--resource-dir",
    default="/pyscenic/resources",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Directory with pySCENIC resource files (TF lists, ranking databases, motif annotations)",
)
@click.option(
    "--tf-file",
    default=None,
    type=click.Path(dir_okay=False),
    help="Path to TF list file (defaults to resource_dir/allTFs_hg38.txt)",
)
@click.option(
    "--n-cores",
    default=16,
    show_default=True,
    type=int,
    help="Number of worker cores for each step",
)
@click.option(
    "--seed",
    default=42,
    show_default=True,
    type=int,
    help="Random seed passed to pySCENIC",
)
@click.option(
    "--max-cells",
    default=None,
    type=int,
    help="Export only the first N cells (optional; defaults to all cells)",
)
@click.option(
    "--log-level",
    default="INFO",
    show_default=True,
    type=click.Choice(
        ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], case_sensitive=False
    ),
    help="Logging level",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing outputs instead of skipping",
)
@click.option(
    "--motif-f",
    default=None,
    type=click.Path(dir_okay=False, exists=True),
    help="Motif annotations file; defaults to resource_dir/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
)
@click.option(
    "--ranking-files",
    "ranking_files",
    multiple=True,
    type=click.Path(dir_okay=False, exists=True),
    help="Ranking database files (at least one required; can pass multiple)",
)
@click.option(
    "--skip-ctx",
    is_flag=True,
    help="Skip ctx (motif enrichment) step",
)
@click.option(
    "--skip-grn",
    is_flag=True,
    help="Skip GRN (GRNBoost2) step; requires existing adjacency.tsv",
)
@click.option(
    "--grn-method",
    default="grnboost2",
    show_default=True,
    type=click.Choice(["grnboost2", "regdiff"], case_sensitive=False),
    help="GRN inference method: grnboost2 (via Singularity) or regdiff (RegDiffusion)",
)
@click.option(
    "--regdiff-percentile",
    default=50,
    show_default=True,
    type=int,
    help="RegDiffusion: percentile threshold for GRN edge weights (only used with --grn-method=regdiff)",
)
@click.option(
    "--skip-aucell",
    is_flag=True,
    help="Skip aucell (AUCell activity scoring) step",
)
def main(
    anndata_path,
    results_dir,
    resource_dir,
    tf_file,
    n_cores,
    seed,
    max_cells,
    log_level,
    overwrite,
    motif_f,
    ranking_files,
    skip_ctx,
    skip_grn,
    grn_method,
    regdiff_percentile,
    skip_aucell,
):
    """Run pySCENIC workflow (GRNBoost2 or RegDiffusion + ctx + AUCell) from an AnnData input."""

    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s %(levelname)s %(message)s",
    )

    # Resolve real (non-symlink) absolute paths
    anndata_path = _real(anndata_path)
    results_dir = _real(results_dir)
    resource_dir = _real(resource_dir)
    tf_file = _real(tf_file) if tf_file else tf_file
    motif_f = _real(motif_f) if motif_f else motif_f
    ranking_files = (
        tuple(_real(rf) for rf in ranking_files) if ranking_files else ranking_files
    )

    os.makedirs(results_dir, exist_ok=True)
    logging.info("Results directory: %s", results_dir)
    logging.info("[SETUP] anndata=%s", anndata_path)
    logging.info("[SETUP] results_dir=%s", results_dir)
    logging.info("[SETUP] resource_dir=%s", resource_dir)
    logging.info(
        "[SETUP] cores=%s seed=%s max_cells=%s overwrite=%s",
        n_cores,
        seed,
        max_cells,
        overwrite,
    )

    # Validate required files before running
    base_checks = [
        anndata_path,
        tf_file or os.path.join(resource_dir, "allTFs_hg38.txt"),
    ]

    if not skip_ctx:
        ranking_defaults = (
            _real(
                os.path.join(
                    resource_dir,
                    "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                )
            ),
            _real(
                os.path.join(
                    resource_dir,
                    "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                )
            ),
        )
        motif_default = _real(
            os.path.join(resource_dir, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
        )
        motif_resolved_check = motif_f or motif_default
        if ranking_files:
            ranking_files_resolved_check = ranking_files
        else:
            ranking_files_resolved_check = ranking_defaults
        base_checks.extend([motif_resolved_check, *ranking_files_resolved_check])

    _ensure_files_exist(base_checks)

    if not skip_grn:
        logging.info("[STEP] Export expression matrix")
        # Determine ranking files to use
        if ranking_files:
            ranking_files_for_filter = ranking_files
        else:
            ranking_files_for_filter = (
                _real(
                    os.path.join(
                        resource_dir,
                        "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                    )
                ),
                _real(
                    os.path.join(
                        resource_dir,
                        "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                    )
                ),
            )
        expression_path = write_expression_tsv(
            anndata_path, results_dir, max_cells, overwrite, ranking_files_for_filter
        )

        if grn_method.lower() == "grnboost2":
            logging.info("[STEP] Running GRNBoost2")
            adjacency_path = run_grnboost2(
                results_dir,
                resource_dir,
                tf_file,
                expression_path,
                n_cores,
                seed,
                overwrite,
            )
            logging.info("[STEP] GRNBoost2 complete")
        elif grn_method.lower() == "regdiff":
            logging.info("[STEP] Running RegDiffusion")
            adjacency_path = run_regdiff(
                results_dir,
                expression_path,
                n_cores,
                regdiff_percentile,
                overwrite,
                tf_file,
                resource_dir,
            )
            logging.info("[STEP] RegDiffusion complete")
        else:
            raise ValueError(f"Unknown GRN method: {grn_method}")
    else:
        logging.info("[STEP] Skipping GRN step (--skip-grn enabled)")
        adjacency_path = os.path.join(results_dir, "adjacency.tsv")
        expression_path = os.path.join(results_dir, "expression.tsv")
        if not os.path.exists(adjacency_path):
            raise FileNotFoundError(
                f"--skip-grn requires existing adjacency file at {adjacency_path}"
            )
        if not os.path.exists(expression_path):
            raise FileNotFoundError(
                f"--skip-grn requires existing expression matrix at {expression_path}"
            )
        logging.info("Found existing adjacency: %s", adjacency_path)
        logging.info("Found existing expression: %s", expression_path)

    if not skip_ctx:
        logging.info("[STEP] ctx (motif enrichment)")
        motif_resolved = motif_f or _real(
            os.path.join(resource_dir, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
        )
        if ranking_files:
            ranking_files_resolved = ranking_files
        else:
            ranking_files_resolved = (
                _real(
                    os.path.join(
                        resource_dir,
                        "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                    )
                ),
                _real(
                    os.path.join(
                        resource_dir,
                        "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                    )
                ),
            )
        # Validate that n_cores >= number of ranking files for ctx step
        num_ranking_files = len(ranking_files_resolved)
        if n_cores < num_ranking_files:
            raise ValueError(
                f"For ctx step, --n-cores ({n_cores}) must be >= number of ranking files ({num_ranking_files})"
            )
        regulons_path = run_ctx(
            results_dir,
            resource_dir,
            adjacency_path,
            expression_path,
            motif_resolved,
            ranking_files_resolved,
            n_cores,
            overwrite,
        )
    else:
        regulons_path = os.path.join(results_dir, "regulons.csv")

    if not skip_aucell:
        logging.info("[STEP] aucell (AUCell activity scoring)")
        if not os.path.exists(regulons_path):
            raise FileNotFoundError(
                f"aucell requires regulons file at {regulons_path}. "
                "Either run without --skip-ctx or provide existing regulons."
            )
        run_aucell(
            results_dir,
            expression_path,
            regulons_path,
            n_cores,
            overwrite,
        )
        logging.info("[STEP] aucell complete")
    else:
        logging.info("[STEP] Skipping aucell step (--skip-aucell enabled)")


if __name__ == "__main__":
    main()
