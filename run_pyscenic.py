import logging
import os
import subprocess

import click
import pandas as pd
import scanpy as sc


def _ensure_files_exist(paths):
    missing = [p for p in paths if p and not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(f"Missing required files: {missing}")


def _real(path):
    return os.path.realpath(path) if path else path


def write_expression_tsv(anndata_path, results_dir, max_cells, overwrite, ranking_files, resource_dir):
    """Load AnnData, filter to genes in ranking databases, export cells x genes TSV; return path. Skip if present unless overwrite.
    
    Args:
        ranking_files: Tuple of ranking database files to filter genes against.
        resource_dir: Path to resource directory (for default ranking files if none provided).
    """
    expr_t_path = os.path.join(results_dir, "expression_t.tsv")
    if os.path.exists(expr_t_path) and not overwrite:
        logging.info(
            "Found existing expression matrix at %s; skipping export (use --overwrite to regenerate)",
            expr_t_path,
        )
        return expr_t_path

    logging.info("Reading AnnData from %s", anndata_path)
    adata = sc.read_h5ad(anndata_path)
    logging.info("Initial gene count: %d", adata.n_vars)

    # Load ranking databases and filter to genes present in all of them
    logging.info("Loading ranking databases to identify common genes")
    db_genes_list = []
    for ranking_file in ranking_files:
        df = pd.read_feather(ranking_file)
        # Get all column names except the last one (which is typically a ranking/metadata column)
        db_genes = set(df.columns[:-1])
        db_genes_list.append(db_genes)
        logging.info("Database %s has %d genes", os.path.basename(ranking_file), len(db_genes))
    
    # Find intersection of all database genes
    if db_genes_list:
        common_genes = set.intersection(*db_genes_list)
        logging.info("Common genes across all databases: %d", len(common_genes))
        
        # Subset AnnData to genes present in all databases
        adata = adata[:, adata.var_names.isin(common_genes)]
        logging.info("Filtered AnnData to genes in all databases: %d genes", adata.n_vars)

    cell_slice = slice(None) if max_cells is None else slice(0, max_cells)
    expr = adata.X[cell_slice, :]
    expr_dense = expr.toarray() if hasattr(expr, "toarray") else expr

    obs_names = adata.obs_names[cell_slice]
    df_cells_x_genes = pd.DataFrame(
        expr_dense, index=obs_names, columns=adata.var_names
    )
    df_cells_x_genes.to_csv(expr_t_path, sep="\t", index=True)

    logging.info("Wrote expression matrix to %s", expr_t_path)
    return expr_t_path


def run_grnboost2(
    sif_img,
    results_dir,
    resource_dir,
    tf_file,
    expression_t_path,
    n_cores,
    seed,
    overwrite,
):
    """Run pySCENIC GRNBoost2 inside Singularity; skip if output exists unless overwrite."""
    tf_resolved = tf_file or os.path.join(resource_dir, "allTFs_hg38.txt")
    adjacency_path = os.path.join(results_dir, "adjacency.tsv")
    if os.path.exists(adjacency_path) and not overwrite:
        logging.info(
            "Found existing adjacency at %s; skipping grn (use --overwrite to rerun)",
            adjacency_path,
        )
        return adjacency_path

    cmd = [
        "singularity",
        "run",
        "-B",
        f"{results_dir}:{results_dir}",
        "-B",
        f"{resource_dir}:{resource_dir}",
        sif_img,
        "pyscenic",
        "grn",
        "--num_workers",
        str(n_cores),
        "--seed",
        str(seed),
        "--method",
        "grnboost2",
        "--output",
        adjacency_path,
        expression_t_path,
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

    logging.info("[GRN] Completed; adjacency at %s", adjacency_path)
    return adjacency_path

def run_regdiff():
    """Placeholder for future regdiff implementation."""
    pass

def run_ctx(
    sif_img,
    results_dir,
    resource_dir,
    adjacency_path,
    expression_t_path,
    motif_f,
    ranking_files,
    n_cores,
    overwrite,
):
    """Run pySCENIC ctx (motif enrichment) inside Singularity; skip if output exists unless overwrite.

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

    cmd = [
        "singularity",
        "run",
        "-B",
        f"{results_dir}:{results_dir}",
        "-B",
        f"{resource_dir}:{resource_dir}",
        sif_img,
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
            expression_t_path,
            "--mode",
            "custom_multiprocessing",
            "--output",
            regulons_path,
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

    logging.info("[CTX] Completed; regulons at %s", regulons_path)
    return regulons_path


@click.command()
@click.option(
    "--anndata-path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to input .h5ad file",
)
@click.option(
    "--sif-img",
    default="/data1/shahs3/users/zatzmanm/work/images/pyscenic/aertslab-pyscenic-0.12.1.sif",
    show_default=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to pySCENIC Singularity image",
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
    default="/data1/shahs3/users/zatzmanm/work/images/pyscenic/pyscenic_resources",
    show_default=True,
    type=click.Path(file_okay=False, exists=True),
    help="Directory with pySCENIC resource files",
)
@click.option(
    "--tf-file",
    required=True,
    type=click.Path(dir_okay=False, exists=True),
    help="Path to TF list file",
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
def main(
    anndata_path,
    sif_img,
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
):
    """Run pySCENIC GRNBoost2 from an AnnData input."""

    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s %(levelname)s %(message)s",
    )

    # Resolve real (non-symlink) absolute paths for Singularity binds
    anndata_path = _real(anndata_path)
    sif_img = _real(sif_img)
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
    logging.info("[SETUP] sif_img=%s", sif_img)
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
        sif_img,
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
        expression_t_path = write_expression_tsv(
            anndata_path, results_dir, max_cells, overwrite, ranking_files_for_filter, resource_dir
        )
        logging.info("[STEP] Running GRNBoost2")
        adjacency_path = run_grnboost2(
            sif_img,
            results_dir,
            resource_dir,
            tf_file,
            expression_t_path,
            n_cores,
            seed,
            overwrite,
        )
        logging.info("[STEP] GRNBoost2 complete")
    else:
        logging.info("[STEP] Skipping GRN step (--skip-grn enabled)")
        adjacency_path = os.path.join(results_dir, "adjacency.tsv")
        expression_t_path = os.path.join(results_dir, "expression_t.tsv")
        if not os.path.exists(adjacency_path):
            raise FileNotFoundError(
                f"--skip-grn requires existing adjacency file at {adjacency_path}"
            )
        if not os.path.exists(expression_t_path):
            raise FileNotFoundError(
                f"--skip-grn requires existing expression matrix at {expression_t_path}"
            )
        logging.info("Found existing adjacency: %s", adjacency_path)
        logging.info("Found existing expression: %s", expression_t_path)

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
                        "hg38-tss-centered-5kb-10species.mc9nr.genes_vs_motifs.rankings.feather",
                    )
                ),
                _real(
                    os.path.join(
                        resource_dir,
                        "hg38-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather",
                    )
                ),
            )
        # Validate that n_cores >= number of ranking files for ctx step
        num_ranking_files = len(ranking_files_resolved)
        if n_cores < num_ranking_files:
            raise ValueError(
                f"For ctx step, --n-cores ({n_cores}) must be >= number of ranking files ({num_ranking_files})"
            )
        run_ctx(
            sif_img,
            results_dir,
            resource_dir,
            adjacency_path,
            expression_t_path,
            motif_resolved,
            ranking_files_resolved,
            n_cores,
            overwrite,
        )


if __name__ == "__main__":
    main()
