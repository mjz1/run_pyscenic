#!/usr/bin/env python3
import ast
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile

import click
import numpy as np
import pandas as pd
import scanpy as sc

from binarize import binarize
from subsample import subsample

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


def _resolve_geosketch_embedding(
    adata,
    results_dir,
    overwrite,
    embedding_key=None,
    batch_correction_method="none",
    batch_key=None,
    n_pca_components=50,
    seed=42,
):
    """Resolve the embedding matrix used by geosketch subsampling.

    Precedence:
      1. ``adata.obsm[embedding_key]`` if *embedding_key* is provided.
      2. In-pipeline batch-corrected embedding (harmony / scanorama).
      3. Naive PCA fallback.

    The computed embedding is cached to ``{results_dir}/subsample_embedding.npy``
    (respects *overwrite*).  adata.X is never modified.

    Returns:
        (embedding, embedding_meta) where *embedding* is an (n_cells, n_dims)
        numpy array and *embedding_meta* is a dict of metadata fields.
    """
    cache_path = os.path.join(results_dir, "subsample_embedding.npy")
    meta = {}

    # ----- Case 1: pre-existing key in obsm -----
    if embedding_key is not None:
        if embedding_key not in adata.obsm:
            raise ValueError(
                f"--subsample-embedding-key '{embedding_key}' not found in "
                f"adata.obsm. Available keys: {list(adata.obsm.keys())}"
            )
        emb = np.asarray(adata.obsm[embedding_key])
        meta["embedding_source"] = "obsm"
        meta["embedding_key"] = embedding_key
        meta["embedding_dimensions"] = emb.shape[1]
        if emb.shape[1] < 5:
            logging.warning(
                "Embedding '%s' has only %d dimensions (< 5); geosketch quality may suffer",
                embedding_key,
                emb.shape[1],
            )
        if emb.shape[1] > 100:
            logging.warning(
                "Embedding '%s' has %d dimensions (> 100); consider reducing dimensionality",
                embedding_key,
                emb.shape[1],
            )
        # Cache even user-supplied embeddings for reproducibility
        _save_embedding(emb, cache_path, results_dir, overwrite)
        return emb, meta

    # ----- Case 2 / 3: compute embedding -----
    if os.path.exists(cache_path) and not overwrite:
        logging.info(
            "Loading cached embedding from %s (use --overwrite to recompute)",
            cache_path,
        )
        emb = np.load(cache_path)
        meta["embedding_source"] = "cached"
        meta["embedding_dimensions"] = emb.shape[1]
        return emb, meta

    # Work on a copy so adata.X (raw counts) is never touched
    import scanpy as sc

    adata_tmp = adata.copy()
    sc.pp.normalize_total(adata_tmp)
    sc.pp.log1p(adata_tmp)
    sc.pp.pca(adata_tmp, n_comps=n_pca_components, random_state=seed)
    meta["n_pca_components"] = n_pca_components

    bc = batch_correction_method.lower()
    if bc == "none":
        emb = adata_tmp.obsm["X_pca"]
        meta["embedding_source"] = "pca"
        meta["batch_correction_method"] = "none"
    elif bc == "harmony":
        if batch_key is None or batch_key not in adata.obs.columns:
            raise ValueError(
                f"--subsample-batch-key '{batch_key}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        try:
            import harmonypy
        except ImportError:
            raise ImportError(
                "harmonypy is not installed. Install it with: pip install harmonypy"
            )
        logging.info(
            "[Embedding] Running Harmony batch correction on batch_key='%s'",
            batch_key,
        )
        ho = harmonypy.run_harmony(
            adata_tmp.obsm["X_pca"],
            adata_tmp.obs,
            batch_key,
            random_state=seed,
        )
        emb = ho.Z_corr.T  # harmony returns (n_dims, n_cells)
        meta["embedding_source"] = "harmony"
        meta["batch_correction_method"] = "harmony"
        meta["batch_key"] = batch_key
    elif bc == "scanorama":
        if batch_key is None or batch_key not in adata.obs.columns:
            raise ValueError(
                f"--subsample-batch-key '{batch_key}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        try:
            import scanorama
        except ImportError:
            raise ImportError(
                "scanorama is not installed. Install it with: pip install scanorama"
            )
        logging.info(
            "[Embedding] Running Scanorama batch correction on batch_key='%s'",
            batch_key,
        )
        # scanorama.integrate_scanpy expects a list of AnnData per batch
        batches = adata_tmp.obs[batch_key].unique()
        adata_list = [adata_tmp[adata_tmp.obs[batch_key] == b].copy() for b in batches]
        scanorama.integrate_scanpy(adata_list)
        # Reassemble in original cell order
        emb = np.zeros((adata_tmp.n_obs, adata_list[0].obsm["X_scanorama"].shape[1]))
        for b, ad_b in zip(batches, adata_list):
            mask = adata_tmp.obs[batch_key] == b
            emb[mask.values] = ad_b.obsm["X_scanorama"]
        meta["embedding_source"] = "scanorama"
        meta["batch_correction_method"] = "scanorama"
        meta["batch_key"] = batch_key
    else:
        raise ValueError(f"Unknown batch correction method: {batch_correction_method}")

    meta["embedding_dimensions"] = emb.shape[1]
    if emb.shape[1] < 5:
        logging.warning(
            "Computed embedding has only %d dimensions (< 5); geosketch quality may suffer",
            emb.shape[1],
        )
    if emb.shape[1] > 100:
        logging.warning(
            "Computed embedding has %d dimensions (> 100); consider reducing dimensionality",
            emb.shape[1],
        )

    _save_embedding(emb, cache_path, results_dir, overwrite)
    return emb, meta


def _save_embedding(emb, cache_path, results_dir, overwrite):
    """Save embedding array to cache_path using atomic write."""
    if os.path.exists(cache_path) and not overwrite:
        return
    with tempfile.NamedTemporaryFile(
        suffix=".npy", dir=results_dir, delete=False
    ) as tmp:
        tmp_path = tmp.name
    try:
        np.save(tmp_path, emb)
        _atomic_move(tmp_path, cache_path)
        logging.info("Cached embedding to %s", cache_path)
    except Exception as e:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise e


def _load_and_filter_adata(anndata_path, ranking_files):
    """Load AnnData, validate/swap to integer counts, filter genes to ranking DB intersection, remove zero genes.

    Args:
        anndata_path: Path to input .h5ad file.
        ranking_files: Tuple/iterable of ranking database feather files.

    Returns:
        Filtered AnnData object (all cells, filtered genes).
    """
    logging.info("Reading AnnData from %s", anndata_path)
    adata = sc.read_h5ad(anndata_path)
    logging.info("Initial gene count: %d", adata.n_vars)

    # Check if the main matrix contains integer counts
    test_data = adata.X[: min(100, adata.n_obs), : min(100, adata.n_vars)]
    if hasattr(test_data, "toarray"):
        test_data = test_data.toarray()

    is_integer = np.allclose(test_data, np.round(test_data))

    if not is_integer:
        logging.warning(
            "Main matrix (adata.X) does not appear to contain integer counts. "
            "This may be normalized data, which is not suitable for GRN inference."
        )
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
        db_genes = set(df.columns[:-1])
        db_genes_list.append(db_genes)
        logging.info(
            "Database %s has %d genes", os.path.basename(ranking_file), len(db_genes)
        )

    if db_genes_list:
        common_genes = set.intersection(*db_genes_list)
        logging.info("Common genes across all databases: %d", len(common_genes))
        adata = adata[:, adata.var_names.isin(common_genes)]
        logging.info(
            "Filtered AnnData to genes in all databases: %d genes", adata.n_vars
        )

    # Remove non-expressed genes
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

    return adata


def _write_expression_tsv(adata, output_path, overwrite):
    """Write an AnnData slice to a cells x genes TSV with atomic write.

    Args:
        adata: AnnData object (or slice) to export.
        output_path: Destination TSV file path.
        overwrite: Whether to overwrite an existing file.

    Returns:
        output_path.
    """
    if os.path.exists(output_path) and not overwrite:
        logging.info(
            "Found existing expression matrix at %s; skipping export (use --overwrite to regenerate)",
            output_path,
        )
        return output_path

    expr = adata.X
    expr_dense = expr.toarray() if hasattr(expr, "toarray") else expr

    df_cells_x_genes = pd.DataFrame(
        expr_dense, index=adata.obs_names, columns=adata.var_names
    )

    results_dir = os.path.dirname(output_path)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", dir=results_dir, delete=False
    ) as tmp_file:
        tmp_path = tmp_file.name

    try:
        df_cells_x_genes.to_csv(tmp_path, sep="\t", index=True)
        _atomic_move(tmp_path, output_path)
        logging.info("Wrote expression matrix to %s", output_path)
    except Exception as e:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise e

    return output_path


def write_expression_tsv(anndata_path, results_dir, overwrite, ranking_files):
    """Load AnnData, filter to genes in ranking databases, export cells x genes TSV; return path. Skip if present unless overwrite.

    This is kept for backward compatibility.  Internally delegates to
    ``_load_and_filter_adata`` + ``_write_expression_tsv``.

    Args:
        anndata_path: Path to input .h5ad file
        results_dir: Directory to write output files
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

    adata = _load_and_filter_adata(anndata_path, ranking_files)
    return _write_expression_tsv(adata, expr_path, overwrite)


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


def run_binarize(
    results_dir,
    auc_mtx_path,
    n_cores,
    overwrite,
):
    """Binarize AUCell matrix using adaptive thresholds; skip if output exists unless overwrite.

    Args:
        results_dir: Directory to save output files.
        auc_mtx_path: Path to AUC matrix CSV (from aucell step).
        n_cores: Number of workers for binarization.
        overwrite: If True, regenerate binarized matrix even if it exists.

    Returns:
        Tuple of (binarized_matrix_path, thresholds_path).
    """
    bin_mtx_path = os.path.join(results_dir, "auc_mtx_binarized.csv")
    thresholds_path = os.path.join(results_dir, "auc_binarization_thresholds.csv")

    if (
        os.path.exists(bin_mtx_path)
        and os.path.exists(thresholds_path)
        and not overwrite
    ):
        logging.info(
            "Found existing binarized matrix at %s; skipping binarization (use --overwrite to rerun)",
            bin_mtx_path,
        )
        return bin_mtx_path, thresholds_path

    logging.info("[Binarize] Loading AUC matrix from %s", auc_mtx_path)
    auc_mtx = pd.read_csv(auc_mtx_path, index_col=0)

    logging.info(
        "[Binarize] Computing adaptive thresholds for %d regulons across %d cells",
        auc_mtx.shape[1],
        auc_mtx.shape[0],
    )
    bin_mtx, thresholds = binarize(auc_mtx, num_workers=n_cores)

    # Use temporary files for output to avoid partial writes
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".csv", dir=results_dir, delete=False
    ) as tmp_bin:
        tmp_bin_path = tmp_bin.name

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".csv", dir=results_dir, delete=False
    ) as tmp_thresh:
        tmp_thresh_path = tmp_thresh.name

    try:
        bin_mtx.to_csv(tmp_bin_path)
        thresholds.to_csv(tmp_thresh_path)

        _atomic_move(tmp_bin_path, bin_mtx_path)
        _atomic_move(tmp_thresh_path, thresholds_path)

        logging.info("[Binarize] Completed; binarized matrix at %s", bin_mtx_path)
        logging.info("[Binarize] Thresholds saved to %s", thresholds_path)
        return bin_mtx_path, thresholds_path
    except Exception as e:
        # Clean up temporary files on error
        if os.path.exists(tmp_bin_path):
            os.remove(tmp_bin_path)
        if os.path.exists(tmp_thresh_path):
            os.remove(tmp_thresh_path)
        raise e


def flatten_regulons(regulons_path):
    """Parse a pySCENIC regulons CSV into a flat DataFrame with one row per TF-motif-target.

    All regulon-level metadata is preserved alongside each target gene, producing
    a complete unfiltered table suitable for downstream filtering or gene-set analyses.

    Args:
        regulons_path: Path to the regulons CSV produced by ``pyscenic ctx``.

    Returns:
        pandas.DataFrame with columns: TF, MotifID, AUC, NES, MotifSimilarityQvalue,
        OrthologousIdentity, Annotation, Context, RankAtMax, n_targets, target, weight.
    """
    df = pd.read_csv(regulons_path, header=[0, 1], index_col=[0, 1])

    rows = []
    for (tf, motif_id), row in df.iterrows():
        target_genes = ast.literal_eval(row[("Enrichment", "TargetGenes")])
        n_targets = len(target_genes)
        base = {
            "TF": tf,
            "MotifID": motif_id,
            "AUC": row[("Enrichment", "AUC")],
            "NES": row[("Enrichment", "NES")],
            "MotifSimilarityQvalue": row[("Enrichment", "MotifSimilarityQvalue")],
            "OrthologousIdentity": row[("Enrichment", "OrthologousIdentity")],
            "Annotation": row[("Enrichment", "Annotation")],
            "Context": row[("Enrichment", "Context")],
            "RankAtMax": row[("Enrichment", "RankAtMax")],
            "n_targets": n_targets,
        }
        for gene, weight in target_genes:
            rows.append({**base, "target": gene, "weight": weight})

    return pd.DataFrame(rows)


def filter_regulons(
    flat_df,
    min_nes=3.0,
    require_activating=True,
    require_direct_annotation=True,
    min_targets=10,
):
    """Filter a flattened regulons table to high-confidence entries.

    Designed to operate on the output of :func:`flatten_regulons`.  All
    parameters have sensible defaults matching standard pySCENIC practice.

    Args:
        flat_df: DataFrame produced by ``flatten_regulons``.
        min_nes: Minimum NES score (default 3.0).
        require_activating: Keep only rows whose Context contains 'activating' (default True).
        require_direct_annotation: Keep only directly annotated or orthologous-directly-annotated
            TF-motif links (default True).
        min_targets: Minimum number of targets per TF-motif pair (default 10).

    Returns:
        Filtered copy of the DataFrame.
    """
    df = flat_df.copy()

    df = df[df["NES"] >= min_nes]

    if require_activating:
        df = df[df["Context"].str.contains("activating", na=False)]

    if require_direct_annotation:
        is_direct = df["Annotation"] == "gene is directly annotated"
        is_ortho_direct = (
            df["Annotation"].str.startswith("gene is orthologous to", na=False)
            & df["Annotation"].str.endswith(
                "which is directly annotated for motif", na=False
            )
        ) | df["Annotation"].str.contains(
            r"which is directly annotated$", regex=True, na=False
        )
        df = df[is_direct | is_ortho_direct]

    if min_targets > 0:
        df = df[df["n_targets"] >= min_targets]

    return df.reset_index(drop=True)


def deduplicate_regulons(df):
    """Collapse a regulons table to one row per unique TF-target pair.

    When a TF-target pair appears under multiple motifs, the row with the
    highest NES is retained.  This produces a clean gene-set table suitable
    for downstream scoring (e.g. AUCell, ssGSEA, UCell).

    Args:
        df: DataFrame produced by :func:`flatten_regulons` or
            :func:`filter_regulons`.

    Returns:
        Deduplicated DataFrame with the same columns, sorted by TF then target.
    """
    return (
        df.sort_values("NES", ascending=False)
        .drop_duplicates(subset=["TF", "target"])
        .sort_values(["TF", "target"])
        .reset_index(drop=True)
    )


def run_flatten_regulons(results_dir, regulons_path, overwrite):
    """Flatten regulons CSV and write unfiltered, filtered, and deduplicated tables.

    Writes four files to *results_dir*:
      - ``regulons_flat.tsv``  — complete unfiltered table (one row per TF-motif-target)
      - ``regulons_flat_filtered.tsv`` — filtered with default quality thresholds
      - ``regulons_flat_dedup.tsv`` — deduplicated unfiltered (one row per TF-target, best NES)
      - ``regulons_flat_filtered_dedup.tsv`` — deduplicated filtered

    Args:
        results_dir: Directory to write output files.
        regulons_path: Path to the regulons CSV from ctx step.
        overwrite: If True, regenerate even if outputs exist.

    Returns:
        Tuple of (flat_path, filtered_path, dedup_path, filtered_dedup_path).
    """
    flat_path = os.path.join(results_dir, "regulons_flat.tsv")
    filtered_path = os.path.join(results_dir, "regulons_flat_filtered.tsv")
    dedup_path = os.path.join(results_dir, "regulons_flat_dedup.tsv")
    filtered_dedup_path = os.path.join(results_dir, "regulons_flat_filtered_dedup.tsv")

    all_paths = (flat_path, filtered_path, dedup_path, filtered_dedup_path)

    if all(os.path.exists(p) for p in all_paths) and not overwrite:
        logging.info(
            "Found existing flattened regulons at %s; skipping (use --overwrite to rerun)",
            flat_path,
        )
        return all_paths

    logging.info("[Flatten] Parsing regulons from %s", regulons_path)
    flat_df = flatten_regulons(regulons_path)
    logging.info(
        "[Flatten] Unfiltered: %d rows (%d unique TFs)",
        len(flat_df),
        flat_df["TF"].nunique(),
    )

    filtered_df = filter_regulons(flat_df)
    logging.info(
        "[Flatten] Filtered: %d rows (%d unique TFs)",
        len(filtered_df),
        filtered_df["TF"].nunique(),
    )

    dedup_df = deduplicate_regulons(flat_df)
    logging.info(
        "[Flatten] Deduplicated (unfiltered): %d unique TF-target pairs (%d unique TFs)",
        len(dedup_df),
        dedup_df["TF"].nunique(),
    )

    filtered_dedup_df = deduplicate_regulons(filtered_df)
    logging.info(
        "[Flatten] Deduplicated (filtered): %d unique TF-target pairs (%d unique TFs)",
        len(filtered_dedup_df),
        filtered_dedup_df["TF"].nunique(),
    )

    # Atomic writes
    tmp_paths = []
    try:
        for df, dest in (
            (flat_df, flat_path),
            (filtered_df, filtered_path),
            (dedup_df, dedup_path),
            (filtered_dedup_df, filtered_dedup_path),
        ):
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tsv", dir=results_dir, delete=False
            ) as tmp:
                tmp_paths.append(tmp.name)
            df.to_csv(tmp_paths[-1], sep="\t", index=False)

        for tmp_p, dest_p in zip(tmp_paths, all_paths):
            _atomic_move(tmp_p, dest_p)

        logging.info("[Flatten] Unfiltered table: %s", flat_path)
        logging.info("[Flatten] Filtered table: %s", filtered_path)
        logging.info("[Flatten] Deduplicated (unfiltered) table: %s", dedup_path)
        logging.info("[Flatten] Deduplicated (filtered) table: %s", filtered_dedup_path)
        return all_paths
    except Exception as e:
        for p in tmp_paths:
            if os.path.exists(p):
                os.remove(p)
        raise e


def run_aucell(
    results_dir,
    expression_path,
    regulons_path,
    n_cores,
    rank_threshold,
    auc_threshold,
    nes_threshold,
    overwrite,
):
    """Run pySCENIC aucell (AUCell activity scoring); skip if output exists unless overwrite.

    Args:
        results_dir: Directory to save output AUC matrix file.
        expression_path: Path to expression matrix TSV (cells x genes).
        regulons_path: Path to regulons CSV (from ctx step).
        n_cores: Number of workers for AUCell computation.
        rank_threshold: Rank threshold for motif enrichment (passed to aucell).
        auc_threshold: AUC threshold for motif enrichment (passed to aucell).
        nes_threshold: NES threshold for motif enrichment (passed to aucell).
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
            "--rank_threshold",
            str(rank_threshold),
            "--auc_threshold",
            str(auc_threshold),
            "--nes_threshold",
            str(nes_threshold),
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
    "--subsample-method",
    default="none",
    show_default=True,
    type=click.Choice(
        ["none", "random", "stratified", "geosketch"], case_sensitive=False
    ),
    help="Subsampling strategy for GRN/ctx steps (none = use all cells)",
)
@click.option(
    "--subsample-n",
    default=None,
    type=int,
    help="Target number of cells for subsampling (required if method != none)",
)
@click.option(
    "--subsample-annotation",
    default=None,
    type=str,
    help="adata.obs column for stratified sampling (required if method=stratified)",
)
@click.option(
    "--subsample-min-per-group",
    default=50,
    show_default=True,
    type=int,
    help="Minimum cells per group for stratified sampling",
)
@click.option(
    "--subsample-embedding-key",
    default=None,
    type=str,
    help="adata.obsm key with pre-computed embedding for geosketch (overrides batch correction)",
)
@click.option(
    "--subsample-batch-correction-method",
    default="none",
    show_default=True,
    type=click.Choice(["none", "harmony", "scanorama"], case_sensitive=False),
    help="Batch correction for geosketch embedding (ignored if --subsample-embedding-key is set)",
)
@click.option(
    "--subsample-batch-key",
    default=None,
    type=str,
    help="adata.obs column with batch labels (required if batch-correction-method != none)",
)
@click.option(
    "--subsample-n-pca-components",
    default=50,
    show_default=True,
    type=int,
    help="PCA components before batch correction for geosketch embedding",
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
@click.option(
    "--rank-threshold",
    default=5000,
    show_default=True,
    type=int,
    help="AUCell: rank threshold for motif enrichment",
)
@click.option(
    "--auc-threshold",
    default=0.05,
    show_default=True,
    type=float,
    help="AUCell: AUC threshold (fraction of ranked genes)",
)
@click.option(
    "--nes-threshold",
    default=3.0,
    show_default=True,
    type=float,
    help="AUCell: NES threshold for enriched features",
)
@click.option(
    "--skip-binarize",
    is_flag=True,
    help="Skip binarization of AUCell matrix",
)
@click.option(
    "--skip-flatten",
    is_flag=True,
    help="Skip flattening regulons to per-target TSV tables",
)
def main(
    anndata_path,
    results_dir,
    resource_dir,
    tf_file,
    n_cores,
    seed,
    subsample_method,
    subsample_n,
    subsample_annotation,
    subsample_min_per_group,
    subsample_embedding_key,
    subsample_batch_correction_method,
    subsample_batch_key,
    subsample_n_pca_components,
    log_level,
    overwrite,
    motif_f,
    ranking_files,
    skip_ctx,
    skip_grn,
    grn_method,
    regdiff_percentile,
    skip_aucell,
    rank_threshold,
    auc_threshold,
    nes_threshold,
    skip_binarize,
    skip_flatten,
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
        "[SETUP] cores=%s seed=%s overwrite=%s",
        n_cores,
        seed,
        overwrite,
    )

    # Validate subsampling options
    subsampling_enabled = subsample_method.lower() != "none"
    if subsampling_enabled:
        if subsample_n is None:
            raise click.BadParameter(
                "--subsample-n is required when --subsample-method is not 'none'"
            )
        if subsample_method.lower() == "stratified" and subsample_annotation is None:
            raise click.BadParameter(
                "--subsample-annotation is required when --subsample-method is 'stratified'"
            )
        if subsample_n < 10000:
            logging.warning(
                "subsample_n=%d is below 10000; GRN inference quality may be impacted",
                subsample_n,
            )
        # Validate geosketch-specific batch correction options
        bc_method = subsample_batch_correction_method.lower()
        if subsample_method.lower() == "geosketch":
            if subsample_embedding_key and bc_method != "none":
                logging.warning(
                    "--subsample-embedding-key is set; ignoring "
                    "--subsample-batch-correction-method=%s",
                    subsample_batch_correction_method,
                )
            if bc_method != "none" and not subsample_embedding_key:
                if subsample_batch_key is None:
                    raise click.BadParameter(
                        "--subsample-batch-key is required when "
                        "--subsample-batch-correction-method is not 'none'"
                    )
        else:
            # Not geosketch — warn if batch correction options are set
            if bc_method != "none":
                logging.warning(
                    "--subsample-batch-correction-method is ignored when "
                    "--subsample-method is not 'geosketch'"
                )
        logging.info(
            "[SETUP] subsampling: method=%s n=%s annotation=%s seed=%s",
            subsample_method,
            subsample_n,
            subsample_annotation,
            seed,
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

        if subsampling_enabled:
            # Load and filter AnnData once, then write two matrices
            adata = _load_and_filter_adata(anndata_path, ranking_files_for_filter)

            # Check if subsampling is meaningful
            if subsample_n >= adata.n_obs:
                logging.warning(
                    "subsample_n (%d) >= total cells after filtering (%d); "
                    "skipping subsampling and using all cells",
                    subsample_n,
                    adata.n_obs,
                )
                subsampling_enabled = False
                expression_path = _write_expression_tsv(
                    adata,
                    os.path.join(results_dir, "expression.tsv"),
                    overwrite,
                )
                expression_full_path = expression_path
            else:
                # Resolve embedding for geosketch (if applicable)
                embedding = None
                embedding_meta = {}
                if subsample_method.lower() == "geosketch":
                    embedding, embedding_meta = _resolve_geosketch_embedding(
                        adata,
                        results_dir,
                        overwrite,
                        embedding_key=subsample_embedding_key,
                        batch_correction_method=subsample_batch_correction_method,
                        batch_key=subsample_batch_key,
                        n_pca_components=subsample_n_pca_components,
                        seed=seed,
                    )

                # Run subsampling
                logging.info(
                    "[Subsample] Running %s subsampling: %d -> %d cells",
                    subsample_method,
                    adata.n_obs,
                    subsample_n,
                )
                indices = subsample(
                    adata,
                    method=subsample_method.lower(),
                    n=subsample_n,
                    seed=seed,
                    annotation_column=subsample_annotation,
                    min_per_group=subsample_min_per_group,
                    embedding=embedding,
                )
                actual_n = len(indices)
                logging.info(
                    "[Subsample] Selected %d cells (target: %d)",
                    actual_n,
                    subsample_n,
                )

                # Build metadata
                meta = {
                    "method": subsample_method.lower(),
                    "target_n": subsample_n,
                    "actual_n": actual_n,
                    "seed": seed,
                    "total_cells_original": adata.n_obs,
                    "total_genes": adata.n_vars,
                }
                # Include embedding metadata for geosketch
                if embedding_meta:
                    meta.update(embedding_meta)
                if subsample_method.lower() == "stratified" and subsample_annotation:
                    meta["annotation_column"] = subsample_annotation
                    labels = adata.obs[subsample_annotation].values
                    unique_labels, orig_counts = np.unique(labels, return_counts=True)
                    meta["per_group_counts"] = {
                        str(k): int(v) for k, v in zip(unique_labels, orig_counts)
                    }
                    sub_labels = labels[indices]
                    sub_unique, sub_counts = np.unique(sub_labels, return_counts=True)
                    meta["per_group_subsampled"] = {
                        str(k): int(v) for k, v in zip(sub_unique, sub_counts)
                    }
                    # Log per-group counts
                    for grp in unique_labels:
                        orig_c = int(orig_counts[np.where(unique_labels == grp)[0][0]])
                        sub_c = meta["per_group_subsampled"].get(str(grp), 0)
                        logging.info(
                            "[Subsample]   %s: %d -> %d cells", grp, orig_c, sub_c
                        )
                    # Warn about small groups
                    for grp, cnt in meta["per_group_subsampled"].items():
                        if cnt < subsample_min_per_group:
                            logging.warning(
                                "[Subsample] Group '%s' has %d cells after sampling "
                                "(< min_per_group=%d)",
                                grp,
                                cnt,
                                subsample_min_per_group,
                            )

                # Save metadata JSON
                meta_path = os.path.join(results_dir, "subsample_metadata.json")
                with tempfile.NamedTemporaryFile(
                    mode="w", suffix=".json", dir=results_dir, delete=False
                ) as tmp_meta:
                    tmp_meta_path = tmp_meta.name
                try:
                    json.dump(meta, open(tmp_meta_path, "w"), indent=2)
                    _atomic_move(tmp_meta_path, meta_path)
                    logging.info("[Subsample] Metadata saved to %s", meta_path)
                except Exception as e:
                    if os.path.exists(tmp_meta_path):
                        os.remove(tmp_meta_path)
                    raise e

                # Save selected barcodes
                barcodes_path = os.path.join(results_dir, "subsampled_barcodes.txt")
                selected_barcodes = adata.obs_names[indices]
                with tempfile.NamedTemporaryFile(
                    mode="w", suffix=".txt", dir=results_dir, delete=False
                ) as tmp_bc:
                    tmp_bc_path = tmp_bc.name
                try:
                    with open(tmp_bc_path, "w") as f:
                        f.write("\n".join(selected_barcodes) + "\n")
                    _atomic_move(tmp_bc_path, barcodes_path)
                    logging.info("[Subsample] Barcodes saved to %s", barcodes_path)
                except Exception as e:
                    if os.path.exists(tmp_bc_path):
                        os.remove(tmp_bc_path)
                    raise e

                # Write subsampled expression matrix (for GRN + ctx)
                expression_path = _write_expression_tsv(
                    adata[indices, :],
                    os.path.join(results_dir, "expression_subsampled.tsv"),
                    overwrite,
                )
                # Write full expression matrix (for AUCell)
                expression_full_path = _write_expression_tsv(
                    adata,
                    os.path.join(results_dir, "expression_full.tsv"),
                    overwrite,
                )
                # Assert identical gene columns
                sub_cols = pd.read_csv(
                    expression_path, sep="\t", index_col=0, nrows=0
                ).columns.tolist()
                full_cols = pd.read_csv(
                    expression_full_path, sep="\t", index_col=0, nrows=0
                ).columns.tolist()
                assert (
                    sub_cols == full_cols
                ), "Gene columns differ between subsampled and full matrices"
        else:
            # No subsampling — single expression.tsv
            expression_path = write_expression_tsv(
                anndata_path, results_dir, overwrite, ranking_files_for_filter
            )
            expression_full_path = expression_path

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
        # When skipping GRN with subsampling, expect the full matrix for aucell
        if subsampling_enabled:
            expression_path = os.path.join(results_dir, "expression_subsampled.tsv")
            expression_full_path = os.path.join(results_dir, "expression_full.tsv")
        else:
            expression_path = os.path.join(results_dir, "expression.tsv")
            expression_full_path = expression_path
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
        # ctx uses subsampled expression (or single expression.tsv if no subsampling)
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

    if not skip_flatten:
        logging.info("[STEP] flatten regulons to per-target tables")
        if not os.path.exists(regulons_path):
            raise FileNotFoundError(
                f"flatten requires regulons file at {regulons_path}. "
                "Either run without --skip-ctx or provide existing regulons.csv."
            )
        flat_path, filtered_path, dedup_path, filtered_dedup_path = run_flatten_regulons(
            results_dir,
            regulons_path,
            overwrite,
        )
        logging.info("[STEP] flatten complete")
    else:
        logging.info("[STEP] Skipping flatten step (--skip-flatten enabled)")

    # AUCell uses the FULL expression matrix when subsampling is enabled
    aucell_expression = expression_full_path if subsampling_enabled else expression_path

    if not skip_aucell:
        logging.info("[STEP] aucell (AUCell activity scoring)")
        if not os.path.exists(regulons_path):
            raise FileNotFoundError(
                f"aucell requires regulons file at {regulons_path}. "
                "Either run without --skip-ctx or provide existing regulons."
            )
        auc_mtx_path = run_aucell(
            results_dir,
            aucell_expression,
            regulons_path,
            n_cores,
            rank_threshold,
            auc_threshold,
            nes_threshold,
            overwrite,
        )
        logging.info("[STEP] aucell complete")
    else:
        logging.info("[STEP] Skipping aucell step (--skip-aucell enabled)")
        auc_mtx_path = os.path.join(results_dir, "auc_mtx.csv")

    if not skip_binarize:
        logging.info("[STEP] binarize (AUCell matrix binarization)")
        if not os.path.exists(auc_mtx_path):
            raise FileNotFoundError(
                f"binarization requires AUC matrix file at {auc_mtx_path}. "
                "Either run without --skip-aucell or provide existing auc_mtx.csv."
            )
        bin_mtx_path, thresholds_path = run_binarize(
            results_dir,
            auc_mtx_path,
            n_cores,
            overwrite,
        )
        logging.info("[STEP] binarization complete")
    else:
        logging.info("[STEP] Skipping binarization step (--skip-binarize enabled)")


if __name__ == "__main__":
    main()
