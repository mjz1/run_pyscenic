"""Integration tests for the pySCENIC pipeline."""

import os

import pandas as pd
import pytest


def test_pipeline_expression_export(synthetic_adata, adata_h5ad_file, mock_ranking_file, tmp_results_dir):
    """Test the first step of the pipeline: expression matrix export."""
    from run_pyscenic import write_expression_tsv
    
    expr_path = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        max_cells=None,
        overwrite=False,
        ranking_files=(mock_ranking_file,),
    )
    
    assert os.path.exists(expr_path)
    expr_df = pd.read_csv(expr_path, sep='\t', index_col=0)
    assert expr_df.shape[0] > 0 and expr_df.shape[1] > 0


def test_pipeline_with_binarization(mock_auc_matrix, tmp_results_dir):
    """Test binarization as a final pipeline step."""
    from run_pyscenic import run_binarize
    
    bin_mtx_path, thresholds_path = run_binarize(
        results_dir=tmp_results_dir,
        auc_mtx_path=mock_auc_matrix,
        n_cores=2,
        overwrite=False,
    )
    
    # Verify outputs
    bin_mtx = pd.read_csv(bin_mtx_path, index_col=0)
    thresholds = pd.read_csv(thresholds_path, index_col=0)
    
    # Binarized matrix should be 0/1
    assert set(bin_mtx.values.flatten()).issubset({0, 1})
    # Thresholds should have one per regulon
    assert len(thresholds) == bin_mtx.shape[1]


def test_atomic_writes(synthetic_adata, adata_h5ad_file, mock_ranking_file, tmp_results_dir):
    """Test that atomic writes prevent partial files."""
    from run_pyscenic import write_expression_tsv
    import time
    
    # Write expression matrix
    expr_path = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        max_cells=None,
        overwrite=False,
        ranking_files=(mock_ranking_file,),
    )
    
    # Verify the file is complete (can be read)
    expr_df = pd.read_csv(expr_path, sep='\t', index_col=0)
    assert expr_df.shape[0] > 0
    
    # Overwrite should produce a valid file too
    time.sleep(0.05)
    expr_path_2 = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        max_cells=None,
        overwrite=True,
        ranking_files=(mock_ranking_file,),
    )
    
    expr_df_2 = pd.read_csv(expr_path_2, sep='\t', index_col=0)
    assert expr_df_2.shape[0] > 0
