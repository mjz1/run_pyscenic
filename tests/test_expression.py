"""Tests for expression matrix export functionality."""

import os

import numpy as np
import pandas as pd
import pytest

from run_pyscenic import write_expression_tsv


def test_write_expression_tsv_basic(synthetic_adata, adata_h5ad_file, mock_ranking_file, tmp_results_dir):
    """Test basic expression matrix export."""
    expr_path = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        max_cells=None,
        overwrite=False,
        ranking_files=(mock_ranking_file,),
    )
    
    assert os.path.exists(expr_path)
    assert expr_path.endswith('expression.tsv')
    
    # Verify the exported file is a valid TSV
    expr_df = pd.read_csv(expr_path, sep='\t', index_col=0)
    assert expr_df.shape[0] > 0  # Has cells
    assert expr_df.shape[1] > 0  # Has genes
    assert (expr_df >= 0).all().all()  # All non-negative


def test_write_expression_tsv_max_cells(synthetic_adata, adata_h5ad_file, mock_ranking_file, tmp_results_dir):
    """Test expression matrix export with max_cells limit."""
    max_cells = 50
    expr_path = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        max_cells=max_cells,
        overwrite=False,
        ranking_files=(mock_ranking_file,),
    )
    
    expr_df = pd.read_csv(expr_path, sep='\t', index_col=0)
    assert expr_df.shape[0] == max_cells  # Should have exactly max_cells rows


def test_write_expression_tsv_skip_when_exists(synthetic_adata, adata_h5ad_file, mock_ranking_file, tmp_results_dir):
    """Test that existing expression matrix is skipped unless overwrite=True."""
    expr_path = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        max_cells=None,
        overwrite=False,
        ranking_files=(mock_ranking_file,),
    )
    
    # Get the modification time
    original_mtime = os.path.getmtime(expr_path)
    
    # Run again without overwrite
    expr_path_2 = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        max_cells=None,
        overwrite=False,
        ranking_files=(mock_ranking_file,),
    )
    
    # File should not have been modified
    assert os.path.getmtime(expr_path_2) == original_mtime


def test_write_expression_tsv_overwrite(synthetic_adata, adata_h5ad_file, mock_ranking_file, tmp_results_dir):
    """Test that expression matrix is regenerated with overwrite=True."""
    import time
    
    expr_path = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        max_cells=None,
        overwrite=False,
        ranking_files=(mock_ranking_file,),
    )
    
    original_mtime = os.path.getmtime(expr_path)
    
    # Wait a bit to ensure mtime would change
    time.sleep(0.1)
    
    # Run again with overwrite
    expr_path_2 = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        max_cells=None,
        overwrite=True,
        ranking_files=(mock_ranking_file,),
    )
    
    # File should have been modified
    assert os.path.getmtime(expr_path_2) > original_mtime
