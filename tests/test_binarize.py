"""Tests for AUCell matrix binarization functionality."""

import os

import numpy as np
import pandas as pd
import pytest

from run_pyscenic import run_binarize


def test_run_binarize_basic(mock_auc_matrix, tmp_results_dir):
    """Test basic binarization of AUCell matrix."""
    bin_mtx_path, thresholds_path = run_binarize(
        results_dir=tmp_results_dir,
        auc_mtx_path=mock_auc_matrix,
        n_cores=2,
        overwrite=False,
    )
    
    assert os.path.exists(bin_mtx_path)
    assert os.path.exists(thresholds_path)
    
    # Check binarized matrix
    bin_mtx = pd.read_csv(bin_mtx_path, index_col=0)
    assert set(bin_mtx.values.flatten()).issubset({0, 1})  # Only 0s and 1s
    
    # Check thresholds
    thresholds = pd.read_csv(thresholds_path, index_col=0)
    assert len(thresholds) > 0  # Has thresholds


def test_run_binarize_output_names(mock_auc_matrix, tmp_results_dir):
    """Test that binarization produces correctly named output files."""
    bin_mtx_path, thresholds_path = run_binarize(
        results_dir=tmp_results_dir,
        auc_mtx_path=mock_auc_matrix,
        n_cores=2,
        overwrite=False,
    )
    
    assert bin_mtx_path.endswith('auc_mtx_binarized.csv')
    assert thresholds_path.endswith('auc_binarization_thresholds.csv')


def test_run_binarize_skip_when_exists(mock_auc_matrix, tmp_results_dir):
    """Test that existing binarized matrix is skipped unless overwrite=True."""
    bin_mtx_path, thresholds_path = run_binarize(
        results_dir=tmp_results_dir,
        auc_mtx_path=mock_auc_matrix,
        n_cores=2,
        overwrite=False,
    )
    
    original_bin_mtime = os.path.getmtime(bin_mtx_path)
    
    # Run again without overwrite
    bin_mtx_path_2, thresholds_path_2 = run_binarize(
        results_dir=tmp_results_dir,
        auc_mtx_path=mock_auc_matrix,
        n_cores=2,
        overwrite=False,
    )
    
    # Files should not have been modified
    assert os.path.getmtime(bin_mtx_path_2) == original_bin_mtime


def test_run_binarize_shapes_match(mock_auc_matrix, tmp_results_dir):
    """Test that binarized matrix has same shape as input."""
    auc_mtx = pd.read_csv(mock_auc_matrix, index_col=0)
    
    bin_mtx_path, thresholds_path = run_binarize(
        results_dir=tmp_results_dir,
        auc_mtx_path=mock_auc_matrix,
        n_cores=2,
        overwrite=False,
    )
    
    bin_mtx = pd.read_csv(bin_mtx_path, index_col=0)
    assert auc_mtx.shape == bin_mtx.shape
    assert list(auc_mtx.index) == list(bin_mtx.index)
    assert list(auc_mtx.columns) == list(bin_mtx.columns)
