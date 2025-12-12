"""Pytest fixtures and configuration for run_pyscenic tests."""

import os
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def tmp_results_dir():
    """Create a temporary results directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def synthetic_adata():
    """Create a small synthetic AnnData object for testing.
    
    Returns an AnnData object with:
    - 100 cells x 200 genes
    - Integer count matrix
    - Basic cell and gene annotations
    """
    np.random.seed(42)
    n_obs, n_vars = 100, 200
    
    # Create synthetic count data
    X = np.random.poisson(lam=5, size=(n_obs, n_vars))
    
    # Create obs (cell metadata)
    obs = pd.DataFrame({
        'cell_id': [f'cell_{i}' for i in range(n_obs)],
        'n_counts': X.sum(axis=1),
    }, index=[f'cell_{i}' for i in range(n_obs)])
    
    # Create var (gene metadata)
    var = pd.DataFrame({
        'gene_name': [f'gene_{i}' for i in range(n_vars)],
        'n_cells': (X > 0).sum(axis=0),
    }, index=[f'gene_{i}' for i in range(n_vars)])
    
    adata = ad.AnnData(X=X, obs=obs, var=var)
    return adata


@pytest.fixture
def adata_h5ad_file(synthetic_adata, tmp_results_dir):
    """Save synthetic AnnData to an h5ad file and return the path."""
    h5ad_path = os.path.join(tmp_results_dir, 'test_data.h5ad')
    synthetic_adata.write_h5ad(h5ad_path)
    yield h5ad_path


@pytest.fixture
def mock_ranking_file(tmp_results_dir):
    """Create a mock ranking database (feather format) for testing."""
    # Create a mock ranking file with ~150 genes (subset of test genes)
    genes = [f'gene_{i}' for i in range(150)]
    rankings_data = {gene: np.random.random(100) for gene in genes}
    rankings_df = pd.DataFrame(rankings_data)
    
    ranking_path = os.path.join(tmp_results_dir, 'mock_rankings.feather')
    rankings_df.to_feather(ranking_path)
    yield ranking_path


@pytest.fixture
def mock_tf_file(tmp_results_dir):
    """Create a mock TF list file for testing."""
    tfs = ['gene_10', 'gene_20', 'gene_30', 'gene_40', 'gene_50']
    tf_df = pd.DataFrame({'TF': tfs})
    
    tf_path = os.path.join(tmp_results_dir, 'tfs.txt')
    tf_df.to_csv(tf_path, header=False, index=False)
    yield tf_path


@pytest.fixture
def mock_motif_file(tmp_results_dir):
    """Create a mock motif annotations file for testing."""
    motif_data = {
        'motif_id': ['MOTIF_1', 'MOTIF_2', 'MOTIF_3'],
        'gene_id': ['gene_10', 'gene_20', 'gene_30'],
        'orthologous_gene_id': ['ENSG000001', 'ENSG000002', 'ENSG000003'],
    }
    motif_df = pd.DataFrame(motif_data)
    
    motif_path = os.path.join(tmp_results_dir, 'motifs.tbl')
    motif_df.to_csv(motif_path, sep='\t', index=False)
    yield motif_path


@pytest.fixture
def mock_adjacency_file(tmp_results_dir):
    """Create a mock adjacency (GRN) file for testing."""
    adjacency_data = {
        'TF': ['gene_10', 'gene_20', 'gene_30'] * 5,
        'target': [f'gene_{50 + i}' for i in range(15)],
        'importance': np.random.random(15),
    }
    adj_df = pd.DataFrame(adjacency_data)
    
    adj_path = os.path.join(tmp_results_dir, 'adjacency.tsv')
    adj_df.to_csv(adj_path, sep='\t', index=False)
    yield adj_path


@pytest.fixture
def mock_regulons_file(tmp_results_dir):
    """Create a mock regulons file (from ctx step) for testing."""
    # Simplified regulons CSV with minimal columns
    regulons_data = {
        'TF': ['gene_10', 'gene_20', 'gene_30'],
        'gene_id': ['gene_100', 'gene_120', 'gene_140'],
        'motif_id': ['MOTIF_1', 'MOTIF_2', 'MOTIF_3'],
    }
    regulons_df = pd.DataFrame(regulons_data)
    
    regulons_path = os.path.join(tmp_results_dir, 'regulons.csv')
    regulons_df.to_csv(regulons_path, index=False)
    yield regulons_path


@pytest.fixture
def mock_auc_matrix(tmp_results_dir):
    """Create a mock AUC matrix (from aucell step) for testing."""
    n_cells, n_regulons = 50, 3
    auc_data = np.random.random((n_cells, n_regulons))
    
    auc_df = pd.DataFrame(
        auc_data,
        index=[f'cell_{i}' for i in range(n_cells)],
        columns=['gene_10_regulon', 'gene_20_regulon', 'gene_30_regulon']
    )
    
    auc_path = os.path.join(tmp_results_dir, 'auc_mtx.csv')
    auc_df.to_csv(auc_path)
    yield auc_path
