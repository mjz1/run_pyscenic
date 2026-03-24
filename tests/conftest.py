"""Pytest fixtures and configuration for run_pyscenic tests."""

import os
import tempfile

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
    obs = pd.DataFrame(
        {
            "cell_id": [f"cell_{i}" for i in range(n_obs)],
            "n_counts": X.sum(axis=1),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )

    # Create var (gene metadata)
    var = pd.DataFrame(
        {
            "gene_name": [f"gene_{i}" for i in range(n_vars)],
            "n_cells": (X > 0).sum(axis=0),
        },
        index=[f"gene_{i}" for i in range(n_vars)],
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)
    return adata


@pytest.fixture
def adata_h5ad_file(synthetic_adata, tmp_results_dir):
    """Save synthetic AnnData to an h5ad file and return the path."""
    h5ad_path = os.path.join(tmp_results_dir, "test_data.h5ad")
    synthetic_adata.write_h5ad(h5ad_path)
    yield h5ad_path


@pytest.fixture
def mock_ranking_file(tmp_results_dir):
    """Create a mock ranking database (feather format) for testing."""
    # Create a mock ranking file with ~150 genes (subset of test genes)
    genes = [f"gene_{i}" for i in range(150)]
    rankings_data = {gene: np.random.random(100) for gene in genes}
    rankings_df = pd.DataFrame(rankings_data)

    ranking_path = os.path.join(tmp_results_dir, "mock_rankings.feather")
    rankings_df.to_feather(ranking_path)
    yield ranking_path


@pytest.fixture
def mock_tf_file(tmp_results_dir):
    """Create a mock TF list file for testing."""
    tfs = ["gene_10", "gene_20", "gene_30", "gene_40", "gene_50"]
    tf_df = pd.DataFrame({"TF": tfs})

    tf_path = os.path.join(tmp_results_dir, "tfs.txt")
    tf_df.to_csv(tf_path, header=False, index=False)
    yield tf_path


@pytest.fixture
def mock_motif_file(tmp_results_dir):
    """Create a mock motif annotations file for testing."""
    motif_data = {
        "motif_id": ["MOTIF_1", "MOTIF_2", "MOTIF_3"],
        "gene_id": ["gene_10", "gene_20", "gene_30"],
        "orthologous_gene_id": ["ENSG000001", "ENSG000002", "ENSG000003"],
    }
    motif_df = pd.DataFrame(motif_data)

    motif_path = os.path.join(tmp_results_dir, "motifs.tbl")
    motif_df.to_csv(motif_path, sep="\t", index=False)
    yield motif_path


@pytest.fixture
def mock_adjacency_file(tmp_results_dir):
    """Create a mock adjacency (GRN) file for testing."""
    adjacency_data = {
        "TF": ["gene_10", "gene_20", "gene_30"] * 5,
        "target": [f"gene_{50 + i}" for i in range(15)],
        "importance": np.random.random(15),
    }
    adj_df = pd.DataFrame(adjacency_data)

    adj_path = os.path.join(tmp_results_dir, "adjacency.tsv")
    adj_df.to_csv(adj_path, sep="\t", index=False)
    yield adj_path


@pytest.fixture
def mock_regulons_file(tmp_results_dir):
    """Create a mock regulons file (from ctx step) for testing."""
    # Simplified regulons CSV with minimal columns
    regulons_data = {
        "TF": ["gene_10", "gene_20", "gene_30"],
        "gene_id": ["gene_100", "gene_120", "gene_140"],
        "motif_id": ["MOTIF_1", "MOTIF_2", "MOTIF_3"],
    }
    regulons_df = pd.DataFrame(regulons_data)

    regulons_path = os.path.join(tmp_results_dir, "regulons.csv")
    regulons_df.to_csv(regulons_path, index=False)
    yield regulons_path


@pytest.fixture
def mock_auc_matrix(tmp_results_dir):
    """Create a mock AUC matrix (from aucell step) for testing."""
    n_cells, n_regulons = 50, 3
    auc_data = np.random.random((n_cells, n_regulons))

    auc_df = pd.DataFrame(
        auc_data,
        index=[f"cell_{i}" for i in range(n_cells)],
        columns=["gene_10_regulon", "gene_20_regulon", "gene_30_regulon"],
    )

    auc_path = os.path.join(tmp_results_dir, "auc_mtx.csv")
    auc_df.to_csv(auc_path)
    yield auc_path


@pytest.fixture
def mock_pyscenic_regulons_csv(tmp_results_dir):
    """Create a mock regulons CSV matching the real pySCENIC ctx multi-header format.

    Includes rows with varying NES, Annotation, and Context values so that
    filter_regulons tests can exercise each filter independently.
    """
    rows = [
        # Row 0: passes all default filters
        {
            "TF": "TF_A",
            "MotifID": "motif_1",
            "AUC": 0.06,
            "NES": 4.5,
            "MotifSimilarityQvalue": 1e-6,
            "OrthologousIdentity": 1.0,
            "Annotation": "gene is directly annotated",
            "Context": "frozenset({'weight>75.0%', 'activating', 'hg38_db.rankings'})",
            "TargetGenes": "[('geneA', 3.5), ('geneB', 2.1), ('geneC', 1.8), "
            "('geneD', 4.0), ('geneE', 2.5), ('geneF', 3.0), "
            "('geneG', 1.2), ('geneH', 2.8), ('geneI', 3.3), "
            "('geneJ', 4.1), ('geneK', 1.5)]",
            "RankAtMax": 500,
        },
        # Row 1: NES below threshold
        {
            "TF": "TF_B",
            "MotifID": "motif_2",
            "AUC": 0.04,
            "NES": 2.5,
            "MotifSimilarityQvalue": 0.01,
            "OrthologousIdentity": 0.9,
            "Annotation": "gene is directly annotated",
            "Context": "frozenset({'weight>75.0%', 'activating', 'hg38_db.rankings'})",
            "TargetGenes": "[('geneX', 1.0), ('geneY', 2.0), ('geneZ', 3.0), "
            "('geneW', 1.5), ('geneV', 2.5), ('geneU', 3.5), "
            "('geneT', 1.1), ('geneS', 2.2), ('geneR', 3.3), "
            "('geneQ', 4.4), ('geneP', 5.5)]",
            "RankAtMax": 600,
        },
        # Row 2: not activating context
        {
            "TF": "TF_C",
            "MotifID": "motif_3",
            "AUC": 0.05,
            "NES": 3.5,
            "MotifSimilarityQvalue": 1e-4,
            "OrthologousIdentity": 1.0,
            "Annotation": "gene is directly annotated",
            "Context": "frozenset({'weight>75.0%', 'repressing', 'hg38_db.rankings'})",
            "TargetGenes": "[('geneM', 2.0), ('geneN', 3.0)]",
            "RankAtMax": 700,
        },
        # Row 3: annotation is inferred only (not direct)
        {
            "TF": "TF_D",
            "MotifID": "motif_4",
            "AUC": 0.05,
            "NES": 3.2,
            "MotifSimilarityQvalue": 0.001,
            "OrthologousIdentity": 0.8,
            "Annotation": "gene is inferred by gene family",
            "Context": "frozenset({'weight>75.0%', 'activating', 'hg38_db.rankings'})",
            "TargetGenes": "[('gene1', 1.0), ('gene2', 2.0), ('gene3', 3.0), "
            "('gene4', 1.5), ('gene5', 2.5), ('gene6', 3.5), "
            "('gene7', 1.1), ('gene8', 2.2), ('gene9', 3.3), "
            "('gene10', 4.4), ('gene11', 5.5)]",
            "RankAtMax": 800,
        },
        # Row 4: too few targets (< 10)
        {
            "TF": "TF_E",
            "MotifID": "motif_5",
            "AUC": 0.07,
            "NES": 5.0,
            "MotifSimilarityQvalue": 1e-8,
            "OrthologousIdentity": 1.0,
            "Annotation": "gene is directly annotated",
            "Context": "frozenset({'weight>75.0%', 'activating', 'hg38_db.rankings'})",
            "TargetGenes": "[('gAA', 2.0), ('gBB', 3.0), ('gCC', 4.0)]",
            "RankAtMax": 300,
        },
        # Row 5: orthologous + directly annotated (should pass)
        {
            "TF": "TF_F",
            "MotifID": "motif_6",
            "AUC": 0.055,
            "NES": 3.8,
            "MotifSimilarityQvalue": 5e-6,
            "OrthologousIdentity": 0.95,
            "Annotation": "gene is orthologous to ENSMUSG00000040310 in M. musculus "
            "(identity = 95%) which is directly annotated for motif",
            "Context": "frozenset({'weight>75.0%', 'activating', 'hg38_db.rankings'})",
            "TargetGenes": "[('oA', 1.0), ('oB', 2.0), ('oC', 3.0), "
            "('oD', 1.5), ('oE', 2.5), ('oF', 3.5), "
            "('oG', 1.1), ('oH', 2.2), ('oI', 3.3), "
            "('oJ', 4.4)]",
            "RankAtMax": 450,
        },
    ]

    # Build the multi-header CSV matching pySCENIC ctx output format
    enrichment_cols = [
        "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity",
        "Annotation", "Context", "TargetGenes", "RankAtMax",
    ]
    header_row0 = ",," + ",".join(["Enrichment"] * len(enrichment_cols))
    header_row1 = ",," + ",".join(enrichment_cols)
    header_row2 = "TF,MotifID," + "," * (len(enrichment_cols) - 1)

    lines = [header_row0, header_row1, header_row2]
    for r in rows:
        vals = [str(r[c]) for c in enrichment_cols]
        # Quote fields that contain commas
        quoted = []
        for v in vals:
            if "," in v or "'" in v:
                quoted.append('"' + v.replace('"', '""') + '"')
            else:
                quoted.append(v)
        lines.append(f"{r['TF']},{r['MotifID']}," + ",".join(quoted))

    csv_path = os.path.join(tmp_results_dir, "regulons.csv")
    with open(csv_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    yield csv_path
