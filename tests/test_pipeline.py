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


def test_aucell_thresholds_passed(tmp_results_dir, monkeypatch):
    """Ensure aucell thresholds are passed through to pyscenic command."""
    import run_pyscenic

    expression_path = os.path.join(tmp_results_dir, "expression.tsv")
    regulons_path = os.path.join(tmp_results_dir, "regulons.csv")
    with open(expression_path, "w", encoding="utf-8") as expr_file:
        expr_file.write("cell\tgene_1\ncell_1\t1\n")
    with open(regulons_path, "w", encoding="utf-8") as reg_file:
        reg_file.write("TF,gene_id,motif_id\nTF1,gene_1,M1\n")

    captured = {}

    def fake_run(cmd, capture_output=True, text=True):
        captured["cmd"] = cmd

        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

        return Result()

    monkeypatch.setattr(run_pyscenic.subprocess, "run", fake_run)

    rank_threshold = 1234
    auc_threshold = 0.12
    nes_threshold = 4.5

    auc_mtx_path = run_pyscenic.run_aucell(
        results_dir=tmp_results_dir,
        expression_path=expression_path,
        regulons_path=regulons_path,
        n_cores=2,
        rank_threshold=rank_threshold,
        auc_threshold=auc_threshold,
        nes_threshold=nes_threshold,
        overwrite=True,
    )

    cmd = captured["cmd"]
    assert "--rank_threshold" in cmd
    assert cmd[cmd.index("--rank_threshold") + 1] == str(rank_threshold)
    assert "--auc_threshold" in cmd
    assert cmd[cmd.index("--auc_threshold") + 1] == str(auc_threshold)
    assert "--nes_threshold" in cmd
    assert cmd[cmd.index("--nes_threshold") + 1] == str(nes_threshold)
    assert os.path.exists(auc_mtx_path)
