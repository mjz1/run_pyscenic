"""Integration tests for the pySCENIC pipeline."""

import os

import pandas as pd


def test_pipeline_expression_export(
    synthetic_adata, adata_h5ad_file, mock_ranking_file, tmp_results_dir
):
    """Test the first step of the pipeline: expression matrix export."""
    from run_pyscenic import write_expression_tsv

    expr_path = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        overwrite=False,
        ranking_files=(mock_ranking_file,),
    )

    assert os.path.exists(expr_path)
    expr_df = pd.read_csv(expr_path, sep="\t", index_col=0)
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


def test_atomic_writes(
    synthetic_adata, adata_h5ad_file, mock_ranking_file, tmp_results_dir
):
    """Test that atomic writes prevent partial files."""
    from run_pyscenic import write_expression_tsv
    import time

    # Write expression matrix
    expr_path = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        overwrite=False,
        ranking_files=(mock_ranking_file,),
    )

    # Verify the file is complete (can be read)
    expr_df = pd.read_csv(expr_path, sep="\t", index_col=0)
    assert expr_df.shape[0] > 0

    # Overwrite should produce a valid file too
    time.sleep(0.05)
    expr_path_2 = write_expression_tsv(
        anndata_path=adata_h5ad_file,
        results_dir=tmp_results_dir,
        overwrite=True,
        ranking_files=(mock_ranking_file,),
    )

    expr_df_2 = pd.read_csv(expr_path_2, sep="\t", index_col=0)
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


def test_flatten_regulons_unfiltered(mock_pyscenic_regulons_csv):
    """Test that flatten_regulons produces all rows with correct columns."""
    from run_pyscenic import flatten_regulons

    df = flatten_regulons(mock_pyscenic_regulons_csv)

    expected_cols = {
        "TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue",
        "OrthologousIdentity", "Annotation", "Context", "RankAtMax",
        "n_targets", "target", "weight",
    }
    assert set(df.columns) == expected_cols

    # No NaN in key columns
    assert df["TF"].notna().all()
    assert df["target"].notna().all()
    assert df["weight"].notna().all()

    # Total rows should equal sum of all target lists:
    # Row0=11, Row1=11, Row2=2, Row3=11, Row4=3, Row5=10 => 48
    assert len(df) == 48

    # n_targets should match per-group count
    for (tf, motif), group in df.groupby(["TF", "MotifID"]):
        assert group["n_targets"].iloc[0] == len(group)

    # Check all 6 TFs present
    assert df["TF"].nunique() == 6


def test_filter_regulons_defaults(mock_pyscenic_regulons_csv):
    """Test that filter_regulons with defaults keeps only high-quality entries."""
    from run_pyscenic import flatten_regulons, filter_regulons

    flat_df = flatten_regulons(mock_pyscenic_regulons_csv)
    filtered = filter_regulons(flat_df)

    # Should keep TF_A (passes all) and TF_F (orthologous + direct)
    # TF_B: NES 2.5 < 3.0 => removed
    # TF_C: repressing context => removed
    # TF_D: inferred annotation => removed
    # TF_E: only 3 targets < 10 => removed
    assert set(filtered["TF"].unique()) == {"TF_A", "TF_F"}
    assert len(filtered) == 11 + 10  # 11 targets for TF_A, 10 for TF_F


def test_filter_regulons_custom_params(mock_pyscenic_regulons_csv):
    """Test filter_regulons with relaxed parameters."""
    from run_pyscenic import flatten_regulons, filter_regulons

    flat_df = flatten_regulons(mock_pyscenic_regulons_csv)

    # Relax all filters
    filtered = filter_regulons(
        flat_df,
        min_nes=0.0,
        require_activating=False,
        require_direct_annotation=False,
        min_targets=0,
    )
    # With no filters, should return everything
    assert len(filtered) == len(flat_df)

    # Only relax NES
    filtered_nes = filter_regulons(flat_df, min_nes=2.0)
    assert "TF_B" in filtered_nes["TF"].values


def test_run_flatten_regulons_writes_files(mock_pyscenic_regulons_csv, tmp_results_dir):
    """Test the full run_flatten_regulons pipeline step."""
    from run_pyscenic import run_flatten_regulons

    flat_path, filtered_path = run_flatten_regulons(
        results_dir=tmp_results_dir,
        regulons_path=mock_pyscenic_regulons_csv,
        overwrite=False,
    )

    assert os.path.exists(flat_path)
    assert os.path.exists(filtered_path)

    flat_df = pd.read_csv(flat_path, sep="\t")
    filtered_df = pd.read_csv(filtered_path, sep="\t")

    assert len(flat_df) == 48
    assert len(filtered_df) < len(flat_df)
    assert set(filtered_df["TF"].unique()) == {"TF_A", "TF_F"}


def test_run_flatten_regulons_skip_existing(mock_pyscenic_regulons_csv, tmp_results_dir):
    """Test that run_flatten_regulons skips when outputs exist and overwrite=False."""
    from run_pyscenic import run_flatten_regulons

    # First run
    run_flatten_regulons(tmp_results_dir, mock_pyscenic_regulons_csv, overwrite=False)

    # Modify the flat file to detect if it gets overwritten
    flat_path = os.path.join(tmp_results_dir, "regulons_flat.tsv")
    mtime_before = os.path.getmtime(flat_path)

    import time
    time.sleep(0.05)

    # Second run — should skip
    run_flatten_regulons(tmp_results_dir, mock_pyscenic_regulons_csv, overwrite=False)
    assert os.path.getmtime(flat_path) == mtime_before

    # With overwrite — should regenerate
    run_flatten_regulons(tmp_results_dir, mock_pyscenic_regulons_csv, overwrite=True)
    assert os.path.getmtime(flat_path) > mtime_before
