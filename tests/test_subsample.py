"""Tests for the subsampling module."""

import numpy as np
import pandas as pd
import pytest

import anndata as ad

from subsample import (
    subsample,
    subsample_random,
    subsample_stratified,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def simple_adata():
    """AnnData with 500 cells, 100 genes, 3 cell types."""
    np.random.seed(0)
    n_obs, n_vars = 500, 100
    X = np.random.poisson(lam=3, size=(n_obs, n_vars)).astype(np.float32)

    labels = np.array(["A"] * 250 + ["B"] * 150 + ["C"] * 100)
    obs = pd.DataFrame({"cell_type": labels}, index=[f"cell_{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])

    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def single_group_adata():
    """AnnData where all cells belong to one group."""
    np.random.seed(1)
    n_obs, n_vars = 200, 50
    X = np.random.poisson(lam=2, size=(n_obs, n_vars)).astype(np.float32)
    obs = pd.DataFrame(
        {"cell_type": ["only_group"] * n_obs},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    return ad.AnnData(X=X, obs=obs, var=var)


# ---------------------------------------------------------------------------
# Random subsampling
# ---------------------------------------------------------------------------


class TestSubsampleRandom:
    def test_returns_correct_count(self):
        indices = subsample_random(1000, 200, seed=42)
        assert len(indices) == 200

    def test_indices_within_bounds(self):
        n_total = 500
        indices = subsample_random(n_total, 100, seed=7)
        assert indices.min() >= 0
        assert indices.max() < n_total

    def test_indices_sorted(self):
        indices = subsample_random(1000, 300, seed=99)
        assert np.all(indices[:-1] <= indices[1:])

    def test_no_duplicates(self):
        indices = subsample_random(500, 200, seed=5)
        assert len(set(indices)) == 200

    def test_seed_reproducibility(self):
        a = subsample_random(1000, 200, seed=123)
        b = subsample_random(1000, 200, seed=123)
        np.testing.assert_array_equal(a, b)

    def test_different_seeds_differ(self):
        a = subsample_random(1000, 200, seed=1)
        b = subsample_random(1000, 200, seed=2)
        assert not np.array_equal(a, b)


# ---------------------------------------------------------------------------
# Stratified subsampling
# ---------------------------------------------------------------------------


class TestSubsampleStratified:
    def test_returns_correct_approximate_count(self, simple_adata):
        n = 100
        labels = simple_adata.obs["cell_type"].values
        indices = subsample_stratified(labels, n, seed=42)
        assert len(indices) == n

    def test_indices_within_bounds(self, simple_adata):
        labels = simple_adata.obs["cell_type"].values
        indices = subsample_stratified(labels, 100, seed=42)
        assert indices.min() >= 0
        assert indices.max() < simple_adata.n_obs

    def test_indices_sorted(self, simple_adata):
        labels = simple_adata.obs["cell_type"].values
        indices = subsample_stratified(labels, 100, seed=42)
        assert np.all(indices[:-1] <= indices[1:])

    def test_preserves_proportions(self, simple_adata):
        """Sampled group proportions should be close to original proportions."""
        labels = simple_adata.obs["cell_type"].values
        n = 200
        indices = subsample_stratified(labels, n, seed=42, min_per_group=1)

        orig_counts = pd.Series(labels).value_counts(normalize=True).sort_index()
        sub_labels = labels[indices]
        sub_counts = pd.Series(sub_labels).value_counts(normalize=True).sort_index()

        for grp in orig_counts.index:
            assert grp in sub_counts.index
            assert abs(orig_counts[grp] - sub_counts[grp]) < 0.1  # within 10%

    def test_min_per_group_respected(self):
        """Small groups should get at least min_per_group cells (if available)."""
        labels = np.array(["big"] * 900 + ["tiny"] * 10)
        indices = subsample_stratified(labels, 100, seed=42, min_per_group=10)
        sub_labels = labels[indices]
        tiny_count = (sub_labels == "tiny").sum()
        assert tiny_count == 10  # tiny group has exactly 10 cells total

    def test_single_group(self, single_group_adata):
        labels = single_group_adata.obs["cell_type"].values
        indices = subsample_stratified(labels, 50, seed=42)
        assert len(indices) == 50
        assert indices.min() >= 0
        assert indices.max() < single_group_adata.n_obs

    def test_seed_reproducibility(self, simple_adata):
        labels = simple_adata.obs["cell_type"].values
        a = subsample_stratified(labels, 100, seed=42)
        b = subsample_stratified(labels, 100, seed=42)
        np.testing.assert_array_equal(a, b)


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


class TestSubsampleDispatcher:
    def test_random_dispatch(self, simple_adata):
        indices = subsample(simple_adata, method="random", n=50, seed=42)
        assert len(indices) == 50
        assert indices.min() >= 0
        assert indices.max() < simple_adata.n_obs

    def test_stratified_dispatch(self, simple_adata):
        indices = subsample(
            simple_adata,
            method="stratified",
            n=50,
            seed=42,
            annotation_column="cell_type",
        )
        assert len(indices) == 50

    def test_stratified_missing_column_raises(self, simple_adata):
        with pytest.raises(ValueError, match="not found in adata.obs"):
            subsample(
                simple_adata,
                method="stratified",
                n=50,
                annotation_column="nonexistent",
            )

    def test_stratified_no_column_raises(self, simple_adata):
        with pytest.raises(ValueError, match="annotation_column is required"):
            subsample(simple_adata, method="stratified", n=50)

    def test_unknown_method_raises(self, simple_adata):
        with pytest.raises(ValueError, match="Unknown subsampling method"):
            subsample(simple_adata, method="bogus", n=50)

    def test_geosketch_import_error(self, simple_adata, monkeypatch):
        """If geosketch is not installed, a clear ImportError is raised."""
        import builtins

        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "geosketch":
                raise ImportError("No module named 'geosketch'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)
        with pytest.raises(ImportError, match="geosketch is not installed"):
            subsample(simple_adata, method="geosketch", n=50)

    def test_geosketch_with_precomputed_embedding(self, simple_adata, monkeypatch):
        """Dispatcher passes embedding through to subsample_geosketch."""
        captured = {}

        def fake_geosketch(adata, n, seed=42, embedding=None):
            captured["embedding"] = embedding
            return np.arange(n)

        monkeypatch.setattr("subsample.subsample_geosketch", fake_geosketch)
        emb = np.random.default_rng(0).standard_normal((simple_adata.n_obs, 20))
        subsample(simple_adata, method="geosketch", n=50, embedding=emb)
        np.testing.assert_array_equal(captured["embedding"], emb)


# ---------------------------------------------------------------------------
# _resolve_geosketch_embedding
# ---------------------------------------------------------------------------


class TestResolveGeosketchEmbedding:
    """Tests for the embedding resolution logic in run_pyscenic."""

    @pytest.fixture
    def adata_with_obsm(self, simple_adata):
        """AnnData with a pre-computed embedding in obsm."""
        emb = np.random.default_rng(99).standard_normal((simple_adata.n_obs, 30))
        simple_adata.obsm["X_pca_custom"] = emb
        return simple_adata

    @pytest.fixture
    def adata_with_batch(self, simple_adata):
        """AnnData with a batch column in obs."""
        simple_adata.obs["batch"] = ["batch_0"] * 250 + ["batch_1"] * 250
        return simple_adata

    def test_obsm_key(self, adata_with_obsm, tmp_results_dir):
        from run_pyscenic import _resolve_geosketch_embedding

        emb, meta = _resolve_geosketch_embedding(
            adata_with_obsm,
            tmp_results_dir,
            overwrite=False,
            embedding_key="X_pca_custom",
        )
        assert emb.shape == (adata_with_obsm.n_obs, 30)
        assert meta["embedding_source"] == "obsm"
        assert meta["embedding_key"] == "X_pca_custom"
        assert meta["embedding_dimensions"] == 30

    def test_obsm_key_missing_raises(self, simple_adata, tmp_results_dir):
        from run_pyscenic import _resolve_geosketch_embedding

        with pytest.raises(ValueError, match="not found in adata.obsm"):
            _resolve_geosketch_embedding(
                simple_adata,
                tmp_results_dir,
                overwrite=False,
                embedding_key="X_nonexistent",
            )

    def test_naive_pca_fallback(self, simple_adata, tmp_results_dir):
        from run_pyscenic import _resolve_geosketch_embedding

        emb, meta = _resolve_geosketch_embedding(
            simple_adata,
            tmp_results_dir,
            overwrite=False,
            n_pca_components=20,
            seed=42,
        )
        assert emb.shape == (simple_adata.n_obs, 20)
        assert meta["embedding_source"] == "pca"
        assert meta["n_pca_components"] == 20
        assert meta["batch_correction_method"] == "none"

    def test_adata_X_unchanged_after_pca(self, simple_adata, tmp_results_dir):
        """adata.X must not be modified by PCA computation."""
        from run_pyscenic import _resolve_geosketch_embedding

        X_orig = simple_adata.X.copy()
        _resolve_geosketch_embedding(
            simple_adata,
            tmp_results_dir,
            overwrite=False,
            seed=42,
        )
        np.testing.assert_array_equal(simple_adata.X, X_orig)

    def test_embedding_cached(self, simple_adata, tmp_results_dir):
        import os
        from run_pyscenic import _resolve_geosketch_embedding

        _resolve_geosketch_embedding(
            simple_adata, tmp_results_dir, overwrite=False, seed=42
        )
        cache_path = os.path.join(tmp_results_dir, "subsample_embedding.npy")
        assert os.path.exists(cache_path)

    def test_cached_embedding_loaded(self, simple_adata, tmp_results_dir):
        from run_pyscenic import _resolve_geosketch_embedding

        # First call: compute and cache
        emb1, _ = _resolve_geosketch_embedding(
            simple_adata, tmp_results_dir, overwrite=False, seed=42
        )
        # Second call: should load from cache
        emb2, meta2 = _resolve_geosketch_embedding(
            simple_adata, tmp_results_dir, overwrite=False, seed=42
        )
        assert meta2["embedding_source"] == "cached"
        np.testing.assert_array_equal(emb1, emb2)

    def test_overwrite_recomputes(self, simple_adata, tmp_results_dir):
        from run_pyscenic import _resolve_geosketch_embedding

        emb1, _ = _resolve_geosketch_embedding(
            simple_adata, tmp_results_dir, overwrite=False, seed=42
        )
        emb2, meta2 = _resolve_geosketch_embedding(
            simple_adata, tmp_results_dir, overwrite=True, seed=42
        )
        assert meta2["embedding_source"] == "pca"  # recomputed, not cached
        np.testing.assert_array_equal(emb1, emb2)

    def test_dimension_warning_low(self, simple_adata, tmp_results_dir, caplog):
        from run_pyscenic import _resolve_geosketch_embedding
        import logging

        with caplog.at_level(logging.WARNING):
            _resolve_geosketch_embedding(
                simple_adata,
                tmp_results_dir,
                overwrite=True,
                n_pca_components=3,
                seed=42,
            )
        assert any("< 5" in msg for msg in caplog.messages)

    def test_harmony_batch_correction(
        self, adata_with_batch, tmp_results_dir, monkeypatch
    ):
        """Test Harmony path with a mock."""
        from run_pyscenic import _resolve_geosketch_embedding

        class FakeHarmonyResult:
            def __init__(self, Z):
                self.Z_corr = Z.T  # harmony returns (n_dims, n_cells)

        def fake_run_harmony(X_pca, obs, batch_key, random_state=0):
            return FakeHarmonyResult(X_pca)

        import types

        fake_module = types.ModuleType("harmonypy")
        fake_module.run_harmony = fake_run_harmony
        monkeypatch.setitem(__import__("sys").modules, "harmonypy", fake_module)

        emb, meta = _resolve_geosketch_embedding(
            adata_with_batch,
            tmp_results_dir,
            overwrite=True,
            batch_correction_method="harmony",
            batch_key="batch",
            n_pca_components=20,
            seed=42,
        )
        assert meta["embedding_source"] == "harmony"
        assert meta["batch_correction_method"] == "harmony"
        assert meta["batch_key"] == "batch"
        assert emb.shape[0] == adata_with_batch.n_obs

    def test_scanorama_batch_correction(
        self, adata_with_batch, tmp_results_dir, monkeypatch
    ):
        """Test Scanorama path with a mock."""
        from run_pyscenic import _resolve_geosketch_embedding

        def fake_integrate_scanpy(adata_list):
            for ad_b in adata_list:
                ad_b.obsm["X_scanorama"] = np.random.default_rng(0).standard_normal(
                    (ad_b.n_obs, 15)
                )

        import types

        fake_module = types.ModuleType("scanorama")
        fake_module.integrate_scanpy = fake_integrate_scanpy
        monkeypatch.setitem(__import__("sys").modules, "scanorama", fake_module)

        emb, meta = _resolve_geosketch_embedding(
            adata_with_batch,
            tmp_results_dir,
            overwrite=True,
            batch_correction_method="scanorama",
            batch_key="batch",
            n_pca_components=20,
            seed=42,
        )
        assert meta["embedding_source"] == "scanorama"
        assert meta["batch_correction_method"] == "scanorama"
        assert meta["batch_key"] == "batch"
        assert emb.shape == (adata_with_batch.n_obs, 15)

    def test_harmony_missing_raises(
        self, adata_with_batch, tmp_results_dir, monkeypatch
    ):
        """Missing harmonypy should raise clear ImportError."""
        from run_pyscenic import _resolve_geosketch_embedding
        import sys

        monkeypatch.delitem(sys.modules, "harmonypy", raising=False)

        import builtins

        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "harmonypy":
                raise ImportError("No module named 'harmonypy'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)
        with pytest.raises(ImportError, match="harmonypy is not installed"):
            _resolve_geosketch_embedding(
                adata_with_batch,
                tmp_results_dir,
                overwrite=True,
                batch_correction_method="harmony",
                batch_key="batch",
            )

    def test_scanorama_missing_raises(
        self, adata_with_batch, tmp_results_dir, monkeypatch
    ):
        """Missing scanorama should raise clear ImportError."""
        from run_pyscenic import _resolve_geosketch_embedding
        import sys

        monkeypatch.delitem(sys.modules, "scanorama", raising=False)

        import builtins

        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "scanorama":
                raise ImportError("No module named 'scanorama'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)
        with pytest.raises(ImportError, match="scanorama is not installed"):
            _resolve_geosketch_embedding(
                adata_with_batch,
                tmp_results_dir,
                overwrite=True,
                batch_correction_method="scanorama",
                batch_key="batch",
            )

    def test_batch_key_missing_raises(self, simple_adata, tmp_results_dir):
        from run_pyscenic import _resolve_geosketch_embedding

        with pytest.raises(ValueError, match="not found in adata.obs"):
            _resolve_geosketch_embedding(
                simple_adata,
                tmp_results_dir,
                overwrite=True,
                batch_correction_method="harmony",
                batch_key="nonexistent_column",
            )

    def test_unknown_method_raises(self, simple_adata, tmp_results_dir):
        from run_pyscenic import _resolve_geosketch_embedding

        with pytest.raises(ValueError, match="Unknown batch correction method"):
            _resolve_geosketch_embedding(
                simple_adata,
                tmp_results_dir,
                overwrite=True,
                batch_correction_method="bogus",
            )


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_subsample_n_equals_total(self):
        """When n == n_total, all indices should be returned."""
        indices = subsample_random(100, 100, seed=42)
        assert len(indices) == 100
        np.testing.assert_array_equal(indices, np.arange(100))

    def test_stratified_n_exceeds_group(self):
        """When n > total in one group, that group contributes all its cells."""
        labels = np.array(["A"] * 5 + ["B"] * 995)
        indices = subsample_stratified(labels, 500, seed=42, min_per_group=5)
        sub_labels = labels[indices]
        # Group A has only 5 cells; all should be selected
        assert (sub_labels == "A").sum() == 5

    def test_stratified_large_min_per_group(self):
        """min_per_group larger than some groups should still work."""
        labels = np.array(["X"] * 3 + ["Y"] * 997)
        indices = subsample_stratified(labels, 100, seed=42, min_per_group=50)
        sub_labels = labels[indices]
        # X only has 3 cells, can't meet min_per_group, takes all 3
        assert (sub_labels == "X").sum() == 3
