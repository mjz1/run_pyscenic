"""Subsampling strategies for single-cell expression matrices.

Provides random, stratified, and geosketch subsampling methods that return
sorted integer index arrays suitable for slicing AnnData objects.
"""

import numpy as np


def subsample(
    adata,
    method,
    n,
    seed=42,
    annotation_column=None,
    min_per_group=50,
    embedding=None,
):
    """Dispatch to the requested subsampling strategy.

    Args:
        adata: AnnData object (used for cell count, obs, and X).
        method: One of "random", "stratified", "geosketch".
        n: Target number of cells to select.
        seed: Random seed for reproducibility.
        annotation_column: obs column for stratified sampling.
        min_per_group: Minimum cells per group for stratified sampling.
        embedding: Optional pre-computed embedding matrix (n_cells x n_dims)
            for geosketch. If None, a naive PCA is computed.

    Returns:
        Sorted 1-D numpy int array of selected cell indices.
    """
    if method == "random":
        return subsample_random(adata.n_obs, n, seed=seed)
    elif method == "stratified":
        if annotation_column is None:
            raise ValueError("annotation_column is required for stratified subsampling")
        if annotation_column not in adata.obs.columns:
            raise ValueError(
                f"annotation_column '{annotation_column}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        labels = adata.obs[annotation_column].values
        return subsample_stratified(labels, n, seed=seed, min_per_group=min_per_group)
    elif method == "geosketch":
        return subsample_geosketch(adata, n, seed=seed, embedding=embedding)
    else:
        raise ValueError(f"Unknown subsampling method: {method}")


def subsample_random(n_total, n, seed=42):
    """Uniform random subsampling.

    Args:
        n_total: Total number of cells.
        n: Target number of cells to select.
        seed: Random seed.

    Returns:
        Sorted 1-D numpy int array of selected cell indices.
    """
    rng = np.random.default_rng(seed)
    indices = rng.choice(n_total, size=n, replace=False)
    return np.sort(indices)


def subsample_stratified(labels, n, seed=42, min_per_group=50):
    """Stratified subsampling that preserves group proportions.

    For each group, the number of sampled cells is proportional to the group's
    share of the total, but at least ``min_per_group`` cells are taken (when
    available). If the proportional allocation for a group is less than
    ``min_per_group``, the group receives ``min(group_size, min_per_group)``
    cells. Remaining budget is distributed proportionally among larger groups.

    Args:
        labels: 1-D array-like of group labels (length = total cells).
        n: Target number of cells to select.
        seed: Random seed.
        min_per_group: Minimum cells per group.

    Returns:
        Sorted 1-D numpy int array of selected cell indices.
    """
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    unique_groups, group_counts = np.unique(labels, return_counts=True)

    # Compute proportional allocations, enforce min_per_group
    proportions = group_counts / group_counts.sum()
    raw_alloc = proportions * n

    alloc = {}
    locked_budget = 0
    unlocked_groups = []

    for grp, cnt, raw in zip(unique_groups, group_counts, raw_alloc):
        grp_key = grp
        floor = min(cnt, min_per_group)
        if raw < floor:
            # Give this group its floor (or full group if smaller)
            alloc[grp_key] = min(cnt, floor)
            locked_budget += alloc[grp_key]
        else:
            unlocked_groups.append((grp_key, cnt, raw))

    # Distribute remaining budget proportionally among unlocked groups
    remaining = n - locked_budget
    if unlocked_groups:
        unlocked_raw_sum = sum(r for _, _, r in unlocked_groups)
        for grp_key, cnt, raw in unlocked_groups:
            share = int(round(remaining * (raw / unlocked_raw_sum)))
            alloc[grp_key] = min(cnt, share)

    # Adjust total to hit exactly n (rounding may cause drift)
    total_alloc = sum(alloc.values())
    diff = n - total_alloc
    if diff != 0:
        # Add/remove from largest groups first
        sorted_groups = sorted(alloc, key=lambda g: alloc[g], reverse=True)
        for grp_key in sorted_groups:
            if diff == 0:
                break
            grp_cnt = group_counts[np.where(unique_groups == grp_key)[0][0]]
            if diff > 0 and alloc[grp_key] < grp_cnt:
                add = min(diff, grp_cnt - alloc[grp_key])
                alloc[grp_key] += add
                diff -= add
            elif diff < 0 and alloc[grp_key] > 0:
                remove = min(-diff, alloc[grp_key])
                alloc[grp_key] -= remove
                diff += remove

    # Sample indices for each group
    selected = []
    for grp in unique_groups:
        grp_indices = np.where(labels == grp)[0]
        grp_n = alloc[grp]
        if grp_n >= len(grp_indices):
            chosen = grp_indices
        else:
            chosen = rng.choice(grp_indices, size=grp_n, replace=False)
        selected.append(chosen)

    indices = np.concatenate(selected)
    return np.sort(indices)


def subsample_geosketch(adata, n, seed=42, embedding=None):
    """Geometric sketching via the geosketch package.

    If ``embedding`` is provided it is used directly; otherwise a naive PCA
    is computed on a copy of the data.

    Args:
        adata: AnnData object (only used when embedding is None).
        n: Target number of cells to select.
        seed: Random seed.
        embedding: Optional pre-computed embedding array (n_cells x n_dims).

    Returns:
        Sorted 1-D numpy int array of selected cell indices.
    """
    try:
        from geosketch import gs
    except ImportError:
        raise ImportError(
            "geosketch is not installed. Install it with: pip install geosketch"
        )

    if embedding is not None:
        X_emb = np.asarray(embedding)
    else:
        import scanpy as sc

        # Work on a copy to avoid mutating the caller's object
        adata_tmp = adata.copy()
        sc.pp.pca(adata_tmp, random_state=seed)
        X_emb = adata_tmp.obsm["X_pca"]

    sketch_indices = gs(X_emb, n, replace=False, seed=seed)
    indices = np.sort(np.asarray(sketch_indices).flatten().astype(int))
    return indices
