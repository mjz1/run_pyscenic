# This code is adapted from pyscenic source files and modified such that it works with the run_pyscenic venv (with numpy>2)
import collections
from typing import Mapping, Optional

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats
from scipy.optimize import minimize_scalar
from sklearn import mixture


def _gcm_(cdf, idxs):
    work_cdf = cdf
    work_idxs = idxs
    gcm = [work_cdf[0]]
    touchpoints = [0]
    while len(work_cdf) > 1:
        distances = work_idxs[1:] - work_idxs[0]
        slopes = (work_cdf[1:] - work_cdf[0]) / distances
        minslope = slopes.min()
        minslope_idx = np.where(slopes == minslope)[0][0] + 1
        gcm.extend(work_cdf[0] + distances[:minslope_idx] * minslope)
        touchpoints.append(touchpoints[-1] + minslope_idx)
        work_cdf = work_cdf[minslope_idx:]
        work_idxs = work_idxs[minslope_idx:]
    return np.array(np.array(gcm)), np.array(touchpoints)


def _lcm_(cdf, idxs):
    g, t = _gcm_(1 - cdf[::-1], idxs.max() - idxs[::-1])
    return 1 - g[::-1], len(cdf) - 1 - t[::-1]


def _touch_diffs_(part1, part2, touchpoints):
    diff = np.abs((part2[touchpoints] - part1[touchpoints]))
    return diff.max(), diff


def diptst(dat, is_hist=False, numt=1000):
    """diptest with pval"""
    # sample dip
    d, (_, idxs, left, _, right, _) = dip_fn(dat, is_hist)

    # simulate from null uniform
    unifs = np.random.uniform(size=numt * idxs.shape[0]).reshape([numt, idxs.shape[0]])
    unif_dips = np.apply_along_axis(dip_fn, 1, unifs, is_hist, True)

    # count dips greater or equal to d, add 1/1 to prevent a pvalue of 0
    pval = (
        None
        if unif_dips.sum() == 0
        else (np.less(d, unif_dips).sum() + 1) / (float(numt) + 1)
    )

    return (d, pval, (len(left) - 1, len(idxs) - len(right)))  # dip, pvalue  # indices


def dip_fn(dat, is_hist=False, just_dip=False):
    """
    Compute the Hartigans' dip statistic either for a histogram of
    samples (with equidistant bins) or for a set of samples.
    """
    if is_hist:
        histogram = dat
        idxs = np.arange(len(histogram))
    else:
        counts = collections.Counter(dat)
        idxs = np.sort(list(counts.keys()))
        histogram = np.array([counts[i] for i in idxs])

    # check for case 1<N<4 or all identical values
    if len(idxs) <= 4 or idxs[0] == idxs[-1]:
        left = []
        right = [1]
        d = 0.0
        return d if just_dip else (d, (None, idxs, left, None, right, None))

    cdf = np.cumsum(histogram, dtype=float)
    cdf /= cdf[-1]

    work_idxs = idxs
    work_histogram = np.asarray(histogram, dtype=float) / np.sum(histogram)
    work_cdf = cdf

    D = 0
    left = [0]
    right = [1]

    while True:
        left_part, left_touchpoints = _gcm_(work_cdf - work_histogram, work_idxs)
        right_part, right_touchpoints = _lcm_(work_cdf, work_idxs)

        d_left, left_diffs = _touch_diffs_(left_part, right_part, left_touchpoints)
        d_right, right_diffs = _touch_diffs_(left_part, right_part, right_touchpoints)

        if d_right > d_left:
            xr = right_touchpoints[d_right == right_diffs][-1]
            xl = left_touchpoints[left_touchpoints <= xr][-1]
            d = d_right
        else:
            xl = left_touchpoints[d_left == left_diffs][0]
            xr = right_touchpoints[right_touchpoints >= xl][0]
            d = d_left

        left_diff = np.abs(left_part[: xl + 1] - work_cdf[: xl + 1]).max()
        right_diff = np.abs(right_part[xr:] - work_cdf[xr:] + work_histogram[xr:]).max()

        if d <= D or xr == 0 or xl == len(work_cdf):
            the_dip = max(
                np.abs(cdf[: len(left)] - left).max(),
                np.abs(cdf[-len(right) - 1 : -1] - right).max(),
            )
            if just_dip:
                return the_dip / 2
            else:
                return the_dip / 2, (cdf, idxs, left, left_part, right, right_part)
        else:
            D = max(D, left_diff, right_diff)

        work_cdf = work_cdf[xl : xr + 1]
        work_idxs = work_idxs[xl : xr + 1]
        work_histogram = work_histogram[xl : xr + 1]

        left[len(left) :] = left_part[1 : xl + 1]
        right[:0] = right_part[xr:-1]


def derive_threshold(
    auc_mtx: pd.DataFrame, regulon_name: str, seed=None, method: str = "hdt"
) -> float:
    """
    Derive threshold on the AUC values of the given regulon to binarize the cells in two clusters: "on" versus "off"
    state of the regulator.

    :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
    :param regulon_name: the name of the regulon for which to predict the threshold.
    :param method: The method to use to decide if the distribution of AUC values for the given regulon is not unimodel.
        Can be either Hartigan's Dip Test (HDT) or Bayesian Information Content (BIC). The former method performs better
        but takes considerable more time to execute (40min for 350 regulons). The BIC compares the BIC for two Gaussian
        Mixture Models: single versus two components.
    :return: The threshold on the AUC values.
    """
    assert auc_mtx is not None and not auc_mtx.empty
    assert regulon_name in auc_mtx.columns
    assert method in {"hdt", "bic"}

    data = auc_mtx[regulon_name].values

    if seed:
        np.random.seed(seed=seed)

    def isbimodal(data, method):
        if method == "hdt":
            # Use Hartigan's dip statistic to decide if distribution deviates from unimodality.
            _, pval, _ = diptst(np.sort(data))
            return (pval is not None) and (pval <= 0.05)
        else:
            # Compare Bayesian Information Content of two Gaussian Mixture Models.
            X = data.reshape(-1, 1)
            gmm2 = mixture.GaussianMixture(
                n_components=2, covariance_type="full", random_state=seed
            ).fit(X)
            gmm1 = mixture.GaussianMixture(
                n_components=1, covariance_type="full", random_state=seed
            ).fit(X)
            return gmm2.bic(X) <= gmm1.bic(X)

    if not isbimodal(data, method):
        # For a unimodal distribution the threshold is set as mean plus two standard deviations.
        return data.mean() + 2.0 * data.std()
    else:
        # Fit a two component Gaussian Mixture model on the AUC distribution using an Expectation-Maximization algorithm
        # to identify the peaks in the distribution.
        gmm2 = mixture.GaussianMixture(
            n_components=2, covariance_type="full", random_state=seed
        ).fit(data.reshape(-1, 1))
        # For a bimodal distribution the threshold is defined as the "trough" in between the two peaks.
        # This is solved as a minimization problem on the kernel smoothed density.
        return minimize_scalar(
            fun=stats.gaussian_kde(data), bounds=sorted(gmm2.means_), method="bounded"
        ).x[0]


def binarize(
    auc_mtx: pd.DataFrame,
    threshold_overides: Optional[Mapping[str, float]] = None,
    seed=None,
    num_workers=1,
) -> (pd.DataFrame, pd.Series):
    """
    "Binarize" the supplied AUC matrix, i.e. decide if for each cells in the matrix a regulon is active or not based
    on the bimodal distribution of the AUC values for that regulon.

    :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
    :param threshold_overides: A dictionary that maps name of regulons to manually set thresholds.
    :return: A "binarized" dataframe and a series containing the AUC threshold used for each regulon.
    """

    def derive_thresholds(auc_mtx, seed=seed):
        thrs = Parallel(n_jobs=num_workers, backend="loky", verbose=10)(
            delayed(derive_threshold)(auc_mtx, c, seed) for c in auc_mtx.columns
        )
        return pd.Series(index=auc_mtx.columns, data=thrs)

    thresholds = derive_thresholds(auc_mtx)
    if threshold_overides is not None:
        thresholds[list(threshold_overides.keys())] = list(threshold_overides.values())
    return (auc_mtx > thresholds).astype(int), thresholds
