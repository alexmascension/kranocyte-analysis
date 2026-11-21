import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture

# ---------- helpers ----------
def _prep_Z_from_df(df):
    nf = np.asarray(df["nuclear_RNA_frac"].values, float)
    umi = np.asarray(df["umi"].values, float)

    # mask of valid entries
    ok = np.isfinite(nf) & np.isfinite(umi) & (nf >= 0) & (umi >= 0)
    X = np.c_[np.log10(umi[ok] + 1.0), nf[ok]]  # (logUMI, NF)
    scaler = StandardScaler().fit(X)
    Z = scaler.transform(X)
    return Z, X, scaler, ok

def _pick_empty_by_center(centers_X):
    # "emptiest" center = lowest logUMI and lowest NF
    ranks = np.argsort(centers_X, axis=0)  # per column
    score = ranks[:, 0] + ranks[:, 1]
    return int(np.argmin(score))

# ---------- Stage 1: empties with KMeans ----------
def kmeans_empty(df, n_init=50, random_state=0, nf_rescue=0.05, umi_rescue=1000):
    Z, X, scaler, ok = _prep_Z_from_df(df)

    km = KMeans(n_clusters=2, n_init=n_init, random_state=random_state).fit(Z)
    centers_X = scaler.inverse_transform(km.cluster_centers_)
    empty_k = _pick_empty_by_center(centers_X)

    # labels for the valid rows (ok) only
    lab_valid = np.where(km.labels_ == empty_k, "empty_droplet", "cell")

    # conservative rescue
    nf_m = X[:, 1]
    umi_m = (10**X[:, 0]) - 1.0
    rescue = (nf_m >= nf_rescue) & (umi_m >= umi_rescue)
    lab_valid[rescue] = "cell"

    # build the full column, aligned by index
    labels = np.full(df.shape[0], "unassigned", dtype=object)
    labels[np.where(ok)[0]] = lab_valid
    df["empty_droplet"] = labels

    return centers_X  # centers only

# ---------- Stage 2: damaged vs cell (within 'cells') ----------
def split_cells_into_damaged(df, method="kmeans", random_state=0, nf_damaged_threshold=0.40):
    # subset of cells
    sub = df[df["empty_droplet"] == "cell"].copy()
    if sub.empty:
        df[f"damaged_cell_{method}"] = "unassigned"
        return None, np.nan

    Z, X, scaler, ok = _prep_Z_from_df(sub)
    sub_valid_index = sub.index[ok]  # original (valid) indices

    if method == "kmeans":
        model = KMeans(n_clusters=2, n_init=50, random_state=random_state).fit(Z)
        centers_X = scaler.inverse_transform(model.cluster_centers_)
        comp = model.labels_
    elif method == "gmm":
        model = GaussianMixture(n_components=2, covariance_type="tied",
                                random_state=random_state).fit(Z)
        centers_X = scaler.inverse_transform(model.means_)
        comp = model.predict(Z)
    else:
        raise ValueError("method must be 'kmeans' or 'gmm'")

    # cluster with the highest NF = 'damaged' candidate
    dmg_k = int(np.argmax(centers_X[:, 1]))
    nf_mean_dmg = float(centers_X[dmg_k, 1])

    labels_cells = np.where(comp == dmg_k, "damaged_cell", "cell")
    if nf_mean_dmg < nf_damaged_threshold:
        labels_cells[:] = "cell"  # does not exceed the absolute threshold

    # initialize the column and assign ALIGNED by index
    col = f"damaged_cell_{method}"
    if col not in df.columns:
        df[col] = "unassigned"

    df.loc[sub_valid_index, col] = pd.Series(labels_cells, index=sub_valid_index, dtype=object)

    return centers_X, nf_mean_dmg

# ---------- Pipeline ----------
def classify_empty_and_damaged(df, method_cells="kmeans",
                               nf_rescue=0.05, umi_rescue=1000,
                               nf_damaged_threshold=0.40, random_state=0):
    empty_centers = kmeans_empty(df, nf_rescue=nf_rescue,
                                 umi_rescue=umi_rescue, random_state=random_state)
    cells_centers, nf_mean_dmg = split_cells_into_damaged(
        df, method=method_cells, random_state=random_state,
        nf_damaged_threshold=nf_damaged_threshold
    )

    df[f"classification_{method_cells}"] = df["empty_droplet"]
    df.loc[df["empty_droplet"] == "cell", f"classification_{method_cells}"] = df.loc[df["empty_droplet"] == "cell", f"damaged_cell_{method_cells}"]

    info = {
        "empty_centers_(logUMI,NF)": empty_centers,
        "cells_centers_(logUMI,NF)": cells_centers,
        "nf_mean_damaged": nf_mean_dmg,
        "nf_damaged_threshold": nf_damaged_threshold
    }
    return info
