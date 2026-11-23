"""
Normalisation of the concatenated, filtered object.

scran size factors are a PER-CELL scale, so computing them over the full set
of samples gives more stable pooling and puts all samples on the same scale.
Note, however, that scran pooling requires a prior clustering computed on the
ALREADY merged object, and that clustering must exist before normalising;
hence the preliminary step.

Workflow:

    adata = concat_samples(adatas)
    adata = adata[adata.obs["qc_pass"]].copy()        # final filtering

    preliminary_clusters(adata)                       # for the pooling
    scran_size_factors(adata, cluster_key="scran_clusters")
    apply_size_factors(adata)                         # writes the layers
    compare_normalizations(adata)                     # inspect before choosing

    adata.X = adata.layers["norm_scran_log1p"].copy()

WHICH MATRIX IS NORMALISED
--------------------------
The counts layer passed in (by default 'counts_cellbender'). scran expects
counts, not already normalised values. CellBender output with the MCKP
estimator consists of integer counts, so it is suitable. To normalise the raw
counts instead, pass counts_layer='counts_raw', but then the decontamination
does not propagate downstream.

MEMORY
------
The transfer to R is the bottleneck: rpy2 needs a dense genes x cells
matrix in float64. With 40k genes and 20k cells that is ~6.4 GB, which is
exactly what crashes the kernel. For this reason ONLY a subset of expressed
genes is sent. This does not change the result: size factors are per cell and
scran already discards low-expression genes internally via min.mean.
"""

from __future__ import annotations

import warnings
from typing import Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp




def collect_ambient_profiles(
    adatas: Mapping[str, ad.AnnData],
    key: str = "ambient_expression",
) -> pd.DataFrame:
    """
    Per-sample soup profiles in a genes x samples DataFrame.
 
    They must be extracted BEFORE concatenating: they differ between samples,
    so merge='same' drops them. And they are worth keeping: comparing the soup
    across timepoints indicates how much tissue was damaged in each one, which
    in an injury model is information, not noise.
    """
    cols = {}
    for name, a in adatas.items():
        if key not in a.var.columns:
            warnings.warn(f"{name}: var['{key}'] not found")
            continue
        cols[name] = pd.Series(a.var[key].to_numpy(dtype=float), index=a.var_names)
    if not cols:
        raise ValueError(f"No sample has var['{key}']")
    return pd.DataFrame(cols)
 
 

def concat_samples(
    adatas: Mapping[str, ad.AnnData],
    sample_key: str = "gsm",
    index_unique: str = "-",
    drop_obsm: Sequence[str] = ("X_umap_qc",),
    prefix_obs: Sequence[str] = ("qc_leiden",),
    verbose: bool = True,
) -> ad.AnnData:
    """
    Concatenates while preserving var and uns, which ad.concat drops by default.
 
    Three things it does that are easily forgotten when done by hand:
 
    merge='same' / uns_merge='same'
        Without them the whole of var (entrez_id, seqname, gene_biotype...) is
        lost, together with uns['annotation'] holding the GTF provenance. The
        ad.concat defaults are None: discard.
    drop_obsm
        Each sample's X_umap_qc lives in its own space. Concatenating those
        coordinates yields a meaningless embedding that looks valid. Dropped.
    prefix_obs
        Cluster '0' of one sample has nothing to do with cluster '0' of another.
        They are prefixed with the sample name so that they are not merged.
    """
    names = list(adatas.keys())
    var_sets = [tuple(adatas[n].var_names) for n in names]
    if len(set(var_sets)) != 1:
        warnings.warn(
            "The samples do not share the same gene set; with join='outer' "
            "missing genes are filled with zeros. Same GTF for all of them?"
        )
 
    prepared = []
    for n in names:
        a = adatas[n].copy()
        a.obs[sample_key] = n
        for k in drop_obsm:
            a.obsm.pop(k, None)
        for c in prefix_obs:
            if c in a.obs.columns:
                a.obs[c] = pd.Categorical(f"{n}_" + a.obs[c].astype(str))
        prepared.append(a)
 
    out = ad.concat(
        prepared, label=sample_key, keys=names, index_unique=index_unique,
        join="outer", merge="same", uns_merge="same",
    )
 
    # soup profiles do not survive merge='same': they are stored separately
    try:
        out.varm["ambient_expression"] = collect_ambient_profiles(adatas).reindex(
            out.var_names).to_numpy(dtype=np.float32)
        out.uns.setdefault("qc", {})["ambient_samples"] = names
    except ValueError:
        pass
 
    if verbose:
        print(f"[concat] {out.n_obs} cells x {out.n_vars} genes from {len(names)} samples")
        print(out.obs[sample_key].value_counts().to_string())
        kept = [c for c in ("entrez_id", "seqname", "gene_biotype", "gene_symbol")
                if c in out.var.columns]
        print(f"[concat] var retained: {kept if kept else 'NOTHING — check merge='}")
        if not kept:
            warnings.warn("No annotation column in var has survived")
        print(f"[concat] uns retained: {sorted(out.uns)}")
        if "ambient_expression" in out.varm:
            print(f"[concat] soup profiles in varm['ambient_expression'] "
                  f"({out.varm['ambient_expression'].shape})")
    return out
 


# --------------------------------------------------------------------------
# 1. preliminary clustering for scran pooling
# --------------------------------------------------------------------------
 
def preliminary_clusters(
    adata: ad.AnnData,
    counts_layer: str = "counts_cellbender",
    key_added: str = "scran_clusters",
    resolution: float = 0.5,
    n_comps: int = 30,
    n_neighbors: int = 15,
    min_cells_per_cluster: int = 102,
    seed: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Quick, disposable clustering, used only so that scran pools similar
    cells together.
 
    This is not the final clustering: it is computed on a crude normalisation
    (CPM + log1p) precisely because the good one does not exist yet. It is an
    unavoidable circular dependency and is adequately resolved this way; scran
    only needs approximately homogeneous groups, not a correct annotation.
 
    min_cells_per_cluster
        computeSumFactors fails if a cluster has fewer cells than the largest
        pool size. In scranPY the default is sizes=arange(21,102,5), i.e.
        pools of up to 101 cells: hence the minimum here is 102, not 100.
        Smaller clusters are merged into a single one, with a warning.
    """
    tmp = ad.AnnData(
        X=adata.layers[counts_layer].copy() if counts_layer in adata.layers else adata.X.copy(),
        obs=adata.obs[[]].copy(), var=adata.var[[]].copy(),
    )
    sc.pp.normalize_total(tmp)
    sc.pp.log1p(tmp)
    sc.pp.pca(tmp, n_comps=min(n_comps, min(tmp.shape) - 1), random_state=seed)
    sc.pp.neighbors(tmp, n_neighbors=n_neighbors, random_state=seed)
    sc.tl.leiden(tmp, resolution=resolution, key_added="cl", flavor="igraph",
                 n_iterations=2, directed=False)
    pcs = tmp.obsm["X_pca"].copy()
 
    cl = tmp.obs["cl"].astype(str)
    del tmp
 
    sizes = cl.value_counts()
    small = sizes[sizes < min_cells_per_cluster].index.tolist()
    if small:
        cl = cl.where(~cl.isin(small), "small_merged")
        n_small = int(sizes[small].sum())
        if verbose:
            print(f"[scran] {len(small)} clusters with <{min_cells_per_cluster} "
                  f"cells merged into 'small_merged' ({n_small} cells)")
 
        # Merging the small clusters together is not enough if there was only ONE:
        # it stays just as small and the pooling remains degraded. In that case it
        # is absorbed into the nearest large cluster in PCA space.
        if n_small < min_cells_per_cluster:
            big = cl.value_counts().drop("small_merged", errors="ignore")
            if big.empty:
                raise ValueError(
                    f"All clusters have <{min_cells_per_cluster} cells. "
                    f"Lower 'resolution'."
                )
            centroids = pd.DataFrame(pcs).groupby(cl.values).mean()
            target = (centroids.drop(index="small_merged")
                      .sub(centroids.loc["small_merged"])
                      .pow(2).sum(axis=1).idxmin())
            cl = cl.replace("small_merged", target)
            if verbose:
                print(f"[scran] 'small_merged' ({n_small}) was still below the "
                      f"minimum: absorbed into cluster '{target}' "
                      f"({int(big[target])} cells)")
 
    adata.obs[key_added] = pd.Categorical(cl.values)
 
    if verbose:
        vc = adata.obs[key_added].value_counts()
        print(f"[scran] {len(vc)} clusters for the pooling; "
              f"sizes {int(vc.min())}-{int(vc.max())}")
        if vc.min() < min_cells_per_cluster:
            warnings.warn(
                f"There is still a cluster of {int(vc.min())} cells. Lower "
                f"'resolution' or raise 'min_cells_per_cluster'."
            )
    return adata
 
 
# --------------------------------------------------------------------------
# 2. size factors with scran (via rpy2)
# --------------------------------------------------------------------------
 
def scran_size_factors(
    adata: ad.AnnData,
    cluster_key: str = "scran_clusters",
    counts_layer: str = "counts_cellbender",
    key_added: str = "size_factors",
    backend: str = "scranPY",
    min_cells_expressed: int = 20,
    min_mean: Optional[float] = 0.1,
    max_genes: Optional[int] = 8000,
    max_cluster_size: int = 3000,
    algorithm: str = "CVXPY",
    seed: int = 0,
    do_plot: bool = False,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Deconvolution of Lun et al. (computeSumFactors), via two routes.
 
    backend
        'scranPY'  -> sfortma2/scranPY, which implements the full algorithm in
                      Python (pooling + linear system + QR/CVXPY). No R.
        'r'        -> Bioconductor scran via rpy2.
 
        Not to be confused with 'scranpy' from BiocPy/libscran, a different
        package with an almost identical name: it does NOT provide
        computeSumFactors and normalises by library size (that is what
        library_size_factors() is for).
 
    In both cases the bottleneck is the same: both routes need the DENSE
    matrix, so the gene subset is not optional on large
    datasets.
 
    min_cells_expressed / max_genes
        Only genes expressed in at least N cells enter the computation, and
        at most the 'max_genes' with the highest mean expression.
 
        CAVEAT, because this was previously overstated: this is NOT
        exactly identical to using all genes. The deconvolution divides
        by library size and takes medians of ratios, so the gene set
        enters the computation. It is a good approximation because the
        excluded genes contribute a tiny fraction of the counts and would
        have been discarded by min_mean anyway, but it is an approximation,
        not an identity. With max_genes=None all genes are used, memory permitting.
 
    Stores the size factors in obs[key_added] and checks that they are positive:
    scran can yield factors <= 0 for very low-quality cells, and that
    produces silent infinities on division.
    """
    if backend not in ("scranPY", "r"):
        raise ValueError("backend must be 'scranPY' or 'r'")
 
    if cluster_key not in adata.obs:
        raise ValueError(
            f"obs['{cluster_key}'] not found; run preliminary_clusters() first"
        )
 
    X = adata.layers[counts_layer] if counts_layer in adata.layers else adata.X
    X = sp.csc_matrix(X)
 
    n_cells_expr = np.asarray((X > 0).sum(axis=0)).ravel()
    keep = n_cells_expr >= min_cells_expressed
    if max_genes is not None and keep.sum() > max_genes:
        mean_expr = np.asarray(X.mean(axis=0)).ravel()
        mean_expr[~keep] = -1
        keep = np.zeros_like(keep)
        keep[np.argsort(mean_expr)[::-1][:max_genes]] = True
 
    if keep.sum() < 500:
        raise ValueError(
            f"Only {int(keep.sum())} genes pass the filter; scran needs more. "
            f"Lower min_cells_expressed."
        )
 
    if verbose:
        gb = int(keep.sum()) * adata.n_obs * 8 / 1e9
        print(f"[scran] {int(keep.sum())} genes x {adata.n_obs} cells, "
              f"dense: {gb:.2f} GB")
        if gb > 8:
            warnings.warn(f"{gb:.1f} GB when dense; lower max_genes if it crashes")
 
    clusters = adata.obs[cluster_key].astype(str).to_numpy()
 
    # ---------------------------------------------------------------- scranPY
    if backend == "scranPY":
        try:
            from scranPY import compute_sum_factors
        except ImportError as e:
            raise ImportError(
                "pip install git+https://github.com/sfortma2/scranPY.git"
            ) from e
 
        # scranPY requires a dense X and reads the clusters by column NAME
        tmp = ad.AnnData(
            X=np.asarray(X[:, keep].todense(), dtype=np.float64),
            obs=pd.DataFrame({cluster_key: adata.obs[cluster_key].astype("category").values},
                             index=adata.obs_names.copy()),
        )
        compute_sum_factors(
            adata=tmp,
            clusters=cluster_key,
            min_mean=min_mean,
            max_size=max_cluster_size,
            algorithm=algorithm,
            normalize_counts=False,     # do NOT modify X: apply_size_factors does that
            log1p=False,
            plotting=do_plot,
            stopwatch=verbose,
        )
        sf = tmp.obs["size_factors"].to_numpy(dtype=np.float64)
        del tmp
 
        validate_size_factors(sf, verbose=verbose)
        adata.obs[key_added] = sf
        adata.uns.setdefault("normalization", {})["scran"] = {
            "backend": "scranPY", "cluster_key": cluster_key,
            "counts_layer": counts_layer, "n_genes_used": int(keep.sum()),
            "min_mean": min_mean, "algorithm": algorithm,
            "max_cluster_size": max_cluster_size,
        }
        return adata
 
    # -------------------------------------------------------------- rpy2 / R
    mat = np.asarray(X[:, keep].todense(), dtype=np.float64).T   # R: genes x cells
    try:
        import anndata2ri  # noqa: F401
        from rpy2.robjects import numpy2ri, r
        from rpy2.robjects.conversion import localconverter
        import rpy2.robjects as ro
    except ImportError as e:
        raise ImportError(
            "rpy2 (and R with scran) is required. "
            "conda install -c conda-forge rpy2 bioconductor-scran"
        ) from e
 
    with localconverter(ro.default_converter + numpy2ri.converter):
        ro.globalenv["mat"] = mat
        ro.globalenv["clusters"] = ro.StrVector(clusters)
        ro.globalenv["min_mean"] = float(min_mean)
        ro.globalenv["seed"] = int(seed)
        sf = r(
            """
            suppressPackageStartupMessages({library(scran); library(BiocParallel)})
            set.seed(seed)
            sce <- SingleCellExperiment(list(counts = mat))
            sce <- computeSumFactors(sce, clusters = clusters,
                                     min.mean = min_mean,
                                     BPPARAM = MulticoreParam())
            sizeFactors(sce)
            """
        )
        sf = np.asarray(sf, dtype=np.float64).ravel()
 
    if sf.shape[0] != adata.n_obs:
        raise ValueError(
            f"scran returned {sf.shape[0]} size factors for {adata.n_obs} cells"
        )
 
    validate_size_factors(sf, verbose=verbose)
    adata.obs[key_added] = sf
    adata.uns.setdefault("normalization", {})["scran"] = {
        "backend": "r", "cluster_key": cluster_key,
        "counts_layer": counts_layer, "n_genes_sent": int(keep.sum()),
        "min_mean": min_mean, "min_cells_expressed": min_cells_expressed,
    }
    return adata
 
 
def validate_size_factors(sf: np.ndarray, low_warn: float = 0.1,
                          verbose: bool = True) -> None:
    """
    scran can return factors <= 0 for very low-quality cells. Dividing by
    them gives infinities or negative values that later surface as NaN midway
    through the pipeline, far from here. The check is done on the spot.
    """
    bad = ~np.isfinite(sf) | (sf <= 0)
    if bad.any():
        raise ValueError(
            f"{int(bad.sum())} size factors are non-positive or non-finite. This usually "
            f"indicates cells with almost no content that have survived "
            f"QC: inspect them and filter them out before normalising."
        )
    q = np.quantile(sf, [0, .01, .25, .5, .75, .99, 1])
    if verbose:
        print(f"[scran] size factors  min={q[0]:.3f}  1%={q[1]:.3f}  "
              f"median={q[3]:.3f}  99%={q[5]:.3f}  max={q[6]:.3f}")
 
    # The LOW tail is the dangerous one, and the effect is asymmetric: normalising
    # divides, so a size factor of 0.02 multiplies that cell's counts by 50. Those
    # cells end up with huge, noisy values and can form their own spurious
    # cluster. The high tail only compresses, which is benign.
    n_low = int((sf < low_warn).sum())
    if n_low:
        warnings.warn(
            f"{n_low} cells with size factor < {low_warn}: on normalisation their "
            f"counts are amplified up to {1 / max(q[0], 1e-9):.0f}x. These are usually "
            f"low-content cells that passed QC. Inspect them with "
            f"size_factor_report() before continuing."
        )
 
 
def library_size_factors(
    adata: ad.AnnData,
    counts_layer: str = "counts_cellbender",
    key_added: str = "size_factors",
    block_key: Optional[str] = None,
    use_scranpy: bool = True,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Library-size size factors, centred. An alternative to scran WITHOUT R.
 
    BEWARE of the usual confusion: `scranpy` does NOT implement computeSumFactors.
    Its RNA normalisation uses the column sums ("If None, this defaults
    to the column sums of the count matrix"), and its only factor computation
    of its own, compute_clrm1_factors, is for ADT. Pooling with the deconvolution
    of Lun et al. (what distinguishes scran from a library-size
    normalisation) is only available in the R package.
 
    So this is library size, whatever it is called. With use_scranpy=True
    the centring is done by scranpy (mode='lowest', as in its workflow); with
    False it is centred on the mean, which is what scran does in R. The difference
    is a global scale factor, irrelevant after the log.
 
    block_key
        With scranpy, centres per block (e.g. 'gsm') instead of globally.
        Use only with a clear reason: centring per sample removes REAL
        differences in RNA content between samples, which in a time course
        may be biology.
    """
    X = adata.layers[counts_layer] if counts_layer in adata.layers else adata.X
    sf = np.asarray(sp.csr_matrix(X).sum(axis=1)).ravel().astype(np.float64)
 
    if (sf <= 0).any():
        n = int((sf <= 0).sum())
        raise ValueError(
            f"{n} cells with 0 counts in '{counts_layer}'. Filter them out before "
            f"normalising: their size factor is undefined."
        )
 
    how = "mean"
    if use_scranpy:
        try:
            import scranpy
            block = (adata.obs[block_key].astype(str).to_numpy()
                     if block_key is not None else None)
            sf = scranpy.center_size_factors(
                sf, block=block, mode="per-block" if block is not None else "lowest")
            how = f"scranpy ({'per-block' if block is not None else 'lowest'})"
        except ImportError:
            warnings.warn("scranpy is not installed; centring on the mean")
            sf = sf / sf.mean()
    else:
        sf = sf / sf.mean()
 
    validate_size_factors(sf, verbose=verbose)
    adata.obs[key_added] = sf
    adata.uns.setdefault("normalization", {})["library_size"] = {
        "counts_layer": counts_layer, "centering": how, "block_key": block_key,
    }
    if verbose:
        print(f"[libsize] library-size size factors, centred by {how}")
        print("[libsize] this is NOT computeSumFactors: no pooling and no deconvolution")
    return adata
 
 
def size_factor_report(
    adata: ad.AnnData,
    key: str = "size_factors",
    low: float = 0.1,
    cols: Optional[Sequence[str]] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Compares cells with a low size factor against the rest.
 
    A small size factor is not a problem of the normaliser: it is a cell
    with very little content that has survived QC. Dividing by it amplifies
    its counts, and those cells end up dominating the variance and
    forming clusters that do not exist.
 
    If the 'sf_low' row has far fewer counts and genes than 'rest', go back
    to QC and set a min_log1p_counts / min_log1p_genes. That is cleaner than
    patching the size factors afterwards.
    """
    if key not in adata.obs:
        raise ValueError(f"obs['{key}'] not found")
    if cols is None:
        cols = ["n_counts_raw", "n_counts_cellbender", "log1p_total_counts",
                "log1p_n_genes_by_counts", "pct_counts_mt", "nuclear_frac",
                "malat1_log1p", "cellbender_removed_frac", key]
    cols = [c for c in cols if c in adata.obs.columns]
 
    grp = np.where(adata.obs[key].to_numpy() < low, "sf_low", "rest")
    out = adata.obs[cols].groupby(grp, observed=True).median().round(3)
    out.insert(0, "n", pd.Series(grp).value_counts())
 
    if verbose:
        print(out.to_string())
        n_low = int((grp == "sf_low").sum())
        if n_low:
            print(f"\n[sf] {n_low} cells ({100 * n_low / adata.n_obs:.2f}%) "
                  f"with size factor < {low}")
    return out
 
 
# --------------------------------------------------------------------------
# 3. apply and compare
# --------------------------------------------------------------------------
 
def apply_size_factors(
    adata: ad.AnnData,
    counts_layer: str = "counts_cellbender",
    key: str = "size_factors",
    target_sum: float = 1e6,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Stores three normalisations as layers, without touching X. The choice is made later.
 
        norm_log1p        log1p of the unscaled counts (control)
        norm_scran_log1p  log(1 + x / size_factor)
        norm_cpm_log1p    log1p of CPM, the usual reference
 
    X is deliberately left as it is, so that the choice is an explicit step
    by the user (`adata.X = adata.layers[...]`) and not a side effect.
    """
    X = adata.layers[counts_layer] if counts_layer in adata.layers else adata.X
    X = sp.csr_matrix(X, dtype=np.float32)
 
    adata.layers["norm_log1p"] = X.copy()
    adata.layers["norm_log1p"].data = np.log1p(adata.layers["norm_log1p"].data)
 
    if key in adata.obs:
        sf = adata.obs[key].to_numpy(dtype=np.float64)
        validate_size_factors(sf, verbose=False)
        scaled = sp.diags(1.0 / sf) @ X
        scaled = sp.csr_matrix(scaled, dtype=np.float32)
        scaled.data = np.log1p(scaled.data)
        adata.layers["norm_scran_log1p"] = scaled
    elif verbose:
        print(f"[norm] obs['{key}'] not found: the scran layer is not generated")
 
    cpm = sc.pp.normalize_total(
        ad.AnnData(X=X.copy()), target_sum=target_sum, inplace=False)["X"]
    cpm = sp.csr_matrix(cpm, dtype=np.float32)
    cpm.data = np.log1p(cpm.data)
    adata.layers["norm_cpm_log1p"] = cpm
 
    if verbose:
        print(f"[norm] layers: {[k for k in adata.layers if k.startswith('norm_')]}")
        print("[norm] X has NOT been modified; choose with adata.X = adata.layers['...']")
    return adata
 
 
def compare_normalizations(
    adata: ad.AnnData,
    layers: Optional[Sequence[str]] = None,
    group_key: Optional[str] = None,
    n_comps_pca: int = 20,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Compares the normalisations on the only thing that matters here: whether
    sequencing depth still dominates afterwards.
 
    'pc1_vs_depth' is the Spearman correlation between PC1 and the raw counts.
    That is the real failure mode: if depth dominates, it comes out as the first
    principal component and structures the entire UMAP. Close to 0 = good.
 
    BEWARE of the temptation to measure it on the SUM of the normalised row: that
    sum depends mostly on how many genes are detected, which correlates with
    depth even when the normalisation is perfect. It is always ~1 and does not
    discriminate anything. That is why it is measured on PC1.
 
    With group_key (e.g. 'gsm'), the BETWEEN-sample spread of the normalised sum
    is added: if one sample sits systematically higher, the samples are not on
    the same scale and this will be confounded with a batch effect.
    """
    if layers is None:
        layers = [k for k in adata.layers if k.startswith("norm_")]
    if not layers:
        raise ValueError("No normalisation layers found; run apply_size_factors()")
 
    from scipy.stats import spearmanr
 
    depth = (adata.obs["n_counts_raw"].to_numpy()
             if "n_counts_raw" in adata.obs
             else np.asarray(adata.X.sum(axis=1)).ravel())
 
    rows = []
    for lay in layers:
        s = np.asarray(adata.layers[lay].sum(axis=1)).ravel()
 
        # PC1 on this normalisation: the failure mode that matters
        tmp = ad.AnnData(X=sp.csr_matrix(adata.layers[lay], dtype=np.float32))
        sc.pp.pca(tmp, n_comps=min(n_comps_pca, min(tmp.shape) - 1),
                  random_state=0)
        pc1 = tmp.obsm["X_pca"][:, 0]
        var1 = float(tmp.uns["pca"]["variance_ratio"][0])
        del tmp
 
        row = {
            "layer": lay,
            "cv_sum": round(float(s.std() / s.mean()), 4),
            "pc1_vs_depth": round(abs(float(spearmanr(pc1, depth).statistic)), 3),
            "pc1_var_ratio": round(var1, 3),
        }
        if group_key is not None and group_key in adata.obs:
            per = pd.Series(s).groupby(adata.obs[group_key].to_numpy()).mean()
            row["spread_across_samples"] = round(
                float(per.max() / per.min()) if per.min() > 0 else np.nan, 3)
        rows.append(row)
 
    out = pd.DataFrame(rows).set_index("layer")
    if verbose:
        print(out.to_string())
        best = out["pc1_vs_depth"].idxmin()
        print(f"\n[norm] PC1 least contaminated by depth: {best} "
              f"(|rho| = {out.loc[best, 'pc1_vs_depth']:.3f})")
        if out.loc[best, "pc1_vs_depth"] > 0.5:
            print("[norm] WARNING: even the best one leaves PC1 strongly tied to "
                  "depth. Check whether very low-quality cells remain.")
        print("[norm] This is a necessary criterion, not a sufficient one: it shows that "
              "depth no longer dominates, not that the biology is correctly scaled.")
    return out