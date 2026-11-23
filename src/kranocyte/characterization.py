"""
Characterization: subpopulations, markers and annotation.

Reprocessing a subpopulation is NOT a rerun of the whole pipeline. Normalization
has already been done and is global (size factors were computed over all
cells, which is correct: recomputing them per subpopulation would place them on
different scales). What must be redone is gene selection and the embedding,
because the genes that separate tenocytes from FAPs are not the ones that
separate tenocytes FROM EACH OTHER.

    adata_sub = subset_and_reprocess(adata, "major_population", "Tenocyte")
    check_marker_dict(adata, dict_markers_teno)     # BEFORE annotating
    ...

DOUBLE DIPPING
--------------
The p-values from rank_genes_groups on Leiden clusters computed from the same
data are not valid as evidence: the clusters were defined to maximize
precisely the differences that are tested afterwards. They serve to RANK
genes, not to claim that a group is significantly different. If the existence
of a population is the result, external evidence is required: known markers,
another modality, another dataset, or experimental validation.
"""

from __future__ import annotations

import warnings
from typing import Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


# --------------------------------------------------------------------------
# 1. subset and reprocessing
# --------------------------------------------------------------------------

def subset_and_reprocess(
    adata: ad.AnnData,
    key: str,
    value,
    batch_key: str = "gsm",
    counts_layer: str = "counts_cellbender",
    min_cells: int = 10,
    n_comps: int = 15,
    n_neighbors: Optional[int] = None,
    metric: str = "correlation",
    n_features: Optional[int] = None,
    integrate: bool = True,
    umap_min_dist: float = 0.1,
    seed: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Extracts a subpopulation and recomputes variable genes, PCA, harmony,
    neighbors and UMAP on it. Does NOT renormalize.

    min_cells
        Gene filter by the number of CELLS expressing the gene, on the counts
        layer. The previous version used filter_genes(min_counts=20) on X,
        which is now log-normalized: there "counts" are sums of log values and
        the threshold means nothing comparable across subpopulations.

    integrate
        True -> neighbors are built on X_harmony. False -> on X_pca.
        Both embeddings are computed and stored regardless, together with both
        graphs ('neighbors_pca' / 'neighbors_harmony'), to allow comparison.

        With a single sample in the subset, harmony has nothing to correct
        and is skipped automatically.
    """
    if type(value) == str:
        mask = adata.obs[key] == value
    else:
        mask = adata.obs[key].isin(value)


    n = int(mask.sum())
    if n < 50:
        raise ValueError(f"Only {n} cells with {key}=={value!r}; too few to reprocess")

    sub = adata[mask].copy()

    # --- gene filtering on COUNTS, not on the normalized matrix
    if counts_layer in sub.layers:
        n_cells_expr = np.asarray((sub.layers[counts_layer] > 0).sum(axis=0)).ravel()
        src = counts_layer
    else:
        n_cells_expr = np.asarray((sub.X > 0).sum(axis=0)).ravel()
        src = "X"
        warnings.warn(
            f"No layer '{counts_layer}'; filtering genes by cells with a nonzero "
            f"value in X. Valid, but record it."
        )
    keep = n_cells_expr >= min_cells
    if keep.sum() < 200:
        raise ValueError(f"Only {int(keep.sum())} genes pass min_cells={min_cells}")
    sub = sub[:, keep].copy()

    if verbose:
        print(f"[sub] {key}=={value!r}: {n} cells, {int(keep.sum())} genes "
              f"(min_cells={min_cells} on '{src}')")

    reprocess_in_place(
        sub, batch_key=batch_key, n_comps=n_comps, n_neighbors=n_neighbors,
        metric=metric, n_features=n_features, integrate=integrate,
        umap_min_dist=umap_min_dist, seed=seed, verbose=verbose,
    )

    sub.uns.setdefault("processing", {})["subset"] = {
        "key": key, "value": str(value), "n_cells": n,
        "n_genes": int(keep.sum()), "min_cells": min_cells,
    }
    return sub


def reprocess_in_place(
    adata: ad.AnnData,
    *,
    batch_key: str = "gsm",
    n_comps: int = 30,
    n_neighbors: Optional[int] = None,
    metric: str = "cosine",
    n_features: Optional[int] = None,
    integrate: bool = True,
    umap_min_dist: float = 0.1,
    seed: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Recomputes variable genes, PCA, harmony, neighbors and UMAP ON THE OBJECT
    PASSED IN. Does not renormalize and does not modify X.

    Truly IN PLACE: modifies `adata` and returns the same object, so it does
    not matter whether the result is assigned. This is what makes it safe to
    use as `reproc(a, integrate=True)` in a cell, where discarding the return
    value of a function that returns a copy is a silent error: it does not
    fail, it simply does nothing and work continues on the old object.

    Does NOT FILTER GENES, precisely for that reason: removing columns requires
    building a new object and would break the contract. For the filter by
    number of cells, use subset_and_reprocess(), which does return a new object.

    Two rounds, as in subset_and_reprocess: first a scaffold with all genes
    to obtain a graph, then triku on that graph and the final embedding. The
    genes that separate tenocytes from FAPs are not the ones that separate
    tenocytes from each other, so the selection has to be redone in each
    subset.

    integrate
        True  -> neighbors are built on X_harmony
        False -> on X_pca
        Both embeddings and both graphs ('neighbors_pca' / 'neighbors_harmony')
        are computed and stored regardless, to allow comparison. With a single
        batch_key category harmony has nothing to correct and is skipped
        automatically.
    """
    from kranocyte.processing import build_embeddings, select_hvgs_triku

    n_batches = adata.obs[batch_key].nunique() if batch_key in adata.obs else 1
    do_harmony = integrate and n_batches > 1
    if integrate and not do_harmony and verbose:
        print(f"[reproc] a single category of '{batch_key}': harmony does not "
              f"apply, neighbors are built on X_pca")

    if n_neighbors is None:
        n_neighbors = max(10, int(0.5 * adata.n_obs ** 0.5))


    neighbors_on = "harmony" if do_harmony else "pca"

    # --- round 1: scaffold with all genes
    adata.var["hvg_all"] = True
    build_embeddings(adata, "hvg_all", batch_key=batch_key, n_comps=n_comps,
                     n_neighbors=n_neighbors, metric=metric,
                     neighbors_on=neighbors_on, skip_harmony=not do_harmony,
                     seed=seed, verbose=False)

    # --- round 2: triku on that graph, and final embedding
    select_hvgs_triku(adata, n_features=n_features, verbose=verbose)
    build_embeddings(adata, "hvg_triku", batch_key=batch_key, n_comps=n_comps,
                     n_neighbors=n_neighbors, metric=metric,
                     neighbors_on=neighbors_on, skip_harmony=not do_harmony,
                     seed=seed, verbose=verbose)

    sc.tl.umap(adata, min_dist=umap_min_dist, random_state=seed)

    adata.uns.setdefault("processing", {})["reprocess"] = {
        "batch_key": batch_key, "n_comps": n_comps, "n_neighbors": n_neighbors,
        "metric": metric, "integrate": bool(do_harmony),
        "n_batches": int(n_batches), "neighbors_on": neighbors_on,
        "n_hvg_triku": int(adata.var["hvg_triku"].sum())
        if "hvg_triku" in adata.var else None,
    }
    return adata


# --------------------------------------------------------------------------
# 2. check the marker dictionary BEFORE using it
# --------------------------------------------------------------------------

def check_marker_dict(
    adata: ad.AnnData,
    markers: Mapping[str, Sequence[str]],
    ambient_key: str = "ambient_expression",
    ambient_top: int = 50,
    verbose: bool = True,
):
    """
    Checks a marker dictionary against the object, before annotating.

    Three problems that go unnoticed because nothing warns about them:

    1. GENES THAT DO NOT EXIST. With NCBI RefSeq annotation many Ensembl/MGI
       symbols are absent or named differently. assign_cats simply ignores
       them, so a population may be scored with half of its markers and its
       score drops without any indication why.
    2. DUPLICATES within a list: that gene counts twice in the score.
    3. MARKERS THAT ARE SOUP. If a "marker" is among the genes with the highest
       ambient expression, the cluster expressing it may simply be more
       contaminated.

    Returns (resumen, detalle_ausentes).
    """
    var_names = set(adata.var_names)

    soup = set()
    if ambient_key in adata.varm:
        amb = pd.Series(np.asarray(adata.varm[ambient_key]).mean(axis=1),
                        index=adata.var_names)
        soup = set(amb.sort_values(ascending=False).head(ambient_top).index)
    elif ambient_key in adata.var:
        soup = set(adata.var[ambient_key].astype(float)
                   .sort_values(ascending=False).head(ambient_top).index)

    rows, faltan = [], {}
    for pop, genes in markers.items():
        genes = list(genes)
        presentes = [g for g in genes if g in var_names]
        ausentes = [g for g in genes if g not in var_names]
        dups = sorted({g for g in genes if genes.count(g) > 1})
        en_soup = sorted(set(presentes) & soup)
        rows.append({
            "population": pop,
            "n_markers": len(genes),
            "n_unique": len(set(genes)),
            "present": len(set(presentes)),
            "missing": len(set(ausentes)),
            "pct_present": round(100 * len(set(presentes)) / max(len(set(genes)), 1), 1),
            "duplicates": ", ".join(dups),
            "en_soup": ", ".join(en_soup),
        })
        if ausentes:
            faltan[pop] = sorted(set(ausentes))

    resumen = pd.DataFrame(rows).set_index("population")

    if verbose:
        print(resumen.to_string())
        if faltan:
            print("\ngenes NOT found in var_names:")
            for pop, gs in faltan.items():
                print(f"  {pop}: {gs}")
            print("\n[markers] assign_cats ignores them silently. Check whether "
                  "this is a nomenclature issue (RefSeq vs Ensembl) before "
                  "accepting the scores.")
        flojas = resumen.index[resumen["pct_present"] < 70].tolist()
        if flojas:
            print(f"\n[markers] WARNING: {flojas} have less than 70% of their "
                  f"markers present; their score is not comparable to the others.")
        consoup = resumen.index[resumen["en_soup"] != ""].tolist()
        if consoup:
            print(f"[markers] {consoup} include ambient genes. A cluster "
                  f"may score high there only because it is more contaminated.")
    return resumen, faltan


# --------------------------------------------------------------------------
# 3. composition by sample / condition
# --------------------------------------------------------------------------

def population_composition(
    adata: ad.AnnData,
    population_key: str,
    group_key: str,
    normalize: str = "index",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Proportion of each population per sample or condition.

    This IS read from the labels, not from the embedding, so it is
    independent of whether integration was applied. It is the correct route
    to compare composition across conditions even if annotation used X_harmony.

    Interpret with care: in a time course the proportion is compositional
    (sums to 1), so if one population grows ALL the others decrease in
    proportion even if their absolute number does not change. Check the counts too.
    """
    cont = pd.crosstab(adata.obs[group_key], adata.obs[population_key])
    prop = pd.crosstab(adata.obs[group_key], adata.obs[population_key],
                       normalize=normalize)
    if verbose:
        print("counts:")
        print(cont.to_string())
        print("\nproportions:")
        print((100 * prop).round(2).to_string())
        print("\n[comp] These are proportions: each row sums to 100. A relative "
              "increase can arise because ANOTHER population has decreased.")
    return prop