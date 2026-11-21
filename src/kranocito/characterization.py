"""
Caracterización: subpoblaciones, marcadores y anotación.

Reprocesar una subpoblación NO es repetir el pipeline entero. La normalización
ya está hecha y es global (los size factors se calcularon sobre todas las
células, que es lo correcto: rehacerlos por subpoblación las pondría en escalas
distintas). Lo que sí hay que rehacer es la selección de genes y el embedding,
porque los genes que separan tenocitos de FAPs no son los que separan tenocitos
ENTRE SÍ.

    adata_sub = subset_and_reprocess(adata, "major_population", "Tenocyte")
    check_marker_dict(adata, dict_markers_teno)     # ANTES de anotar
    ...

DOBLE INMERSIÓN
---------------
Los p-valores de rank_genes_groups sobre clusters de Leiden calculados con los
mismos datos no son válidos como evidencia: los clusters se han definido para
maximizar precisamente las diferencias que después se testean. Sirven para
ORDENAR genes, no para afirmar que un grupo es significativamente distinto. Si
la existencia de una población es el resultado, hace falta evidencia externa:
marcadores conocidos, otra modalidad, otro dataset, o validación experimental.
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
# 1. subconjunto y reprocesado
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
    Extrae una subpoblación y recalcula genes variables, PCA, harmony, vecinos
    y UMAP sobre ella. NO renormaliza.

    min_cells
        Filtro de genes por número de CÉLULAS que lo expresan, sobre la capa de
        counts. Tu versión anterior usaba filter_genes(min_counts=20) sobre X,
        que ahora está log-normalizada: ahí "counts" son sumas de valores log y
        el umbral no significa nada comparable entre subpoblaciones.

    integrate
        True -> los vecinos se construyen sobre X_harmony. False -> sobre X_pca.
        Ambos embeddings se calculan y se guardan igualmente, junto con los dos
        grafos ('neighbors_pca' / 'neighbors_harmony'), para poder comparar.

        Con una sola muestra en el subconjunto, harmony no tiene nada que
        corregir y se salta automáticamente.
    """
    if type(value) == str:
        mask = adata.obs[key] == value
    else:
        mask = adata.obs[key].isin(value)


    n = int(mask.sum())
    if n < 50:
        raise ValueError(f"Solo {n} células con {key}=={value!r}; muy pocas para reprocesar")

    sub = adata[mask].copy()

    # --- filtrado de genes sobre COUNTS, no sobre la matriz normalizada
    if counts_layer in sub.layers:
        n_cells_expr = np.asarray((sub.layers[counts_layer] > 0).sum(axis=0)).ravel()
        src = counts_layer
    else:
        n_cells_expr = np.asarray((sub.X > 0).sum(axis=0)).ravel()
        src = "X"
        warnings.warn(
            f"No hay capa '{counts_layer}'; filtro genes por células con valor "
            f"distinto de 0 en X. Válido, pero deja constancia."
        )
    keep = n_cells_expr >= min_cells
    if keep.sum() < 200:
        raise ValueError(f"Solo {int(keep.sum())} genes pasan min_cells={min_cells}")
    sub = sub[:, keep].copy()

    if verbose:
        print(f"[sub] {key}=={value!r}: {n} células, {int(keep.sum())} genes "
              f"(min_cells={min_cells} sobre '{src}')")

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
    Recalcula genes variables, PCA, harmony, vecinos y UMAP SOBRE EL OBJETO QUE
    SE LE PASA. No renormaliza y no toca X.

    IN PLACE de verdad: modifica `adata` y devuelve el mismo objeto, así que da
    igual si asignas el resultado o no. Es lo que la hace segura de usar como
    `reproc(a, integrate=True)` en una celda, donde descartar el retorno de una
    función que devuelve copia es un error silencioso: no falla, simplemente no
    hace nada y sigues trabajando con el objeto viejo.

    NO FILTRA GENES, precisamente por eso: quitar columnas obliga a construir
    un objeto nuevo y rompería el contrato. Si quieres el filtro por número de
    células, usa subset_and_reprocess(), que sí devuelve un objeto nuevo.

    Dos rondas, como en subset_and_reprocess: primero un andamio con todos los
    genes para tener un grafo, después triku sobre ese grafo y el embedding
    definitivo. Los genes que separan tenocitos de FAPs no son los que separan
    tenocitos entre sí, así que la selección hay que rehacerla en cada
    subconjunto.

    integrate
        True  -> los vecinos se construyen sobre X_harmony
        False -> sobre X_pca
        Los dos embeddings y los dos grafos ('neighbors_pca' / 'neighbors_harmony')
        se calculan y guardan igualmente, para poder comparar. Con una sola
        categoría de batch_key harmony no tiene nada que corregir y se salta
        automáticamente.
    """
    from pyfuncs.processing import build_embeddings, select_hvgs_triku

    n_batches = adata.obs[batch_key].nunique() if batch_key in adata.obs else 1
    do_harmony = integrate and n_batches > 1
    if integrate and not do_harmony and verbose:
        print(f"[reproc] una sola categoría de '{batch_key}': harmony no "
              f"aplica, los vecinos van sobre X_pca")

    if n_neighbors is None:
        n_neighbors = max(10, int(0.5 * adata.n_obs ** 0.5))


    neighbors_on = "harmony" if do_harmony else "pca"

    # --- ronda 1: andamio con todos los genes
    adata.var["hvg_all"] = True
    build_embeddings(adata, "hvg_all", batch_key=batch_key, n_comps=n_comps,
                     n_neighbors=n_neighbors, metric=metric,
                     neighbors_on=neighbors_on, skip_harmony=not do_harmony,
                     seed=seed, verbose=False)

    # --- ronda 2: triku sobre ese grafo, y embedding definitivo
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
# 2. revisar el diccionario de marcadores ANTES de usarlo
# --------------------------------------------------------------------------

def check_marker_dict(
    adata: ad.AnnData,
    markers: Mapping[str, Sequence[str]],
    ambient_key: str = "ambient_expression",
    ambient_top: int = 50,
    verbose: bool = True,
):
    """
    Comprueba un diccionario de marcadores contra el objeto, antes de anotar.

    Tres problemas que pasan desapercibidos porque nadie avisa:

    1. GENES QUE NO EXISTEN. Con anotación NCBI RefSeq muchos símbolos de
       Ensembl/MGI no están, o se llaman distinto. assign_cats simplemente los
       ignora, así que una población puede estar puntuándose con la mitad de
       sus marcadores y su score baja sin que sepas por qué.
    2. DUPLICADOS dentro de una lista: ese gen pesa doble en el score.
    3. MARCADORES QUE SON SOUP. Si un "marcador" está entre los genes de mayor
       ambient, el cluster que lo exprese puede estar simplemente más
       contaminado.

    Devuelve (resumen, detalle_ausentes).
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
            "poblacion": pop,
            "n_marcadores": len(genes),
            "n_unicos": len(set(genes)),
            "presentes": len(set(presentes)),
            "ausentes": len(set(ausentes)),
            "pct_presentes": round(100 * len(set(presentes)) / max(len(set(genes)), 1), 1),
            "duplicados": ", ".join(dups),
            "en_soup": ", ".join(en_soup),
        })
        if ausentes:
            faltan[pop] = sorted(set(ausentes))

    resumen = pd.DataFrame(rows).set_index("poblacion")

    if verbose:
        print(resumen.to_string())
        if faltan:
            print("\ngenes NO encontrados en var_names:")
            for pop, gs in faltan.items():
                print(f"  {pop}: {gs}")
            print("\n[markers] assign_cats los ignora en silencio. Comprueba si "
                  "es un problema de nomenclatura (RefSeq vs Ensembl) antes de "
                  "aceptar los scores.")
        flojas = resumen.index[resumen["pct_presentes"] < 70].tolist()
        if flojas:
            print(f"\n[markers] AVISO: {flojas} tienen menos del 70% de sus "
                  f"marcadores presentes; su score no es comparable al de las demás.")
        consoup = resumen.index[resumen["en_soup"] != ""].tolist()
        if consoup:
            print(f"[markers] {consoup} incluyen genes del ambient. Un cluster "
                  f"puede puntuar alto ahí solo por estar más contaminado.")
    return resumen, faltan


# --------------------------------------------------------------------------
# 3. composición por muestra / condición
# --------------------------------------------------------------------------

def population_composition(
    adata: ad.AnnData,
    population_key: str,
    group_key: str,
    normalize: str = "index",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Proporción de cada población por muestra o condición.

    Esto SÍ se lee sobre las etiquetas, no sobre el embedding, así que es
    independiente de si integraste o no. Es la vía correcta para comparar
    composición entre condiciones aunque hayas anotado sobre X_harmony.

    Ojo con la interpretación: en un timecourse la proporción es composicional
    (suma 1), así que si una población crece TODAS las demás bajan en
    proporción aunque su número absoluto no cambie. Mira también los conteos.
    """
    cont = pd.crosstab(adata.obs[group_key], adata.obs[population_key])
    prop = pd.crosstab(adata.obs[group_key], adata.obs[population_key],
                       normalize=normalize)
    if verbose:
        print("conteos:")
        print(cont.to_string())
        print("\nproporciones:")
        print((100 * prop).round(2).to_string())
        print("\n[comp] Son proporciones: suman 100 por fila. Un aumento "
              "relativo puede venir de que OTRA población haya bajado.")
    return prop