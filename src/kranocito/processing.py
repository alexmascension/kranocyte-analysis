"""
Selección de genes, PCA e integración, sobre el objeto ya normalizado.

Orden de las llamadas (dos rondas, como pide triku):

    HVGs preliminares (scanpy)  ->  PCA  ->  harmony  ->  neighbors  ->  triku
                                                                          |
                                     HVGs definitivos <-------------------+
                                                |
                              PCA + harmony sobre el set definitivo
                                                |
                                   X_pca                X_harmony
                             (biología temporal)         (anotar)

La primera ronda existe solo para tener un grafo de vecinos razonable: triku
puntúa genes según cómo se distribuye su expresión sobre ese grafo, así que
necesita uno construido con algo mejor que todos los genes. Los HVGs
preliminares son un andamio, no un resultado.

LAS DOS VÍAS
------------
Ambas se calculan sobre EL MISMO set de HVGs, para que sean comparables. La
diferencia es solo si se aplica harmony:

    X_pca       sin corregir. Para composición y dinámica entre condiciones.
    X_harmony   corregido. Para clusterizar y anotar tipos celulares.

Cuando muestra y condición coinciden (una réplica por timepoint), harmony
alinea los timepoints entre sí y puede fusionar estados que solo existen en
uno. harmony_merge_report() busca exactamente eso.
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
# 1. HVGs preliminares (andamio)
# --------------------------------------------------------------------------

def select_hvgs_preliminary(
    adata: ad.AnnData,
    n_top_genes: int = 3000,
    flavor: str = "seurat",
    batch_key: Optional[str] = None,
    key_added: str = "hvg_prelim",
    verbose: bool = True,
) -> ad.AnnData:
    """
    HVGs de scanpy sobre X ya normalizada y log-transformada.

    Solo sirven para construir el grafo que después usa triku. No los uses
    como set final ni los reportes como resultado.

    batch_key
        Con batch_key, scanpy selecciona por muestra y combina, lo que evita
        que un gen variable en una sola muestra domine. Útil, pero ojo si
        muestra == condición: entonces penaliza genes que varían ENTRE
        condiciones, que es justo tu señal. Con pocas muestras yo lo dejaría
        en None y confiaría en la segunda ronda.
    """
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, flavor=flavor,
        batch_key=batch_key, inplace=True,
    )
    adata.var[key_added] = adata.var["highly_variable"].to_numpy()
    adata.uns.setdefault("processing", {})["hvg_prelim"] = {
        "n_top_genes": n_top_genes, "flavor": flavor, "batch_key": batch_key,
        "n_selected": int(adata.var[key_added].sum()),
    }
    if verbose:
        print(f"[hvg] preliminares: {int(adata.var[key_added].sum())} genes "
              f"(flavor={flavor}, batch_key={batch_key})")
    return adata


# --------------------------------------------------------------------------
# 2. PCA + harmony + vecinos
# --------------------------------------------------------------------------

def build_embeddings(
    adata: ad.AnnData,
    hvg_key: str,
    batch_key: str = "gsm",
    n_comps: int = 50,
    n_neighbors: Optional[int] = None,
    metric: str = "correlation",
    pca_key: str = "X_pca",
    harmony_key: str = "X_harmony",
    neighbors_on: str = "harmony",
    skip_harmony: bool = False,
    max_iter_harmony: int = 30,
    deterministic_harmony: bool = True,
    seed: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Calcula PCA y su versión corregida con harmony sobre el mismo set de genes,
    y el grafo de vecinos sobre el que se indique.

    neighbors_on
        'harmony' o 'pca'. Determina qué representación alimenta el grafo (y
        por tanto a triku y al clustering). Se guardan AMBOS embeddings pase lo
        que pase; esto solo elige cuál se usa ahora.

    Guarda los grafos con neighbors_key propio ('neighbors_pca' /
    'neighbors_harmony') para poder clusterizar en los dos espacios sin
    pisarse, que es lo que necesita harmony_merge_report().

    deterministic_harmony
        Limita harmony a un hilo para que X_harmony sea reproducible. Ver el
        comentario largo en el cuerpo: con random_state fijo pero varios hilos
        sigue saliendo distinto en cada ejecución. Ponlo a False solo si
        prefieres velocidad y te da igual la reproducibilidad exacta.
    """
    if hvg_key not in adata.var:
        raise ValueError(f"No hay var['{hvg_key}']")
    n_hvg = int(adata.var[hvg_key].sum())
    if n_hvg < 50:
        raise ValueError(f"Solo {n_hvg} genes en '{hvg_key}'")

    adata.var["highly_variable"] = adata.var[hvg_key].to_numpy()  # scanpy lee ésta
    sc.pp.pca(adata, n_comps=min(n_comps, n_hvg - 1, adata.n_obs - 1),
              use_highly_variable=True, random_state=seed)
    if pca_key != "X_pca":
        adata.obsm[pca_key] = adata.obsm["X_pca"].copy()

    if skip_harmony:
        # sin batch que corregir (p.ej. una sola muestra en el subconjunto):
        # X_harmony = X_pca para que el resto del código no tenga que ramificar
        adata.obsm[harmony_key] = adata.obsm["X_pca"].copy()
        neighbors_on = "pca"
    else:
      # harmonypy directo, no sce.pp.harmony_integrate: el wrapper de scanpy
      # transpone Z_corr, y harmonypy >=2.0 ya lo devuelve como (células x PCs).
      # Con esa combinación el wrapper falla con un error de shape confuso.
      import harmonypy

      # DETERMINISMO. harmony inicializa con sklearn.KMeans, y sklearn.KMeans
      # NO es reproducible bit a bit cuando corre multihilo: las reducciones en
      # coma flotante de OpenMP se acumulan en orden distinto según cómo se
      # repartan los hilos. Es un comportamiento documentado de sklearn, no un
      # bug de harmonypy. Con random_state fijo la inicialización es la misma,
      # pero esas diferencias de último bit se amplifican a lo largo de las
      # iteraciones de harmony y salen X_harmony distintos en cada ejecución —
      # y con ellos vecinos y UMAP distintos.
      #
      # Se limita a un hilo SOLO durante esta llamada, con threadpoolctl (viene
      # con sklearn). Es preferible a poner OMP_NUM_THREADS=1 global: no
      # ralentiza el resto del pipeline y no obliga a tocar el entorno antes de
      # los imports.
      def _run():
          try:
              return harmonypy.run_harmony(
                  adata.obsm["X_pca"], adata.obs, [batch_key],
                  max_iter_harmony=max_iter_harmony, random_state=seed,
              )
          except TypeError:
              warnings.warn(
                  "Tu harmonypy no acepta random_state; se siembra el RNG global "
                  "de numpy como paliativo. Comprueba con dos ejecuciones que "
                  f"obsm['{harmony_key}'] sale idéntico."
              )
              np.random.seed(seed)
              return harmonypy.run_harmony(
                  adata.obsm["X_pca"], adata.obs, [batch_key],
                  max_iter_harmony=max_iter_harmony,
              )

      np.random.seed(seed)
      if deterministic_harmony:
          try:
              from threadpoolctl import threadpool_limits
          except ImportError:
              warnings.warn(
                  "Sin threadpoolctl no puedo forzar un solo hilo en harmony; "
                  "X_harmony puede variar entre ejecuciones. pip install "
                  "threadpoolctl, o exporta OMP_NUM_THREADS=1 antes de importar."
              )
              ho = _run()
          else:
              with threadpool_limits(limits=1):
                  ho = _run()
      else:
          ho = _run()
      Z = np.asarray(ho.Z_corr)
      n_pc = adata.obsm["X_pca"].shape[1]
      if Z.shape == (adata.n_obs, n_pc):
          pass
      elif Z.shape == (n_pc, adata.n_obs):
          Z = Z.T
      else:
          raise ValueError(
              f"harmonypy ha devuelto Z_corr con shape {Z.shape}; se esperaba "
              f"({adata.n_obs}, {n_pc}) o su transpuesta"
          )
      adata.obsm[harmony_key] = np.ascontiguousarray(Z, dtype=np.float32)

    if n_neighbors is None:
        n_neighbors = max(10, int(0.5 * adata.n_obs ** 0.5))

    for name, rep in (("pca", pca_key), ("harmony", harmony_key)):
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=rep,
                        metric=metric, random_state=seed,
                        key_added=f"neighbors_{name}")

    # el grafo "activo" (.obsp por defecto) es el que usarán triku y leiden
    src = f"neighbors_{neighbors_on}"
    adata.uns["neighbors"] = dict(adata.uns[src])
    adata.obsp["distances"] = adata.obsp[adata.uns[src]["distances_key"]]
    adata.obsp["connectivities"] = adata.obsp[adata.uns[src]["connectivities_key"]]
    adata.uns["neighbors"]["distances_key"] = "distances"
    adata.uns["neighbors"]["connectivities_key"] = "connectivities"

    adata.uns.setdefault("processing", {})["embeddings"] = {
        "hvg_key": hvg_key, "n_hvg": n_hvg, "batch_key": batch_key,
        "n_comps": int(adata.obsm["X_pca"].shape[1]), "n_neighbors": n_neighbors,
        "metric": metric, "neighbors_on": neighbors_on,
        "seed": seed, "deterministic_harmony": bool(deterministic_harmony),
    }
    if verbose:
        print(f"[emb] PCA sobre {n_hvg} HVGs -> {adata.obsm['X_pca'].shape[1]} PCs")
        print(f"[emb] harmony por '{batch_key}' -> obsm['{harmony_key}']")
        print(f"[emb] vecinos (k={n_neighbors}, {metric}) en 'neighbors_pca' y "
              f"'neighbors_harmony'; activo = {neighbors_on}")
    return adata


# --------------------------------------------------------------------------
# 3. HVGs definitivos con triku
# --------------------------------------------------------------------------

def select_hvgs_triku(
    adata: ad.AnnData,
    n_features: Optional[int] = None,
    key_added: str = "hvg_triku",
    use_raw: bool = False,
    random_state: int = 0,
    n_procs: Optional[int] = 1,
    verbose: bool = True,
    **kwargs,
) -> ad.AnnData:
    """
    triku sobre el grafo de vecinos activo. Segunda ronda.

    triku puntúa cada gen por cómo se concentra su expresión en el grafo de
    vecinos, no por la relación media-varianza. Por eso necesita el grafo, y
    por eso hay dos rondas.

    random_state
        NO ES OPCIONAL AUNQUE LO PAREZCA. triku compara la distribución
        observada contra un nulo que construye ALEATORIZANDO la expresión
        (apply_background_correction). Sin semilla, cada ejecución selecciona un
        conjunto de genes distinto, y a partir de ahí cambian PCA, harmony,
        vecinos, clusters y UMAP. Es la causa clásica de "el UMAP me sale
        distinto cada vez que corro el notebook".
    n_procs
        1 por defecto, también por reproducibilidad: con varios procesos el
        orden de reducción puede variar. Súbelo si prefieres velocidad a
        determinismo exacto.
    """
    try:
        import triku as tk
    except ImportError as e:
        raise ImportError("pip install triku") from e

    if "neighbors" not in adata.uns:
        raise ValueError("No hay grafo de vecinos; ejecuta build_embeddings() antes")

    tk_kwargs = dict(n_features=n_features, use_raw=use_raw, **kwargs)
    if n_procs is not None:
        tk_kwargs["n_procs"] = n_procs
    try:
        tk.tl.triku(adata, random_state=random_state, **tk_kwargs)
    except TypeError:
        warnings.warn(
            "Tu versión de triku no acepta random_state y/o n_procs. Se siembra "
            "el RNG global de numpy como paliativo, pero NO garantiza "
            "reproducibilidad: comprueba con dos ejecuciones que var['%s'] "
            "sale igual." % key_added
        )
        np.random.seed(random_state)
        tk.tl.triku(adata, n_features=n_features, use_raw=use_raw, **kwargs)

    adata.var[key_added] = adata.var["highly_variable"].to_numpy()
    adata.uns.setdefault("processing", {})["hvg_triku"] = {
        "n_features": n_features, "n_selected": int(adata.var[key_added].sum()),
        "random_state": random_state, "n_procs": n_procs,
    }
    if verbose:
        print(f"[hvg] triku: {int(adata.var[key_added].sum())} genes")
    return adata


# --------------------------------------------------------------------------
# 4. ¿se cuelan los genes del soup entre los HVGs?
# --------------------------------------------------------------------------

def soup_vs_hvg_report(
    adata: ad.AnnData,
    hvg_keys: Sequence[str] = ("hvg_prelim", "hvg_triku"),
    top_soup: int = 50,
    ambient_key: str = "ambient_expression",
    verbose: bool = True,
):
    """
    Comprueba si los genes del ambient acaban seleccionados como variables.

    La pregunta es legítima porque tanto scanpy como triku corrigen por nivel
    de expresión: si están bien calibrados, un gen muy expresado NO debería
    salir seleccionado solo por serlo. Así que esto es un test, no una excusa
    para excluirlos de antemano.

    Devuelve (detalle, resumen):
      detalle : por gen del soup, su rango de ambient y si está en cada set
      resumen : por set de HVGs, cuántos genes del soup contiene, cuántos se
                esperarían por azar, y el p-valor hipergeométrico

    Cómo leerlo: si 'observados' ~ 'esperados' y p > 0.05, los selectores están
    haciendo su trabajo y NO hace falta excluir nada. Si hay enriquecimiento
    claro, entonces sí: quítalos del set antes de PCA, no después.
    """
    from scipy.stats import hypergeom

    # el perfil de soup puede venir de var (una muestra) o de varm (concatenado)
    if ambient_key in adata.varm:
        amb = pd.Series(np.asarray(adata.varm[ambient_key]).mean(axis=1),
                        index=adata.var_names, name=ambient_key)
        src = f"varm['{ambient_key}'] (media de {adata.varm[ambient_key].shape[1]} muestras)"
    elif ambient_key in adata.var:
        amb = adata.var[ambient_key].astype(float)
        src = f"var['{ambient_key}']"
    else:
        raise ValueError(f"No encuentro '{ambient_key}' ni en var ni en varm")

    amb = amb.sort_values(ascending=False)
    soup_genes = amb.head(top_soup).index

    detalle = pd.DataFrame({ambient_key: amb.head(top_soup)})
    detalle.insert(0, "rango_soup", np.arange(1, len(soup_genes) + 1))
    for k in hvg_keys:
        if k in adata.var:
            detalle[k] = adata.var.loc[soup_genes, k].to_numpy().astype(bool)

    rows = []
    M = adata.n_vars
    for k in hvg_keys:
        if k not in adata.var:
            warnings.warn(f"No hay var['{k}']")
            continue
        sel = adata.var[k].to_numpy().astype(bool)
        N = int(sel.sum())
        obs = int(adata.var.loc[soup_genes, k].to_numpy().astype(bool).sum())
        esp = top_soup * N / M
        # P(X >= obs) bajo selección aleatoria
        p = float(hypergeom.sf(obs - 1, M, top_soup, N)) if obs > 0 else 1.0
        rows.append({
            "hvg_set": k, "n_hvg": N, "n_genes": M,
            "soup_observados": obs, "soup_esperados": round(esp, 2),
            "enriquecimiento": round(obs / esp, 2) if esp > 0 else np.nan,
            "p_hipergeom": p,
        })
    resumen = pd.DataFrame(rows).set_index("hvg_set")

    if verbose:
        print(f"perfil de soup desde {src}; top {top_soup} genes\n")
        print(resumen.to_string())
        print()
        for k in resumen.index:
            o, e = resumen.loc[k, "soup_observados"], resumen.loc[k, "soup_esperados"]
            p = resumen.loc[k, "p_hipergeom"]
            if p < 0.05 and o > e:
                print(f"[soup] {k}: ENRIQUECIDO ({o} vs {e:.1f} esperados, "
                      f"p={p:.2g}). Considera excluirlos antes del PCA.")
            else:
                print(f"[soup] {k}: sin enriquecimiento ({o} vs {e:.1f}, "
                      f"p={p:.2g}). El selector corrige bien por expresión; "
                      f"no hace falta excluir nada.")
        sel_cols = [k for k in hvg_keys if k in detalle.columns]
        if sel_cols:
            hit = detalle[detalle[sel_cols].any(axis=1)]
            if len(hit):
                print(f"\ngenes del soup seleccionados en algún set:")
                print(hit.head(25).to_string())
    return detalle, resumen


def exclude_genes_from_hvg(
    adata: ad.AnnData,
    genes: Sequence[str],
    hvg_key: str,
    key_added: Optional[str] = None,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Quita genes concretos de un set de HVGs. Úsalo SOLO si soup_vs_hvg_report()
    ha mostrado enriquecimiento: excluir a ciegas es tirar biología por si
    acaso.
    """
    key_added = key_added or hvg_key
    present = [g for g in genes if g in adata.var_names]
    missing = set(genes) - set(present)
    if missing:
        warnings.warn(f"{len(missing)} genes no están en var_names: {sorted(missing)[:5]}")

    sel = adata.var[hvg_key].to_numpy().astype(bool).copy()
    idx = adata.var_names.get_indexer(present)
    removed = int(sel[idx].sum())
    sel[idx] = False
    adata.var[key_added] = sel

    if verbose:
        print(f"[hvg] {removed} de {len(present)} genes excluidos de '{hvg_key}' "
              f"-> '{key_added}' ({int(sel.sum())} genes)")
    return adata


# --------------------------------------------------------------------------
# 5. ¿qué fusiona harmony?
# --------------------------------------------------------------------------

def harmony_merge_report(
    adata: ad.AnnData,
    condition_key: str,
    resolution: float = 1.0,
    dominance: float = 0.8,
    seed: int = 0,
    verbose: bool = True,
):
    """
    Clusteriza en los dos espacios y busca los clusters que existen sin
    corregir pero se disuelven al integrar.

    El fallo de harmony cuando muestra == condición no se ve como "la señal
    desaparece del UMAP". Se ve en la GRANULARIDAD: un estado que solo aparece
    en un timepoint forma su cluster en X_pca, y harmony lo empuja hacia el
    estado equivalente de los demás hasta que leiden los fusiona. El estado no
    queda atenuado: deja de existir como categoría, y entonces ningún análisis
    posterior puede encontrarlo.

    dominance
        Un cluster se considera "específico de condición" si esa fracción de
        sus células viene de una sola condición.

    REQUISITO: condition_key debe tener AL MENOS DOS valores. Con uno solo todos
    los clusters salen 100% puros trivialmente y la tabla es ruido. Ese caso, de
    hecho, significa que batch y condición no están confundidos y que puedes usar
    X_harmony sin reservas: el dilema de las dos vías no existe ahí.

    Devuelve (sospechosos, contingencia). 'sospechosos' lista los clusters
    específicos de condición en X_pca junto a en cuántos clusters de X_harmony
    se reparten: si se reparte en varios mixtos, harmony lo ha disuelto.
    """
    for nk in ("neighbors_pca", "neighbors_harmony"):
        if nk not in adata.uns:
            raise ValueError(f"Falta uns['{nk}']; ejecuta build_embeddings() antes")
    if condition_key not in adata.obs:
        raise ValueError(f"No hay obs['{condition_key}']")

    n_cond = adata.obs[condition_key].nunique()
    if n_cond < 2:
        raise ValueError(
            f"obs['{condition_key}'] tiene un solo valor "
            f"({adata.obs[condition_key].iloc[0]!r}). Con una sola condición TODOS "
            f"los clusters son 100% puros por definición y el informe no dice "
            f"nada.\n"
            f"  Y es buena noticia: si batch y condición no están confundidos, "
            f"harmony hace corrección de batch legítima y no hay dilema de dos "
            f"vías. Usa X_harmony sin más.\n"
            f"  Si lo que quieres es ver poblaciones exclusivas de UNA MUESTRA, "
            f"pásale condition_key='<tu columna de muestra>'."
        )
    if verbose and n_cond == adata.obs.groupby(condition_key, observed=True).ngroups:
        pass

    for name in ("pca", "harmony"):
        sc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_{name}",
                     neighbors_key=f"neighbors_{name}", flavor="igraph",
                     n_iterations=2, directed=False, random_state=seed)

    comp = pd.crosstab(adata.obs["leiden_pca"], adata.obs[condition_key],
                       normalize="index")
    especificos = comp.index[comp.max(axis=1) >= dominance].tolist()

    cont = pd.crosstab(adata.obs["leiden_pca"], adata.obs["leiden_harmony"])

    rows = []
    for c in especificos:
        fila = cont.loc[c]
        fila = fila[fila > 0].sort_values(ascending=False)
        n = int(fila.sum())
        principal = fila.index[0]
        # ¿el cluster de harmony al que va a parar es mixto?
        mezcla = comp.loc[c].idxmax()
        harm_comp = pd.crosstab(adata.obs["leiden_harmony"],
                                adata.obs[condition_key], normalize="index")
        rows.append({
            "leiden_pca": c, "n_celulas": n,
            "condicion": mezcla,
            "pureza_pca": round(float(comp.loc[c].max()), 3),
            "n_clusters_harmony": int(len(fila)),
            "harmony_principal": principal,
            "pureza_harmony_destino": round(float(harm_comp.loc[principal].max()), 3),
        })
    sospechosos = pd.DataFrame(rows)

    if verbose:
        n_pca = adata.obs["leiden_pca"].nunique()
        n_h = adata.obs["leiden_harmony"].nunique()
        print(f"[harmony] {n_pca} clusters en X_pca, {n_h} en X_harmony "
              f"(resolution={resolution})")
        if sospechosos.empty:
            print(f"[harmony] ningún cluster de X_pca supera {dominance:.0%} de "
                  f"una sola condición: no hay candidatos a estado específico.")
        else:
            print(f"[harmony] {len(sospechosos)} clusters específicos de condición "
                  f"en X_pca:\n")
            print(sospechosos.to_string(index=False))
            disueltos = sospechosos[sospechosos["pureza_harmony_destino"] < dominance]
            if len(disueltos):
                print(f"\n[harmony] AVISO: {len(disueltos)} de ellos acaban en un "
                      f"cluster de harmony MIXTO. Mira sus DEGs antes de aceptar "
                      f"la versión integrada: puede ser un artefacto de batch que "
                      f"harmony ha hecho bien en absorber, o el estado transitorio "
                      f"que buscas.")
    return sospechosos, cont