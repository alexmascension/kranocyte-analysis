"""
Normalización sobre el objeto concatenado y filtrado.

Los size factors de scran son una escala POR CÉLULA, así que calcularlos sobre
el conjunto de muestras da un pooling más estable y deja a todas en la misma
escala. Eso sí: el pooling de scran necesita un clustering previo hecho sobre
el objeto YA unido, y ese clustering tiene que existir antes de normalizar —
de ahí el paso preliminar.

Flujo:

    adata = concat_samples(adatas)
    adata = adata[adata.obs["qc_pass"]].copy()        # filtrado definitivo

    preliminary_clusters(adata)                       # para el pooling
    scran_size_factors(adata, cluster_key="scran_clusters")
    apply_size_factors(adata)                         # deja las layers
    compare_normalizations(adata)                     # mira antes de elegir

    adata.X = adata.layers["norm_scran_log1p"].copy()

QUÉ MATRIZ SE NORMALIZA
-----------------------
La capa de counts que le pases (por defecto 'counts_cellbender'). scran espera
counts, no algo ya normalizado. La salida de CellBender con el estimador MCKP
son counts enteros, así que vale. Si prefieres normalizar sobre los crudos,
pásale counts_layer='counts_raw' — pero entonces la descontaminación no llega
a downstream.

MEMORIA
-------
El paso a R es el cuello de botella: rpy2 necesita una matriz densa y
genes x células en float64. Con 40k genes y 20k células son ~6.4 GB, que es
justo lo que hace petar el kernel. Por eso se envía SOLO un subconjunto de
genes expresados. No cambia el resultado: los size factors son por célula y
scran ya descarta internamente los genes de baja expresión con min.mean.
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
    Perfiles de soup por muestra en un DataFrame genes x muestras.
 
    Hay que sacarlos ANTES de concatenar: difieren entre muestras, así que
    merge='same' los tira. Y merece la pena conservarlos: comparar el soup
    entre timepoints te dice cuánto tejido se ha destrozado en cada uno, que
    en un modelo de lesión es información, no ruido.
    """
    cols = {}
    for name, a in adatas.items():
        if key not in a.var.columns:
            warnings.warn(f"{name}: no tiene var['{key}']")
            continue
        cols[name] = pd.Series(a.var[key].to_numpy(dtype=float), index=a.var_names)
    if not cols:
        raise ValueError(f"Ninguna muestra tiene var['{key}']")
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
    Concatena preservando var y uns, que es lo que ad.concat tira por defecto.
 
    Tres cosas que hace y que a mano se olvidan:
 
    merge='same' / uns_merge='same'
        Sin ellos pierdes var entera (entrez_id, seqname, gene_biotype...) y
        uns['annotation'] con la procedencia del GTF. Los defaults de ad.concat
        son None: descartar.
    drop_obsm
        X_umap_qc de cada muestra vive en un espacio propio. Concatenar esas
        coordenadas produce un embedding sin sentido que parece válido. Fuera.
    prefix_obs
        El cluster '0' de una muestra no tiene nada que ver con el '0' de otra.
        Se prefijan con el nombre de la muestra para que no se fusionen.
    """
    names = list(adatas.keys())
    var_sets = [tuple(adatas[n].var_names) for n in names]
    if len(set(var_sets)) != 1:
        warnings.warn(
            "Las muestras no comparten el mismo conjunto de genes; con join='outer' "
            "los ausentes se rellenan con ceros. ¿Mismo GTF en todas?"
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
 
    # los perfiles de soup no sobreviven a merge='same': se guardan aparte
    try:
        out.varm["ambient_expression"] = collect_ambient_profiles(adatas).reindex(
            out.var_names).to_numpy(dtype=np.float32)
        out.uns.setdefault("qc", {})["ambient_samples"] = names
    except ValueError:
        pass
 
    if verbose:
        print(f"[concat] {out.n_obs} células x {out.n_vars} genes de {len(names)} muestras")
        print(out.obs[sample_key].value_counts().to_string())
        kept = [c for c in ("entrez_id", "seqname", "gene_biotype", "gene_symbol")
                if c in out.var.columns]
        print(f"[concat] var conservado: {kept if kept else 'NADA — revisa merge='}")
        if not kept:
            warnings.warn("No ha sobrevivido ninguna columna de anotación en var")
        print(f"[concat] uns conservado: {sorted(out.uns)}")
        if "ambient_expression" in out.varm:
            print(f"[concat] perfiles de soup en varm['ambient_expression'] "
                  f"({out.varm['ambient_expression'].shape})")
    return out
 


# --------------------------------------------------------------------------
# 1. clustering preliminar para el pooling de scran
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
    Clustering rápido y desechable, solo para que scran agrupe células
    parecidas al hacer el pooling.
 
    No es tu clustering final: se calcula sobre una normalización burda
    (CPM + log1p) precisamente porque la buena todavía no existe. Es una
    dependencia circular inevitable y está bien resuelta así — scran solo
    necesita grupos aproximadamente homogéneos, no una anotación correcta.
 
    min_cells_per_cluster
        computeSumFactors falla si un cluster tiene menos células que el mayor
        tamaño de pool. En scranPY el default es sizes=arange(21,102,5), o sea
        pools de hasta 101 células: por eso el mínimo aquí es 102, no 100.
        Los clusters más pequeños se fusionan en uno solo, y se avisa.
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
            print(f"[scran] {len(small)} clusters con <{min_cells_per_cluster} "
                  f"células fusionados en 'small_merged' ({n_small} células)")
 
        # Fusionar los pequeños entre sí no basta si solo había UNO: se queda
        # igual de pequeño y el pooling sigue degradado. En ese caso se absorbe
        # en el cluster grande más cercano en espacio PCA.
        if n_small < min_cells_per_cluster:
            big = cl.value_counts().drop("small_merged", errors="ignore")
            if big.empty:
                raise ValueError(
                    f"Todos los clusters tienen <{min_cells_per_cluster} células. "
                    f"Baja 'resolution'."
                )
            centroids = pd.DataFrame(pcs).groupby(cl.values).mean()
            target = (centroids.drop(index="small_merged")
                      .sub(centroids.loc["small_merged"])
                      .pow(2).sum(axis=1).idxmin())
            cl = cl.replace("small_merged", target)
            if verbose:
                print(f"[scran] 'small_merged' ({n_small}) seguía por debajo del "
                      f"mínimo: absorbido en el cluster '{target}' "
                      f"({int(big[target])} células)")
 
    adata.obs[key_added] = pd.Categorical(cl.values)
 
    if verbose:
        vc = adata.obs[key_added].value_counts()
        print(f"[scran] {len(vc)} clusters para el pooling; "
              f"tamaños {int(vc.min())}-{int(vc.max())}")
        if vc.min() < min_cells_per_cluster:
            warnings.warn(
                f"Sigue habiendo un cluster de {int(vc.min())} células. Baja "
                f"'resolution' o sube 'min_cells_per_cluster'."
            )
    return adata
 
 
# --------------------------------------------------------------------------
# 2. size factors con scran (vía rpy2)
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
    Deconvolución de Lun et al. (computeSumFactors), por dos vías.
 
    backend
        'scranPY'  -> sfortma2/scranPY, que implementa el algoritmo completo en
                      Python (pooling + sistema lineal + QR/CVXPY). Sin R.
        'r'        -> scran de Bioconductor vía rpy2.
 
        No confundir con 'scranpy' de BiocPy/libscran, que es otro paquete con
        nombre casi idéntico: ése NO trae computeSumFactors y normaliza por
        tamaño de librería (para eso está library_size_factors()).
 
    En ambos casos el cuello de botella es el mismo: las dos vías necesitan la
    matriz DENSA, así que el subconjunto de genes no es opcional en datasets
    grandes.
 
    min_cells_expressed / max_genes
        Solo entran al cálculo los genes expresados en al menos N células, y
        como mucho los 'max_genes' de mayor expresión media.
 
        MATIZ, porque antes te lo vendí demasiado fuerte: esto NO es
        exactamente idéntico a usar todos los genes. La deconvolución divide
        por el tamaño de librería y toma medianas de ratios, así que el
        conjunto de genes entra en el cálculo. Es una buena aproximación
        porque los genes excluidos aportan una fracción ínfima de los counts y
        los habría descartado igualmente min_mean — pero es aproximación, no
        identidad. Con max_genes=None se usan todos, si te cabe en memoria.
 
    Devuelve los size factors en obs[key_added] y valida que sean positivos:
    scran puede dar factores <= 0 en células de muy baja calidad, y eso
    produce infinitos silenciosos al dividir.
    """
    if backend not in ("scranPY", "r"):
        raise ValueError("backend debe ser 'scranPY' o 'r'")
 
    if cluster_key not in adata.obs:
        raise ValueError(
            f"No hay obs['{cluster_key}']; ejecuta preliminary_clusters() antes"
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
            f"Solo {int(keep.sum())} genes pasan el filtro; scran necesita más. "
            f"Baja min_cells_expressed."
        )
 
    if verbose:
        gb = int(keep.sum()) * adata.n_obs * 8 / 1e9
        print(f"[scran] {int(keep.sum())} genes x {adata.n_obs} células, "
              f"densa: {gb:.2f} GB")
        if gb > 8:
            warnings.warn(f"{gb:.1f} GB en denso; baja max_genes si peta")
 
    clusters = adata.obs[cluster_key].astype(str).to_numpy()
 
    # ---------------------------------------------------------------- scranPY
    if backend == "scranPY":
        try:
            from scranPY import compute_sum_factors
        except ImportError as e:
            raise ImportError(
                "pip install git+https://github.com/sfortma2/scranPY.git"
            ) from e
 
        # scranPY exige X densa y lee los clusters por NOMBRE de columna
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
            normalize_counts=False,     # NO tocar X: eso lo hace apply_size_factors
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
    mat = np.asarray(X[:, keep].todense(), dtype=np.float64).T   # R: genes x células
    try:
        import anndata2ri  # noqa: F401
        from rpy2.robjects import numpy2ri, r
        from rpy2.robjects.conversion import localconverter
        import rpy2.robjects as ro
    except ImportError as e:
        raise ImportError(
            "Hace falta rpy2 (y R con scran). "
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
            f"scran devolvió {sf.shape[0]} size factors para {adata.n_obs} células"
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
    scran puede devolver factores <= 0 en células de muy baja calidad. Dividir
    por ellos da infinitos o negativos que luego aparecen como NaN a mitad del
    pipeline, muy lejos de aquí. Se comprueba en el sitio.
    """
    bad = ~np.isfinite(sf) | (sf <= 0)
    if bad.any():
        raise ValueError(
            f"{int(bad.sum())} size factors no positivos o no finitos. Suele "
            f"significar células con contenido casi nulo que han sobrevivido al "
            f"QC: revísalas y fíltralas antes de normalizar."
        )
    q = np.quantile(sf, [0, .01, .25, .5, .75, .99, 1])
    if verbose:
        print(f"[scran] size factors  min={q[0]:.3f}  1%={q[1]:.3f}  "
              f"mediana={q[3]:.3f}  99%={q[5]:.3f}  max={q[6]:.3f}")
 
    # La cola BAJA es la peligrosa y es asimétrica: normalizar divide, así que
    # un size factor de 0.02 multiplica los counts de esa célula por 50. Esas
    # células salen con valores enormes y ruidosos y pueden formar su propio
    # cluster espurio. La cola alta solo comprime, que es benigno.
    n_low = int((sf < low_warn).sum())
    if n_low:
        warnings.warn(
            f"{n_low} células con size factor < {low_warn}: al normalizar sus "
            f"counts se amplifican hasta {1 / max(q[0], 1e-9):.0f}x. Suelen ser "
            f"células de bajo contenido que han pasado el QC. Míralas con "
            f"size_factor_report() antes de seguir."
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
    Size factors por tamaño de librería, centrados. Alternativa a scran SIN R.
 
    OJO con la confusión habitual: `scranpy` NO implementa computeSumFactors.
    Su normalización de RNA usa las sumas por columna ("If None, this defaults
    to the column sums of the count matrix"), y su único cálculo propio de
    factores, compute_clrm1_factors, es para ADT. El pooling con deconvolución
    de Lun et al. — lo que distingue a scran de una normalización por library
    size — solo está en el paquete de R.
 
    Así que esto es library size, se llame como se llame. Con use_scranpy=True
    el centrado lo hace scranpy (mode='lowest', igual que su workflow); con
    False se centra a la media, que es lo que hace scran en R. La diferencia
    es un factor global de escala, irrelevante tras el log.
 
    block_key
        Con scranpy, centra por bloque (p.ej. 'gsm') en vez de globalmente.
        Úsalo solo si sabes por qué: centrar por muestra elimina diferencias
        REALES de contenido de RNA entre muestras, que en un timecourse pueden
        ser biología.
    """
    X = adata.layers[counts_layer] if counts_layer in adata.layers else adata.X
    sf = np.asarray(sp.csr_matrix(X).sum(axis=1)).ravel().astype(np.float64)
 
    if (sf <= 0).any():
        n = int((sf <= 0).sum())
        raise ValueError(
            f"{n} células con 0 counts en '{counts_layer}'. Fíltralas antes de "
            f"normalizar: su size factor no está definido."
        )
 
    how = "media"
    if use_scranpy:
        try:
            import scranpy
            block = (adata.obs[block_key].astype(str).to_numpy()
                     if block_key is not None else None)
            sf = scranpy.center_size_factors(
                sf, block=block, mode="per-block" if block is not None else "lowest")
            how = f"scranpy ({'per-block' if block is not None else 'lowest'})"
        except ImportError:
            warnings.warn("scranpy no está instalado; centro a la media")
            sf = sf / sf.mean()
    else:
        sf = sf / sf.mean()
 
    validate_size_factors(sf, verbose=verbose)
    adata.obs[key_added] = sf
    adata.uns.setdefault("normalization", {})["library_size"] = {
        "counts_layer": counts_layer, "centering": how, "block_key": block_key,
    }
    if verbose:
        print(f"[libsize] size factors por library size, centrados por {how}")
        print("[libsize] NO es computeSumFactors: sin pooling ni deconvolución")
    return adata
 
 
def size_factor_report(
    adata: ad.AnnData,
    key: str = "size_factors",
    low: float = 0.1,
    cols: Optional[Sequence[str]] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Compara las células de size factor bajo con el resto.
 
    Un size factor pequeño no es un problema del normalizador: es una célula
    con muy poco contenido que ha sobrevivido al QC. Al dividir por él sus
    counts se amplifican, y esas células acaban dominando la varianza y
    formando clusters que no existen.
 
    Si la fila 'sf_bajo' tiene muchos menos counts y genes que 'resto', vuelve
    al QC y ponles un min_log1p_counts / min_log1p_genes. Es más limpio que
    parchear los size factors después.
    """
    if key not in adata.obs:
        raise ValueError(f"No hay obs['{key}']")
    if cols is None:
        cols = ["n_counts_raw", "n_counts_cellbender", "log1p_total_counts",
                "log1p_n_genes_by_counts", "pct_counts_mt", "nuclear_frac",
                "malat1_log1p", "cellbender_removed_frac", key]
    cols = [c for c in cols if c in adata.obs.columns]
 
    grp = np.where(adata.obs[key].to_numpy() < low, "sf_bajo", "resto")
    out = adata.obs[cols].groupby(grp, observed=True).median().round(3)
    out.insert(0, "n", pd.Series(grp).value_counts())
 
    if verbose:
        print(out.to_string())
        n_low = int((grp == "sf_bajo").sum())
        if n_low:
            print(f"\n[sf] {n_low} células ({100 * n_low / adata.n_obs:.2f}%) "
                  f"con size factor < {low}")
    return out
 
 
# --------------------------------------------------------------------------
# 3. aplicar y comparar
# --------------------------------------------------------------------------
 
def apply_size_factors(
    adata: ad.AnnData,
    counts_layer: str = "counts_cellbender",
    key: str = "size_factors",
    target_sum: float = 1e6,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Deja tres normalizaciones como layers, sin tocar X. Eliges después.
 
        norm_log1p        log1p de los counts sin escalar (control)
        norm_scran_log1p  log(1 + x / size_factor)
        norm_cpm_log1p    log1p de CPM, la referencia habitual
 
    X se queda como está a propósito: que la elección sea un paso explícito
    tuyo (`adata.X = adata.layers[...]`) y no un efecto lateral.
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
        print(f"[norm] sin obs['{key}']: no se genera la capa de scran")
 
    cpm = sc.pp.normalize_total(
        ad.AnnData(X=X.copy()), target_sum=target_sum, inplace=False)["X"]
    cpm = sp.csr_matrix(cpm, dtype=np.float32)
    cpm.data = np.log1p(cpm.data)
    adata.layers["norm_cpm_log1p"] = cpm
 
    if verbose:
        print(f"[norm] layers: {[k for k in adata.layers if k.startswith('norm_')]}")
        print("[norm] X NO se ha tocado; elige con adata.X = adata.layers['...']")
    return adata
 
 
def compare_normalizations(
    adata: ad.AnnData,
    layers: Optional[Sequence[str]] = None,
    group_key: Optional[str] = None,
    n_comps_pca: int = 20,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Compara las normalizaciones por lo único que importa aquí: si la
    profundidad de secuenciación sigue dominando después.
 
    'pc1_vs_depth' es la correlación de Spearman entre PC1 y los counts crudos.
    Ése es el modo de fallo real: si la profundidad domina, sale como primer
    componente principal y te estructura el UMAP entero. Cerca de 0 = bien.
 
    OJO con la tentación de medirlo sobre la SUMA de la fila normalizada: esa
    suma depende sobre todo de cuántos genes se detectan, que correlaciona con
    la profundidad aunque la normalización sea perfecta. Da ~1 siempre y no
    distingue nada. Por eso se mide sobre PC1.
 
    Con group_key (p.ej. 'gsm') añade la dispersión ENTRE muestras de la suma
    normalizada: si una muestra queda sistemáticamente por encima, no están en
    la misma escala y eso se confundirá con efecto de batch.
    """
    if layers is None:
        layers = [k for k in adata.layers if k.startswith("norm_")]
    if not layers:
        raise ValueError("No hay layers de normalización; ejecuta apply_size_factors()")
 
    from scipy.stats import spearmanr
 
    depth = (adata.obs["n_counts_raw"].to_numpy()
             if "n_counts_raw" in adata.obs
             else np.asarray(adata.X.sum(axis=1)).ravel())
 
    rows = []
    for lay in layers:
        s = np.asarray(adata.layers[lay].sum(axis=1)).ravel()
 
        # PC1 sobre esta normalización: el modo de fallo que importa
        tmp = ad.AnnData(X=sp.csr_matrix(adata.layers[lay], dtype=np.float32))
        sc.pp.pca(tmp, n_comps=min(n_comps_pca, min(tmp.shape) - 1),
                  random_state=0)
        pc1 = tmp.obsm["X_pca"][:, 0]
        var1 = float(tmp.uns["pca"]["variance_ratio"][0])
        del tmp
 
        row = {
            "layer": lay,
            "cv_suma": round(float(s.std() / s.mean()), 4),
            "pc1_vs_depth": round(abs(float(spearmanr(pc1, depth).statistic)), 3),
            "pc1_var_ratio": round(var1, 3),
        }
        if group_key is not None and group_key in adata.obs:
            per = pd.Series(s).groupby(adata.obs[group_key].to_numpy()).mean()
            row["spread_entre_muestras"] = round(
                float(per.max() / per.min()) if per.min() > 0 else np.nan, 3)
        rows.append(row)
 
    out = pd.DataFrame(rows).set_index("layer")
    if verbose:
        print(out.to_string())
        best = out["pc1_vs_depth"].idxmin()
        print(f"\n[norm] PC1 menos contaminado por la profundidad: {best} "
              f"(|rho| = {out.loc[best, 'pc1_vs_depth']:.3f})")
        if out.loc[best, "pc1_vs_depth"] > 0.5:
            print("[norm] AVISO: incluso la mejor deja PC1 muy ligado a la "
                  "profundidad. Revisa si quedan células de muy baja calidad.")
        print("[norm] Es un criterio necesario, no suficiente: mide que la "
              "profundidad ya no domina, no que la biología esté bien escalada.")
    return out