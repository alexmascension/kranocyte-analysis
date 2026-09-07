"""
QC por muestra, previo a la concatenación.

Filosofía: CALCULAR -> INSPECCIONAR -> DECIDIR -> APLICAR, en ese orden y en
pasos separados. Ninguna función de este módulo filtra células por su cuenta
ni elige umbrales por ti. Los flags se añaden como columnas booleanas y el
filtrado real se hace al final, cuando ya has visto todas las muestras.

    for gsm in samples:
        a = ad.read_h5ad(...)
        compute_qc_metrics(a, mt_contig=MT_CONTIG)   # métricas, sin flags
        add_droplet_qc(a, **params[gsm])             # DropletQC, sin flags
        flag_doublets(a)                             # por muestra, nunca concatenado
        qc_embedding(a)                              # UMAP desechable para mirar
        adatas[gsm] = a

    plot_qc_overview(adatas[gsm])                    # decides umbrales mirando
    nf_band_report(adatas[gsm], band=0.07)           # ¿qué hay en la banda plana?
    ambient_top_genes(adatas[gsm])                   # ¿de qué está hecho el soup?

    apply_qc_flags(adatas[gsm], **thresholds[gsm])
    qc_summary(adatas)                               # cuánto pierde cada muestra

AVISO SOBRE SINGLE-CELL vs SINGLE-NUCLEI
----------------------------------------
Varios umbrales de aquí significan cosas OPUESTAS según la preparación:

    métrica          single-cell                  single-nuclei
    ---------------------------------------------------------------------
    nuclear_frac     alta = célula dañada         alta = normal (es un núcleo)
                     baja = debris/soup           -> la rama de "damaged"
                                                     de DropletQC no aplica
    pct_counts_mt    alto = célula dañada         alto = contaminación
                                                     citoplasmática, núcleo
                                                     mal aislado
    Malat1           bajo = debris sin núcleo     alto siempre (es nuclear)

Este módulo NO decide por ti: te obliga a pasar los umbrales explícitamente y
te recuerda el problema por pantalla. La decisión es tuya, muestra a muestra.

NOTA SOBRE LA FRACCIÓN NUCLEAR
------------------------------
En single-cell es el ÚNICO criterio de este módulo que separa célula real de
debris anucleado, porque es el único que mide algo que el debris no puede
tener. Ni cell_probability de CellBender ni los counts totales lo hacen:
CellBender clasifica "por encima del fondo ambiente", no "es una célula", y un
fragmento de citoplasma con miles de UMIs le parece célula con p ~ 1.
"""

from __future__ import annotations

import os
import warnings
from typing import Callable, Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


PREP_WARNING = (
    "[qc] Recuerda: nuclear_frac, pct_counts_mt y Malat1 se interpretan al "
    "revés en single-nuclei que en single-cell. Revisa los umbrales para ESTA "
    "preparación antes de aplicarlos (ver docstring del módulo)."
)

# Contigs mitocondriales habituales (accesión RefSeq del genoma mitocondrial)
MT_CONTIG_MOUSE_REFSEQ = "NC_005089.1"
MT_CONTIG_HUMAN_REFSEQ = "NC_012920.1"


# --------------------------------------------------------------------------
# 1. métricas (no marcan nada)
# --------------------------------------------------------------------------

def compute_qc_metrics(
    adata: ad.AnnData,
    mt_contig: Optional[str] = None,
    mt_prefix: str = "mt-",
    malat_gene: str = "Malat1",
    counts_layer_for_umi: str = "counts_raw",
    verbose: bool = True,
) -> ad.AnnData:
    """
    Añade a obs todas las métricas de QC. No crea ningún flag.

    mt_contig
        Accesión del contig mitocondrial, p.ej. 'NC_005089.1' en ratón RefSeq.
        Es la vía preferente: selecciona por CONTIG usando var['seqname'], que
        viene de annotate_var_from_gtf(). Independiente de la nomenclatura.
    mt_prefix
        Fallback por nombre cuando no hay 'seqname' o no se da mt_contig. OJO:
        el prefijo 'mt-' es convención de Ensembl/GENCODE. En NCBI RefSeq los
        mitocondriales se llaman ND1..ND6, CYTB, COX1..3, ATP6/8, TrnX, Rnr1/2
        y este prefijo NO coge ninguno.
    counts_layer_for_umi
        Capa usada para el eje de UMIs de DropletQC. Por defecto 'counts_raw':
        los umbrales del paper de DropletQC (umi_rescue=3000 y similares) están
        definidos sobre counts crudos, no descontaminados.
    """
    # ---- mitocondriales: por contig si se puede, por nombre si no
    if mt_contig is not None and "seqname" in adata.var.columns:
        adata.var["mt"] = (adata.var["seqname"] == mt_contig).to_numpy()
        mt_how = f"contig {mt_contig}"
    else:
        if mt_contig is not None:
            warnings.warn(
                f"Se ha pedido mt_contig='{mt_contig}' pero no hay var['seqname']. "
                f"¿Corriste annotate_var_from_gtf()? Uso el prefijo '{mt_prefix}'."
            )
        elif "seqname" in adata.var.columns:
            warnings.warn(
                f"Hay var['seqname'] pero no has dado mt_contig: uso el prefijo "
                f"'{mt_prefix}', que falla con anotaciones NCBI RefSeq. "
                f"Considera mt_contig='{MT_CONTIG_MOUSE_REFSEQ}' (ratón)."
            )
        adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
        mt_how = f"prefijo '{mt_prefix}'"

    n_mt = int(adata.var["mt"].sum())
    if n_mt == 0:
        warnings.warn(
            f"0 genes mitocondriales seleccionados por {mt_how}. pct_counts_mt "
            f"saldrá 0 para todo y el criterio quedará inerte."
        )
    elif n_mt > 100:
        warnings.warn(f"{n_mt} genes marcados como mitocondriales; parecen demasiados")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True
    )

    # counts por capa, para diagnóstico y plots
    for layer, name in (("counts_raw", "n_counts_raw"),
                        ("counts_cellbender", "n_counts_cellbender")):
        if layer in adata.layers:
            adata.obs[name] = np.asarray(adata.layers[layer].sum(axis=1)).ravel()

    if "n_counts_raw" in adata.obs and "n_counts_cellbender" in adata.obs:
        raw = adata.obs["n_counts_raw"].to_numpy()
        cb = adata.obs["n_counts_cellbender"].to_numpy()
        with np.errstate(invalid="ignore", divide="ignore"):
            adata.obs["cellbender_removed_frac"] = np.where(raw > 0, 1 - cb / raw, np.nan)

    # ---- velocyto: fracción nuclear
    if all(k in adata.layers for k in ("spliced", "unspliced")):
        for layer in ("spliced", "unspliced", "ambiguous"):
            if layer in adata.layers:
                v = np.asarray(adata.layers[layer].sum(axis=1)).ravel()
                adata.obs[f"n_{layer}"] = v
                adata.obs[f"n_{layer}_log1p"] = np.log1p(v)

        s = adata.obs["n_spliced"].to_numpy()
        u = adata.obs["n_unspliced"].to_numpy()
        denom = s + u
        with np.errstate(invalid="ignore", divide="ignore"):
            adata.obs["nuclear_frac"] = np.where(denom > 0, u / denom, np.nan)

        # cobertura de Velocyto sobre GeneFull: si cae mucho en una subpoblación,
        # sospecha de alineamiento antes que de biología
        if "n_counts_raw" in adata.obs:
            velo = s + u + adata.obs.get("n_ambiguous", 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                adata.obs["velocyto_coverage"] = np.where(
                    adata.obs["n_counts_raw"] > 0, velo / adata.obs["n_counts_raw"], np.nan
                )
    else:
        warnings.warn("No hay layers spliced/unspliced: sin nuclear_frac ni DropletQC")

    # ---- UMIs para el eje de DropletQC
    if counts_layer_for_umi in adata.layers:
        umi = np.asarray(adata.layers[counts_layer_for_umi].sum(axis=1)).ravel()
        src = counts_layer_for_umi
    else:
        umi = np.asarray(adata.X.sum(axis=1)).ravel()
        src = "X"
        warnings.warn(
            f"No existe la capa '{counts_layer_for_umi}'; uso X para el eje de "
            f"UMIs. Si X es CellBender, los umbrales de DropletQC no significan "
            f"lo mismo que en el paper."
        )
    adata.obs["umi_dropletqc"] = umi
    adata.obs["log10_umi"] = np.log10(umi + 1.0)
    adata.uns.setdefault("qc", {})["umi_source"] = src

    # ---- Malat1
    if malat_gene in adata.var_names:
        x = adata[:, malat_gene].X
        x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
        adata.obs["malat1_log1p"] = np.log1p(x)
    else:
        warnings.warn(f"'{malat_gene}' no está en var_names; sin métrica de Malat1")

    if verbose:
        print(PREP_WARNING)
        print(f"[qc] {adata.n_obs} droplets | UMIs de '{src}' | "
              f"{n_mt} genes mt por {mt_how}")
    return adata


# --------------------------------------------------------------------------
# 2. DropletQC (no marca nada: solo guarda la clasificación)
# --------------------------------------------------------------------------

def add_droplet_qc(
    adata: ad.AnnData,
    *,
    nf_rescue: float,
    umi_rescue: float,
    nf_damaged_threshold: float,
    method: str = "kmeans",
    prep: str = "cell",
    classifier: Optional[Callable] = None,
    random_state: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Ejecuta classify_empty_and_damaged sobre (nuclear_frac, umi) y guarda el
    resultado en obs['droplet_class']. NO marca células para eliminar.

    Todos los umbrales son obligatorios y por muestra: no hay defaults.

    umi_rescue
        Rescata como célula lo que tenga nf baja pero muchos UMIs. Cuidado:
        el debris de mayor tamaño (p.ej. fragmentos de sarcoplasma en músculo)
        tiene MÁS UMIs que una célula real, así que esta regla puede rescatar
        justo lo peor. Usa nf_band_report() para ver qué está rescatando antes
        de fiarte del valor por defecto del paper.
    prep
        'cell' o 'nucleus'. NO cambia el cálculo — solo el aviso que imprime.
    """
    if prep not in ("cell", "nucleus"):
        raise ValueError("prep debe ser 'cell' o 'nucleus'")
    if "nuclear_frac" not in adata.obs:
        raise ValueError("Falta nuclear_frac; ejecuta compute_qc_metrics() antes")

    if classifier is None:
        from pyfuncs.dropletQC import classify_empty_and_damaged as classifier

    nf_df = pd.DataFrame({
        "nuclear_RNA_frac": adata.obs["nuclear_frac"].to_numpy(),
        "umi": adata.obs["umi_dropletqc"].to_numpy(),
        "logumi": adata.obs["log10_umi"].to_numpy(),
    }, index=adata.obs_names)

    classifier(
        nf_df, method_cells=method, nf_rescue=nf_rescue, umi_rescue=umi_rescue,
        nf_damaged_threshold=nf_damaged_threshold, random_state=random_state,
    )

    col = f"classification_{method}"
    if col not in nf_df.columns:
        raise ValueError(
            f"El clasificador no ha creado '{col}'. Columnas: {list(nf_df.columns)}"
        )

    adata.obs["droplet_class"] = pd.Categorical(nf_df[col].values)
    adata.uns.setdefault("qc", {})["dropletqc"] = {
        "method": method, "nf_rescue": nf_rescue, "umi_rescue": umi_rescue,
        "nf_damaged_threshold": nf_damaged_threshold, "prep": prep,
    }

    if verbose:
        print(f"[dropletqc] {dict(adata.obs['droplet_class'].value_counts())}")
        if prep == "nucleus":
            print(
                "[dropletqc] PREPARACIÓN = NUCLEUS: la clase 'damaged_cell' se "
                "ha calculado igualmente pero NO es interpretable aquí — en "
                "single-nuclei la fracción nuclear es alta para todo. Decide a "
                "mano si la ignoras (usa pct_counts_mt alto como proxy de "
                "núcleo mal aislado) antes de excluir nada."
            )
    return adata


def nf_band_report(
    adata: ad.AnnData,
    band: float,
    cols: Optional[Sequence[str]] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Qué hay dentro de la banda plana de nuclear_frac, y qué está rescatando
    umi_rescue.

    Crea obs['nf_group'] con tres grupos:
        celula      fuera de la banda, clasificada como célula
        debris      dentro de la banda, clasificada como no-célula
        rescatado   dentro de la banda pero clasificada como CÉLULA
                    (esto es lo que rescata umi_rescue)
        otro        fuera de la banda pero clasificada como no-célula

    Cómo leer la tabla de medianas: si 'rescatado' es idéntico a 'debris' en
    todo salvo en UMIs, umi_rescue está partiendo una población homogénea por
    la única variable que no discrimina -> súbelo o desactívalo. Si 'rescatado'
    se parece a 'celula' en Malat1 y complejidad, es un tipo celular real con
    poco contenido intrónico y el rescate está haciendo su trabajo.
    """
    for req in ("nuclear_frac", "droplet_class"):
        if req not in adata.obs:
            raise ValueError(f"Falta obs['{req}']")

    inband = adata.obs["nuclear_frac"] < band
    iscell = adata.obs["droplet_class"] == "cell"

    if verbose:
        print(pd.crosstab(inband, adata.obs["droplet_class"],
                          rownames=[f"nuclear_frac < {band}"],
                          colnames=["droplet_class"]))

    adata.obs["nf_group"] = pd.Categorical(np.select(
        [inband & iscell, inband & ~iscell, ~inband & iscell],
        ["rescatado", "debris", "celula"], default="otro"))

    if cols is None:
        cols = ["umi_dropletqc", "nuclear_frac", "n_unspliced_log1p",
                "malat1_log1p", "log1p_n_genes_by_counts", "pct_counts_mt",
                "velocyto_coverage", "cell_probability"]
    cols = [c for c in cols if c in adata.obs.columns]

    out = adata.obs.groupby("nf_group", observed=True)[cols].median().round(3)
    out.insert(0, "n", adata.obs["nf_group"].value_counts())

    if verbose:
        print()
        print(out.to_string())
        n_resc = int((adata.obs["nf_group"] == "rescatado").sum())
        if n_resc:
            print(f"\n[nf_band] {n_resc} droplets rescatados por umi_rescue. "
                  f"Compara su fila con 'debris' antes de aceptarlos.")
    return out


def ambient_top_genes(adata: ad.AnnData, n: int = 30, verbose: bool = True) -> pd.DataFrame:
    """
    De qué está hecho el soup, según CellBender.

    Es un dato sobre el TEJIDO, no sobre la muestra: en músculo esperas
    transcritos de fibra (Ckm, Acta1, troponinas, glucolíticos) porque las
    fibras son sincitios enormes que se destrozan al disociar. Los genes que
    salgan aquí arriba son sospechosos en dos sitios:
      - anotación: un cluster con esa firma puede ser soup, no biología
      - velocity: el soup es 100% spliced, así que su U/S está sesgado a la
        baja y dará velocidad negativa espuria -> sácalos de velocity_genes
    """
    if "ambient_expression" not in adata.var.columns:
        raise ValueError("No hay var['ambient_expression']; ¿cargaste la salida de CellBender?")

    cols = ["ambient_expression"]
    for extra in ("ambient_fraction", "cellbender_analyzed", "gene_biotype"):
        if extra in adata.var.columns:
            cols.append(extra)

    out = adata.var[cols].sort_values("ambient_expression", ascending=False).head(n)
    if verbose:
        print(out.to_string())
    return out


# --------------------------------------------------------------------------
# 3. doublets (por muestra, nunca sobre el objeto concatenado)
# --------------------------------------------------------------------------

def flag_doublets(
    adata: ad.AnnData,
    use_doubletdetection: bool = True,
    use_scrublet: bool = True,
    seed: int = 0,
    dd_n_iters: int = 10,
    dd_p_thresh: float = 1e-16,
    dd_voter_thresh: float = 0.5,
    dd_n_jobs: int = 1,
    expected_doublet_rate: float = 0.05,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Scores de doublets. Se ejecuta sobre UNA muestra.

    Un doublet se forma cuando dos células entran en la misma gota, en la misma
    carrera de 10x: células de muestras distintas no pueden formar uno. Correr
    esto sobre el objeto concatenado simula híbridos que no existen.

    Ninguno de los dos detectores necesita clusters externos: DoubletDetection
    clusteriza por dentro en cada iteración, y Scrublet construye su propio
    grafo KNN. El 'batch_key' de scrublet es para MUESTRAS; pasarle clusters lo
    deja casi ciego, porque los doublets que importan mezclan tipos celulares.

    Consejo: córrelo DESPUÉS de haber mirado el debris. Un dataset con 40% de
    fragmentos anucleados le da a los detectores un espacio de expresión que no
    es el de las células, y los scores salen menos informativos.
    """
    if use_doubletdetection:
        import doubletdetection
        # dd_n_jobs=1 por defecto: DoubletDetection clusteriza por dentro
        # (louvain/phenograph) y con varios procesos el resultado puede variar
        # entre ejecuciones aunque pases random_state. Si el número de doublets
        # cambia de una corrida a otra, cambian las células que pasan el QC y
        # con ellas TODO lo de aguas abajo, UMAP incluido. Súbelo solo si has
        # comprobado que tu versión es reproducible en paralelo.
        np.random.seed(seed)
        clf = doubletdetection.BoostClassifier(
            n_iters=dd_n_iters, clustering_algorithm="louvain",
            standard_scaling=True, pseudocount=0.1, n_jobs=dd_n_jobs,
            random_state=seed,
        )
        calls = clf.fit(adata.X).predict(p_thresh=dd_p_thresh, voter_thresh=dd_voter_thresh)
        score = clf.doublet_score()
        adata.obs["dd_score"] = np.asarray(getattr(score, "data", score)).ravel()
        adata.obs["dd_doublet"] = np.asarray(calls).ravel() == 1
        adata.uns.setdefault("qc", {})["doubletdetection"] = {
            "n_iters": dd_n_iters, "p_thresh": dd_p_thresh,
            "voter_thresh": dd_voter_thresh, "n_jobs": dd_n_jobs, "seed": seed,
        }
        if verbose:
            print(f"[doublets] DoubletDetection: {int(adata.obs['dd_doublet'].sum())} "
                  f"de {adata.n_obs}")

    if use_scrublet:
        import scanpy.external as sce
        res = sce.pp.scrublet(              # sin batch_key: ya estamos por muestra
            adata.copy(), expected_doublet_rate=expected_doublet_rate,
            random_state=seed, knn_dist_metric="cosine", log_transform=True, copy=True,
        )
        adata.obs["scrublet_score"] = res.obs["doublet_score"].values
        adata.obs["scrublet_doublet"] = res.obs["predicted_doublet"].values
        adata.uns.setdefault("qc", {})["scrublet"] = dict(
            res.uns.get("scrublet", {}).get("parameters", {})
        )
        del res
        if verbose:
            print(f"[doublets] Scrublet: {int(adata.obs['scrublet_doublet'].sum())} "
                  f"de {adata.n_obs}")
    return adata


# --------------------------------------------------------------------------
# 4. embedding desechable para inspeccionar
# --------------------------------------------------------------------------

def qc_embedding(
    adata: ad.AnnData,
    resolution: float = 0.5,
    n_comps: int = 15,
    n_neighbors: int = 15,
    seed: int = 0,
) -> ad.AnnData:
    """
    Clustering y UMAP provisionales, SOLO para mirar dónde caen los flags.

    Se guardan como 'qc_leiden' y 'X_umap_qc' a propósito: no como 'X_umap'.
    Así no se cuelan en las figuras finales cuando concatenes, que es un error
    silencioso y clásico. Todo esto se recalcula tras la concatenación.
    """
    tmp = adata.copy()
    if "counts_cellbender" in tmp.layers:
        tmp.X = tmp.layers["counts_cellbender"].copy()
    sc.pp.normalize_total(tmp)
    sc.pp.log1p(tmp)
    sc.pp.pca(tmp, n_comps=min(n_comps, min(tmp.shape) - 1), random_state=seed)
    sc.pp.neighbors(tmp, n_neighbors=n_neighbors, random_state=seed)
    sc.tl.leiden(tmp, resolution=resolution, key_added="qc_leiden")
    sc.tl.umap(tmp, min_dist=0.1, random_state=seed)

    adata.obs["qc_leiden"] = tmp.obs["qc_leiden"].values
    adata.obsm["X_umap_qc"] = tmp.obsm["X_umap"]
    del tmp
    return adata



# --------------------------------------------------------------------------
# 5. aplicar umbrales (obligatorios y explícitos)
# --------------------------------------------------------------------------

def apply_qc_flags(
    adata: ad.AnnData,
    min_nuclear_frac: float | None = None,
    max_nuclear_frac: Optional[float] | None = None,
    min_pct_mt: Optional[float] | None = None,
    max_pct_mt: Optional[float] | None = None,
    min_log1p_genes: Optional[float] | None = None,
    max_log1p_genes: Optional[float] | None = None,
    min_log1p_counts: Optional[float] | None = None,
    max_log1p_counts: Optional[float] | None = None,
    min_log10_umi: Optional[float] | None = None,
    max_log10_umi: Optional[float] | None = None,
    min_log1p_malat1: Optional[float] | None = None,
    max_log1p_malat1: Optional[float] | None = None,
    max_cellbender_removed: Optional[float] | None = None,
    min_n_unspliced_log1p: Optional[float] | None = None,
    max_n_unspliced_log1p: Optional[float] | None = None,
    min_n_spliced_log1p: Optional[float] | None = None,
    max_n_spliced_log1p: Optional[float] | None = None,
    droplet_exclude: Sequence[str] = (),
    doublet_cols: Sequence[str] = (),
    verbose: bool = True,
) -> ad.AnnData:
    """
    Crea los flags booleanos y la columna 'qc_pass'. NO elimina células.

    Todos los umbrales son keyword-only y obligatorios: pásalos como None si
    quieres desactivar un criterio, pero tienes que escribirlo. No hay defaults
    para que ningún umbral se aplique sin que lo hayas decidido para ESTA
    muestra y ESTA preparación.

    min_nuclear_frac
        En single-cell es el criterio principal: separa célula real de debris
        anucleado, y es el único que lo hace. Un corte directo aquí es más
        transparente y reproducible que depender de las reglas internas de
        DropletQC (que pueden rescatar debris grande vía umi_rescue).
        En single-nuclei NO tiene el mismo sentido: todo es nuclear.
    max_nuclear_frac
        Células dañadas en single-cell (perdieron citoplasma). En single-nuclei
        déjalo en None: marcaría toda la muestra.
    min_counts
        Sobre 'umi_dropletqc' (counts crudos), no sobre X.
    max_cellbender_removed
        Fracción de counts que CellBender le quitó a la célula. Una célula a la
        que se le quita el 40% es mayoritariamente soup.
    droplet_exclude
        Clases de DropletQC a marcar, p.ej. ('empty_droplet',). Vacío para no
        usar el criterio (razonable si ya cortas por min_nuclear_frac).
    doublet_cols
        Columnas booleanas de doublets a incluir, p.ej. ('dd_doublet',).
    """
    flags = {}

    def _need(col):
        if col not in adata.obs:
            raise ValueError(f"No hay obs['{col}']; ¿corriste compute_qc_metrics()?")
        return adata.obs[col].to_numpy()

    if min_nuclear_frac is not None:
        flags["flag_low_nf"] = _need("nuclear_frac") < min_nuclear_frac
    if max_nuclear_frac is not None:
        flags["flag_high_nf"] = _need("nuclear_frac") > max_nuclear_frac
    if max_pct_mt is not None:
        flags["flag_high_mt"] = _need("pct_counts_mt") > max_pct_mt
    if min_pct_mt is not None:
        flags["flag_low_mt"] = _need("pct_counts_mt") < min_pct_mt
    if min_log1p_genes is not None:
        flags["flag_low_genes"] = _need("log1p_n_genes_by_counts") < min_log1p_genes
    if max_log1p_genes is not None:
        flags["flag_high_genes"] = _need("log1p_n_genes_by_counts") > max_log1p_genes
    if min_log1p_counts is not None:
        flags["flag_low_counts"] = _need("log1p_total_counts") < min_log1p_counts
    if max_log1p_counts is not None:
        flags["flag_high_counts"] = _need("log1p_total_counts") > max_log1p_counts
    if min_log10_umi is not None:
        flags["flag_low_umi"] = _need("log10_umi") < min_log10_umi
    if max_log10_umi is not None:
        flags["flag_high_umi"] = _need("log10_umi") > max_log10_umi
    if min_log1p_malat1 is not None:
        flags["flag_low_malat1"] = _need("malat1_log1p") < min_log1p_malat1
    if max_log1p_malat1 is not None:
        flags["flag_high_malat1"] = _need("malat1_log1p") > max_log1p_malat1

    if min_n_unspliced_log1p is not None:
        flags["flag_min_n_unspliced_log1p"] = _need("log10_umi") < min_n_unspliced_log1p
    if max_n_unspliced_log1p is not None:
        flags["flag_max_n_unspliced_log1p"] = _need("log10_umi") > max_n_unspliced_log1p
    if min_n_spliced_log1p is not None:
        flags["flag_min_n_spliced_log1p"] = _need("malat1_log1p") < min_n_spliced_log1p
    if max_n_spliced_log1p is not None:
        flags["flag_max_n_spliced_log1p"] = _need("malat1_log1p") > max_n_spliced_log1p

    if max_cellbender_removed is not None:
        flags["flag_high_ambient"] = _need("cellbender_removed_frac") > max_cellbender_removed
    if len(droplet_exclude):
        if "droplet_class" not in adata.obs:
            raise ValueError("No hay 'droplet_class'; ¿corriste add_droplet_qc?")
        flags["flag_droplet"] = adata.obs["droplet_class"].isin(droplet_exclude).to_numpy()
    for col in doublet_cols:
        if col not in adata.obs:
            raise ValueError(f"No hay '{col}' en obs")
        flags[f"flag_{col}"] = adata.obs[col].to_numpy().astype(bool)

    checks = [
        ("flag_low_nf",       "nuclear_frac",            min_nuclear_frac,   "lt"),
        ("flag_high_nf",      "nuclear_frac",            max_nuclear_frac,   "gt"),
        ("flag_low_mt",       "pct_counts_mt",           min_pct_mt,         "lt"),
        ("flag_high_mt",      "pct_counts_mt",           max_pct_mt,         "gt"),
        ("flag_low_genes",    "log1p_n_genes_by_counts", min_log1p_genes,    "lt"),
        ("flag_high_genes",   "log1p_n_genes_by_counts", max_log1p_genes,    "gt"),
        ("flag_low_counts",   "log1p_total_counts",      min_log1p_counts,   "lt"),
        ("flag_high_counts",  "log1p_total_counts",      max_log1p_counts,   "gt"),
        ("flag_low_umi",      "log10_umi",               min_log10_umi,      "lt"),
        ("flag_high_umi",     "log10_umi",               max_log10_umi,      "gt"),
        ("flag_low_malat1",   "malat1_log1p",            min_log1p_malat1,   "lt"),
        ("flag_high_malat1",  "malat1_log1p",            max_log1p_malat1,   "gt"),

        ("flag_min_n_unspliced_log1p",   "n_unspliced_log1p",            min_n_unspliced_log1p,   "lt"),
        ("flag_max_n_unspliced_log1p",  "n_unspliced_log1p",            max_n_unspliced_log1p,   "gt"),
        ("flag_min_n_spliced_log1p",   "n_spliced_log1p",            min_n_spliced_log1p,   "lt"),
        ("flag_max_n_spliced_log1p",  "n_spliced_log1p",            max_n_spliced_log1p,   "gt"),

        ("flag_high_ambient", "cellbender_removed_frac", max_cellbender_removed, "gt"),
    ]
    for name, col, thr_val, op in checks:
        if thr_val is None:
            continue
        v = _need(col)
        flags[name] = v < thr_val if op == "lt" else v > thr_val


    if not flags:
        raise ValueError("No has activado ningún criterio")

    fail = np.zeros(adata.n_obs, dtype=bool)
    for k, v in flags.items():
        v = np.asarray(v)
        v = np.where(np.isnan(v.astype(float)), False, v) if v.dtype.kind == "f" else v
        adata.obs[k] = v.astype(bool)
        fail |= adata.obs[k].to_numpy()
    adata.obs["qc_pass"] = ~fail

    adata.uns.setdefault("qc", {})["thresholds"] = {
        "min_nuclear_frac": min_nuclear_frac, "max_nuclear_frac": max_nuclear_frac,
        "min_pct_mt": min_pct_mt, "max_pct_mt": max_pct_mt, 
        "min_log1p_genes": min_log1p_genes, "max_log1p_genes": max_log1p_genes,
        "min_log1p_counts": min_log1p_counts, "max_log1p_counts": max_log1p_counts, 
        "min_log1p_malat1": min_log1p_malat1, "max_log1p_malat1": max_log1p_malat1,
        "max_cellbender_removed": max_cellbender_removed,
        "droplet_exclude": list(droplet_exclude), "doublet_cols": list(doublet_cols),
    }


    if verbose:
        print(f"[qc] {int(adata.obs['qc_pass'].sum())}/{adata.n_obs} pasan "
              f"({100 * adata.obs['qc_pass'].mean():.1f}%)")
    return adata



# --------------------------------------------------------------------------
# 6. resumen entre muestras
# --------------------------------------------------------------------------

def qc_summary(
    adatas: Mapping[str, ad.AnnData],
    spread_warning: float = 15.0,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Tabla por muestra: cuántas células marca cada criterio, cuántas marca EN
    EXCLUSIVA, y cuántas sobreviven.

    La columna '_solo' es la que interesa: si un criterio marca 2000 células
    pero solo 30 en exclusiva, es redundante con los demás y su umbral casi da
    igual. Si marca 2000 y 1800 en exclusiva, ese umbral está decidiendo tu
    dataset él solo.

    spread_warning
        Puntos porcentuales de diferencia en pct_pass entre la mejor y la peor
        muestra a partir de los cuales avisa. En un timecourse, una retención
        muy desigual significa pérdida sesgada por condición: o la muestra es
        mala, o el umbral está midiendo biología en vez de calidad.
    """
    rows = []
    for name, a in adatas.items():
        flag_cols = [c for c in a.obs.columns if c.startswith("flag_")]
        row = {"sample": name, "n_cells": a.n_obs}
        for c in flag_cols:
            v = a.obs[c].to_numpy().astype(bool)
            others = np.zeros(a.n_obs, dtype=bool)
            for o in flag_cols:
                if o != c:
                    others |= a.obs[o].to_numpy().astype(bool)
            short = c.replace("flag_", "")
            row[short] = int(v.sum())
            row[short + "_solo"] = int((v & ~others).sum())
        if "qc_pass" in a.obs:
            row["n_pass"] = int(a.obs["qc_pass"].sum())
            row["pct_pass"] = round(100 * float(a.obs["qc_pass"].mean()), 1)
        rows.append(row)

    df = pd.DataFrame(rows).set_index("sample")

    if verbose and "pct_pass" in df.columns and len(df) > 1:
        spread = float(df["pct_pass"].max() - df["pct_pass"].min())
        print(f"[qc_summary] retención entre {df['pct_pass'].min()}% "
              f"({df['pct_pass'].idxmin()}) y {df['pct_pass'].max()}% "
              f"({df['pct_pass'].idxmax()}) — spread {spread:.1f} puntos")
        if spread > spread_warning:
            print(
                f"[qc_summary] AVISO: spread > {spread_warning} puntos. Si las "
                f"muestras son condiciones distintas, revisa si el umbral está "
                f"midiendo biología. Mira la columna '_solo' para saber cuál."
            )
    return df


# --------------------------------------------------------------------------
# 7. figura de inspección
# --------------------------------------------------------------------------

def plot_qc_overview(adata: ad.AnnData, title: str = "", figsize=(13, 12)):
    """
    Los nueve paneles para elegir umbrales en una muestra. Dibuja las líneas
    de los umbrales que ya estén en uns['qc']['thresholds'].

    No modifica adata: todo se calcula sobre una copia de obs.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    thr = adata.uns.get("qc", {}).get("thresholds", {})
    obs = adata.obs.copy()                      # <- sin efectos secundarios
    hue = obs["droplet_class"] if "droplet_class" in obs else None

    fig, axs = plt.subplots(3, 3, figsize=figsize, constrained_layout=True)

    def hline(ax, key, **kw):
        if thr.get(key) is not None:
            ax.axhline(thr[key], c="r", ls="--", **kw)

    # --- fila 0: mitocondrial y distribuciones
    sns.scatterplot(x=obs["log1p_total_counts"], y=obs["pct_counts_mt"],
                    hue=hue, s=3, alpha=.4, ax=axs[0, 0], legend=False)
    hline(axs[0, 0], "max_pct_mt")
    axs[0, 0].set(title="mitocondrial")

    sns.violinplot(obs, x="droplet_class", y="pct_counts_mt", ax=axs[0, 1])
    hline(axs[0, 1], "max_pct_mt")
    axs[0, 1].set(title="mitocondrial por clase")

    sns.violinplot(obs, x="droplet_class", y="log1p_n_genes_by_counts", ax=axs[0, 2])
    hline(axs[0, 2], "min_log1p_genes")          # <- era max_pct_mt (bug)
    axs[0, 2].set(title="complejidad por clase")

    # --- fila 1: DropletQC, complejidad, Malat1
    sns.scatterplot(x=obs["log10_umi"], y=obs["nuclear_frac"], hue=hue,
                    s=3, alpha=.4, ax=axs[1, 0], legend="brief")
    hline(axs[1, 0], "min_nuclear_frac")
    hline(axs[1, 0], "max_nuclear_frac")
    axs[1, 0].set(xlabel="log10 UMI", ylabel="nuclear fraction", title="DropletQC")

    sns.scatterplot(x=obs["log1p_total_counts"], y=obs["log1p_n_genes_by_counts"],
                    hue=hue, s=3, alpha=.4, ax=axs[1, 1], legend=False)
    hline(axs[1, 1], "min_log1p_genes")
    axs[1, 1].set(title="complejidad")

    if "malat1_log1p" in obs:
        sns.scatterplot(x=obs["nuclear_frac"], y=obs["malat1_log1p"], hue=hue,
                        s=3, alpha=.4, ax=axs[1, 2], legend=False)
        hline(axs[1, 2], "min_malat1")
        axs[1, 2].set(title="Malat1")

    # --- fila 2: ambient, doublets, spliced/unspliced
    if "cellbender_removed_frac" in obs:
        sns.kdeplot(obs, x="cellbender_removed_frac", hue=hue, ax=axs[2, 0])
        if thr.get("max_cellbender_removed") is not None:
            axs[2, 0].axvline(thr["max_cellbender_removed"], c="r", ls="--")
        axs[2, 0].set(title="fracción eliminada por CellBender")

    score = next((c for c in ("scrublet_score", "dd_score") if c in obs), None)
    if score:
        sns.scatterplot(x=obs["log1p_n_genes_by_counts"], y=obs[score], hue=hue,
                        s=3, alpha=.4, ax=axs[2, 1], legend=False)
        axs[2, 1].set(title=score)

    if "n_spliced_log1p" in obs and "n_unspliced_log1p" in obs:
        sns.scatterplot(x=obs["n_spliced_log1p"], y=obs["n_unspliced_log1p"],
                        hue=hue, s=3, alpha=.4, ax=axs[2, 2], legend=False)
        axs[2, 2].set(title="spliced vs unspliced")

    fig.suptitle(title or adata.uns.get("qc", {}).get("sample", ""))
    return fig, axs


# --------------------------------------------------------------------------
# 8. filtrado de debris y ensamblaje
# --------------------------------------------------------------------------

def filter_debris(
    adata: ad.AnnData,
    min_nuclear_frac: float,
    sample: Optional[str] = None,
    report_cols: Optional[Sequence[str]] = None,
    verbose: bool = True,
) -> ad.AnnData:
    """
    ÚNICO filtro duro del módulo, y a propósito: el debris anucleado no aporta
    nada aguas abajo y distorsiona todo lo que venga después (detección de
    doublets, HVGs, PCA, vecinos). Todo lo demás sigue siendo flags.

    Nunca elimina en silencio: imprime cuánto se va y el perfil medio de lo
    eliminado, para que puedas comprobar que es lo que crees que es. Deja
    constancia en uns['qc']['debris_filter'].

    En single-nuclei este filtro NO aplica tal cual: la fracción nuclear es
    alta para todo y un mínimo bajo no elimina nada. Ahí el proxy es otro
    (counts bajos, mitocondrial alto).
    """
    if "nuclear_frac" not in adata.obs:
        raise ValueError("Falta nuclear_frac; ejecuta compute_qc_metrics() antes")

    nf = adata.obs["nuclear_frac"].to_numpy()
    keep = np.where(np.isnan(nf), False, nf >= min_nuclear_frac)
    n_before, n_keep = adata.n_obs, int(keep.sum())
    n_nan = int(np.isnan(nf).sum())

    if n_keep == 0:
        raise ValueError(
            f"min_nuclear_frac={min_nuclear_frac} elimina TODAS las células. "
            f"¿Es una muestra de single-nuclei, o el umbral está invertido?"
        )

    if report_cols is None:
        report_cols = ["umi_dropletqc", "nuclear_frac", "malat1_log1p",
                       "log1p_n_genes_by_counts", "pct_counts_mt",
                       "n_unspliced_log1p", "cell_probability"]
    report_cols = [c for c in report_cols if c in adata.obs.columns]

    if verbose:
        tag = f"{sample}: " if sample else ""
        pct = 100 * n_keep / n_before
        print(f"[filter] {tag}{n_before} -> {n_keep} droplets ({pct:.1f}% retenidos); "
              f"{n_before - n_keep} eliminados por nuclear_frac < {min_nuclear_frac}"
              + (f" ({n_nan} con nuclear_frac indefinida)" if n_nan else ""))
        prof = pd.DataFrame({
            "eliminado": adata.obs.loc[~keep, report_cols].median(),
            "conservado": adata.obs.loc[keep, report_cols].median(),
        }).round(3)
        print(prof.to_string())
        if pct < 40:
            print(f"[filter] AVISO: se va más del 60% de los droplets. Comprueba "
                  f"que el perfil de lo eliminado es realmente debris.")

    out = adata[keep].copy()
    out.uns.setdefault("qc", {})["debris_filter"] = {
        "min_nuclear_frac": min_nuclear_frac, "n_before": n_before,
        "n_after": n_keep, "n_removed": n_before - n_keep,
        "n_nan_nuclear_frac": n_nan,
    }
    return out


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


def compare_ambient_profiles(
    profiles: pd.DataFrame,
    method: str = "spearman",
    top: int = 10,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Correlación del perfil de soup entre muestras, y los genes que lo dominan.

    Si todas las muestras correlacionan alto, el soup es el mismo tejido
    destrozado en distinta cantidad. Si alguna se descuelga, esa muestra tiene
    una composición de contaminación distinta — y conviene saberlo antes de
    interpretar diferencias de expresión entre condiciones.
    """
    corr = profiles.corr(method=method)
    if verbose:
        print(f"correlación ({method}) del perfil ambiente entre muestras:")
        print(corr.round(3).to_string())
        print(f"\ntop {top} genes del soup por muestra:")
        print(pd.DataFrame({
            s: profiles[s].sort_values(ascending=False).head(top).index
            for s in profiles.columns
        }).to_string(index=False))
        off = corr.where(~np.eye(len(corr), dtype=bool))
        if off.min().min() < 0.5:
            # la descolgada es la de menor correlación MEDIA con el resto, no la
            # del par más bajo: la correlación es simétrica y el par acusa a dos
            worst = off.mean().idxmin()
            print(f"\n[ambient] AVISO: '{worst}' es la muestra con el soup más "
                  f"distinto (correlación media {off.mean().min():.3f} con el resto). "
                  f"Su contaminación no es del mismo tipo que la de las demás.")
    return corr


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


def resolve_cellbender_h5(
    cb_dir: str,
    fpr: Optional[str] = None,
    filtered: bool = False,
    stem: str = "adata_raw_cellbender",
    strict: bool = True,
) -> str:
    """
    Localiza la salida de CellBender independientemente de cómo se nombró.

    CellBender cambia el patrón de fichero según cuántos --fpr se le pasen:
        --fpr 0.01         ->  <stem>.h5              y  <stem>_filtered.h5
        --fpr 0.01 0.05    ->  <stem>_FPR_0.01.h5     y  <stem>_FPR_0.01_filtered.h5

    Así que muestras corridas en momentos distintos acaban con nombres
    distintos, y una guarda tipo os.path.exists() falla en silencio contra la
    mezcla. Esta función acepta ambos.

    OJO con '_filtered': NO significa "células limpias", significa
    cell_probability > 0.5. En tejidos con mucho debris anucleado esa
    probabilidad es ~1 también para los fragmentos, así que el fichero
    filtrado los trae igual. El filtro real es filter_debris().
    """
    import glob

    suffix = "_filtered" if filtered else ""
    cands = []
    if fpr is not None:
        cands.append(f"{cb_dir}/{stem}_FPR_{fpr}{suffix}.h5")
    cands.append(f"{cb_dir}/{stem}{suffix}.h5")

    for c in cands:
        if os.path.exists(c):
            return c

    # último recurso: cualquier FPR presente, para poder avisar de cuál hay
    found = sorted(glob.glob(f"{cb_dir}/{stem}_FPR_*{suffix}.h5"))
    if not strict and found:
        warnings.warn(f"No encuentro FPR={fpr}; uso {os.path.basename(found[0])}")
        return found[0]

    raise FileNotFoundError(
        f"No encuentro la salida de CellBender en {cb_dir}.\n"
        f"  probado: {[os.path.basename(c) for c in cands]}\n"
        f"  presentes: {[os.path.basename(f) for f in sorted(glob.glob(f'{cb_dir}/{stem}*.h5'))]}"
    )


# --------------------------------------------------------------------------
# 9. elegir el FPR de CellBender con un criterio medible
# --------------------------------------------------------------------------

def fpr_tradeoff(
    adatas_fpr: Mapping[str, ad.AnnData],
    group_key: str,
    negative: Mapping[str, Sequence[str]],
    positive: Mapping[str, Sequence[str]],
    baseline: Optional[str] = None,
    layer: str = "counts_cellbender",
    verbose: bool = True,
):
    """
    Las dos curvas que hacen falta para elegir FPR sin caer en la circularidad.

    No existe un FPR "correcto" recuperable de los datos: es una decisión sobre
    qué error prefieres. Y las métricas internas al uso (silhouette,
    especificidad de marcadores, "qué limpios quedan los clusters") premian la
    homogeneidad, así que mejoran de forma casi monótona cuanto más
    descontamines — te empujan a sobre-corregir. Son circulares.

    Lo único con verdad conocida son los genes que una población NO PUEDE
    expresar: Ckm en un macrófago es soup al 100%, el valor correcto es cero.
    Eso da dos magnitudes medibles y opuestas:

        especificidad ganada  cuánto baja el gen imposible donde no debe estar
        señal perdida         cuánto baja ESE MISMO gen donde sí debe estar

    El criterio es el balance entre ambas, nunca una sola. Si al subir el FPR
    el gen cae un 90% en la población negativa y un 5% en la positiva, la
    corrección funciona. Si cae un 60% y un 40%, no has ganado especificidad:
    has aplanado el eje entero.

    Parameters
    ----------
    adatas_fpr
        {fpr: AnnData}, mismas células y genes en todos. p.ej.
        {"0.0": a0, "0.01": a1, "0.05": a5}
    group_key
        Columna de obs con las poblaciones. Basta una anotación gruesa: no
        necesitas la definitiva para saber qué es un macrófago.
    negative
        {población: [genes que NO puede expresar]}. El valor correcto es 0.
    positive
        {población: [genes que SÍ debe expresar]}. Usa los MISMOS genes que en
        'negative': comparar el mismo gen en los dos sitios es lo que hace que
        el balance signifique algo.
    baseline
        FPR de referencia. Por defecto el más bajo, idealmente '0.0'. Sin una
        referencia sin corregir no puedes medir "cuánto ha bajado".

    Returns
    -------
    (detalle, resumen)
        detalle : una fila por (fpr, tipo, población, gen)
        resumen : una fila por fpr, con el balance agregado
    """
    fprs = sorted(adatas_fpr, key=float)
    if baseline is None:
        baseline = fprs[0]
    if baseline not in adatas_fpr:
        raise ValueError(f"baseline '{baseline}' no está en adatas_fpr")
    if float(baseline) > 0 and verbose:
        print(f"[fpr] AVISO: la referencia es {baseline}, no 0.0. Las caídas se "
              f"miden respecto a una matriz ya corregida, así que subestiman "
              f"tanto la especificidad ganada como la señal perdida.")

    ref = adatas_fpr[baseline]
    for f, a in adatas_fpr.items():
        if not a.obs_names.equals(ref.obs_names):
            raise ValueError(f"FPR {f}: las células no coinciden con la referencia")
        if not a.var_names.equals(ref.var_names):
            raise ValueError(f"FPR {f}: los genes no coinciden con la referencia")
    if group_key not in ref.obs:
        raise ValueError(f"No hay obs['{group_key}'] en la referencia")

    groups = ref.obs[group_key].astype(str).to_numpy()   # mismas células: vale para todos
    gene_pos = {g: i for i, g in enumerate(ref.var_names)}

    rows = []
    for f in fprs:
        a = adatas_fpr[f]
        X = a.layers[layer] if layer and layer in a.layers else a.X
        X = sp.csc_matrix(X)
        tot = np.asarray(X.sum(axis=1)).ravel()
        tot[tot == 0] = np.nan                       # evita dividir por cero

        for kind, spec in (("negativo", negative), ("positivo", positive)):
            for pop, genes in spec.items():
                mask = groups == pop
                if not mask.any():
                    warnings.warn(f"'{pop}' no aparece en obs['{group_key}']")
                    continue
                for g in genes:
                    if g not in gene_pos:
                        warnings.warn(f"'{g}' no está en var_names")
                        continue
                    col = np.asarray(X[:, gene_pos[g]].todense()).ravel()[mask]
                    cpm = np.nanmean(col / tot[mask] * 1e6)
                    rows.append({
                        "fpr": f, "tipo": kind, "poblacion": pop, "gen": g,
                        "n_celulas": int(mask.sum()),
                        "pct_expresando": round(100 * float((col > 0).mean()), 2),
                        "mean_cpm": round(float(cpm), 2),
                        "mean_counts": round(float(col.mean()), 3),
                    })

    detalle = pd.DataFrame(rows)
    if detalle.empty:
        raise ValueError("Ninguna combinación población/gen ha producido datos")

    # caída relativa a la referencia, por (tipo, población, gen)
    idx = pd.MultiIndex.from_frame(detalle[["tipo", "poblacion", "gen"]])
    for col_src, col_dst in (("mean_cpm", "frac_eliminada"),
                             ("mean_counts", "frac_eliminada_counts")):
        base = (detalle[detalle["fpr"] == baseline]
                .set_index(["tipo", "poblacion", "gen"])[col_src])
        ref_vals = base.reindex(idx).to_numpy()
        with np.errstate(invalid="ignore", divide="ignore"):
            detalle[col_dst] = np.where(
                ref_vals > 0, 1 - detalle[col_src].to_numpy() / ref_vals, np.nan)
        detalle[col_dst] = detalle[col_dst].round(3)

    # resumen por fpr
    res = []
    for f in fprs:
        d = detalle[detalle["fpr"] == f]
        neg, pos = d[d["tipo"] == "negativo"], d[d["tipo"] == "positivo"]
        gain = float(neg["frac_eliminada"].mean())
        loss = float(pos["frac_eliminada"].mean())
        loss_counts = float(pos["frac_eliminada_counts"].mean())
        res.append({
            "fpr": f,
            "n_neg": len(neg), "n_pos": len(pos),
            "especificidad_ganada": round(gain, 3),
            "senal_perdida": round(loss, 3),
            "balance": round(gain - loss, 3),
            # en counts crudos: la normalización por célula puede enmascarar
            # una eliminación uniforme, así que esta columna es el control
            "senal_perdida_counts": round(loss_counts, 3),
            "especif_ganada_counts": round(float(neg["frac_eliminada_counts"].mean()), 3),
            "neg_pct_expresando": round(float(neg["pct_expresando"].mean()), 2),
            "neg_cpm_residual": round(float(neg["mean_cpm"].mean()), 2),
        })
    resumen = pd.DataFrame(res).set_index("fpr")

    if verbose:
        print(resumen.to_string())
        print()
        cand = resumen.drop(index=baseline, errors="ignore")
        if len(cand):
            best = cand["balance"].idxmax()
            print(f"[fpr] mejor balance: {best} "
                  f"(especificidad +{cand.loc[best, 'especificidad_ganada']:.2f}, "
                  f"señal -{cand.loc[best, 'senal_perdida']:.2f})")
            if cand["balance"].max() <= 0:
                print("[fpr] Ningún FPR gana más especificidad de la señal que "
                      "pierde: quédate con el más conservador.")
            elif cand.loc[best, "senal_perdida_counts"] > 0.25:
                print(f"[fpr] AVISO: en counts crudos se pierde el "
                      f"{100*cand.loc[best,'senal_perdida_counts']:.0f}% de la señal "
                      f"real de la población positiva. En CPM apenas se nota "
                      f"({100*cand.loc[best,'senal_perdida']:.0f}%) porque la "
                      f"normalización por célula compensa una eliminación "
                      f"uniforme. Para velocity y cinética cuantitativa, que "
                      f"trabajan sobre counts, esa pérdida SÍ importa.")
        print("\n[fpr] Recuerda: esto ordena los FPR, no dicta el valor. Si tu "
              "conclusión biológica no cambia entre ellos, la elección no "
              "importa y se dice. Si cambia, el resultado depende de un "
              "parámetro de molestia y eso es en sí mismo el hallazgo.")
    return detalle, resumen


def plot_fpr_tradeoff(resumen: pd.DataFrame, figsize=(10, 4)):
    """Las dos curvas y el balance, en función del FPR."""
    import matplotlib.pyplot as plt

    x = [float(f) for f in resumen.index]
    fig, axs = plt.subplots(1, 2, figsize=figsize, constrained_layout=True)

    axs[0].plot(x, resumen["especificidad_ganada"], "o-", label="especificidad ganada")
    axs[0].plot(x, resumen["senal_perdida"], "s-", label="señal perdida")
    axs[0].set(xlabel="FPR", ylabel="fracción eliminada vs referencia",
               title="las dos curvas")
    axs[0].legend(); axs[0].grid(alpha=.3)

    axs[1].plot(x, resumen["balance"], "o-", c="k")
    axs[1].axhline(0, c="r", ls="--", lw=1)
    axs[1].set(xlabel="FPR", ylabel="ganada - perdida", title="balance")
    axs[1].grid(alpha=.3)
    return fig, axs


# --------------------------------------------------------------------------
# 10. umbrales por muestra a partir de números, no de figuras
# --------------------------------------------------------------------------
#
# plot_qc_overview sirve para ver la ESTRUCTURA de una muestra. No sirve para
# fijar umbrales de varias muestras a la vez: las densidades están escaladas
# por panel, los ejes cambian entre figuras, y lo que parece "un poco más a la
# izquierda" pueden ser 0.05 o 0.4 de nuclear_frac. Estas funciones sacan los
# mismos datos en números para que el umbral sea reproducible y quede escrito.


def _kde_modes(x: np.ndarray, n_grid: int = 512, bw: str = "scott"):
    """
    Modas y valles de una distribución 1D por KDE.

    Devuelve (grid, dens, modas, valles) con modas/valles ordenados por
    posición. Se usa para responder dos preguntas concretas:
      - ¿hay dos nubes de nuclear_frac o una? (bimodalidad)
      - ¿dónde está el mínimo entre ellas? (ahí va el umbral, no en un número
        redondo elegido a ojo)
    """
    from scipy.stats import gaussian_kde

    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 50 or np.allclose(x, x[0]):
        return None, None, np.array([]), np.array([])

    lo, hi = np.nanmin(x), np.nanmax(x)
    pad = 0.02 * (hi - lo)
    grid = np.linspace(lo - pad, hi + pad, n_grid)
    try:
        dens = gaussian_kde(x, bw_method=bw)(grid)
    except Exception:
        return None, None, np.array([]), np.array([])

    up = np.r_[False, dens[1:] > dens[:-1]]
    down = np.r_[dens[:-1] > dens[1:], False]
    modas = grid[up & down]
    alturas = dens[up & down]
    # descarta modas que no llegan al 10% de la principal: ruido del KDE
    if modas.size:
        keep = alturas >= 0.10 * alturas.max()
        modas, alturas = modas[keep], alturas[keep]

    valles = []
    for a, b in zip(modas[:-1], modas[1:]):
        m = (grid >= a) & (grid <= b)
        valles.append(grid[m][np.argmin(dens[m])])
    return grid, dens, modas, np.asarray(valles)


def _auc(x_pos: np.ndarray, x_neg: np.ndarray) -> float:
    """
    AUC = P(un elemento de pos > uno de neg). 0.5 = indistinguibles.

    Es la medida honesta de "¿separa esta variable las dos clases?": no
    depende de la escala ni asume normalidad, y a diferencia de comparar
    medianas no se deja engañar por distribuciones bimodales solapadas.
    """
    from scipy.stats import rankdata

    x_pos = np.asarray(x_pos, float); x_pos = x_pos[np.isfinite(x_pos)]
    x_neg = np.asarray(x_neg, float); x_neg = x_neg[np.isfinite(x_neg)]
    n1, n2 = x_pos.size, x_neg.size
    if n1 == 0 or n2 == 0:
        return np.nan
    r = rankdata(np.concatenate([x_pos, x_neg]))
    return float((r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n2))


_TR_COLS = (
    "nuclear_frac", "log1p_total_counts", "log1p_n_genes_by_counts",
    "pct_counts_mt", "malat1_log1p", "cellbender_removed_frac",
    "n_spliced_log1p", "n_unspliced_log1p", "cell_probability",
)


def threshold_report(
    adatas,
    cols: Optional[Sequence[str]] = None,
    class_key: str = "droplet_class",
    debris_nf: float = 0.05,
    verbose: bool = True,
) -> dict:
    """
    Los números que hay detrás de plot_qc_overview, para una o varias muestras.

    adatas
        Un AnnData, o un dict {nombre: AnnData}. Si le pasas un objeto ya
        concatenado, sepáralo antes por muestra: los umbrales son por muestra
        y una tabla conjunta los promedia justo donde más difieren.

    Devuelve un dict con tres tablas:

    'medianas'   mediana de cada métrica por muestra y droplet_class.
    'separacion' AUC(cell vs empty_droplet) por métrica. Esta es la tabla que
                 decide si te puedes fiar de droplet_class en esa muestra:
                   ~0.5  -> la clase no discrimina en esa variable
                   >0.9  -> separa bien
                   <0.1  -> separa bien PERO al revés (las "cells" tienen
                            MENOS de esa métrica que los "empty"); si pasa en
                            complejidad o Malat1, la clasificación está
                            invertida y no debes excluir por ella.
    'modas'      modas y valles de nuclear_frac y de cellbender_removed_frac.
                 Dos modas en nuclear_frac por encima de debris_nf = dos
                 poblaciones celulares reales (p.ej. mono- vs multinucleadas,
                 o núcleos frente a células), NO debris: ahí max_nuclear_frac
                 debe quedarse en None hasta saber qué son.
    """
    if isinstance(adatas, ad.AnnData):
        adatas = {"muestra": adatas}
    cols = list(cols) if cols is not None else list(_TR_COLS)

    med_rows, sep_rows, mod_rows = [], [], []

    for name, a in adatas.items():
        obs = a.obs
        pres = [c for c in cols if c in obs.columns]

        if class_key in obs.columns:
            g = obs.groupby(class_key, observed=True)[pres].median()
            for cls, row in g.iterrows():
                d = {"muestra": name, "clase": str(cls),
                     "n": int((obs[class_key] == cls).sum())}
                d.update(row.round(3).to_dict())
                med_rows.append(d)

            es_cell = (obs[class_key] == "cell").to_numpy()
            es_empty = (obs[class_key] == "empty_droplet").to_numpy()
            if es_cell.any() and es_empty.any():
                d = {"muestra": name,
                     "n_cell": int(es_cell.sum()),
                     "n_empty": int(es_empty.sum())}
                for c in pres:
                    v = obs[c].to_numpy(dtype=float)
                    d[c] = round(_auc(v[es_cell], v[es_empty]), 3)
                sep_rows.append(d)
        else:
            d = {"muestra": name, "clase": "(sin droplet_class)", "n": a.n_obs}
            d.update(obs[pres].median().round(3).to_dict())
            med_rows.append(d)

        for c in ("nuclear_frac", "cellbender_removed_frac"):
            if c not in obs.columns:
                continue
            x = obs[c].to_numpy(dtype=float)
            if c == "nuclear_frac":
                x_sub = x[x >= debris_nf]          # sin el pico de debris
                etiqueta = f"{c} (>= {debris_nf})"
            else:
                x_sub, etiqueta = x, c
            _, _, modas, valles = _kde_modes(x_sub)
            mod_rows.append({
                "muestra": name, "variable": etiqueta,
                "n_modas": int(modas.size),
                "modas": ", ".join(f"{m:.3f}" for m in modas),
                "valles": ", ".join(f"{v:.3f}" for v in valles),
                "frac_debajo_debris_nf": (round(float((x < debris_nf).mean()), 3)
                                          if c == "nuclear_frac" else np.nan),
            })

    out = {
        "medianas": pd.DataFrame(med_rows),
        "separacion": pd.DataFrame(sep_rows),
        "modas": pd.DataFrame(mod_rows),
    }

    if verbose:
        for k, df in out.items():
            if df.empty:
                continue
            print(f"--- {k} ---")
            print(df.to_string(index=False))
            print()
        sep = out["separacion"]
        if not sep.empty:
            print("[thr] Cómo leer 'separacion': AUC de cell frente a "
                  "empty_droplet.\n"
                  "      0.5 = droplet_class no distingue nada en esa métrica.\n"
                  "      <0.5 en log1p_n_genes_by_counts o malat1_log1p = la\n"
                  "      clasificación está INVERTIDA en esa muestra; usa\n"
                  "      droplet_exclude=() y umbrales directos.")
            for _, r in sep.iterrows():
                malos = [c for c in ("malat1_log1p", "log1p_n_genes_by_counts")
                         if c in r.index and pd.notna(r[c]) and r[c] < 0.60]
                if malos:
                    print(f"      -> {r['muestra']}: {malos} con AUC "
                          f"{[r[c] for c in malos]}; no te fíes de la clase aquí.")
    return out


def suggest_thresholds(
    adata: ad.AnnData,
    class_key: str = "droplet_class",
    debris_nf: float = 0.05,
    mad_k: float = 3.0,
    verbose: bool = True,
) -> dict:
    """
    Propuesta de umbrales para apply_qc_flags a partir de la distribución.

    NO es automático en el sentido de "acéptalo y sigue". Es un punto de
    partida trazable: cada valor sale de una regla explícita sobre los datos y
    se imprime de dónde. Cámbialo si conoces la biología; pero entonces el
    cambio es una decisión tuya y queda escrita, que es distinto de haber
    puesto un número redondo porque en la figura parecía bien.

    Reglas:
      min_nuclear_frac      valle del KDE entre el pico de debris y la primera
                            nube celular; si no hay dos modas, debris_nf.
      max_nuclear_frac      None SIEMPRE. Una nuclear_frac alta no es un
                            defecto: es contenido intrónico, y separar dos
                            nubes celulares por ahí es tirar un tipo celular.
      min_log1p_counts      percentil 1 de la población por encima de
                            min_nuclear_frac.
      min_log1p_genes       ídem.
      max_pct_mt            mediana + mad_k*MAD de esa población, acotado a
                            [10, 25]. Si la distribución es bimodal con un modo
                            cerca de 100, se avisa: eso no es "mito alto", son
                            droplets que solo tienen mitocondrial.
      min_log1p_malat1      solo si Malat1 separa (AUC >= 0.75 frente a los de
                            nuclear_frac baja); si no, None.
      max_cellbender_removed  valle si hay segundo modo, si no percentil 99.
    """
    obs = adata.obs
    razones = {}

    def _q(col, q):
        return float(np.nanquantile(obs.loc[nucl, col].to_numpy(float), q))

    # --- min_nuclear_frac
    if "nuclear_frac" not in obs:
        raise ValueError("Falta nuclear_frac; ejecuta compute_qc_metrics() antes")
    x = obs["nuclear_frac"].to_numpy(float)
    _, _, modas, valles = _kde_modes(x)
    if valles.size and modas.size >= 2 and modas[0] < debris_nf * 2:
        min_nf = float(round(valles[0], 3))
        razones["min_nuclear_frac"] = (
            f"valle del KDE entre la moda de debris ({modas[0]:.3f}) y la "
            f"siguiente ({modas[1]:.3f})")
    else:
        min_nf = debris_nf
        razones["min_nuclear_frac"] = (
            f"sin bimodalidad clara en nuclear_frac; se usa debris_nf={debris_nf}")

    nucl = obs["nuclear_frac"] >= min_nf
    if nucl.sum() < 100:
        raise ValueError(f"Solo {int(nucl.sum())} droplets por encima de "
                         f"min_nuclear_frac={min_nf}; revisa el umbral")

    out = {"min_nuclear_frac": min_nf, "max_nuclear_frac": None}
    razones["max_nuclear_frac"] = ("None por norma: nuclear_frac alta es "
                                   "contenido intrónico, no un defecto")

    # --- profundidad y complejidad
    for k, col in (("min_log1p_counts", "log1p_total_counts"),
                   ("min_log1p_genes", "log1p_n_genes_by_counts")):
        if col in obs:
            out[k] = round(_q(col, 0.01), 2)
            razones[k] = f"percentil 1 de {col} en droplets con nuclear_frac >= {min_nf}"

    # --- mitocondrial
    if "pct_counts_mt" in obs:
        v = obs.loc[nucl, "pct_counts_mt"].to_numpy(float)
        v = v[np.isfinite(v)]
        med = float(np.median(v))
        mad = float(np.median(np.abs(v - med))) * 1.4826
        prop = med + mad_k * mad
        out["max_pct_mt"] = float(np.clip(round(prop, 1), 10, 25))
        razones["max_pct_mt"] = (f"mediana {med:.1f} + {mad_k}*MAD {mad:.1f} = "
                                 f"{prop:.1f}, acotado a [10, 25]")
        _, _, m_mt, _ = _kde_modes(obs["pct_counts_mt"].to_numpy(float))
        if m_mt.size >= 2 and m_mt.max() > 80:
            razones["max_pct_mt"] += (
                f"  AVISO: hay una moda en pct_counts_mt={m_mt.max():.0f}. "
                f"Eso no son células con estrés, son droplets cuyo transcriptoma "
                f"es casi solo mitocondrial: citoplasma perdido. Míralos aparte.")

    # --- Malat1
    if "malat1_log1p" in obs:
        v = obs["malat1_log1p"].to_numpy(float)
        auc = _auc(v[nucl.to_numpy()], v[~nucl.to_numpy()])
        if np.isfinite(auc) and auc >= 0.75:
            out["min_log1p_malat1"] = round(_q("malat1_log1p", 0.01), 2)
            razones["min_log1p_malat1"] = (
                f"Malat1 separa (AUC={auc:.2f} frente a nuclear_frac baja); "
                f"percentil 1 de la población retenida")
        else:
            out["min_log1p_malat1"] = None
            razones["min_log1p_malat1"] = (
                f"Malat1 NO separa aquí (AUC={auc:.2f}); filtrar por él quitaría "
                f"células reales. None.")

    # --- ambient
    if "cellbender_removed_frac" in obs:
        v = obs["cellbender_removed_frac"].to_numpy(float)
        _, _, m_cb, v_cb = _kde_modes(v)
        if v_cb.size and m_cb.size >= 2:
            out["max_cellbender_removed"] = float(round(v_cb[-1], 2))
            razones["max_cellbender_removed"] = (
                f"valle tras la moda principal; hay {m_cb.size} modas "
                f"({', '.join(f'{m:.2f}' for m in m_cb)}). Un segundo modo alto "
                f"es un grupo de droplets del que CellBender ha quitado mucho: "
                f"míralos antes de tirarlos, pueden ser un tipo celular de baja "
                f"complejidad, no basura.")
        else:
            out["max_cellbender_removed"] = float(round(np.nanquantile(v, 0.99), 2))
            razones["max_cellbender_removed"] = "percentil 99 (distribución unimodal)"

    # --- ¿se puede usar droplet_class?
    fiable = None
    if class_key in obs.columns:
        c = (obs[class_key] == "cell").to_numpy()
        e = (obs[class_key] == "empty_droplet").to_numpy()
        if c.any() and e.any():
            aucs = {}
            for col in ("nuclear_frac", "malat1_log1p", "log1p_n_genes_by_counts"):
                if col in obs:
                    v = obs[col].to_numpy(float)
                    aucs[col] = _auc(v[c], v[e])
            fiable = all(a >= 0.75 for a in aucs.values() if np.isfinite(a))
            razones["droplet_exclude"] = (
                f"AUC cell vs empty: "
                + ", ".join(f"{k}={v:.2f}" for k, v in aucs.items())
                + ("  -> la clase separa, puedes usar droplet_exclude"
                   if fiable else
                   "  -> la clase NO separa; droplet_exclude=() y filtra por "
                     "umbrales directos"))
    out["_droplet_class_fiable"] = fiable

    if verbose:
        print("umbrales propuestos:")
        for k, v in out.items():
            if k.startswith("_"):
                continue
            print(f"  {k:24s} = {v!r}")
        print("\npor qué:")
        for k, r in razones.items():
            print(f"  {k}: {r}")
        print("\n[thr] Son un punto de partida derivado de la distribución de "
              "ESTA muestra. Antes de aceptarlos, pasa apply_qc_flags y mira "
              "qué se lleva por delante cada bandera en qc_summary.")
    out["_razones"] = razones
    return out








