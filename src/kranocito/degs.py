"""
Expresión diferencial entre CONDICIONES, dentro de cada población.

    res   = deg_by_group(adata_FAP, "cell_type", "condition",
                         reference="non_DEN", test="DEN_2d")
    sets  = deg_sets(res, n_top=100, direction="up")
    plot_upset(sets)
    lfc_concordance(res)
    core_vs_specific(sets)

LA ADVERTENCIA QUE NO PUEDE FALTAR
----------------------------------
Comparar condiciones tratando las CÉLULAS como réplicas es el caso donde la
pseudorreplicación hace más daño. Con un ratón por condición, el test contesta
"¿difieren estas células?" y no "¿difieren estas condiciones?", y no hay
corrección posible: la n verdadera es 1 vs 1. Squair et al. (Nat Commun 2021)
mostraron que en este escenario los p-valores están inflados por órdenes de
magnitud.

Si tienes varios ratones por condición, la respuesta correcta es PSEUDOBULK:
sumar counts por (muestra x población) y hacer el test entre muestras. Si no
los tienes, se puede seguir — pero entonces esto ORDENA genes candidatos y no
demuestra diferencias, y así hay que escribirlo.

Y EL CONFUSOR ESPECÍFICO DE ESTE DATASET
----------------------------------------
La profundidad por célula difiere hasta 5x entre timepoints (D2 ~70.000 counts
medianos frente a ~15.000 en NI). La profundidad y la condición son la misma
variable en este diseño, y la tasa de detección de un gen depende fuertemente
de la profundidad. Resultado: genes poco expresados aparecerán "subidos" en la
condición secuenciada más hondo sin que haya cambiado nada.

Por eso deg_by_group() exige pts=True y devuelve la diferencia de detección, y
por eso existe depth_report(): míralo ANTES de interpretar la lista.
"""

from __future__ import annotations

import warnings
from typing import Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


# --------------------------------------------------------------------------
# 0. el confusor, antes que nada
# --------------------------------------------------------------------------

def depth_report(adata: ad.AnnData, condition_key: str,
                 group_key: Optional[str] = None,
                 counts_key: str = "total_counts",
                 verbose: bool = True) -> pd.DataFrame:
    """
    Profundidad y número de células por condición (y población).

    Si las dos condiciones que vas a comparar difieren mucho en counts por
    célula, cualquier lista de DEGs va a estar contaminada por eso: los genes
    poco expresados se detectan más donde hay más profundidad. Mira el ratio
    antes de mirar los genes.
    """
    cols = [condition_key] + ([group_key] if group_key else [])
    g = adata.obs.groupby(cols, observed=True)
    out = pd.DataFrame({
        "n_celulas": g.size(),
        "counts_mediana": g[counts_key].median() if counts_key in adata.obs else np.nan,
        "genes_mediana": (g["n_genes_by_counts"].median()
                          if "n_genes_by_counts" in adata.obs else np.nan),
    })
    if verbose:
        print(out.round(1).to_string())
        if group_key is None and len(out) == 2:
            r = out["counts_mediana"].max() / max(out["counts_mediana"].min(), 1)
            print(f"\n[depth] ratio de profundidad entre condiciones: {r:.1f}x")
            if r > 1.5:
                print("[depth] AVISO: por encima de 1.5x, la tasa de detección "
                      "ya sesga la lista de DEGs. Contrasta con "
                      "deg_by_group(..., downsample=True).")
    return out


# --------------------------------------------------------------------------
# 1. el test, por población
# --------------------------------------------------------------------------

def deg_by_group(
    adata: ad.AnnData,
    group_key: str,
    condition_key: str,
    reference: str,
    test: str,
    groups: Optional[Sequence[str]] = None,
    method: str = "wilcoxon",
    tie_correct: bool = True,
    min_cells: int = 20,
    downsample: bool = False,
    counts_layer: str = "counts_cellbender",
    seed: int = 0,
    verbose: bool = True,
) -> dict:
    """
    Corre el test test-vs-reference DENTRO de cada población. {población: df}.

    El signo queda fijado por construcción: `logfoldchanges` positivo = ARRIBA
    en `test` respecto a `reference`. En tu bucle esto lo decidía el argumento
    `cluster=` del volcano, y elegir el grupo equivocado invierte la lista
    entera sin que nada falle.

    downsample
        Iguala la profundidad entre las dos condiciones submuestreando counts
        antes del test, sobre `counts_layer`. No es gratis —tiras información—
        pero es la única forma de saber si tu lista sobrevive al confusor de
        profundidad. Corre las dos versiones y compara.

    min_cells
        Poblaciones con menos de esto en alguna de las dos condiciones se
        saltan: con 8 células no hay test que valga y la lista sería ruido.
    """
    obs = adata.obs
    if condition_key not in obs or group_key not in obs:
        raise ValueError(f"Faltan obs['{condition_key}'] u obs['{group_key}']")

    if groups is None:
        groups = [g for g in obs[group_key].astype(str).unique()]
    salida, resumen = {}, []

    for g in groups:
        m = (obs[group_key].astype(str) == g) & \
            (obs[condition_key].astype(str).isin([reference, test]))
        sub = adata[m.to_numpy()].copy()
        cond = sub.obs[condition_key].astype(str)
        n_ref, n_test = int((cond == reference).sum()), int((cond == test).sum())
        if min(n_ref, n_test) < min_cells:
            if verbose:
                print(f"[deg] {g}: {n_ref} vs {n_test} células -> me la salto "
                      f"(min_cells={min_cells})")
            continue

        if downsample:
            if counts_layer in sub.layers:
                sub.X = sub.layers[counts_layer].copy()
            objetivo = int(np.percentile(
                np.asarray(sub.X.sum(axis=1)).ravel(), 5))
            sc.pp.downsample_counts(sub, counts_per_cell=objetivo,
                                    random_state=seed)
            sc.pp.normalize_total(sub, target_sum=1e4)
            sc.pp.log1p(sub)
            if verbose:
                print(f"[deg] {g}: submuestreado a {objetivo} counts/célula")

        sub.obs["_cond"] = pd.Categorical(cond, categories=[reference, test])
        kw = dict(groupby="_cond", groups=[test], reference=reference,
                  method=method, pts=True, key_added="deg")
        if method == "wilcoxon":
            kw["tie_correct"] = tie_correct
        sc.tl.rank_genes_groups(sub, **kw)

        df = sc.get.rank_genes_groups_df(sub, group=test, key="deg")
        # diferencia de DETECCIÓN: la señal del confusor de profundidad
        pts = sub.uns["deg"].get("pts")
        if pts is not None:
            df = df.merge(
                pts.rename(columns={reference: "pct_ref", test: "pct_test"})
                   .reset_index().rename(columns={"index": "names"}),
                on="names", how="left")
            df["delta_pct"] = df["pct_test"] - df["pct_ref"]
        df["poblacion"] = g
        salida[g] = df
        resumen.append({"poblacion": g, f"n_{reference}": n_ref,
                        f"n_{test}": n_test, "n_genes": len(df)})
        del sub

    if verbose and resumen:
        print(pd.DataFrame(resumen).set_index("poblacion").to_string())
        print(f"\n[deg] logfoldchanges POSITIVO = arriba en '{test}' "
              f"respecto a '{reference}'.")
    return salida


def depth_bias_report(res: Mapping[str, pd.DataFrame], adata: ad.AnnData,
                      layer: str = "counts_cellbender",
                      verbose: bool = True) -> pd.DataFrame:
    """
    ¿Tiene tu lista de DEGs la firma del sesgo de profundidad?

    La firma es inconfundible y se puede simular: con DOS grupos de expresión
    IDÉNTICA pero secuenciados a 15.000 y 70.000 counts, tras normalizar y
    log1p salen ~7% de genes "significativos", el 98% de ellos AL ALZA en el
    grupo profundo, y el efecto se concentra en los genes poco detectados
    (logFC medio 0.29 y +42 puntos de detección en los de expresión baja,
    frente a 0.08 y +2 puntos en los de expresión alta).

    Esta función reproduce esa tabla sobre TUS resultados. Si ves el mismo
    patrón —asimetría de dirección y sesgo concentrado en lo poco expresado—,
    buena parte de tu lista es profundidad. Si el logFC es plano respecto al
    nivel de expresión, no lo es.
    """
    M = adata.layers[layer] if layer in adata.layers else adata.X
    media = pd.Series(np.asarray(M.mean(axis=0)).ravel(), index=adata.var_names)

    filas = []
    for pop, df in res.items():
        d = df.copy()
        d["media"] = media.reindex(d["names"]).to_numpy()
        d["nivel"] = pd.cut(d["media"], [-.01, .05, .3, 1, 5, 1e9],
                            labels=["~0", "muy bajo", "bajo", "medio", "alto"])
        for niv, sub in d.groupby("nivel", observed=True):
            sig = sub[sub["pvals_adj"] < 0.05]
            filas.append({
                "poblacion": pop, "nivel": str(niv), "n": len(sub),
                "logFC_medio": round(float(sub["logfoldchanges"].mean()), 3),
                "delta_deteccion": (round(float(sub["delta_pct"].mean()), 3)
                                    if "delta_pct" in sub else np.nan),
                "pct_sig": round(100 * len(sig) / max(len(sub), 1), 1),
                "frac_sig_al_alza": (round(float((sig["logfoldchanges"] > 0).mean()), 2)
                                     if len(sig) else np.nan),
            })
    out = pd.DataFrame(filas)
    if verbose:
        print(out.to_string(index=False))
        alza = out.groupby("poblacion")["frac_sig_al_alza"].mean()
        desbal = alza[(alza > 0.85) | (alza < 0.15)]
        if len(desbal):
            print(f"\n[bias] {list(desbal.index)} tienen >85% de los DEGs en "
                  f"UNA sola dirección. Con biología real se espera algo más "
                  f"repartido; esa asimetría es la firma del sesgo técnico.")
        print("[bias] Compara también 'logFC_medio' entre niveles: si crece al "
              "bajar la expresión, es profundidad, no biología.")
    return out


# --------------------------------------------------------------------------
# 2. conjuntos y UpSet
# --------------------------------------------------------------------------

def deg_sets(
    res: Mapping[str, pd.DataFrame],
    n_top: Optional[int] = 100,
    direction: str = "up",
    max_padj: float = 0.05,
    min_lfc: float = 0.5,
    min_delta_pct: Optional[float] = None,
    min_pct_both: Optional[float] = 0.1,
    verbose: bool = True,
) -> dict:
    """
    Conjuntos de genes por población, filtrados y ordenados. {población: set}.

    direction
        'up' o 'down'. SEPARADOS SIEMPRE. Un gen que sube en FAP.1 y baja en
        FAP.3 no es un gen "compartido": es lo contrario. Meterlos en el mismo
        UpSet convierte una divergencia en un solapamiento.

    min_pct_both
        ESTE es el filtro contra el confusor de profundidad: exige que el gen
        se detecte en al menos esta fracción de células EN LAS DOS
        condiciones. El sesgo de profundidad se concentra casi entero en los
        genes poco detectados —ahí la diferencia de profundidad se traduce en
        diferencia de ceros, y el rank-sum lee sobre todo ceros—, mientras que
        los genes bien expresados en ambos grupos apenas se ven afectados.

        (Corrijo lo que dije antes: min_delta_pct NO protege de esto. Con 5x de
        profundidad, el artefacto PRODUCE delta_pct grandes; exigirlo
        seleccionaría artefactos en vez de descartarlos.)

    min_delta_pct
        Diferencia mínima en la fracción de células que expresan el gen. Sirve
        para quedarse con cambios de presencia/ausencia en vez de solo de
        magnitud, pero NO es un control de profundidad. Úsalo junto a
        min_pct_both, nunca en su lugar.
    """
    sets = {}
    for pop, df in res.items():
        d = df[df["pvals_adj"] <= max_padj]
        d = d[d["logfoldchanges"] >= min_lfc] if direction == "up" \
            else d[d["logfoldchanges"] <= -min_lfc]
        if min_pct_both is not None and {"pct_ref", "pct_test"} <= set(d.columns):
            antes = len(d)
            d = d[(d["pct_ref"] >= min_pct_both) & (d["pct_test"] >= min_pct_both)]
            if verbose and antes:
                print(f"[sets] {pop}: {antes} -> {len(d)} tras exigir "
                      f"deteccion >= {min_pct_both} en LAS DOS condiciones")
        if min_delta_pct is not None and "delta_pct" in d:
            antes = len(d)
            d = d[d["delta_pct"].abs() >= min_delta_pct]
            if verbose and antes:
                print(f"[sets] {pop}: {antes} -> {len(d)} tras delta_pct "
                      f">= {min_delta_pct}")
        d = d.reindex(d["logfoldchanges"].abs().sort_values(ascending=False).index)
        genes = d["names"].tolist()
        sets[pop] = set(genes[:n_top] if n_top else genes)

    if verbose:
        print(f"\n[sets] dirección '{direction}': "
              + ", ".join(f"{k}={len(v)}" for k, v in sets.items()))
        tam = [len(v) for v in sets.values()]
        if tam and max(tam) > 3 * max(min(tam), 1):
            print("[sets] Los conjuntos difieren mucho de tamaño. El UpSet "
                  "premia a los grandes: fíjate en las fracciones, no solo en "
                  "las barras.")
    return sets


def plot_upset(sets: Mapping[str, set], min_subset_size: int = 1,
               sort_by: str = "cardinality", show_counts: bool = True,
               figsize=(11, 6), **kwargs):
    """
    UpSet de los conjuntos. Devuelve (fig, dict del upsetplot).

    Cómo NO leerlo: la barra de "solo FAP.3" no significa "gen exclusivo de
    FAP.3". Significa "gen que entró en el top-N de FAP.3 y no en el de los
    demás", que con listas truncadas es muy distinto — puede estar en el
    puesto N+1 de las otras cuatro. Por eso conviene mirar también
    lfc_concordance(), que no depende de ningún corte.
    """
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_contents

    datos = from_contents({k: sorted(v) for k, v in sets.items()})

    def _dibuja(counts):
        fig = plt.figure(figsize=figsize)
        up = UpSet(datos, min_subset_size=min_subset_size, sort_by=sort_by,
                   show_counts=counts, **kwargs)
        ejes = up.plot(fig=fig)
        fig.canvas.draw()      # el fallo de las etiquetas salta AL DIBUJAR,
        return fig, ejes       # no al construir: hay que forzarlo aquí


    try:
        return _dibuja(show_counts)
    except TypeError as e:
        # upsetplot 0.9 + matplotlib recientes revientan al dibujar las
        # etiquetas de conteo. No es tu figura: es la librería.
        if not show_counts:
            raise
        plt.close("all")
        warnings.warn(f"upsetplot falló con show_counts=True ({e}); repito sin "
                      f"las etiquetas. Los tamaños están en core_vs_specific().")
        return _dibuja(False)


# --------------------------------------------------------------------------
# 3. lo que no depende de un corte
# --------------------------------------------------------------------------

def lfc_concordance(
    res: Mapping[str, pd.DataFrame],
    restrict: Optional[Sequence[str]] = None,
    max_padj: Optional[float] = None,
    min_abs_lfc: float = 0.0,
    method: str = "pearson",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Correlación del logFC entre poblaciones, sobre TODOS los genes comunes.

    Es la versión del UpSet que no depende de n_top. El UpSet contesta
    "¿coinciden las listas?", que es una pregunta sobre tu umbral; esto
    contesta "¿es la misma respuesta?", que es la pregunta biológica. Si dos
    poblaciones correlacionan a 0.8 en logFC pero comparten pocos genes en el
    top-100, la respuesta es compartida y el solapamiento pequeño era un
    artefacto del corte.

    restrict
        PÁSALE LA UNIÓN DE TUS CONJUNTOS DE DEGs. Sobre los 25.000 genes, la
        correlación la dominan los miles que no cambian en ninguna población y
        cuyo logFC es ruido: sale baja aunque la respuesta sea idéntica. La
        pregunta útil es "entre los genes que cambian en ALGUNA población,
        ¿cambian igual?".

            union = set().union(*sets.values())
            lfc_concordance(res, restrict=union)

    method
        'pearson' por defecto: aquí la MAGNITUD del cambio es informativa y
        Spearman la tira. Con Spearman, un conjunto de genes que suben todos
        parecido queda dominado por el ruido de los rangos y la correlación
        sale artificialmente baja. Mira las dos si dudas.
    """
    lfc = {}
    for pop, df in res.items():
        d = df
        if max_padj is not None:
            d = d[d["pvals_adj"] <= max_padj]
        if min_abs_lfc:
            d = d[d["logfoldchanges"].abs() >= min_abs_lfc]
        lfc[pop] = d.set_index("names")["logfoldchanges"]

    M = pd.DataFrame(lfc)
    if restrict is not None:
        M = M.loc[M.index.intersection(list(restrict))]
    elif verbose:
        print("[lfc] sin 'restrict': se usan TODOS los genes, y los que no "
              "cambian en ninguna población diluyen la correlación. Pásale la "
              "unión de tus conjuntos de DEGs.")
    corr = M.corr(method=method, min_periods=50)
    if verbose:
        print(f"[lfc] {len(M)} genes; correlación {method} del logFC:")
        print(corr.round(3).to_string())
        vals = corr.where(~np.eye(len(corr), dtype=bool)).stack()
        if len(vals):
            print(f"\n[lfc] media entre pares: {vals.mean():.3f} "
                  f"(min {vals.min():.3f} entre "
                  f"{vals.idxmin()[0]} y {vals.idxmin()[1]})")
            print("[lfc] Alto = respuesta común a la lesión. Bajo = cada "
                  "población responde a lo suyo. Es lo que el UpSet insinúa "
                  "pero no mide.")
    return corr


def query_sets(
    sets: Mapping[str, set],
    include: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
    exact: bool = True,
    n_populations: Optional[int] = None,
    res: Optional[Mapping[str, pd.DataFrame]] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Saca los genes de una celda concreta del UpSet.

        query_sets(up, ["FAP.4"])                      # SOLO en FAP.4
        query_sets(up, list(up))                       # comunes a todas
        query_sets(up, ["FAP.1", "FAP.3"])             # solo en esas dos
        query_sets(up, ["FAP.1", "FAP.3"], exact=False)  # en las dos, y quizá más
        query_sets(up, n_populations=3)                # en exactamente 3, cualesquiera

    LA DISTINCIÓN QUE IMPORTA: `exact`
        True  -> el gen está en las poblaciones de `include` Y EN NINGUNA MÁS.
                 Es lo que dibuja cada barra del UpSet.
        False -> el gen está al menos en las de `include`, pudiendo estar en
                 otras.

        "Genes comunes a FAP.1 y FAP.3" significa cosas distintas según cuál
        elijas, y la diferencia suele ser grande. Con exact=True estás pidiendo
        "compartido por esas dos y por nadie más", que es una afirmación de
        especificidad; con exact=False, solo "presente en ambas". La función
        imprime los dos números para que no los confundas.

    res
        Si le pasas el dict de deg_by_group(), añade el logFC de CADA población
        para los genes encontrados, incluidas aquellas en cuyo conjunto el gen
        no entró. Eso es lo interesante: un gen "específico de FAP.4" cuyo
        logFC en FAP.1 es 0.9 —justo por debajo de tu umbral— no es
        específico, es un artefacto del corte. Si en las demás está en 0.05,
        entonces sí.
    """
    todas = list(sets)
    if n_populations is not None:
        elegidos = [g for g in set().union(*sets.values())
                    if sum(g in v for v in sets.values()) == n_populations]
        etiqueta = f"en exactamente {n_populations} poblaciones"
    else:
        if not include:
            raise ValueError("Da 'include' o 'n_populations'")
        falta = [p for p in include if p not in sets]
        if falta:
            raise ValueError(f"{falta} no están en sets. Disponibles: {todas}")
        dentro = set.intersection(*[sets[p] for p in include])
        fuera = set(exclude) if exclude else (set(todas) - set(include)
                                             if exact else set())
        prohibidos = set().union(*[sets[p] for p in fuera]) if fuera else set()
        elegidos = sorted(dentro - prohibidos)
        etiqueta = (f"solo en {list(include)}" if exact and not exclude
                    else f"en {list(include)}"
                         + (f", excluyendo {list(exclude)}" if exclude else ""))
        if verbose and not exclude:
            otras = set(todas) - set(include)
            en_otras = set().union(*[sets[p] for p in otras]) if otras else set()
            print(f"[query] en {list(include)} (pudiendo estar en otras): "
                  f"{len(dentro)} genes")
            print(f"[query] SOLO en {list(include)}: {len(dentro - en_otras)} genes")

    elegidos = sorted(elegidos)
    out = pd.DataFrame({"gen": elegidos})
    out["poblaciones"] = [";".join(sorted(p for p in todas if g in sets[p]))
                          for g in elegidos]
    out["n_poblaciones"] = out["poblaciones"].str.count(";") + 1

    if res is not None:
        for pop in todas:
            if pop not in res:
                continue
            d = res[pop].set_index("names")
            out[f"lfc_{pop}"] = d["logfoldchanges"].reindex(elegidos).round(2).values
        cols_lfc = [c for c in out.columns if c.startswith("lfc_")]
        if cols_lfc:
            out = out.sort_values(cols_lfc[0], ascending=False)

    if verbose:
        print(f"\n[query] {len(out)} genes {etiqueta}")
        if len(out):
            print(", ".join(out["gen"].tolist()[:60])
                  + (" ..." if len(out) > 60 else ""))
    return out.reset_index(drop=True)


def core_vs_specific(sets: Mapping[str, set], verbose: bool = True) -> pd.DataFrame:
    """
    Reparte los genes en núcleo compartido / intermedios / específicos.

    El núcleo (en todas las poblaciones) es la respuesta genérica a la lesión;
    los específicos son lo que distingue a cada población. Suelen ser dos
    historias distintas y conviene pasarlas a GOEA por separado — un análisis
    conjunto queda dominado por el núcleo, que es el mismo para todos.
    """
    todos = sorted(set().union(*sets.values())) if sets else []
    n = len(sets)
    filas = []
    for g in todos:
        donde = [k for k, v in sets.items() if g in v]
        filas.append({"gen": g, "n_poblaciones": len(donde),
                      "poblaciones": ";".join(sorted(donde))})
    out = pd.DataFrame(filas).sort_values(["n_poblaciones", "gen"],
                                          ascending=[False, True])
    if verbose and len(out):
        print(out["n_poblaciones"].value_counts().sort_index(ascending=False)
              .rename("n_genes").to_string())
        nucleo = out[out["n_poblaciones"] == n]["gen"].tolist()
        print(f"\n[core] núcleo ({n}/{n} poblaciones), {len(nucleo)} genes:")
        print("  " + ", ".join(nucleo[:40]) + (" ..." if len(nucleo) > 40 else ""))
        unicos = out[out["n_poblaciones"] == 1]
        print(f"\n[core] {len(unicos)} genes en una sola población")
        print("[core] Pasa núcleo y específicos a run_goea() por SEPARADO: "
              "juntos, el núcleo domina el resultado y tapa lo propio de cada "
              "población.")
    return out