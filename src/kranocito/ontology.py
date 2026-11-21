"""
Sobre-representación de términos GO (ORA) por población, con gseapy.

    fondo  = build_background(adata, min_cells=10)
    degs   = top_degs(adata, "minor_population", n_top=150)
    res    = run_goea(degs, fondo, organism="Mouse")
    plot_goea(res)

EL FONDO ES LA MITAD DEL RESULTADO
----------------------------------
Un test de sobre-representación pregunta: "de los N genes que PODRÍAN haber
salido, han salido k de este término, ¿es más de lo esperado?". Ese "podrían
haber salido" es el fondo, y elegirlo mal cambia el resultado más que ninguna
otra decisión.

El default de casi todas las herramientas es el genoma entero (~22.000 genes en
ratón). En single-cell eso está mal, y siempre en la misma dirección: tu ensayo
solo puede detectar los genes que se expresan en TU tejido, así que comparar
contra el genoma completo infla sistemáticamente todos los términos propios del
tejido. En un dataset de músculo te saldrá "muscle contraction" enriquecido en
todas las poblaciones, incluidas las que no son musculares.

LA REGLA: el fondo es el conjunto de genes sobre los que CORRIÓ EL TEST de
expresión diferencial. Si rank_genes_groups corrió sobre los 25.000 genes que
quedaron tras filter_genes, el fondo son esos 25.000. Ni más (genes que nunca
pudieron salir) ni menos (genes que sí compitieron).

¿HVGs COMO FONDO? Preguntaste por ellos, y mi respuesta es que no por defecto.
Los HVGs no son "los genes que compitieron": rank_genes_groups no se limita a
ellos. Usarlos como fondo responde a otra pregunta —"entre los genes variables,
¿este término está sobrerrepresentado en los marcadores de esta población?"—
que es legítima y más conservadora, pero no es la que uno cree estar haciendo.
Si la quieres, build_background(use_hvg=True) la construye y deja constancia en
el nombre. Lo que NO puede pasar es mezclar: un fondo de HVGs con DEGs
calculados sobre todos los genes hace que parte del primer plano ni siquiera
esté en el universo, y el test queda mal definido.

DOBLE INMERSIÓN
---------------
Estos DEGs vienen de comparar clusters que se definieron con los mismos datos.
Los p-valores del test de expresión diferencial ya no son válidos como
evidencia, y los de GO heredan el problema entero. Sirven para ORDENAR
hipótesis sobre qué hace cada población, no para afirmar que la hace.

gseapy: OJO CON QUÉ FUNCIÓN
--------------------------
gp.enrichr() usa la API web de Enrichr y IGNORA tu fondo. La que acepta fondo
personalizado es gp.enrich(), que hace el test de Fisher en local. Este módulo
usa gp.enrich() siempre.
"""

from __future__ import annotations

import re
import warnings
from typing import Iterable, Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


# NO se fija la versión aquí a propósito: Enrichr publica GO cada año y
# hardcodear un año significa quedarse con una anotación vieja sin enterarse.
# latest_go_libraries() pregunta a Enrichr cuáles hay AHORA y coge la más nueva.
PREFIJOS_GO = ("GO_Biological_Process", "GO_Molecular_Function",
               "GO_Cellular_Component")

_RE_ANYO = re.compile(r"_(\d{4})$")


def list_libraries(organism: str = "Mouse", pattern: Optional[str] = None,
                   verbose: bool = True) -> list:
    """Librerías disponibles en Enrichr, opcionalmente filtradas por subcadena."""
    import gseapy as gp

    libs = gp.get_library_name(organism=organism)
    if pattern:
        libs = [l for l in libs if pattern.lower() in l.lower()]
    if verbose:
        for l in sorted(libs):
            print(" ", l)
    return sorted(libs)


def latest_go_libraries(organism: str = "Mouse", prefixes: Sequence[str] = PREFIJOS_GO,
                        verbose: bool = True) -> list:
    """
    La versión MÁS RECIENTE de cada rama de GO que Enrichr sirva hoy.

    Preguntar en vez de hardcodear un año tiene dos ventajas: no te quedas con
    una anotación obsoleta, y no revienta cuando Enrichr retira una versión.

    Pero apunta cuál te ha devuelto. GO cambia entre versiones —términos que se
    fusionan, genes que se reanotan— y dos análisis con versiones distintas no
    son comparables. run_goea() guarda los nombres resueltos en el resultado
    justamente por eso.
    """
    import gseapy as gp

    disponibles = gp.get_library_name(organism=organism)
    elegidas = []
    for pref in prefixes:
        cand = [l for l in disponibles if l.startswith(pref)]
        if not cand:
            warnings.warn(f"Enrichr ({organism}) no ofrece ninguna '{pref}*'")
            continue
        def _anyo(l):
            m = _RE_ANYO.search(l)
            return int(m.group(1)) if m else -1
        elegidas.append(max(cand, key=_anyo))
    if verbose:
        print(f"[libs] {organism}: {elegidas}")
    return elegidas


# --------------------------------------------------------------------------
# 1. el fondo
# --------------------------------------------------------------------------

def build_background(
    adata: ad.AnnData,
    min_cells: int = 10,
    group_key: Optional[str] = None,
    layer: Optional[str] = "counts_cellbender",
    use_hvg: bool = False,
    hvg_key: str = "hvg_triku",
    exclude: Optional[Iterable[str]] = None,
    verbose: bool = True,
) -> list:
    """
    Universo de genes para el test. Devuelve una lista de símbolos.

    min_cells
        Un gen detectado en 3 células de 9.000 no tuvo ninguna posibilidad real
        de salir como marcador, así que meterlo en el universo solo diluye.

        POR QUÉ ESTE CRITERIO Y NO EL TOTAL DE COUNTS. Recortar por
        sum(counts) parece equivalente y no lo es: un gen con 100 counts
        totales puede ser 10 counts en cada una de 10 células de una población
        de 12 — es decir, el marcador perfecto de una población rara. El total
        de counts confunde "poco expresado en todas partes" con "muy expresado
        en unas pocas", y solo la segunda te interesa. Contar CÉLULAS separa
        esos dos casos; sumar counts, no.

    group_key
        Si lo das (p.ej. 'minor_population'), un gen entra en el fondo cuando
        se detecta en >= min_cells células DE AL MENOS UNA población, en vez de
        en el objeto entero. Es el criterio correcto para una ORA de marcadores
        de cluster: la pregunta es "¿pudo este gen haber salido como marcador
        de ALGUNA población?", y con poblaciones pequeñas el umbral global
        elimina justo los marcadores de las poblaciones raras. Con kranocitos
        en el objeto, esto no es una sutileza.

    layer
        Se cuenta sobre COUNTS, no sobre la X normalizada, porque "expresado en
        N células" solo significa algo sobre counts. Si la capa no está, cae a
        X con un aviso.

    use_hvg
        Cambia la pregunta (ver la cabecera del módulo). Si lo activas,
        asegúrate de que tus DEGs se calcularon también restringiendo a HVGs.

    exclude
        Genes a sacar del universo Y, por coherencia, del primer plano. Aquí van
        las globinas y lo que sepas que es sopa: si un gen está en el dataset
        por contaminación, no debería contar como "podría haber salido".
    """
    if use_hvg:
        if hvg_key not in adata.var:
            raise ValueError(f"No hay var['{hvg_key}']")
        genes = adata.var_names[adata.var[hvg_key].to_numpy().astype(bool)]
        origen = f"HVGs ({hvg_key})"
    else:
        if layer and layer in adata.layers:
            M = adata.layers[layer]
            fuente = layer
        else:
            M = adata.X
            fuente = "X"
            if layer:
                warnings.warn(
                    f"No hay capa '{layer}'; cuento células con valor distinto "
                    f"de 0 en X. Sobre datos log-normalizados el criterio sigue "
                    f"siendo válido (cero es cero), pero déjalo escrito."
                )
        if group_key is None:
            n_cel = np.asarray((M > 0).sum(axis=0)).ravel()
            keep = n_cel >= min_cells
            origen = f"detectados en >={min_cells} células sobre '{fuente}'"
        else:
            if group_key not in adata.obs:
                raise ValueError(f"No hay obs['{group_key}']")
            keep = np.zeros(adata.n_vars, dtype=bool)
            tam = {}
            for grp in adata.obs[group_key].astype(str).unique():
                m = (adata.obs[group_key].astype(str) == grp).to_numpy()
                n_cel = np.asarray((M[m] > 0).sum(axis=0)).ravel()
                keep |= n_cel >= min_cells
                tam[grp] = int(m.sum())
            origen = (f">={min_cells} células en al menos una de las "
                      f"{len(tam)} categorías de '{group_key}' sobre '{fuente}'")
            if verbose:
                peq = {k: v for k, v in tam.items() if v < 2 * min_cells}
                if peq:
                    print(f"[fondo] OJO: {peq} tienen menos de {2*min_cells} "
                          f"células, así que min_cells={min_cells} es casi todo "
                          f"su tamaño. Baja min_cells o esas poblaciones no "
                          f"aportan ningún gen al universo.")
        genes = adata.var_names[keep]

    genes = pd.Index(genes)
    n0 = len(genes)
    if exclude:
        exclude = set(exclude)
        genes = genes[~genes.isin(exclude)]

    if verbose:
        print(f"[fondo] {len(genes)} genes de {adata.n_vars} ({origen})")
        if exclude:
            print(f"[fondo] {n0 - len(genes)} excluidos explícitamente")
        print("[fondo] Comprueba que este número coincide con los genes sobre "
              "los que corrió rank_genes_groups. Si no coincide, el test está "
              "mal planteado en un sentido o en el otro.")
    return list(genes)


# --------------------------------------------------------------------------
# 2. el primer plano
# --------------------------------------------------------------------------

def top_degs(
    adata: ad.AnnData,
    groupby: str,
    n_top: int = 150,
    key: str = "rank_genes_groups",
    method: str = "wilcoxon",
    tie_correct: bool = True,
    max_padj: float = 0.05,
    min_logfc: float = 0.25,
    min_pct: float = 0.10,
    rank_by: str = "pi",
    exclude: Optional[Iterable[str]] = None,
    recompute: bool = False,
    verbose: bool = True,
) -> dict:
    """
    Top N marcadores por población, con filtros, como dict {población: [genes]}.

    No basta con coger los N primeros por score: entre ellos puede haber genes
    con un efecto minúsculo o expresados en el 2% de las células, y esos
    arrastran términos GO al azar. Los filtros van ANTES del recorte a N.

    tie_correct
        Solo para Wilcoxon, y solo cuando recompute=True. scanpy lo trae en
        False; con una matriz llena de ceros los empates dominan el ranking y
        no corregirlos infla la significación. True por defecto aquí.

    min_pct
        Requiere que rank_genes_groups se haya corrido con pts=True. Si no está
        disponible, el filtro se salta y se avisa.

    rank_by
        Con qué criterio se ordena antes de cortar a n_top:

        'pi'     -log10(padj) * logFC. El pi-value de Xiao et al. (2014):
                 combina evidencia y tamaño de efecto, en vez de dejar que
                 mande solo una de las dos. Es el default.
        'score'  el z-score que devuelve scanpy (su orden nativo).
        'logfc'  solo tamaño de efecto. Útil como control: si el top-N cambia
                 mucho entre 'logfc' y 'pi', tu lista la está decidiendo el
                 p-valor, es decir, el tamaño del cluster.

        AVISO SOBRE 'pi' CON MUCHAS CÉLULAS. Con 9.000 células los p-valores se
        desbordan a 0 y -log10(p) satura; aquí se recortan al menor valor
        positivo observado para que el producto siga siendo finito, pero eso
        significa que decenas de genes comparten el mismo -log10(p) y entre
        ellos ordena el LFC solo. No es un defecto del código, es que el
        p-valor ya no discrimina a ese tamaño muestral.

    Devuelve también, en el dict, cuántos genes ha dejado cada población: si una
    se queda con 12 genes y otra con 150, sus resultados de GO no son
    comparables en potencia.
    """
    if recompute or key not in adata.uns:
        if verbose:
            print(f"[deg] calculando rank_genes_groups ({method}) sobre '{groupby}'")
        kw = dict(groupby=groupby, method=method, key_added=key, pts=True)
        if method == "wilcoxon":
            # Con datos tan dispersos los empates (todos los ceros) están por
            # todas partes. scanpy trae tie_correct=False por defecto y eso
            # infla la significación; aquí se activa.
            kw["tie_correct"] = tie_correct
        sc.tl.rank_genes_groups(adata, **kw)

    grupos = list(adata.uns[key]["names"].dtype.names)
    exclude = set(exclude or ())
    salida, resumen = {}, []

    for g in grupos:
        df = sc.get.rank_genes_groups_df(adata, group=g, key=key)
        n_ini = len(df)
        df = df[df["names"].notna()]

        if "pvals_adj" in df:
            df = df[df["pvals_adj"] <= max_padj]
        if "logfoldchanges" in df:
            df = df[df["logfoldchanges"] >= min_logfc]

        col_pct = next((c for c in ("pct_nz_group", "pts") if c in df.columns), None)
        if col_pct is not None:
            df = df[df[col_pct] >= min_pct]
        elif verbose and min_pct:
            print(f"[deg] sin columna de fracción de células; me salto min_pct "
                  f"(re-córrelo con pts=True si lo quieres)")

        df = df[~df["names"].isin(exclude)]

        # --- ordenación
        if rank_by == "pi":
            col_p = "pvals_adj" if "pvals_adj" in df else "pvals"
            pv = df[col_p].to_numpy(dtype=float)
            positivos = pv[pv > 0]
            suelo = float(positivos.min()) if positivos.size else 1e-300
            pv = np.clip(pv, suelo, 1.0)
            df = df.assign(_pi=-np.log10(pv) * df["logfoldchanges"].to_numpy(float))
            df = df.sort_values("_pi", ascending=False)
        elif rank_by == "logfc":
            df = df.sort_values("logfoldchanges", ascending=False)
        elif rank_by == "score":
            if "scores" in df:
                df = df.sort_values("scores", ascending=False)
        else:
            raise ValueError("rank_by debe ser 'pi', 'score' o 'logfc'")

        genes = df["names"].tolist()[:n_top]
        salida[g] = genes

        fila = {"poblacion": g, "testados": n_ini,
                "tras_filtros": len(df), "usados": len(genes)}
        if verbose and rank_by != "score" and "scores" in df:
            alt = df.sort_values("scores", ascending=False)["names"].tolist()[:n_top]
            inter = len(set(genes) & set(alt))
            fila["jaccard_vs_score"] = round(
                inter / max(len(set(genes) | set(alt)), 1), 2)
        resumen.append(fila)

    tabla = pd.DataFrame(resumen).set_index("poblacion")
    if verbose:
        print(tabla.to_string())
        pocos = tabla.index[tabla["usados"] < 20].tolist()
        if pocos:
            print(f"[deg] AVISO: {pocos} se quedan con menos de 20 genes. Su "
                  f"GOEA tendrá muy poca potencia y la ausencia de términos NO "
                  f"significa ausencia de biología.")
        if "jaccard_vs_score" in tabla:
            j = tabla["jaccard_vs_score"].mean()
            print(f"[deg] Jaccard medio entre el top-{n_top} por '{rank_by}' y "
                  f"por el score de scanpy: {j:.2f}. Si está cerca de 1, la "
                  f"elección del criterio no cambia nada y no merece discusión.")
        desigual = tabla["usados"].max() / max(tabla["usados"].min(), 1)
        if desigual > 3:
            print(f"[deg] Las listas difieren hasta {desigual:.0f}x en tamaño. "
                  f"El número de términos significativos escala con el tamaño "
                  f"de la lista: no compares poblaciones por cuántos términos "
                  f"sacan.")
    return salida


# --------------------------------------------------------------------------
# 3. cobertura: ¿hablan el mismo idioma tus símbolos y la librería?
# --------------------------------------------------------------------------

def _genes_de_libreria(gs: Mapping[str, Sequence[str]]) -> set:
    out = set()
    for v in gs.values():
        out.update(v)
    return out


def coverage_report(background: Sequence[str], gene_sets: Mapping,
                    nombre: str = "", verbose: bool = True) -> dict:
    """
    Qué fracción de tu fondo existe en la librería, tal cual y en mayúsculas.

    Es la comprobación que casi nadie hace y que invalida resultados enteros:
    si tus símbolos son de RefSeq ratón (Col1a1) y la librería es humana
    (COL1A1), el solapamiento es casi nulo, gseapy no falla, y te devuelve una
    tabla vacía o —peor— basada en el puñado de genes que coinciden por azar.
    """
    lib = _genes_de_libreria(gene_sets)
    bg = set(background)
    tal_cual = len(bg & lib) / max(len(bg), 1)
    mayus = len({g.upper() for g in bg} & {g.upper() for g in lib}) / max(len(bg), 1)
    out = {"libreria": nombre, "n_fondo": len(bg), "n_libreria": len(lib),
           "cobertura": round(tal_cual, 3), "cobertura_mayusculas": round(mayus, 3)}
    if verbose:
        print(f"[cov] {nombre}: {tal_cual:.1%} del fondo está en la librería "
              f"({mayus:.1%} ignorando mayúsculas)")
        if tal_cual < 0.30:
            print("[cov] AVISO: cobertura muy baja.")
            if mayus > tal_cual + 0.2:
                print("[cov] Ignorando mayúsculas sube mucho -> tus símbolos y "
                      "los de la librería son de organismos/nomenclaturas "
                      "distintas. Comprueba organism= antes de leer nada.")
    return out


def _resolver_nomenclatura(
    background: Sequence[str],
    gene_dict: Mapping[str, Sequence[str]],
    gene_sets: Mapping[str, Mapping],
    match_case: str = "auto",
    verbose: bool = True,
):
    """
    Pone tus símbolos y los de la librería en el mismo alfabeto.

    POR QUÉ HACE FALTA. Enrichr sirve sus librerías de ratón con los símbolos
    en MAYÚSCULAS (COL1A1), mientras que la nomenclatura MGI/RefSeq de ratón es
    capitalizada (Col1a1). Aunque pidas organism='Mouse', los GENES vienen en
    mayúsculas: el organismo determina QUÉ anotación te dan, no cómo escriben
    los símbolos. Sin normalizar, el solapamiento es del 0.1% y gseapy no
    falla: devuelve tablas casi vacías con los pocos genes que coinciden por
    azar (los que ya eran todo mayúsculas, tipo mitocondriales).

    LO QUE ESTO NO ARREGLA. Poner en mayúsculas empareja los símbolos que solo
    difieren en capitalización, que son la mayoría, pero no los que son nombres
    distintos entre especies: Trp53 no se convierte en TP53. Esos genes se
    quedan fuera y no hay forma de recuperarlos sin una tabla de ortología.

    match_case: 'auto' | 'upper' | 'exact'.
    """
    bgs = list(dict.fromkeys(background))
    lib = set()
    for gs in gene_sets.values():
        lib |= _genes_de_libreria(gs)

    exacta = len(set(bgs) & lib) / max(len(bgs), 1)
    mayus = len({g.upper() for g in bgs} & {g.upper() for g in lib}) / max(len(bgs), 1)

    if match_case == "auto":
        modo = "upper" if mayus > exacta + 0.10 else "exact"
        if verbose:
            print(f"[case] solapamiento exacto {exacta:.1%}, en mayúsculas "
                  f"{mayus:.1%} -> modo '{modo}'")
            if modo == "upper":
                print("[case] Enrichr sirve las librerías de ratón con símbolos "
                      "en MAYÚSCULAS. Se normalizan los dos lados. Ojo: esto NO "
                      "empareja genes con nombre distinto entre especies "
                      "(Trp53 != TP53); esos se pierden.")
    else:
        modo = match_case

    if modo != "upper":
        return bgs, dict(gene_dict), dict(gene_sets), None, modo

    mapa = {}
    colisiones = []
    for g in bgs:
        u = g.upper()
        if u in mapa and mapa[u] != g:
            colisiones.append((mapa[u], g))
        else:
            mapa[u] = g
    if colisiones and verbose:
        print(f"[case] {len(colisiones)} colisiones al pasar a mayúsculas "
              f"(p.ej. {colisiones[:3]}); me quedo con la primera de cada par.")

    bg2 = list(dict.fromkeys(g.upper() for g in bgs))
    gd2 = {k: [g.upper() for g in v] for k, v in gene_dict.items()}
    gs2 = {lib_n: {t: [g.upper() for g in gg] for t, gg in gs.items()}
           for lib_n, gs in gene_sets.items()}
    return bg2, gd2, gs2, mapa, modo


# --------------------------------------------------------------------------
# 4. el test
# --------------------------------------------------------------------------

def run_goea(
    gene_dict: Mapping[str, Sequence[str]],
    background: Sequence[str],
    libraries: Optional[Sequence[str]] = None,
    organism: str = "Mouse",
    gene_sets: Optional[Mapping[str, Mapping]] = None,
    min_set_size: int = 5,
    max_set_size: int = 500,
    match_case: str = "auto",
    global_fdr: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    ORA por población contra un fondo propio, con gp.enrich() (no gp.enrichr()).

    gene_sets
        Opcional, {nombre_libreria: {termino: [genes]}}. Si no lo das se
        descargan con gp.get_library(), lo cual necesita red. Guárdalos y
        pásalos si vas a repetir: así el análisis es reproducible aunque
        Enrichr cambie de versión, que cambia.

    min_set_size / max_set_size
        Se aplican DESPUÉS de intersectar cada término con el fondo, que es
        donde tiene sentido: un término de 300 genes del que solo 4 están en tu
        tejido es, para este test, un término de 4 genes. Los términos enormes
        ('biological process') no discriminan nada y los diminutos dan
        p-valores inestables.

    global_fdr
        gseapy corrige por múltiple test DENTRO de cada llamada. Si corres 8
        poblaciones x 3 librerías, eso son 24 familias de tests y la corrección
        de cada una no sabe de las otras. Con True se añade 'padj_global' con
        Benjamini-Hochberg sobre TODO. Usa esa columna si vas a decir "estos
        términos son los significativos".
    """
    import gseapy as gp
    from statsmodels.stats.multitest import multipletests

    if type(background) is not list:
        bg = list(dict.fromkeys(background))
    else: 
        bg = background
    bg_set = set(bg)

    if gene_sets is None:
        if libraries is None:
            libraries = latest_go_libraries(organism=organism, verbose=verbose)
        gene_sets = {}
        for lib in libraries:
            if verbose:
                print(f"[goea] descargando {lib} ({organism})...")
            gene_sets[lib] = gp.get_library(name=lib, organism=organism)

    # los símbolos de Enrichr y los tuyos pueden estar en alfabetos distintos
    bg, gene_dict, gene_sets, mapa_case, modo_case = _resolver_nomenclatura(
        bg, gene_dict, gene_sets, match_case=match_case, verbose=verbose)
    bg_set = set(bg)

    filas = []
    for lib, gs in gene_sets.items():
        coverage_report(bg, gs, nombre=lib, verbose=verbose)

        # recortar los términos al fondo: el universo del test es el fondo
        gs_bg = {}
        for term, genes in gs.items():
            inter = [g for g in set(genes) if g in bg_set]
            if min_set_size <= len(inter) <= max_set_size:
                gs_bg[term] = inter
        if verbose:
            print(f"[goea] {lib}: {len(gs_bg)} de {len(gs)} términos con "
                  f"{min_set_size}-{max_set_size} genes dentro del fondo")
        if not gs_bg:
            warnings.warn(f"{lib}: ningún término utilizable tras intersectar "
                          f"con el fondo. Revisa la cobertura.")
            continue

        for pop, genes in gene_dict.items():
            genes = [g for g in dict.fromkeys(genes) if g in bg_set]
            if len(genes) < 5:
                warnings.warn(f"{pop}/{lib}: solo {len(genes)} genes del primer "
                              f"plano están en el fondo; me la salto")
                continue
            try:
                enr = gp.enrich(gene_list=genes, gene_sets=gs_bg, background=bg,
                                outdir=None, no_plot=True, verbose=False)
            except Exception as e:
                warnings.warn(f"{pop}/{lib} falló: {e}")
                continue
            if enr is None or enr.results is None or len(enr.results) == 0:
                continue
            df = enr.results.copy()
            df["poblacion"] = pop
            df["libreria"] = lib
            df["n_genes_disponibles"] = len(genes)
            filas.append(df)

    if not filas:
        raise RuntimeError(
            "Ningún test produjo resultados. Mira las líneas [cov] y [case] de "
            "arriba: si el solapamiento sigue siendo bajo incluso en mayúsculas, "
            "el problema no es la capitalización sino que tus símbolos no son de "
            "la especie/nomenclatura de la librería."
        )

    res = pd.concat(filas, ignore_index=True)

    ren = {"Term": "termino", "Overlap": "solapamiento", "P-value": "pval",
           "Adjusted P-value": "padj", "Odds Ratio": "odds_ratio",
           "Combined Score": "score", "Genes": "genes"}
    res = res.rename(columns={k: v for k, v in ren.items() if k in res.columns})

    if mapa_case is not None and "genes" in res:
        # devolver los símbolos a la nomenclatura del usuario
        res["genes"] = res["genes"].astype(str).apply(
            lambda s: ";".join(mapa_case.get(x, x) for x in s.split(";")))

    if global_fdr and "pval" in res:
        res["padj_global"] = multipletests(res["pval"].to_numpy(), method="fdr_bh")[1]

    orden = [c for c in ("poblacion", "libreria", "termino", "solapamiento",
                         "odds_ratio", "pval", "padj", "padj_global", "genes",
                         "n_primer_plano") if c in res.columns]
    res = res[orden + [c for c in res.columns if c not in orden]]
    res = res.sort_values(["poblacion", "padj_global" if global_fdr else "padj"])
    # procedencia: qué versión de GO se ha usado. GO cambia entre versiones y
    # dos análisis con librerías distintas no son comparables.
    res.attrs["librerias"] = list(gene_sets.keys())
    res.attrs["n_fondo"] = len(bg)
    res.attrs["organism"] = organism
    res.attrs["match_case"] = modo_case

    if verbose:
        col = "padj_global" if global_fdr else "padj"
        sig = res[res[col] < 0.05]
        print(f"\n[goea] {len(sig)} términos con {col} < 0.05, de {len(res)} testados")
        print(sig.groupby("poblacion").size().to_string()
              if len(sig) else "  (ninguno)")
        print("\n[goea] Recuerda: los DEGs vienen de clusters definidos con "
              "estos mismos datos. Esto ordena hipótesis, no las confirma.")
    return res


# --------------------------------------------------------------------------
# 5. ¿cambia la conclusión si cambio el fondo?
# --------------------------------------------------------------------------

def background_sensitivity(
    gene_dict: Mapping[str, Sequence[str]],
    backgrounds: Mapping[str, Sequence[str]],
    alpha: float = 0.05,
    col_p: str = "padj_global",
    verbose: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """
    Corre la misma ORA con VARIOS fondos y devuelve qué términos sobreviven a
    todos.

    Discutir cuál es el fondo correcto es interesante, pero la pregunta que de
    verdad importa es si la conclusión depende de esa elección. Si un término
    sale con los tres fondos, la discusión es teórica. Si solo sale con uno,
    tu resultado es un artefacto de una decisión de nuisance, y eso es en sí
    mismo el hallazgo.

        fondos = {
            "detectado_global": build_background(adata, min_cells=10),
            "detectado_por_pop": build_background(adata, min_cells=5,
                                                  group_key="minor_population"),
            "hvg": build_background(adata, use_hvg=True),
        }
        rob = background_sensitivity(degs, fondos, gene_sets=gs)

    Devuelve una fila por (población, término) con el p-valor bajo cada fondo y
    'n_fondos_sig' = en cuántos sale por debajo de alpha.
    """
    trozos = {}
    for nombre, bg in backgrounds.items():
        if verbose:
            print(f"\n===== fondo '{nombre}' ({len(set(bg))} genes) =====")
        r = run_goea(gene_dict, bg, verbose=verbose, **kwargs)
        trozos[nombre] = r.set_index(["poblacion", "libreria", "termino"])[col_p]

    tabla = pd.concat(trozos, axis=1)
    tabla.columns = list(trozos.keys())
    tabla["n_fondos_sig"] = (tabla < alpha).sum(axis=1)
    tabla["n_fondos_testado"] = tabla[list(trozos)].notna().sum(axis=1)
    tabla = tabla.sort_values(["n_fondos_sig", list(trozos)[0]],
                              ascending=[False, True])

    if verbose:
        n = len(backgrounds)
        robustos = tabla[tabla["n_fondos_sig"] == n]
        frag = tabla[(tabla["n_fondos_sig"] > 0) & (tabla["n_fondos_sig"] < n)]
        print(f"\n[sens] {len(robustos)} términos significativos con LOS {n} "
              f"fondos; {len(frag)} solo con algunos.")
        if len(frag):
            print("[sens] Los de la segunda lista dependen de una decisión de "
                  "nuisance. No los reportes sin decir con qué fondo salen.")
        print(tabla.head(20).round(4).to_string())
    return tabla


# --------------------------------------------------------------------------
# 6. redundancia: GO no son términos independientes
# --------------------------------------------------------------------------

def collapse_redundant(
    res: pd.DataFrame,
    col_p: str = "padj_global",
    alpha: float = 0.05,
    max_jaccard: float = 0.5,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Quita términos que son significativos POR LOS MISMOS GENES que otro mejor.

    Salir con 1.500 términos no significa haber encontrado 1.500 cosas. GO es
    un grafo jerárquico: 'muscle cell differentiation' está contenido en
    'muscle organ development', y si tus 20 genes musculares disparan los dos,
    has encontrado UN resultado y lo has contado dos veces. El número de
    términos significativos no es una medida de nada.

    Estrategia (voraz, por población): ordenar por p-valor y quedarse con un
    término solo si el conjunto de genes que lo hace significativo —la columna
    'genes', que es la intersección con tu primer plano, no el término entero—
    solapa menos de max_jaccard con el de todos los ya conservados.

    Es una versión pobre de REVIGO: agrupa por genes compartidos, no por
    similitud semántica en la ontología. A cambio no necesita nada externo y es
    exactamente el criterio que importa aquí, que es "¿son los mismos genes?".
    """
    if col_p not in res:
        col_p = "padj"
    sig = res[res[col_p] < alpha].copy()
    if not len(sig):
        if verbose:
            print(f"[colapso] ningún término con {col_p} < {alpha}")
        return sig

    sig["_set"] = sig["genes"].astype(str).apply(lambda x: frozenset(x.split(";")))

    guardados = []
    for pop, sub in sig.groupby("poblacion", sort=False):
        sub = sub.sort_values(col_p)
        elegidos = []
        for i, fila in sub.iterrows():
            s = fila["_set"]
            if not s:
                continue
            redundante = any(
                len(s & t) / max(len(s | t), 1) > max_jaccard for t in elegidos)
            if not redundante:
                elegidos.append(s)
                guardados.append(i)

    out = sig.loc[guardados].drop(columns=["_set"])
    if verbose:
        antes = sig.groupby("poblacion").size()
        despues = out.groupby("poblacion").size()
        tabla = pd.DataFrame({"significativos": antes, "no_redundantes": despues})
        tabla["reduccion"] = (1 - tabla["no_redundantes"] / tabla["significativos"]).round(2)
        print(tabla.to_string())
        print(f"\n[colapso] {len(sig)} -> {len(out)} términos "
              f"(Jaccard de genes > {max_jaccard} = redundante)")
        print("[colapso] Esto NO es corrección estadística: los p-valores no "
              "cambian. Solo deja de contar el mismo hallazgo varias veces.")
    return out


# --------------------------------------------------------------------------
# 7. dibujo
# --------------------------------------------------------------------------

def plot_goea(res: pd.DataFrame, n_top: int = 6, libreria: Optional[str] = None,
              col_p: str = "padj_global", figsize=(11, 8)):
    """
    Dotplot: términos (filas) x poblaciones (columnas), tamaño = solapamiento,
    color = -log10(p).

    Se muestran los n_top términos de CADA población, no los n_top globales:
    si no, la población con listas más largas se lleva la figura entera.
    """
    import matplotlib.pyplot as plt

    df = res if libreria is None else res[res["libreria"] == libreria]
    if col_p not in df:
        col_p = "padj"

    elegidos = (df.sort_values(col_p).groupby("poblacion").head(n_top)["termino"]
                .drop_duplicates().tolist())
    sub = df[df["termino"].isin(elegidos)].copy()

    sub["k"] = sub["solapamiento"].astype(str).str.split("/").str[0].astype(float)
    sub["mlog"] = -np.log10(sub[col_p].clip(lower=1e-300))

    pops = sorted(sub["poblacion"].unique())
    terms = (sub.groupby("termino")["mlog"].max().sort_values().index.tolist())
    xi = {p: i for i, p in enumerate(pops)}
    yi = {t: i for i, t in enumerate(terms)}

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    sc_ = ax.scatter([xi[p] for p in sub["poblacion"]],
                     [yi[t] for t in sub["termino"]],
                     s=15 * sub["k"], c=sub["mlog"], cmap="magma_r",
                     edgecolors="none")
    ax.set_xticks(range(len(pops)))
    ax.set_xticklabels(pops, rotation=45, ha="right")
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels([t[:70] for t in terms], fontsize=8)
    ax.set(xlabel="", ylabel="")
    ax.grid(alpha=.2)
    fig.colorbar(sc_, ax=ax, label=f"-log10({col_p})", shrink=.5)
    return fig, ax