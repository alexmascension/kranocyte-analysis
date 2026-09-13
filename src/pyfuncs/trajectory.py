"""
Inferencia de trayectorias: RNA velocity (scVelo), grafo dirigido (PAGA) y
embedding PHATE.

Sustituye a velocity_check.py, que se puede borrar.

QUÉ RESPONDE Y QUÉ NO
---------------------
Con scVelo + PAGA puedes preguntar "¿va A hacia B?". No puedes preguntar "¿qué
fracción de A acaba en B?": eso son probabilidades de destino y necesitan
Palantir o CellRank. Tenlo presente al escribir conclusiones.

PHATE está aquí como EMBEDDING, para mirar, no como método de inferencia. Los
grafos son para inferir; los embeddings son para dibujar. Reducir a dos
dimensiones y después medir distancias en ellas invalida lo que midas.

EL PROBLEMA DE FONDO
--------------------
Estos métodos SIEMPRE devuelven una trayectoria. Sobre un cluster homogéneo sin
ninguna transición real, PAGA dibuja un grafo con aristas y scVelo dibuja
flechas coherentes. La figura no es evidencia; la figura es la salida por
defecto del método. Por eso este módulo está construido alrededor de tres
comprobaciones que no dependen de tener etiquetas externas:

  1. ESTABILIDAD. ¿Sobrevive la arista al remuestreo de células?
     -> paga_stability()
  2. COHERENCIA INTERNA. ¿La dirección que da PAGA concuerda con el orden que
     da el pseudotiempo? Son dos lecturas del mismo campo de velocidad y no
     tienen por qué coincidir.
     -> paga_report() la calcula sola
  3. SENSIBILIDAD A LA RAÍZ. Si mover la raíz cambia el orden de poblaciones,
     estás leyendo tu hipótesis, no los datos.
     -> run_velocity(root_key=...) con distintas raíces

Si además tienes una variable ordinal externa —días post-lesión, estadio, dosis,
un score de proliferación—, pásala como `external_key` y se usa como ancla. Es
la evidencia más fuerte disponible, pero es OPCIONAL: el módulo no la necesita
y no está diseñado alrededor de ella.

ORIENTACIÓN DE LA MATRIZ DE TRANSICIONES
----------------------------------------
La convención de si transitions_confidence[i, j] es "de i a j" o "de j a i"
cambia entre versiones de scanpy/scVelo y es una fuente clásica de conclusiones
invertidas. Aquí no se asume: paga_report() contrasta el signo de cada arista
contra el orden del pseudotiempo y avisa si la matriz parece transpuesta.

AVISO
-----
Las llamadas a scVelo de este módulo se escribieron sin poder instalar scvelo
en el entorno donde se desarrolló. La orquestación está probada con mocks; las
llamadas a la librería no. Verifica cada salida antes de creértela.
"""

from __future__ import annotations

import warnings
from typing import Mapping, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


# --------------------------------------------------------------------------
# normalización de las capas de Velocyto
# --------------------------------------------------------------------------

def normalize_velocyto_layers(
    adata: ad.AnnData,
    size_factor_key: str = "size_factors",
    layers: Sequence[str] = ("spliced", "unspliced"),
    keep_raw_suffix: str = "_raw",
    verbose: bool = True,
) -> ad.AnnData:
    """
    Escala spliced/unspliced por los MISMOS size factors que se usaron para X.

    Tus capas de Velocyto son counts crudos a propósito. El flujo documentado
    de scVelo es filter_and_normalize -> moments, y si llamas a moments
    directamente te saltas la normalización. Usar los size factors ya
    calculados mantiene S, U y X en la misma escala, que es lo que quieres.

    OJO, esto NO es corrección de batch. Es normalización por profundidad. Son
    dos cosas distintas que en el 0A quedaron mezcladas dentro de la misma
    función.

    Guarda las originales como '<layer>_raw' para poder volver atrás.
    """
    if size_factor_key not in adata.obs:
        raise ValueError(
            f"No hay obs['{size_factor_key}']; los size factors vienen de "
            f"scran_size_factors() sobre el objeto completo"
        )
    sf = adata.obs[size_factor_key].to_numpy(dtype=np.float64)
    if not np.all(np.isfinite(sf)) or (sf <= 0).any():
        raise ValueError("Hay size factors no positivos o no finitos")

    inv = sp.diags(1.0 / sf)
    for lay in layers:
        if lay not in adata.layers:
            warnings.warn(f"No hay capa '{lay}'")
            continue
        raw_key = f"{lay}{keep_raw_suffix}"
        if raw_key not in adata.layers:
            adata.layers[raw_key] = adata.layers[lay].copy()
        # sparse @ diags mantiene CSR. (Nota: en scipy moderno `csr / sf`
        # tampoco densifica — devuelve COO —, así que la división directa del
        # 0A no era el problema que yo creía; pero @ diags evita el cambio de
        # formato y el coste de reconstruir el COO.)
        adata.layers[lay] = sp.csr_matrix(
            inv @ sp.csr_matrix(adata.layers[raw_key]), dtype=np.float32
        )

    adata.uns.setdefault("velocity", {})["normalization"] = {
        "size_factor_key": size_factor_key, "layers": list(layers),
    }
    if verbose:
        print(f"[velo] {list(layers)} escaladas por '{size_factor_key}'; "
              f"originales en '*{keep_raw_suffix}'")
    return adata


def _harmony_su(
    adata: ad.AnnData,
    batch_key: str,
    n_comps: int = 100,
    seed: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Segunda corrección de batch: sobre las propias matrices S y U (variante
    '2bc' del 0A).

    Descompone M = S + U, corrige M con harmony en espacio PCA, y reparte el M
    corregido entre S y U conservando la fracción spliced R = S/M de cada
    célula y gen.

    Lo que estaba mal en la versión original: Z_corr.T asumía harmonypy 1.x
    (PCs x células). Con 2.0 viene como (células x PCs) y esa línea produce
    basura o revienta. Aquí se comprueba la orientación.

    Se mantiene aquí para poder COMPARARLA, no porque la recomiende: aplicar
    harmony a S y U por separado es corregir dos veces el mismo efecto.
    """
    import harmonypy
    from sklearn.decomposition import PCA

    S = sp.csr_matrix(adata.layers["spliced"], dtype=np.float64)
    U = sp.csr_matrix(adata.layers["unspliced"], dtype=np.float64)
    M = (S + U).tocsr()

    Md = np.asarray(M.todense())
    Sd = np.asarray(S.todense())
    soporte = Md > 0
    with np.errstate(invalid="ignore", divide="ignore"):
        R = np.where(soporte, Sd / Md, 0.0)         # fracción spliced
    del Sd

    M_log = np.log1p(Md)
    del Md

    pca = PCA(n_components=min(n_comps, min(M_log.shape) - 1),
              random_state=seed, svd_solver="arpack")
    Z = pca.fit_transform(M_log)

    ho = harmonypy.run_harmony(Z, adata.obs, [batch_key], max_iter_harmony=30)
    Zc = np.asarray(ho.Z_corr)
    if Zc.shape == (Z.shape[1], adata.n_obs):        # harmonypy 1.x
        Zc = Zc.T
    elif Zc.shape != (adata.n_obs, Z.shape[1]):
        raise ValueError(f"Z_corr con shape inesperado: {Zc.shape}")

    M_log_h = Zc @ pca.components_ + pca.mean_
    M_h = np.expm1(M_log_h)
    np.clip(M_h, 0, None, out=M_h)
    # BUG CORREGIDO. La reconstrucción PCA es DENSA: devuelve un valor
    # positivo en entradas que en el original eran cero. Como ahí R = 0, todo
    # ese valor se iba entero a 'unspliced' (M_h * (1 - 0)), convirtiendo cada
    # cero de la matriz en señal no procesada. Resultado: U/S inflado de forma
    # sistemática en todo el objeto y velocidad sesgada hacia inducción en
    # todas partes, además de perder la dispersidad. Se restringe la
    # reconstrucción al soporte original.
    M_h[~soporte] = 0.0
    del M_log, M_log_h, Zc, Z, soporte

    adata.layers["spliced"] = sp.csr_matrix((M_h * R).astype(np.float32))
    adata.layers["unspliced"] = sp.csr_matrix((M_h * (1.0 - R)).astype(np.float32))
    del M_h, R

    if verbose:
        print(f"[velo] 2bc: S y U reconstruidas tras harmony sobre M = S + U")
    return adata


# --------------------------------------------------------------------------
# una corrida de velocity, sobre el objeto que se le pasa
# --------------------------------------------------------------------------

def run_velocity(
    adata: ad.AnnData,
    use_rep: str = "X_harmony",
    batch_key: Optional[str] = None,
    correct_su: bool = False,
    normalize_su: bool = True,
    ambient_max: Optional[float] = 0.2,
    ambient_key: str = "ambient_fraction",
    n_neighbors: Optional[int] = None,
    mode: str = "dynamical",
    root_key: Optional[str] = None,
    n_jobs: int = 4,
    seed: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    scVelo sobre EL OBJETO QUE SE LE PASA.

    En el 0A esta función escribía sobre el `adata` global en tres de sus
    cuatro llamadas, así que las velocidades se calculaban sobre todas las
    células y el subconjunto se quedaba vacío. Aquí no hay ninguna referencia
    a variables globales.

    correct_su
        False -> '1bc': el batch solo se corrige en la representación que
                 alimenta el grafo de vecinos (use_rep='X_harmony').
        True  -> '2bc': además se corrigen S y U. Para comparar, no por defecto.

    ambient_max
        Excluye del conjunto de genes de velocity los que tengan
        var[ambient_key] por encima de ese valor. El soup es 100% spliced, así
        que en los genes contaminados el U/S observado está sesgado a la baja y
        producen velocidad negativa espuria. None para desactivarlo.

    root_key
        Columna booleana de obs con las células raíz. SIN esto, scVelo elige la
        raíz sola y el signo del pseudotiempo es arbitrario en cada corrida:
        dos ejecuciones del mismo objeto pueden dar pseudotiempos que
        correlacionan a -1 sin que nada haya cambiado. Para COMPARAR corridas
        (variantes, subsets, remuestreos) hay que anclarla y usar la misma en
        todas. Para una corrida suelta da igual.
    """
    import scvelo as scv

    if use_rep not in adata.obsm:
        raise ValueError(f"No hay obsm['{use_rep}']")
    if n_neighbors is None:
        n_neighbors = max(10, int(0.5 * adata.n_obs ** 0.5))

    if normalize_su:
        normalize_velocyto_layers(adata, verbose=verbose)
    if correct_su:
        if batch_key is None:
            raise ValueError("correct_su=True necesita batch_key")
        if adata.obs[batch_key].nunique() < 2:
            raise ValueError(
                f"correct_su=True con una sola categoría de '{batch_key}': "
                f"no hay batch que corregir"
            )
        _harmony_su(adata, batch_key, seed=seed, verbose=verbose)

    n_pcs = adata.obsm[use_rep].shape[1]
    scv.pp.moments(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)
    scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
    scv.tl.velocity(adata, mode=mode)

    # ---- filtrar genes contaminados ANTES de construir el grafo
    if ambient_max is not None:
        if ambient_key not in adata.var:
            warnings.warn(
                f"No hay var['{ambient_key}']; no se filtran genes por ambient. "
                f"Calcúlalo con qc.ambient_fraction_per_gene() sobre el objeto "
                f"concatenado y arrástralo al subconjunto."
            )
        else:
            vg = adata.var["velocity_genes"].to_numpy().astype(bool)
            amb = adata.var[ambient_key].to_numpy(dtype=float)
            drop = vg & (amb > ambient_max)
            adata.var["velocity_genes"] = vg & ~(amb > ambient_max)
            if verbose:
                quitados = adata.var_names[drop].tolist()
                print(f"[velo] {int(drop.sum())} genes de velocity excluidos por "
                      f"{ambient_key} > {ambient_max}: {quitados[:10]}"
                      f"{' ...' if len(quitados) > 10 else ''}")
                print(f"[velo] quedan {int(adata.var['velocity_genes'].sum())} "
                      f"genes de velocity")

    scv.tl.velocity_graph(adata, n_jobs=n_jobs)
    if root_key is not None:
        if root_key not in adata.obs:
            raise ValueError(f"No hay obs['{root_key}']")
        scv.tl.velocity_pseudotime(adata, root_key=root_key)
    else:
        scv.tl.velocity_pseudotime(adata)
        if verbose:
            print("[velo] raíz elegida por scVelo: el SIGNO del pseudotiempo "
                  "es arbitrario. No lo compares entre corridas sin root_key.")
    scv.tl.velocity_confidence(adata)

    adata.uns.setdefault("velocity", {})["run"] = {
        "use_rep": use_rep, "correct_su": bool(correct_su),
        "normalize_su": bool(normalize_su), "ambient_max": ambient_max,
        "mode": mode, "n_neighbors": n_neighbors, "root_key": root_key,
        "n_velocity_genes": int(adata.var["velocity_genes"].sum()),
    }
    if verbose:
        print(f"[velo] confianza media: "
              f"{float(adata.obs['velocity_confidence'].mean()):.3f}")
    return adata


# --------------------------------------------------------------------------
# PAGA: grafo dirigido entre poblaciones
# --------------------------------------------------------------------------

def _parche_scvelo_igraph() -> bool:
    """
    scv.tl.paga revienta con scipy moderno. No es culpa tuya ni del pipeline.

    scvelo/tools/paga.py hace:

        csr_matrix((weights, zip(*edges)), shape=shape)

    Un `zip` es un generador. scipy antiguo lo materializaba por el camino;
    scipy nuevo pide secuencias de verdad y se encuentra con cero arrays de
    índices, de ahí el error:

        ValueError: mismatching number of index arrays for shape;
                    got 0, expected 2

    El parche NO toca la lógica de aristas de scVelo: solo sustituye el
    csr_matrix que ese módulo ve por uno que materializa el generador antes de
    pasárselo a scipy. Idempotente; devuelve True si ha hecho falta parchear.
    """
    import importlib
    import sys

    # OJO: `import scvelo.tools.paga as _pg` NO sirve. El __init__ de
    # scvelo.tools hace `from .paga import paga`, así que el atributo `paga` del
    # paquete apunta a la FUNCIÓN y tapa al submódulo; la forma `import ... as`
    # resuelve por atributo y te devuelve la función. El módulo de verdad está
    # en sys.modules bajo su nombre completo.
    importlib.import_module("scvelo.tools.paga")
    _pg = sys.modules["scvelo.tools.paga"]
    if not hasattr(_pg, "csr_matrix"):
        raise RuntimeError(
            "scvelo.tools.paga no expone csr_matrix: tu versión de scvelo no es "
            "la que este parche espera. Comprueba si scv.tl.paga ya funciona sin él."
        )

    if getattr(_pg.csr_matrix, "_parcheado", False):
        return False

    _csr = sp.csr_matrix

    def csr_matrix(arg1, shape=None, dtype=None, **kw):
        if isinstance(arg1, tuple) and len(arg1) == 2 and not np.isscalar(arg1[0]):
            data, ij = arg1
            if not isinstance(ij, (np.ndarray, list, tuple)):
                arg1 = (np.asarray(data), tuple(np.asarray(x) for x in ij))
        return _csr(arg1, shape=shape, dtype=dtype, **kw)

    csr_matrix._parcheado = True
    _pg.csr_matrix = csr_matrix
    return True


def _cats(adata: ad.AnnData, group_key: str) -> list:
    s = adata.obs[group_key]
    if not isinstance(s.dtype, pd.CategoricalDtype):
        adata.obs[group_key] = s.astype("category")
    elif len(s.cat.categories) != s.nunique():
        # Categorías sin células rompen los índices de PAGA en silencio: la
        # matriz sale con filas vacías y las etiquetas se descolocan.
        adata.obs[group_key] = s.cat.remove_unused_categories()
    return list(adata.obs[group_key].cat.categories)


def _dense(M) -> np.ndarray:
    if M is None:
        return None
    return np.asarray(M.todense()) if sp.issparse(M) else np.asarray(M)


def paga_report(
    adata: ad.AnnData,
    group_key: str,
    use_time_prior: str = "velocity_pseudotime",
    external_key: Optional[str] = None,
    min_cells: int = 20,
    verbose: bool = True,
) -> dict:
    """
    PAGA dirigido por velocidad, devuelto como NÚMEROS.

    Usa scv.tl.paga, no sc.tl.paga(use_rna_velocity=True), que es el camino
    roto/deprecado.

    Devuelve un dict:

      'conectividad'  DataFrame simétrico. Sale del grafo kNN y dice si dos
                      clusters son ADYACENTES. Es una afirmación sobre el
                      embedding, no sobre biología ni sobre dirección. Es la
                      que dibuja PAGA por defecto y la que se malinterpreta
                      como linaje.
      'transiciones'  DataFrame asimétrico. Sale de la velocidad. Esta es la
                      que lleva la flecha.
      'flujo_neto'    T - T.T, antisimétrica. Positivo en [a, b] = el flujo va
                      de a hacia b. Es la lectura robusta: no depende de la
                      raíz ni del signo global del pseudotiempo.
      'aristas'       tabla ordenada por |flujo_neto| con las dos métricas al
                      lado, más el salto de pseudotiempo de cada arista.
      'orden'         pseudotiempo medio por población.
      'coherencia'    fracción de aristas en las que el signo del flujo neto
                      concuerda con el orden del pseudotiempo.

    NO fija un threshold. El aspecto del grafo depende brutalmente de él, y esa
    es una decisión de dibujo, no un resultado: elígelo al plotear y decláralo.

    external_key
        Columna ORDINAL opcional (días, estadio, dosis, un score). Si la das se
        añade su correlación de Spearman con el pseudotiempo. No hace falta.
    """
    import scvelo as scv

    if _parche_scvelo_igraph() and verbose:
        print("[paga] aplicado el parche de compatibilidad scvelo/scipy "
              "(ver _parche_scvelo_igraph)")

    cats = _cats(adata, group_key)
    tam = adata.obs[group_key].value_counts()
    peq = tam[tam < min_cells]
    if len(peq) and verbose:
        print(f"[paga] AVISO: {list(peq.index)} tienen menos de {min_cells} "
              f"células ({dict(peq)}). Sus aristas son ruido; no las leas.")

    if use_time_prior and use_time_prior not in adata.obs:
        raise ValueError(
            f"No hay obs['{use_time_prior}']; ejecuta run_velocity() antes"
        )
    scv.tl.paga(adata, groups=group_key, use_time_prior=use_time_prior)

    P = adata.uns["paga"]
    conect = pd.DataFrame(_dense(P.get("connectivities")), index=cats, columns=cats)
    T = _dense(P.get("transitions_confidence"))
    if T is None:
        raise RuntimeError(
            "scv.tl.paga no ha devuelto 'transitions_confidence'. Sin eso solo "
            "tienes el grafo no dirigido, que no dice nada sobre dirección."
        )
    trans = pd.DataFrame(T, index=cats, columns=cats)
    neto = trans - trans.T

    orden = adata.obs.groupby(group_key, observed=True)[use_time_prior].mean()

    filas = []
    for i, a in enumerate(cats):
        for b in cats[i + 1:]:
            filas.append({
                "a": a, "b": b,
                "flujo_neto": float(neto.loc[a, b]),
                "trans_a_b": float(trans.loc[a, b]),
                "trans_b_a": float(trans.loc[b, a]),
                "conectividad": float(conect.loc[a, b]),
                "delta_pseudotiempo": float(orden[b] - orden[a]),
                "n_a": int(tam[a]), "n_b": int(tam[b]),
            })
    aristas = pd.DataFrame(filas)
    aristas["direccion"] = np.where(aristas["flujo_neto"] > 0,
                                    aristas["a"] + " -> " + aristas["b"],
                                    aristas["b"] + " -> " + aristas["a"])
    aristas = aristas.reindex(
        aristas["flujo_neto"].abs().sort_values(ascending=False).index
    ).reset_index(drop=True)

    # ---- coherencia interna: ¿el flujo va hacia pseudotiempo creciente?
    vivas = aristas[aristas["flujo_neto"].abs() > 1e-8]
    if len(vivas):
        coherencia = float(
            (np.sign(vivas["flujo_neto"]) == np.sign(vivas["delta_pseudotiempo"])).mean()
        )
    else:
        coherencia = np.nan

    out = {"conectividad": conect, "transiciones": trans, "flujo_neto": neto,
           "aristas": aristas, "orden": orden, "coherencia": coherencia}

    if external_key is not None:
        from scipy.stats import spearmanr
        v = adata.obs[external_key]
        v = v.cat.codes.to_numpy(float) if isinstance(v.dtype, pd.CategoricalDtype) \
            else v.to_numpy(float)
        pt = adata.obs[use_time_prior].to_numpy(float)
        m = np.isfinite(v) & np.isfinite(pt)
        out["corr_externa"] = float(spearmanr(pt[m], v[m]).statistic)
        out["orden_externo"] = adata.obs.groupby(external_key, observed=True)[
            use_time_prior].mean()

    adata.uns.setdefault("trajectory", {})["paga"] = {
        "group_key": group_key, "use_time_prior": use_time_prior,
        "coherencia": coherencia, "n_grupos": len(cats),
    }

    if verbose:
        print(f"\n[paga] pseudotiempo medio por población:")
        print(orden.sort_values().round(3).to_string())
        print(f"\n[paga] aristas ordenadas por |flujo neto|:")
        print(aristas.head(15).round(4).to_string(index=False))
        print(f"\n[paga] coherencia interna: {coherencia:.2f} de las aristas "
              f"tienen el flujo hacia pseudotiempo creciente.")
        if np.isfinite(coherencia) and coherencia < 0.5:
            print(
                "[paga] AVISO SERIO: por debajo de 0.5. O la matriz de\n"
                "       transiciones viene con la convención transpuesta en tu\n"
                "       versión de scanpy/scVelo (entonces TODAS tus flechas\n"
                "       están al revés y basta con transponer 'transiciones'),\n"
                "       o el campo de velocidad y el pseudotiempo se contradicen\n"
                "       y no hay trayectoria que leer. Comprueba una arista a\n"
                "       mano contra la figura de velocity_embedding_grid antes\n"
                "       de seguir."
            )
        elif np.isfinite(coherencia) and coherencia < 0.8:
            print("[paga] coherencia mediocre: hay aristas donde PAGA y el "
                  "pseudotiempo dicen cosas distintas. Míralas una a una.")
        if "corr_externa" in out:
            print(f"\n[paga] Spearman(pseudotiempo, {external_key}) = "
                  f"{out['corr_externa']:+.3f}")
            print(out["orden_externo"].round(3).to_string())
        print("\n[paga] 'conectividad' es NO dirigida y sale del grafo kNN: "
              "dice adyacencia, no linaje. La dirección está en 'flujo_neto'.")
    return out


def paga_net_transition(adata: ad.AnnData, group_key: str,
                        pair: Sequence[str]) -> float:
    """
    Flujo neto entre dos poblaciones: T[a,b] - T[b,a].

    Caso particular de paga_report(); se mantiene porque compare_batch_correction
    lo usa como referencia comparable entre corridas (no depende de la raíz).
    """
    import scvelo as scv

    _parche_scvelo_igraph()
    a, b = pair
    if "paga" not in adata.uns or "transitions_confidence" not in adata.uns.get("paga", {}):
        scv.tl.paga(adata, groups=group_key)
    T = _dense(adata.uns["paga"]["transitions_confidence"])
    cats = _cats(adata, group_key)
    if a not in cats or b not in cats:
        raise ValueError(f"{pair} no están ambos en {group_key}: {cats}")
    i, j = cats.index(a), cats.index(b)
    return float(T[i, j] - T[j, i])


def paga_stability(
    adata: ad.AnnData,
    group_key: str,
    n_boot: int = 10,
    frac: float = 0.8,
    use_time_prior: str = "velocity_pseudotime",
    min_cells: int = 20,
    n_jobs: int = 4,
    seed: int = 0,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    ¿Qué aristas de PAGA sobreviven al remuestreo de células?

    Esta es la comprobación que hace que PAGA sea reportable. El grafo siempre
    sale; la pregunta es cuáles de sus aristas son una propiedad de los datos y
    cuáles dependen de las células concretas que tocaron.

    Devuelve, por par de poblaciones: el flujo neto en el objeto completo, la
    mediana y el rango intercuartílico entre repeticiones, y `estabilidad` =
    fracción de repeticiones en las que el signo coincide con el del objeto
    completo. Una arista con estabilidad < 0.8 no aguanta una figura.

    COSTE Y ALCANCE. Se reutilizan los ajustes de recover_dynamics del objeto
    completo y solo se recalculan velocity_graph y paga sobre cada submuestra.
    Es lo que hace esto asumible en tiempo, y significa que lo que se mide es la
    estabilidad DEL GRAFO Y DE LA PARTICIÓN, no la del ajuste cinético. El
    ajuste cinético es otra fuente de variabilidad que esta función no captura.
    """
    import scvelo as scv

    _parche_scvelo_igraph()
    for req in ("velocity", "Ms", "Mu"):
        if req not in adata.layers:
            raise ValueError(
                f"Falta layers['{req}']; ejecuta run_velocity() sobre este "
                f"objeto antes de medir estabilidad"
            )

    cats = _cats(adata, group_key)
    base = paga_report(adata, group_key, use_time_prior=use_time_prior,
                       min_cells=min_cells, verbose=False)
    neto_full = base["flujo_neto"]

    rng = np.random.default_rng(seed)
    n_sub = int(round(frac * adata.n_obs))
    muestras = []

    for k in range(n_boot):
        idx = rng.choice(adata.n_obs, size=n_sub, replace=False)
        sub = adata[np.sort(idx)].copy()
        sub.obs[group_key] = sub.obs[group_key].cat.remove_unused_categories()
        if len(sub.obs[group_key].cat.categories) < len(cats):
            faltan = set(cats) - set(sub.obs[group_key].cat.categories)
            warnings.warn(f"repetición {k}: sin células de {faltan}")
        try:
            scv.tl.velocity_graph(sub, n_jobs=n_jobs)
            scv.tl.paga(sub, groups=group_key, use_time_prior=use_time_prior)
            Tk = _dense(sub.uns["paga"]["transitions_confidence"])
            ck = list(sub.obs[group_key].cat.categories)
            Tk = pd.DataFrame(Tk, index=ck, columns=ck)
            muestras.append(Tk - Tk.T)
        except Exception as e:
            warnings.warn(f"repetición {k} falló: {e}")
        del sub
        if verbose:
            print(f"[estab] {k + 1}/{n_boot}", end="\r")

    if not muestras:
        raise RuntimeError("Ninguna repetición completó")

    filas = []
    for i, a in enumerate(cats):
        for b in cats[i + 1:]:
            vals = np.array([
                float(m.loc[a, b]) for m in muestras
                if a in m.index and b in m.index
            ])
            vals = vals[np.isfinite(vals)]
            ref = float(neto_full.loc[a, b])
            if vals.size == 0:
                est = np.nan
            elif abs(ref) < 1e-8:
                est = np.nan
            else:
                est = float((np.sign(vals) == np.sign(ref)).mean())
            filas.append({
                "a": a, "b": b,
                "flujo_completo": ref,
                "mediana_boot": float(np.median(vals)) if vals.size else np.nan,
                "iqr_boot": float(np.subtract(*np.percentile(vals, [75, 25])))
                            if vals.size > 3 else np.nan,
                "estabilidad": est,
                "n_rep": int(vals.size),
            })

    out = pd.DataFrame(filas)
    out = out.reindex(
        out["flujo_completo"].abs().sort_values(ascending=False).index
    ).reset_index(drop=True)

    if verbose:
        print(" " * 30, end="\r")
        print(out.round(4).to_string(index=False))
        firmes = out[(out["estabilidad"] >= 0.8) & out["estabilidad"].notna()]
        print(f"\n[estab] {len(firmes)} de {len(out)} aristas conservan el signo "
              f"en >=80% de {n_boot} remuestreos al {int(100*frac)}%.")
        if len(firmes):
            print("[estab] firmes: " + ", ".join(
                f"{r.a}-{r.b}" for r in firmes.itertuples()))
        flojas = out[(out["estabilidad"] < 0.8) & out["estabilidad"].notna()]
        if len(flojas):
            print("[estab] NO firmes: " + ", ".join(
                f"{r.a}-{r.b} ({r.estabilidad:.2f})" for r in flojas.itertuples()))
        print("[estab] Recuerda: esto mide estabilidad del grafo y la "
              "partición con los ajustes cinéticos fijos, no del ajuste.")
    return out


# --------------------------------------------------------------------------
# comparación de variantes de corrección de batch
# --------------------------------------------------------------------------

def compare_batch_correction(
    adata: ad.AnnData,
    sample_key: str = "gsm",
    group_key: str = "minor_population",
    pair: Optional[Sequence[str]] = None,
    use_rep_joint: str = "X_harmony",
    ambient_max: Optional[float] = 0.2,
    reprocess_fn=None,
    seed: int = 0,
    verbose: bool = True,
):
    """
    Corre las variantes y las compara numéricamente sobre las MISMAS células.

        per_<muestra>   cada muestra sola, sin integrar. Sin batch que
                        corregir, así que es la referencia.
        1bc             conjunto, harmony solo en la representación del grafo.
        2bc             conjunto, harmony también en S y U.

    Criterio: si las muestras por separado coinciden entre sí, esa es la
    dirección buena. La variante conjunta que la reproduzca es la que usas.
    Si las muestras por separado NO coinciden, la comparación no decide nada
    y el problema es otro.

    reprocess_fn
        Callable(adata_sub, integrate: bool) -> adata_sub, para recalcular
        HVGs/PCA/harmony/vecinos en cada variante. Normalmente
        characterization.subset_and_reprocess envuelto, o processing.build_embeddings.
        Si es None, se asume que el objeto ya trae X_pca/X_harmony.

    Devuelve (pseudotiempos, resumen).
    """
    from scipy.stats import spearmanr

    muestras = list(adata.obs[sample_key].astype(str).unique())
    if len(muestras) < 2:
        raise ValueError(
            f"Solo una categoría en '{sample_key}': no hay nada que comparar"
        )
    if pair is None:
        raise ValueError(
            "pair es obligatorio. El criterio de esta función es '¿reproduce la "
            "variante conjunta lo que dicen las muestras por separado?', y las "
            "corridas por muestra NO comparten ninguna célula: su correlación de "
            "pseudotiempo es NaN por construcción. La única referencia común es "
            "el flujo neto entre dos poblaciones, que además no depende de la "
            "raíz. Pasa pair=('PobA','PobB')."
        )

    variantes = {}
    for m in muestras:
        variantes[f"per_{m}"] = dict(
            mask=(adata.obs[sample_key].astype(str) == m).to_numpy(),
            use_rep="X_pca", correct_su=False, integrate=False)
    variantes["1bc"] = dict(mask=np.ones(adata.n_obs, bool),
                            use_rep=use_rep_joint, correct_su=False, integrate=True)
    variantes["2bc"] = dict(mask=np.ones(adata.n_obs, bool),
                            use_rep=use_rep_joint, correct_su=True, integrate=True)

    pt, conf, trans = {}, {}, {}
    for nombre, cfg in variantes.items():
        if verbose:
            print(f"\n{'='*60}\n{nombre}  ({int(cfg['mask'].sum())} células)\n{'='*60}")
        sub = adata[cfg["mask"]].copy()
        if reprocess_fn is not None:
            sub = reprocess_fn(sub, cfg["integrate"])
        run_velocity(sub, use_rep=cfg["use_rep"], batch_key=sample_key,
                     correct_su=cfg["correct_su"], ambient_max=ambient_max,
                     seed=seed, verbose=verbose)
        pt[nombre] = pd.Series(sub.obs["velocity_pseudotime"].to_numpy(),
                               index=sub.obs_names)
        conf[nombre] = float(sub.obs["velocity_confidence"].mean())
        if pair is not None:
            try:
                trans[nombre] = paga_net_transition(sub, group_key, pair)
            except Exception as e:      # una población puede faltar en un subset
                trans[nombre] = np.nan
                warnings.warn(f"{nombre}: PAGA falló ({e})")
        del sub

    pseudotiempos = pd.DataFrame(pt)

    nombres = list(pseudotiempos.columns)
    corr = pd.DataFrame(np.nan, index=nombres, columns=nombres, dtype=float)
    for i in nombres:
        for j in nombres:
            if i == j:
                corr.loc[i, j] = 1.0
                continue
            # concat con claves nuevas: pseudotiempos[[i, j]] con i == j daría
            # un DataFrame de columnas duplicadas y spearmanr devolvería matriz
            comun = pd.concat([pseudotiempos[i], pseudotiempos[j]],
                              axis=1, keys=["a", "b"]).dropna()
            if len(comun) > 20:
                corr.loc[i, j] = float(spearmanr(comun["a"], comun["b"]).statistic)

    resumen = pd.DataFrame({
        "confianza_media": pd.Series(conf),
        "flujo_neto": pd.Series(trans) if trans else np.nan,
    })
    per = [c for c in nombres if c.startswith("per_")]
    # OJO: el signo de velocity_pseudotime es ARBITRARIO en cada corrida —
    # depende de qué célula acabe elegida como raíz. Promediar correlaciones
    # con signo hace que +0.9 y -0.8 se cancelen y salga ~0, que es justo lo
    # contrario de lo que pasa (las dos son fuertes). Se promedia |r|.
    resumen["corr_con_per_muestra"] = [
        corr.loc[n, per].drop(index=n, errors="ignore").abs().mean()
        for n in nombres
    ]
    resumen["corr_min_abs"] = [
        corr.loc[n, per].drop(index=n, errors="ignore").abs().min()
        for n in nombres
    ]

    if verbose:
        print("\n" + "="*60)
        print("correlación de Spearman del velocity_pseudotime")
        print("="*60)
        print(corr.round(3).to_string())
        print("\n" + resumen.round(3).to_string())
        if len(per) >= 2:
            # OJO: las corridas por muestra no comparten NINGUNA célula, así que
            # la correlación de pseudotiempo entre ellas es NaN por construcción.
            # Para saber si coinciden hay que mirar algo definido a nivel de
            # POBLACIÓN, que sí es común: el flujo neto del par.
            flujos = pd.Series({k: v for k, v in trans.items() if k in per}) \
                if trans else pd.Series(dtype=float)
            if flujos.notna().sum() >= 2:
                de_acuerdo = bool(np.all(np.sign(flujos.dropna()) == np.sign(flujos.dropna().iloc[0]))
                                  and (flujos.dropna().abs() > 1e-6).all())
                print(f"\n[velo] flujo neto por muestra: "
                      f"{flujos.round(4).to_dict()}")
                print(f"[velo] las muestras por separado "
                      f"{'COINCIDEN' if de_acuerdo else 'NO coinciden'} en la dirección.")
                if not de_acuerdo:
                    print("[velo] Sin una referencia consistente esta comparación "
                          "no decide nada: el problema no es la corrección de "
                          "batch. Revisa confianza, número de genes de velocity "
                          "y si el subconjunto tiene estructura suficiente.")
                    return pseudotiempos, resumen
            elif pair is None:
                print("\n[velo] Sin 'pair' no puedo comparar las muestras entre "
                      "sí: no comparten células, así que la correlación de "
                      "pseudotiempo entre ellas es NaN por construcción. Pasa "
                      "pair=('PobA','PobB') para tener la referencia.")
            if True:
                cand = resumen.drop(index=per)["corr_con_per_muestra"]
                mejor = cand.idxmax()
                print(f"[velo] variante más parecida a las muestras sueltas: "
                      f"{mejor} ({cand[mejor]:+.3f})")
                print("[velo] corr_con_per_muestra es |r|: mide si la "
                      "variante conjunta reproduce el ORDEN de las muestras "
                      "sueltas, no su dirección. La dirección la da flujo_neto.")
        if pair is not None and trans:
            print(f"\n[velo] flujo neto {pair[0]} -> {pair[1]} "
                  f"(positivo = de {pair[0]} a {pair[1]}):")
            print(pd.Series(trans).round(4).to_string())
    return pseudotiempos, resumen


# --------------------------------------------------------------------------
# PHATE: embedding para MIRAR
# --------------------------------------------------------------------------

def run_phate(
    adata: ad.AnnData,
    use_rep: str = "X_harmony",
    key_added: str = "X_phate",
    n_components: int = 2,
    knn: Optional[int] = None,
    decay: int = 40,
    t: str | int = "auto",
    n_jobs: int = 1,
    seed: int = 0,
    verbose: bool = True,
    **kwargs,
) -> ad.AnnData:
    """
    PHATE sobre la MISMA representación que alimenta el resto del pipeline.

    PHATE es un EMBEDDING, no un método de inferencia. Sirve para mirar: a
    diferencia de UMAP conserva mejor las distancias de transición y las
    estructuras alargadas, que es justo lo que interesa cuando sospechas un
    continuo. Pero las coordenadas que produce NO deben alimentar ningún
    pseudotiempo ni ninguna distancia: eso sería reducir a dos dimensiones y
    después medir en ellas, que es el error de `ndims_rep=2` en Palantir. La
    inferencia va sobre el grafo o sobre el mapa de difusión completos.

    use_rep
        POR QUÉ ESTO IMPORTA MÁS QUE NADA: `phate.PHATE().fit_transform(adata.X)`
        —que es lo que sale en todos los tutoriales— corre sobre la matriz de
        expresión y por tanto IGNORA harmony. Con varias muestras, el lote
        vuelve a aparecer y te encuentras dos nubes paralelas que parecen
        biología. Aquí se corre sobre `use_rep` (X_harmony por defecto), que ya
        está integrado.

        Como esa representación YA es un espacio reducido, se pasa `n_pca=None`
        para que PHATE no le haga otro PCA encima. Si le pasas una matriz de
        expresión cruda en `use_rep="X"`, entonces sí conviene dejar que haga
        su PCA: pásale n_pca=100 por kwargs.

    knn
        Si no lo das, se toma el mismo k del grafo de vecinos ya construido
        (uns['neighbors']), para que PHATE y el grafo miren el tejido a la
        misma escala. Que dos métodos usen vecindarios de tamaño distinto es
        una fuente clásica de "PHATE dice una cosa y PAGA otra".

    t
        La escala de difusión, y el parámetro que más cambia la figura. 'auto'
        la elige por el codo de la entropía de von Neumann. Queda registrada en
        uns['trajectory']['phate'] porque es una decisión, no una constante.

    n_jobs
        1 por defecto, por reproducibilidad. Ya nos ha pasado con harmony: los
        algoritmos que reducen en paralelo no dan bit a bit lo mismo entre
        ejecuciones. Súbelo si prefieres velocidad y no te importa que el
        embedding se mueva un poco.
    """
    try:
        import phate
    except ImportError as e:
        raise ImportError("pip install phate") from e

    if use_rep == "X":
        X = adata.X
        X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
        n_pca = kwargs.pop("n_pca", 100)
        if verbose:
            print("[phate] corriendo sobre X (expresión). OJO: esto NO está "
                  "integrado; si tienes más de una muestra, el lote va a "
                  "aparecer en el embedding.")
    else:
        if use_rep not in adata.obsm:
            raise ValueError(
                f"No hay obsm['{use_rep}']. Ejecuta build_embeddings() antes, o "
                f"pasa use_rep='X' si de verdad quieres correr sobre expresión."
            )
        X = np.asarray(adata.obsm[use_rep])
        n_pca = kwargs.pop("n_pca", None)      # ya es un espacio reducido

    if knn is None:
        knn = adata.uns.get("neighbors", {}).get("params", {}).get("n_neighbors")
        if knn is None:
            knn = max(5, int(0.5 * adata.n_obs ** 0.5))
        if verbose:
            print(f"[phate] knn={knn} (heredado del grafo de vecinos)")

    op = phate.PHATE(
        n_components=n_components, knn=knn, decay=decay, t=t,
        n_pca=n_pca, n_jobs=n_jobs, random_state=seed,
        verbose=1 if verbose else 0, **kwargs,
    )

    if n_jobs == 1:
        try:
            from threadpoolctl import threadpool_limits
        except ImportError:
            Y = op.fit_transform(X)
        else:
            with threadpool_limits(limits=1):
                Y = op.fit_transform(X)
    else:
        Y = op.fit_transform(X)

    adata.obsm[key_added] = np.ascontiguousarray(Y, dtype=np.float32)

    t_usado = getattr(op, "t", t)
    adata.uns.setdefault("trajectory", {})["phate"] = {
        "use_rep": use_rep, "knn": int(knn), "decay": decay,
        "t": int(t_usado) if isinstance(t_usado, (int, np.integer)) else str(t_usado),
        "n_pca": n_pca, "n_components": n_components, "seed": seed,
        "n_jobs": n_jobs,
    }

    if verbose:
        print(f"[phate] obsm['{key_added}'] listo; t={t_usado}, knn={knn}, "
              f"sobre '{use_rep}'")
        pt = "velocity_pseudotime"
        if pt in adata.obs and n_components >= 1:
            from scipy.stats import spearmanr
            v = adata.obs[pt].to_numpy(float)
            m = np.isfinite(v)
            r = float(spearmanr(Y[m, 0], v[m]).statistic)
            print(f"[phate] Spearman(PHATE1, {pt}) = {r:+.3f}. Es DESCRIPTIVO: "
                  f"dice si el eje 1 del dibujo coincide con el pseudotiempo, "
                  f"no valida ninguno de los dos.")
        print("[phate] Dibuja con sc.pl.embedding(adata, basis='phate', "
              "color=...). No uses estas coordenadas para calcular distancias "
              "ni pseudotiempos.")
    return adata


# --------------------------------------------------------------------------
# dibujo
# --------------------------------------------------------------------------

def plot_velocity(adata: ad.AnnData, basis: str = "umap",
                  color: Optional[str] = None, min_mass: float = 2.0,
                  autoscale: bool = False, **kwargs):
    """
    velocity_embedding_grid con los dos defaults que importan.

    min_mass
        velocity_embedding_grid dibuja una flecha en cada celda de la rejilla,
        también donde apenas hay células. Esas flechas son extrapolación, no
        datos, y son las más largas y llamativas de la figura.
    autoscale
        True normaliza la longitud por figura, así que dos figuras dejan de ser
        comparables en magnitud. Para comparar variantes, False.
    """
    import scvelo as scv
    return scv.pl.velocity_embedding_grid(
        adata, basis=basis, color=color, min_mass=min_mass,
        autoscale=autoscale, **kwargs)