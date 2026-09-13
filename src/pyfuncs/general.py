import scanpy as sc
import scanpy.external as sce
import triku as tk
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import seaborn as sns
import matplotlib as mpl
from matplotlib import font_manager, rcParams

import sys
sys.path.append('..')
from pyfuncs.common_vars import SEED



def preprocessing_adata_sub(adata_sub, integrate = True, k=None, n_comps=15):
    if k is None:
        k = int(0.5 * len(adata_sub) ** 0.5)

    if integrate:
        use_rep = 'X_harmony'
    else:
        use_rep = 'X_pca'

    sc.pp.filter_cells(adata_sub, min_counts=10)
    sc.pp.filter_genes(adata_sub, min_counts=20)
    sc.pp.pca(adata_sub, n_comps=n_comps, random_state=SEED, use_highly_variable=False)
    if integrate:
        sce.pp.harmony_integrate(adata_sub, 'batch', random_state=SEED,
                                    basis='X_pca', adjusted_basis='X_harmony', max_iter_harmony=30, verbose=False)
        
    sc.pp.neighbors(adata_sub, n_neighbors=k, random_state=SEED, metric='correlation', 
                    use_rep=use_rep)
    tk.tl.triku(adata_sub, use_raw=False)

    sc.pp.pca(adata_sub, n_comps=n_comps, random_state=SEED, use_highly_variable=True)
    if integrate:
        sce.pp.harmony_integrate(adata_sub, 'batch', random_state=SEED,
                                basis='X_pca', adjusted_basis='X_harmony', max_iter_harmony=30, verbose=False)
    sc.pp.neighbors(adata_sub, n_neighbors=k, random_state=SEED, metric='correlation', 
                    use_rep=use_rep)



