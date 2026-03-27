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


SEED = 0

magma = [plt.get_cmap('magma')(i) for i in np.linspace(0,1, 80)]
magma[0] = (0.88, 0.88, 0.88, 1)
magma = mpl.colors.LinearSegmentedColormap.from_list("", magma[:65])



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



def set_plotting_style():
    # Set font and plot properties
    sns.set_style("white")

    font_path = "src/fonts/texgyreheros-regular.otf" # Alternative
    font_path = "/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreheros-regular.otf"

    mpl.rcParams.update({'font.size': 22})


    # Path to the Helvetica font file
    custom_font = font_manager.FontProperties(fname=font_path)

    # Add the font to Matplotlib's font manager
    font_manager.fontManager.addfont(font_path)

    # Set the font globally
    rcParams['font.family'] = custom_font.get_name()
    rcParams['font.size'] = 16
    rcParams['figure.dpi'] = 300



def plot_volcano(adata, cluster, pval_threshold=0.0001, lfc_threshold=2, topn=10, bottomn=8, zero_pval='auto', xlim=None, return_df=False, plot_positive_only=True):
    df_pvals = pd.DataFrame({'gene': adata.uns['rank_genes_groups']['names'][cluster], 
                         'adj_pval': adata.uns['rank_genes_groups']['pvals_adj'][cluster], 
                         'logfoldchanges': adata.uns['rank_genes_groups']['logfoldchanges'][cluster],})

    # adjust p values that are zero to a small number
    if zero_pval == 'auto':
        min_nonzero = df_pvals.loc[df_pvals['adj_pval'] > 0, 'adj_pval'].min()
        zero_pval = min_nonzero * 0.1

    df_pvals.loc[df_pvals['adj_pval'] == 0, 'adj_pval'] = zero_pval

    df_pvals['neg_log_pval'] = -np.log10(df_pvals['adj_pval'])

    if plot_positive_only:
        df_pvals = df_pvals[df_pvals['logfoldchanges'] > 0]

    plt.figure(figsize=(5, 4))
    plt.scatter(df_pvals['logfoldchanges'], df_pvals['neg_log_pval'], color='gray', alpha=0.5, s=3)
    
    if xlim is not None:
        df_pvals = df_pvals[df_pvals['logfoldchanges'] <= xlim[1]]
    
    significant = (df_pvals['adj_pval'] < pval_threshold) & (abs(df_pvals['logfoldchanges']) > lfc_threshold)


    plt.scatter(df_pvals['logfoldchanges'][significant], df_pvals['neg_log_pval'][significant], color='red', alpha=0.7, s=3)

    if xlim is not None:
        plt.gca().set_xlim(xlim)
        
    df_pvals['pvalxlfc'] = df_pvals['neg_log_pval'] * df_pvals['logfoldchanges']
    top_genes = df_pvals.nlargest(topn, 'pvalxlfc')
    bottom_genes = df_pvals.nsmallest(bottomn, 'pvalxlfc')

    texts = []
    for _, row in pd.concat([top_genes, bottom_genes]).iterrows():
        texts.append(plt.text(row['logfoldchanges'], row['neg_log_pval'], row['gene'], 
                            fontsize=7, ha='center', color='black', 
                            )) 

    adjust_text(texts, 
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.8)) 

    plt.axhline(y=-np.log10(pval_threshold), color='#232323', linestyle='--', linewidth=0.8)
    plt.axvline(x=lfc_threshold, color='#232323', linestyle='--', linewidth=0.8)

    if not plot_positive_only:
        plt.axvline(x=-lfc_threshold, color='#232323', linestyle='--', linewidth=0.8)

    plt.xlabel('LFC')
    plt.ylabel('-log$_{10}$(Adjusted p-value)')
    plt.title('')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)


    # plt.savefig('../../../figures/4E_volcano_glia.png', dpi=300, bbox_inches='tight')

    plt.show()

    if return_df:
        return df_pvals