import scanpy as sc
import seaborn as sn
from matplotlib import font_manager, rcParams
import seaborn as sns
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
import numpy as np
import pandas as pd
import os
from pathlib import Path
from urllib.request import urlretrieve
import matplotlib as mpl
import matplotlib.pyplot as plt
from adjustText import adjust_text

from .common_vars import BASE_DIR


palette = {
    # ---- Major (igual que antes) ----
    'Endothelial': "#E766E7",
    'FAP':         "#165A81",
    'Immune':      "#6933AF",
    'Kranocyte':   "#C23C64",
    'Satellite':   "#289C79",
    'Tenocyte':    "#D38E28",

    # ---- Minor (más oscuros) ----
    # FAP (7)
    'FAP_A':"#C7CDD3",
    'FAP_AB':"#899199",
    'FAP_AF':"#565D63",

    'FAP_B1':"#6E61BB",
    'FAP_B2':"#3C2F88",

    'FAP_C' :"#6BB6D3",
    'FAP_D' :"#3C85A2",
    'FAP_E' :"#175274",
    'FAP_F' :"#05293C",

    # Tenocyte (4)
    'Teno_A':"#EABD62",
    'Teno_B':"#C6993A",
    'Teno_C':"#6F5214",
    'Teno_D':"#392705",

    # Kranocyte (3)
    'Krano_A':"#E46283",
    'Krano_B':"#9C2B60",
    'Krano_C':"#550D29",

    # Satellite (3)
    'Sat_A1':"#67DAB7",
    'Sat_A2':"#2E9B76",
    'Sat_B' :"#084233",
}





FONT_BASE_URL = (
    "https://raw.githubusercontent.com/alexmascension/myplotfonts/"
    "main/TeXGyreHeros/"
)

FONT_FILES = {
    "regular": "texgyreheros-regular.otf",
    "bold": "texgyreheros-bold.otf",
    "italic": "texgyreheros-italic.otf",
    "bolditalic": "texgyreheros-bolditalic.otf",
}



magma = [plt.get_cmap('magma')(i) for i in np.linspace(0,1, 80)]
magma[0] = (0.88, 0.88, 0.88, 1)
magma = mpl.colors.LinearSegmentedColormap.from_list("", magma[:65])




def set_plotting_style():
    font_dir = Path(BASE_DIR) / ".cache" / "fonts"
    font_dir.mkdir(parents=True, exist_ok=True)

    font_paths = {}

    for style, filename in FONT_FILES.items():
        font_path = font_dir / filename

        if not font_path.exists():
            urlretrieve(FONT_BASE_URL + filename, font_path)

        font_manager.fontManager.addfont(str(font_path))
        font_paths[style] = font_path

    regular_font = font_manager.FontProperties(
        fname=str(font_paths["regular"])
    )
    bold_font = font_manager.FontProperties(
        fname=str(font_paths["bold"])
    )
    italic_font = font_manager.FontProperties(
        fname=str(font_paths["italic"])
    )
    bolditalic_font = font_manager.FontProperties(
        fname=str(font_paths["bolditalic"])
    )

    family = regular_font.get_name()

    sc.set_figure_params(dpi=250)
    sns.set_style("white")

    rcParams.update({
        # Main Matplotlib font family
        "font.family": family,
        "font.sans-serif": [family],
        "font.cursive": [family],

        # MathText
        "mathtext.fontset": "custom",
        "mathtext.rm": family,
        "mathtext.it": f"{family}:italic",
        "mathtext.bf": f"{family}:bold",
        "mathtext.bfit": f"{family}:bold:italic",
        "mathtext.cal": family,
    })

    return {
        "regular": regular_font,
        "bold": bold_font,
        "italic": italic_font,
        "bolditalic": bolditalic_font,
    }


def savefig(fig, filename, fig_dir="../figures/", dpi=600, bbox_inches='tight'):
    filename_safe = filename.replace(' ', '-').replace('/', '|')
    os.makedirs(fig_dir, exist_ok=True)

    fig.savefig(fig_dir + filename_safe + '.png', dpi=dpi, bbox_inches=bbox_inches)
    fig.savefig(fig_dir + filename_safe + '.svg', bbox_inches=bbox_inches, dpi=dpi)
    fig.savefig(fig_dir + filename_safe + '.pdf', bbox_inches=bbox_inches, dpi=dpi)



def plot_cell_stats(adata):
    # Show cell 10X plot (it'll be truncated, but its good enough) and relationship between this and cell probabilities
    fig, axs = plt.subplots(1, 3, figsize=(15,5))

    gene_counts = np.sort(adata.obs['counts_cellbender'])[::-1]
    axs[0].plot(np.log10(np.arange(adata.shape[0]) + 1), np.log10(gene_counts + 1))

    gene_counts_filter = gene_counts[gene_counts > 200]
    axs[0].plot(np.log10(np.arange(len(gene_counts_filter)) + 1), np.log10(gene_counts_filter + 1))
    axs[0].set_xlabel('Cell rank')
    axs[0].set_ylabel('CellBender counts (log1p)')

    axs[1].scatter(adata.obs['cell_probability'], np.log10(adata.obs['counts_cellbender'] + 1), s=2)
    axs[1].set_xlabel('Cell probability')
    axs[1].set_ylabel('CellBender counts (log1p)')

    axs[2].scatter(np.log10(adata.obs['counts_cellbender'] + 1), np.log10(adata.obs['counts_raw'] + 1), s=2)
    axs[2].plot([0, 5], [0, 5], c='#bc0000')
    axs[2].set_xlabel('CellBender counts (cells) - log10')
    axs[2].set_ylabel('Raw counts (cells) -')


    plt.tight_layout()

def plot_gene_stats(adata):
    # Show ambient gene expression
    display(adata.var.sort_values(by='ambient_expression', ascending=False).head(15))

    fig, axs = plt.subplots(1, 3, figsize=(15,5))
    

    axs[0].scatter(np.log10(adata.var['counts_cellbender'] + 1), adata.var['ambient_expression'], s=2)
    axs[0].set_xlabel('CellBender counts (genes)')
    axs[0].set_ylabel('Raw counts (log1p)')

    axs[1].scatter(np.log10(adata.var['counts_cellbender'] + 1), (adata.var['counts_cellbender'] + 1) / (adata.var['counts_raw'] + 1),  s=2)
    axs[1].set_xlabel('CellBender counts (genes)')
    axs[1].set_ylabel('$\\frac{\\text{CellBender counts  + 1}}{\\text{Raw counts + 1}}$ (genes)')

    axs[2].scatter(np.log10(adata.var['counts_cellbender'] + 1), np.log10(adata.var['counts_raw'] + 1), s=2)
    axs[2].plot([0, 6], [0, 6], c='#bc0000')
    axs[2].set_xlabel('CellBender counts (genes)')
    axs[2].set_ylabel('Raw counts (genes)')
    plt.tight_layout()





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