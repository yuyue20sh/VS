"""QC Strategy:

1. Cells with < 500 genes were excluded.
2. Genes expressed in < 0 cells were exclude (not filtering genes before merging).
3. Percentage of counts in mitochondrial genes was calculated in each cell.
4. Doublet scores were calculated using Scrublet.
5. Leiden clustering was done.
6. Cells with percentage of counts in mitochondrial genes greater than 2
   standard deviations from the mean percentage of counts in mitochondrial
   genes of their leiden clusters were removed.
7. The predicted doublets were removed.
8. Leiden clusters with mean doublet score greater than 2 standard deviations
   from the mean doublet score of all cells the were excluded.

"""


import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import logging
from tqdm import tqdm
import warnings


def set_logger(name, log_file, level=logging.INFO, sh_level=logging.DEBUG, fh_level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'):
    """
    Set logger with StreamHandler and FileHandler

    Args:
        name: str, name of the logger
        log_file: str, path to the log file
        level: int, logging level
        sh_level: int, logging level for StreamHandler
        fh_level: int, logging level for FileHandler
        format: str, logging format

    Returns:
        logger: logger object

    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    formatter = logging.Formatter(format)

    sh = logging.StreamHandler()
    sh.setLevel(sh_level)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    fh = logging.FileHandler(log_file)
    fh.setLevel(fh_level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    return logger


def plot_figs(adata, fig_dir, stage):
    """
    Plot figures for quality control

    Args:
        adata: AnnData object, adata to plot
        fig_dir: str, path to the directory to save figures
        stage: str, stage of the quality control, either "before" or "after"

    Returns:
        None

    """
    assert stage in ['before', 'after'], 'stage must be either "before" or "after"'
    sc.pl.violin(
        adata,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'doublet_score'],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    plt.savefig('%s/%s_violin.png' % (fig_dir, stage))
    sc.pl.scatter(
        adata,
        'total_counts', 'n_genes_by_counts',
        color='pct_counts_mt',
        show=False,
    )
    plt.savefig('%s/%s_scatter.png' % (fig_dir, stage))
    sc.pl.umap(
        adata,
        color=['leiden', 'log1p_n_genes_by_counts', 'log1p_total_counts',
            'pct_counts_mt', 'doublet_score', 'predicted_doublet'],
        wspace=0.5,
        ncols=3,
        show=False,
    )
    plt.savefig('%s/%s_umap.png' % (fig_dir, stage))
    plt.close('all')

    return


def qc(adata, logger, fig_dir, min_genes=500, min_cells=0):
    """
    Quality control for raw adata

    Args:
        adata: AnnData object, raw adata
        logger: logger object
        fig_dir: str, path to the directory to save figures
        min_genes: int, minimum number of genes expressed in a cell
        min_cells: int, minimum number of cells expressing a gene

    Returns:
        adata: AnnData object, adata after quality control

    """
    warnings.filterwarnings('ignore', message='are not unique', category=UserWarning)
    Path(fig_dir).mkdir(parents=True, exist_ok=True)

    shape = adata.shape
    logger.info('raw: n_obs x n_vars = %d x %d' % shape)

    # filter cells and genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    logger.info('filter cells: n=%d' % (shape[0] - adata.shape[0]))
    logger.info('filter genes: n=%d' % (shape[1] - adata.shape[1]))
    shape = adata.shape

    # calculate qc metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars='mt', inplace=True, log1p=True)

    # run scrublet
    logger.info('run scrublet... (this may take some time)')
    sc.pp.scrublet(adata)
    logger.info('number of predicted doublets: %d' % sum(adata.obs['predicted_doublet']))

    # run leiden clustering
    logger.info('run leiden clustering...')
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, flavor='igraph')
    adata.X = adata.layers['counts'].copy()
    # free memory
    del adata.layers

    # plot before fileration
    plot_figs(adata, fig_dir, 'before')
    sc.pl.violin(
        adata,
        'doublet_score',
        groupby='leiden',
        stripplot=False,
        inner='box',
        show=False,
    )
    plt.savefig('%s/cluster_doublet.png' % fig_dir)

    # assess qc metrics
    # cells with too much mitochondrial genes
    cells_with_too_much_mt_genes = []
    for c in adata.obs['leiden'].cat.categories:
        pct_counts_mt = adata.obs[adata.obs['leiden'] == c]['pct_counts_mt']
        top = pct_counts_mt.mean() + 2 * pct_counts_mt.std()
        cells = list(pct_counts_mt[pct_counts_mt > top].index)
        cells_with_too_much_mt_genes += cells
    logger.info('cells with too much mt genes: n=%d' % len(cells_with_too_much_mt_genes))

    # clusters with high doublet scores
    clusters_with_high_doublet_score = (adata.obs.groupby('leiden', observed=True)['doublet_score'].mean() >
                                        adata.obs['doublet_score'].mean() + 2 * adata.obs['doublet_score'].std())  # greater than mean + 2 std
    clusters_with_high_doublet_score = list(clusters_with_high_doublet_score[clusters_with_high_doublet_score].index)
    logger.info('clusters with high doublet score: %s' % clusters_with_high_doublet_score)

    # filter cells
    adata = adata[~adata.obs.index.isin(cells_with_too_much_mt_genes)]
    adata = adata[~adata.obs['predicted_doublet']]
    adata = adata[~adata.obs['leiden'].isin(clusters_with_high_doublet_score)]
    logger.info('filter cells: n=%d' % (shape[0] - adata.shape[0]))
    shape = adata.shape
    logger.info('after qc: n_obs x n_vars = %d x %d' % shape)

    # plot after fileration
    plot_figs(adata, fig_dir, 'after')

    # remove unneeded elements
    del adata.uns
    del adata.obsm
    del adata.varm
    del adata.obsp
    adata.obs = adata.obs.drop('leiden', axis=1)

    plt.close('all')

    return adata


if __name__ == "__main__":

    data_dir = './outs/adata/raw/'
    out_dir = './outs/adata/qc/'
    log_file = './logs/qc.log'
    fig_dir = './outs/figs/qc/'


    # load data
    print('loading data...')
    adatas = {}
    for f in tqdm(Path(data_dir).glob('./*.h5ad'), total=len(list(Path(data_dir).glob('./*.h5ad')))):
        id = f.stem
        adatas[id] = sc.read_h5ad(f)

    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    logger = set_logger('qc', log_file)
    for id, adata in adatas.items():
        logger.info('quality control for %s' % id)
        adata = qc(adata, logger, fig_dir='%s/%s/' % (fig_dir, id), min_genes=500, min_cells=0)
        adata.write('%s/%s.h5ad' % (out_dir, id))

    print('Done!')
