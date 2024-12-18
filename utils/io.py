import scanpy as sc
import pandas as pd
import scipy.io as sio

from pathlib import Path
import warnings


def read_mtx(mtx_file):
    """
    Read mtx format data.

    Args:
        mtx_file: str, file path

    Returns:
        adata: AnnData, scRNA data

    """
    adata = sc.read_mtx(mtx_file).T
    
    # get file names
    filename = Path(mtx_file).name
    parent_dir = Path(mtx_file).parent
    barcodes_file = '%s/%s' % (parent_dir, filename.replace('matrix', 'barcodes').replace('mtx', 'tsv'))
    features_file = '%s/%s' % (parent_dir, filename.replace('matrix', 'features').replace('mtx', 'tsv'))
    metadata_file = '%s/%s' % (parent_dir, filename.replace('matrix', 'metadata').replace('mtx', 'tsv'))

    barcodes = pd.read_csv(barcodes_file, sep='\t', header=None)
    features = pd.read_csv(features_file, sep='\t', header=None)

    # add metadata
    adata.obs_names = barcodes[0]
    if features.shape[1] > 1:
        adata.var_names = features[1]
        adata.var['ens_id'] = features[0].to_list()
        adata.var['symbol'] = features[1].to_list()
    else:
        adata.var_names = features[0]
        adata.var['ens_id'] = '-'
        adata.var['symbol'] = features[0].to_list()
    
    if Path(metadata_file).exists():
        metadata = pd.read_csv(metadata_file, sep='\t', index_col=0, header=0)
        adata.obs = metadata.loc[adata.obs_names]

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    return adata


def read_h5(h5_file):
    """
    Read h5 format data.

    Args:
        h5_file: str, file path

    Returns:
        adata: AnnData, scRNA data

    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = sc.read_10x_h5(h5_file)

    # select columns
    adata.var = adata.var.iloc[:, 0].to_frame()
    adata.var.columns = ['ens_id']
    adata.var['symbol'] = adata.var_names

    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    
    return adata


def h5ad2mtx(adata, out_dir):
    """
    Convert anndata object to mtx format data.
    output files: matrix.mtx, barcodes.tsv, features.tsv, metadata.tsv

    Args:
        adata: AnnData object, adata to convert
        out_dir: str, path to the output directory

    Returns:
        None

    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    matrix = adata.X.T
    barcodes = adata.obs.index.to_frame()
    features = adata.var.index.to_frame()
    metadata = adata.obs

    sio.mmwrite('%s/matrix.mtx' % out_dir, matrix)
    barcodes.to_csv('%s/barcodes.tsv' % out_dir, sep='\t', header=False, index=False)
    features.to_csv('%s/features.tsv' % out_dir, sep='\t', header=False, index=False)
    metadata.to_csv('%s/metadata.tsv' % out_dir, sep='\t', header=True, index=True)

    return
