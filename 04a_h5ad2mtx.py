"""Convert h5ad data to mtx format.
"""


import scanpy as sc

import utils.io as io


if __name__ == '__main__':

    data_file = './outs/adata/concatenated.h5ad'
    out_dir = './outs/mtx/concatenated/'

    adata = sc.read_h5ad(data_file)
    io.h5ad2mtx(adata, out_dir)

    print('Done!')
