"""
Functions to compute de adjacency matrix of a given 4D fmri image.
"""
import numpy
import logging
import os
import pickle
from scipy import sparse
from scipy.io import mmwrite
from numpy.linalg import norm
import nibabel as nib
from operator import mul
import itertools as itt
from glob import glob

LOGGER = logging.getLogger(__name__)

def func2adj(func_file, mask_file=None, min_corr=0.7):
    """Computes the adjecency matrix for the functional image in func_file. The
    (i, j) element of adj is the correlation of the time serie in voxel
    node2voxel[i] and the time serie in voxel node2voxel[j]. "mask_file" is the
    mask file for "func_file". "min_corr" is the minimun correlation value that
    will be save to the adjacency (sparse) matrix.

    Parameters
    ----------
    func_file : string.
        Name of the file or folder with the functional image.
    mask_file : string or None.
        Name of the file with the mask to applay to the functional image.
    min_corr : float.
        Minumun correlation coeficient to save. It should be a real number
        greather than 0 an less than 1. The default min_corr is 0.7.

    Returns
    -------
    node2voxel : dictionary
        For each coloumn in the adjacency matrix it gives you the corresponding
        coordinates in the volumen matrix.
    adj : sparse matrix
        Adjacency matrix in spase format.

    """
    LOGGER.info('loading image')
    fimg = get_img(func_file)
    data = fimg.get_data()
    if mask_file != None:
        mask = nib.load(mask_file).get_data()
    else:
        mask = numpy.ones(data.shape[:3])
    if any([m != d for m, d in zip(data.shape, mask.shape)]):
        data = reshape_data2mask(get_img(func_file), mask, mask_file)
    imat = numpy.arange((reduce(mul, data.shape[:3], 1)))
    series = data.reshape(imat.shape[0], data.shape[3])
    imat = imat.reshape(*data.shape[:3])
    fil = (mask.reshape(series.shape[0]) != 0) & (norm(series, axis=1) != 0)
    series = series[fil, :]
    node2voxel = list(itt.product(*[xrange(d) for d in data.shape[:3]]))
    node2voxel = numpy.array(sorted(node2voxel, key=lambda x: imat[x]))
    node2voxel = {i: v for i, v in enumerate(node2voxel[fil])}

    LOGGER.info('number of nodes in scan: %i', series.shape[0])
    series = series - numpy.mean(series, axis=1, keepdims=True)
    series = series / norm(series, axis=1).reshape(series.shape[0], 1)
    shard_size = 1000
    adj = []
    for i in range(0, series.shape[0], shard_size):
        chunk = series[i:i+shard_size, :]
        if i % 10*shard_size == 0:
            LOGGER.info('nodes: %i of %i', shard_size*i, series.shape[0])
        adjnew = numpy.dot(chunk, series.T)
        rows_and_cols = numpy.where(adjnew > min_corr)
        values = [adjnew[i, j] for i, j in zip(*rows_and_cols)]
        adj.append(sparse.coo_matrix((values, rows_and_cols),
                                     shape=adjnew.shape,
                                     dtype=series.dtype))
    adj = sparse.vstack(adj)
    LOGGER.info('Filtering number_of_edges = %i', adj.nnz)
    adj = sparse.triu(adj)
    return adj, node2voxel

def reshape_data2mask(fimg, mask, mask_file):
    """Reshape functional image data to fit into the mask."""
    faff = fimg.get_affine()
    data = fimg.get_data()
    oldsh = data.shape
    LOGGER.info('Fill in out of bounding box')
    mimg = nib.load(mask_file)
    maff = mimg.get_affine()
    temp = numpy.linalg.inv(maff) - numpy.linalg.inv(faff)
    (i0, j0, k0) = temp.round().astype('int32')[0:3, 3]
    LOGGER.info('Shift (%i,%i,%i)', i0, j0, k0)
    newsh = (mask.shape[0], mask.shape[1], mask.shape[2], oldsh[3])
    rs_data = numpy.zeros(newsh)
    rs_data[i0:i0+oldsh[0], j0:j0+oldsh[1], k0:k0+oldsh[2]] = data
    return rs_data

def get_img(func_file):
    """It loads func_file functional image regardless if it es a folder or a
    file.
    """
    if os.path.isfile(func_file):
        LOGGER.info('one file image')
        img = nib.load(func_file)
    elif os.path.isdir(func_file):
        files = sorted(glob(os.path.join(func_file, '*.img')))
        LOGGER.info('multiple images. len(files)=%i', len(files))
        img = nib.funcs.concat_images(files)
    return img

if __name__ == '__main__':
    FUNC_FILE = glob('func*.nii*')[0]
    MASK_FILE = glob('MNI*.nii*')[0]
    ADJ, N2V = func2adj(FUNC_FILE, mask_file=MASK_FILE)
    mmwrite('adjacency.mtx', ADJ, precision=5)
    pickle.dump(N2V, open('node2voxel.pkl', 'w'))
