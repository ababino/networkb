import numpy
import scipy
from numpy.linalg import norm
import nibabel as nib
from operator import mul
import itertools as itt
import logging
import os
from glob import glob

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel('INFO')

def gen_edgelist(func_file, mask_file=None, min_th=0.7):
    """
    Gennerates a list of correlations among voxel's time series of file
    "filename". "mask" is the mask file for "filename. "th" is the minimun
    threshold that will be save to the file.
    """
    LOGGER.info('loading image')
    fimg = get_img(func_file)
    data = fimg.get_data()
    if mask_file != None:
        mask = nib.load(mask_file).get_data()
    else:
        mask = numpy.ones(data.shape[:3])
    if any([m != d for m, d in zip(data.shape, mask.shape)]):
        data = reshape_data2mask(data, mask, mask_file)
    imat = numpy.arange((reduce(mul, data.shape[:3], 1)))
    series = data.reshape(imat.shape[0], data.shape[3])
    imat = imat.reshape(*data.shape[:3])
    fil = (mask.reshape(series.shape[0]) != 0) & (norm(series, axis=1) != 0)
    series = series[fil, :]
    node2voxel = list(itt.product(*[xrange(d) for d in data.shape[:3]]))
    node2voxel = numpy.array(sorted(node2voxel, key=lambda x: imat[x]))
    node2voxel = node2voxel[fil]

    LOGGER.info('number of nodes in scan: %i', series.shape[0])
    series = series - numpy.mean(series, axis=1, keepdims=True)
    series = series / norm(series, axis=1).reshape(series.shape[0], 1)

    shard_size = 1000
    adj = []
    #adjnew = numpy.zeros((shard_size, series.shape[0]), dtype=series.dtype)
    for i, chunk in enumerate(grouper(series, shard_size)):
        if i % 10 == 0:
            LOGGER.info('nodes: %i', shard_size*(i))
        adjnew = numpy.dot(chunk, series.T)
        rows, columns = numpy.where(adjnew > min_th)
        values = [adjnew[i, j] for i, j in zip(rows, columns)]
        #adjnew[adjnew <= min_th] = 0
        #adj.append(scipy.sparse.dok_matrix(adjnew))
        adj.append(scipy.sparse.coo_matrix((values, (rows, columns)), shape=adjnew.shape))
    adj = scipy.sparse.vstack(adj)

    LOGGER.info('Filtering number_of_edges = %i', adj.nnz)
    adj = scipy.sparse.triu(adj)
    adj = scipy.sparse.lil_matrix(adj)
    #LOGGER.info('Saving')
    #scipy.io.mmwrite(listname, S)
    return node2voxel, adj

def reshape_data2mask(fimg, mask, mask_file):
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
    if os.path.isfile(func_file):
        LOGGER.info('one file image')
        img = nib.load(func_file)
    elif os.path.isdir(func_file):
        files = sorted(glob(os.path.join(func_file, '*.img')))
        LOGGER.info('multiple images. len(files)=%i', len(files))
        img = nib.funcs.concat_images(files)
    return img

def grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    I = itt.izip_longest(fillvalue=None, *args)
    I2 = list(itt.imap(lambda v: numpy.array([x for x in v if not x is None]), I))
    return I2
