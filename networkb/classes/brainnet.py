"""
Base class BrainNet
"""
from networkb.algorithms.mass_path_distance import mass_path_distance
from networkb.algorithms.weak_link_distribution import weak_link_distribution
from glob import glob
import json
import gc as garbage
import itertools as itt
import os.path
import os
import networkx as nx
import sys
sys.path.append('/home/andres/08-develop/nibabel')
sys.path.append('/home/andres/08-develop/snap-python/swig')
import snap
import nibabel as nib
import numpy
import scipy
from scipy import spatial
from networkx.readwrite import json_graph
import logging
LOGGER = logging.getLogger('networkb.classes.brainnet')
# pylint: disable-msg=W0311
# pylint: disable-msg=C0103
# pylint: disable-msg=E1101

class BrainNet(object):
  """
  BrainNet is de class representing a brain scan and its network atributtes.
  """
  def __init__(self, directory, name, mask=None, min_th=0.6,
               force_edgelist=False, force_percolation=False,
               correlation='both', func_dir='func'):
    self.dir = directory
    self.subject = os.path.basename(directory)
    self.func_dir = os.path.join(directory, func_dir)
    self.name = name
    if mask != None:
      if os.path.isabs(mask):
        self.mask = mask
      else:
        self.mask = os.path.join(self.func_dir, mask)
    else:
      self.mask = mask
    self.network_dir = os.path.join(self.func_dir, name.split('.')[0] + '_network')

    fh = logging.FileHandler(self.network_dir + '/info.log')
    fh.setLevel(logging.DEBUG)
    LOGGER.addHandler(fh)
    formatter  =  logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    LOGGER.info(self.subject)

    self.edgelist_file = os.path.join(self.network_dir,'edgelist.mtx')
    self.node2voxel_file = os.path.join(self.network_dir,'node2voxel.json')
    self.affine_file = os.path.join(self.network_dir,'affine.npy')
    if not os.path.exists(self.affine_file):
      img = self.get_img()
      self.affine = img.get_affine()
      numpy.save(self.affine_file, self.affine)
    else:
      self.affine = numpy.load(self.affine_file)
    self.min_th = min_th
    if not os.path.exists(self.network_dir):
      os.makedirs(self.network_dir)
    if not os.path.exists(self.edgelist_file) or force_edgelist:
      LOGGER.info('Building edgelist')
      self.__gen_edgelist__()
    self.node2voxel = self.get_node2voxel()
    if (
    not os.path.exists(os.path.join(self.network_dir, 'cluster_dic.json')) or
    force_percolation):
      LOGGER.info('Percolating')
      self.percolation_network(correlation = correlation)

    self.non_file = os.path.join(self.network_dir,'non.json')
    self.cluster_dic_file = os.path.join(self.network_dir,'cluster_dic.json')
    self.gc_file = os.path.join(self.network_dir,'gc.json')
    self.th_file = os.path.join(self.network_dir,'th.json')
    self.cluster_properties_file = os.path.join(self.network_dir,
                                              'cluster_properties.json')
    self.weak_link_distribution_file = os.path.join(self.network_dir,
                                                  'weak_link_distribution.json')
    self.volume_shape_file = os.path.join(self.network_dir,'shape.npy')
    if not os.path.exists(self.volume_shape_file):
      self.volume_shape = self.reshape_data2mask().shape[0:3]
      numpy.save(self.volume_shape_file,self.volume_shape)
    else:
      self.volume_shape = numpy.load(self.volume_shape_file)
    return

  def __check_dir_structure__(self):
    nii_files = glob(self.dir+'*.nii.gz')
    nii_mask = ''
    nii_data = ''
    for nii_file in nii_files:
      if nii_file.find('mask') == -1:
        nii_mask = nii_file
      else:
        nii_data = nii_file
    if len(nii_mask)=='':
      print 'No mask nii.gz file'

    if len(nii_data)=='':
      print 'No mask nii.gz file'

    dat_files = glob(self.dir+'*.dat')
    if len(dat_files)!=1:
      print 'wrong number of .dat files: ' + str(len(dat_files))
    return

#===========================================================================
# Generate edge list
#===========================================================================
  def __gen_edgelist__(self):
    """
    Gennerates a list of correlations among voxel's time series of file
    "filename". "mask" is the mask file for "filename. "th" is the minimun
    threshold that will be save to the file.
    """
    LOGGER.info('loading image')
    M = self.get_mask_data()
    D = self.reshape_data2mask()
    sh = D.shape

    listname = os.path.join(self.network_dir,'edgelist.mtx')

    Dr = numpy.zeros((sh[0]*sh[1]*sh[2],sh[3]))
    node2voxel = {}
    count = 0

    for (i, j, k) in itt.product(range(sh[0]), range(sh[1]), range(sh[2])):
      if M[i, j, k] != 0 and any(D[i, j, k, :] != D[i, j, k, 0]):
        v = D[i, j, k, :] - numpy.mean(D[i, j, k, :])
        node2voxel[count] = (i, j, k)
        Dr[count, :] = v/numpy.linalg.norm(v)
        count = count + 1

    self.count = count
    LOGGER.info('number of nodes in scan: %i', count)
    Dr = Dr[:count, :]
    Drt = Dr.copy()
    Drt = Drt.transpose()

    f = open(os.path.join(self.network_dir, 'node2voxel.json'), 'w')
    json.dump(node2voxel, f)
    f.close()

    n = 1000
    S = self.correlate(Drt, Dr[0:n, :])
    for i, v in enumerate(grouper(Dr[n:, :], n)):
      if i % 10 == 0:
        LOGGER.info('nodes: %i', S.shape[0])
      Snew = self.correlate(Drt, v)
      S = scipy.sparse.vstack([S, Snew])

    LOGGER.info('Filtering number_of_edges = %i', S.nnz)
    S = scipy.sparse.triu(S)
    S = scipy.sparse.lil_matrix(S)
    LOGGER.info('Saving')
    scipy.io.mmwrite(listname, S)

    return

  def __correlate__(self,Drt,v):
    try:
      A = numpy.zeros((v.shape[0],Drt.shape[1]))
      C = numpy.dot(v,Drt,out = A)
    except:
      LOGGER.info('dot')
      C = numpy.dot(v,Drt)      
    C[numpy.absolute(C)<=self.min_th] = 0
    SC = scipy.sparse.coo_matrix(C)
    
    return SC  

  
  def percolation_network(self,correlation = 'both'):
    if correlation not in ['negative','positive','both']:
      raise Exception(
      'correlation must be one of: negative, positive or both')
    G, wdic = self.get_graph(self.min_th, correlation=correlation)

    LOGGER.info('number of nodes in network: %i', G.GetNodes())
    mcs = max(int(0.001*G.GetNodes()),2)
    th = list(numpy.arange(self.min_th, 1,0.001))
    (gc, NON) = self.__percolation__(G, th, mcs, wdic)
    self.save_data(NON,gc,th)  
    return

  def save_data(self, NON, gc, th):
    nondata = json_graph.node_link_data(NON)

    nonfile = open(os.path.join(self.network_dir,'non.json'),'w')
    nonfile.write(json.dumps(nondata))
    nonfile.close()

    gcfile = open(os.path.join(self.network_dir,'gc.json'),'w')
    gcfile.write(json.dumps(gc))
    gcfile.close()

    thfile = open(os.path.join(self.network_dir,'th.json'),'w')
    thfile.write(json.dumps(th))
    thfile.close()
    return

#===========================================================================
# Geters
#===========================================================================
  def get_node2voxel(self):
    node2voxel = json.load(open(self.node2voxel_file))
    return node2voxel

  def get_affine(self):
    img = self.get_img()
    affine = img.get_affine()
    return affine    

  def reshape_data2mask(self):
    M = self.get_mask_data()
    D = self.get_img_data()
    sh = D.shape
    if any([M.shape[i]!=sh[i] for i in [0,1,2]]):
      LOGGER.info('Fill out of bounding box')
      img2 = self.get_mask()      
      A_mask = img2.get_affine()
      A = self.get_affine()
      temp = numpy.linalg.inv(A_mask)-numpy.linalg.inv(A)
      (i0,j0,k0) = temp.round().astype('int32')[0:3,3]
      LOGGER.info('Shift (%i,%i,%i)',i0,j0,k0)
      temp = numpy.zeros((M.shape[0],M.shape[1],M.shape[2],sh[3]))
      temp[i0:i0+sh[0],j0:j0+sh[1],k0:k0+sh[2]] = D
      D = temp
    return D

  def get_mask_data(self):
    if self.mask != None:
      if os.path.isabs(self.mask):
        img2 = nib.load(self.mask)
      else:
        img2 = nib.load(os.path.join(self.func_dir,self.mask))
      M = img2.get_data()
    else:
      D = self.get_img_data()
      sh = D.shape
      M = numpy.ones((sh[0],sh[1],sh[2]))
    return M

  def get_img(self):
    """
    Get 4d functional image.
    """
    img_dir = os.path.join(self.func_dir, self.name)
    if os.path.isfile(img_dir):
      LOGGER.info('one file image')
      img  =  nib.load(img_dir)
    elif os.path.isdir(img_dir):
      files = sorted(glob(os.path.join(img_dir,'*.img')))
      LOGGER.info('multiple images. len(files) = %i', len(files))
      img = nib.funcs.concat_images(files)
    return img

  def get_anatomical(self):
    """
    Get 3d anatomical image.
    """
    ana_dir = os.path.join(self.dir, 'anat/mprage_anonymized_brain.nii.gz')
    if os.path.isfile(ana_dir):
      LOGGER.info('one file image')
      ana =  nib.load(ana_dir)
    elif os.path.isdir(ana_dir):
      files = sorted(glob(os.path.join(ana_dir,'*.img')))
      LOGGER.info('multiple images. len(files) = %i', len(files))
      ana = nib.funcs.concat_images(files)
    return ana


  def get_mask(self):
    if os.path.isabs(self.mask):
      img = nib.load(self.mask)
    else:
      img = nib.load(os.path.join(self.func_dir, self.mask))
    return img

  def get_img_data(self):
    """
    Get 4D data from image.
    """
    img = self.get_img()
    D = img.get_data()
    return D

  def get_gc(self):
    gc = json.load(open(self.gc_file))
    return gc

  def get_th(self):
    th = json.load(open(self.th_file))
    return th

  def get_cluster_dic(self):
    cluster_dic = json.load(open(self.cluster_dic_file))
    return cluster_dic

  def get_non(self):
    NON = json_graph.node_link_graph(json.load(open(self.non_file)))
    return NON

  def get_cluster_properties(self,force = False,N_clus = 3,op = 'first',
                             correlation = 'both'):
    if os.path.exists(self.cluster_properties_file) and not force:
      clusters_properties = json.load(open(self.cluster_properties_file))
    else:
      clusters_properties = mass_path_distance(self, N_clus, op, correlation)
    return clusters_properties

  def get_weak_link_distribution(self, force=False, N_clus=2, mcs=0, n_jumps=1):
    if os.path.exists(self.weak_link_distribution_file) and not force:
      d = json.load(open(self.weak_link_distribution_file))
    else:
      d = weak_link_distribution(self, N_clus=N_clus, mcs=mcs)
    return d

  def get_SubGraph(self, th, nodelist, edge_list=None, correlation='both'):
    if correlation not in ['negative', 'positive', 'both']:
      raise Exception(
      'correlation must be one of: negative, positive or both')
    G = nx.Graph()
    if edge_list is None:
      M = scipy.io.mmread(self.edgelist_file)
      M = scipy.sparse.lil_matrix(M)
      for i, row in enumerate(M.rows):
        for j in row:
          w = M[i, j]
          if correlation == 'both':
            if th <= abs(w) and i in nodelist:
              G.add_edge(i, j)
          elif correlation == 'positive':
            if th <= w and i in nodelist:
              G.add_edge(i, j)
          elif correlation == 'negative':
            if w <= -th and i in nodelist:
              G.add_edge(i, j)
    else:
      for edge in edge_list:
        if edge[0] in nodelist and edge[1] in nodelist:
          G.add_edge(edge[0], edge[1])
    return G

  def get_graph(self, th_low, th_up=1.1, correlation='both'):
    """
    Returns graph with edges in (th_low, th_up) interval.
    """
    if correlation not in ['negative', 'positive', 'both']:
      raise Exception(
      'correlation must be one of: negative, positive or both')
    if correlation == 'both':
      test = lambda th_low, weight, th_up : th_low < abs(weight) <= th_up
    if correlation == 'positive':
      test = lambda th_low, weight, th_up : th_low < weight <= th_up
    if correlation == 'negative':
      test = lambda th_low, weight, th_up : -th_up <= weight < th_low
    graph = snap.TNEANet.New()
    wdic = {}
    """
    if os.path.isfile(self.edgelist_file[:-3]+"graph"):
      FIn  =  snap.TFIn(self.edgelist_file[:-3]+"graph")
      graph = snap.TNEANet.Load(FIn)
    else:
    """
    adj_matrix = scipy.io.mmread(self.edgelist_file)
    adj_matrix = scipy.sparse.lil_matrix(adj_matrix)
    rows, cols = adj_matrix.nonzero()
    for row, col in zip(rows, cols):
        row = int(row)
        col = int(col)
        weight = adj_matrix[row, col]
        if test(th_low, weight, th_up):
            if not graph.IsNode(row):
                graph.AddNode(row)
            if not graph.IsNode(col):
                graph.AddNode(col)
            edge = graph.AddEdge(row, col)
            graph.AddFltAttrDatE(edge, weight, 'float')
    fout = snap.TFOut(self.edgelist_file[:-3]+"graph")
    graph.Save(fout)
    fout.Flush()
    return graph, wdic

  def number_of_nodes(self):
    """
    Number of vocxel with data.
    """
    D = self.get_img_data()
    M = self.get_mask_data()
    sh = D.shape
    count = 0
    for i in range(sh[0]):
      for j in range(sh[1]):
        for k in range(sh[2]):
          if M[i, j, k] != 0 and any(D[i, j, k, :] != 0):
            count = count+1
    return count

  def __prune__(self, G, th, mcs):
    """
    Prunes G and returns a list of clusters biggers than mcs
    wdic = {}
    for edge in G.Edges():
      c = G.GetFltAttrDatE(edge.GetId(),'float')
      if c<th:
        wdic[(edge.GetSrcNId(),edge.GetDstNId())] = c
    """
    r_edges = [(NI.GetId(), Id) for NI in G.Nodes() for
               Id in NI.GetOutEdges() if
               G.GetFltAttrDatE(G.GetEId(NI.GetId(), Id), 'float')<th]
    """
    r_edges = []
    for NI in G.Nodes():
      for Id in NI.GetOutEdges():
        if G.GetFltAttrDatE(G.GetEId(NI.GetId(), Id),'float')<th:
          r_edges.append((NI.GetId(), Id))
    """
    LOGGER.debug('deleting %i edges', len(r_edges))
    for edge in r_edges:#wdic:
        G.DelEdge(edge[0], edge[1])
    LOGGER.debug('deleted edges: %i', len(r_edges))
    CnComV = snap.TCnComV()
    MxWccGraph = snap.GetWccs(G, CnComV)

    cc = [[n for n in clus] for clus in CnComV if clus.Len()>1]

    snap.DelZeroDegNodes_PNEANet(G)
    LOGGER.debug('connected components')
    del r_edges, CnComV, MxWccGraph
    garbage.collect()
    return [G, cc]

  def __percolation_data__(self, gc, cc, nn):
    cc_sizes = [float(len(clus))/nn for clus in cc]
    for j in range(max(len(cc_sizes), len(gc))):
      if j >= len(gc):
        lengs = [len(x) for x in gc]
        gc.append([0 for k in range(max(lengs)-1)])
      if j >= len(cc_sizes):
        gc[j].append(0)
      else:
        gc[j].append(cc_sizes[j])
    return gc

  def __update_non__(self, NON, cc, th, th_old):
    """
    Update netowrk of networks directed graph.
    """
    check = True
    nodes = NON.nodes(data=True)
    if len(NON)==0:
      n = 0
      check = False
    else:
      n = max(NON.nodes())+1
    cc = sorted(cc, key=len, reverse=True)
    for j in range(len(cc)):
      # Ignore small clusters for efficiency
      if len(cc[j])>10:
        NON.add_node(n, th=th, cc=cc[j], order=j)
        for (node, dat) in nodes:
          if dat['th'] == th_old:
            if cc[j][0] in dat['cc']:
                NON.add_edge(node, n)
        if nx.degree(NON, n)==0 and check:
          LOGGER.warn('node in NON without parent. (%i,%i,th = %f,th_old = %f)',
                      n,len(cc[j]),th,th_old)
        n = n + 1
    return NON

  def __percolation__(self, G, th, mcs, wdic):
    """
    Percolation process of G throu the thresholds th. Returns gc: a list of
    sizes of clusters for a given threshold; NON: Network of networks;
    cluster_dic: a dictionary of clusters bigger than mcs for each
    threshold.
    """
    number_of_nodes = self.number_of_nodes()
    #powr of components list
    power = [[]]
    #network of networks
    non = nx.DiGraph()
    for i in range(len(th)):
      LOGGER.debug('threshold: %f', th[i])
      [G, clusters] = self.__prune__(G, th[i], mcs)
      garbage.collect()
      power = self.__percolation_data__(power, clusters, number_of_nodes)
      non = self.__update_non__(non, clusters, th[i], th[i-1])
    return (power, non)

  def nodedistance(self, (n1, n2)):
    """
    node distance in cm. (en general)
    """
    ind1 = self.node2voxel[str(n1)]
    ind2 = self.node2voxel[str(n2)]
    if len(ind1)==3:
      ind1.append(1)
    if len(ind2)==3:
      ind2.append(1)
    v1 = numpy.dot(self.affine, numpy.transpose(ind1))[0:3]
    v2 = numpy.dot(self.affine, numpy.transpose(ind2))[0:3]
    d = spatial.distance.euclidean(v1, v2)
    return d

  def __unicode__(self):
    return unicode(self.name+', '+ self.subject).encode('utf-8')
  def __repr__(self):
    return self.name + ', ' + self.subject

def grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    I = itt.izip_longest(fillvalue=None, *args)
    I2 = list(itt.imap(lambda v: numpy.array([x for x in v if x != None]), I))
    return I2


