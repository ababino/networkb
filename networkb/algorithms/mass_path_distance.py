# -*- coding: utf-8 -*-
"""
Created on Wed May 15 19:23:23 2013

@author: andres
"""
import networkx as nx
import networkb, numpy, json
from scipy import spatial
import logging
logger = logging.getLogger('networkb.algorithms.mass_path_distance')

def mass_path_distance(bn,N_clus,op='first',correlation='both'):
  """
  Returns mass, maximun radio, maximun path length and average path 
  length for all clusters at critical threshold. If op='first' the 
  critical threshold is the nearest to 1 and if op='biggest' the 
  threshold is the one with the biggest jump in the percolation curve.
  """
  N=[]
  lmean=[]
  lmax=[]
  rmax=[]
  voxels=bn.get_node2voxel()
  affine=bn.get_affine()

  jumps=networkb.find_th_jumps(bn,N_clus)
  if op=='first':  
    if len(jumps)==0:
      print 'no jumps'
      return {'N':[],'rmax':[],'lmax':[],'lmean':[]}
    pc=max(jumps)
  elif op=='biggest':
    gap=0
    pc=1
    gc=bn.get_gc()[0]
    th=bn.get_th()
    if len(jumps)==0:
      print 'no jumps'
      return {'N':[],'rmax':[],'lmax':[],'lmean':[]}
    for jump in jumps:
      j=th.index(jump)
      if gc[j-1]-gc[j]>gap:
        pc=jump
        gap=gc[j-1]-gc[j]    
  logger.info('pc=%f',pc)
  G=bn.get_Graph(pc,correlation=correlation)
  logger.info('number of nodes in network at pc: %i',G.number_of_nodes())
  cluster_list=nx.connected_components(G)
  edge_list=G.edges()
  logger.info('number of clusters: %i',len(cluster_list))
  N=[len(x) for x in cluster_list]
  for i,nodelist in enumerate(cluster_list):
    if i<10:
      logger.info('procesing cluster %i of mass %i',i,N[i])
    else:
      logger.debug('procesing cluster %i of mass %i',i,N[i])
    Gsub=bn.get_SubGraph(float(pc),nodelist,edge_list)
    A=numpy.zeros([N[i],3])
    pos=numpy.array([0,0,0,1])
    path_length=[]
    H = Gsub.copy()
    for j,node in enumerate(H):
      temp=[x for x in nx.single_source_shortest_path_length(Gsub, node).values() if x>0]
      path_length.extend(temp)
      pos[0:3]=voxels[str(node)]
      A[j,:]=numpy.dot(affine, numpy.transpose(pos))[0:3]
      Gsub.remove_node(node)
    D=spatial.distance.cdist(A,A)
    rmax.append(numpy.max(D))
    lmax.append(max(path_length))
    lmean.append(numpy.mean(path_length))
    del Gsub, D, A, pos
  json.dump({'N':N,'rmax':rmax,'lmax':lmax,'lmean':lmean},
            open(bn.cluster_properties_file,'w'))
  return {'N':N,'rmax':rmax,'lmax':lmax,'lmean':lmean}
