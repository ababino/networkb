# -*- coding: utf-8 -*-
"""
Created on Wed May 15 19:23:23 2013

@author: andres
"""
import networkx as nx
import networkb, numpy, json
from scipy import spatial
from networkb.algorithms.utils import find_peaks

def mass_path_distance(bn,N_clus,mcs,op='first'):#,jump_size
  N=[]
  lmean=[]
  lmax=[]
  rmax=[]
  voxels=bn.get_node2voxel()
  affine=bn.get_affine()

  #gc=bn.get_gc()
  #th=bn.get_th()
  #pc=str(max(find_peaks(th,gc),key=lambda x:x[1])[0])
  #pc=str(max(find_peaks(th,gc),key=lambda x:x[0])[0])
  #pc=th[th.index(max(networkb.find_th_jumps(bn,N_clus)))+1]
  jumps=networkb.find_th_jumps(bn,N_clus)
  if op=='first':  
    pc=max(jumps)
  elif op=='biggest':
    gap=0
    gc=bn.get_gc()[0]
    th=bn.get_th()
    for jump in jumps:
      j=th.index(jump)
      if gc[j-1]-gc[j]>gap:
        pc=jump
        gap=gc[j-1]-gc[j]    
  print pc
  G=bn.get_Graph(pc)
  print G.number_of_nodes()
  cluster_list=nx.connected_components(G)
  edge_list=G.edges()
  print len(cluster_list)
  #cluster_list=bn.get_cluster_dic()[pc]
  N=[len(x) for x in cluster_list]
  """
  cluster_list=networkb.define_clusters(bn,jump_size,mcs)
  N=[len(x[1]) for x in cluster_list]    
  """
  for i,nodelist in enumerate(cluster_list):
    print str(pc)+'  '+str(i)+', '+str(N[i])
    Gsub=bn.get_SubGraph(float(pc),nodelist,edge_list)
    A=numpy.zeros([N[i],3])
    pos=numpy.array([0,0,0,1])
    path_length=[]
    H = Gsub.copy()
    for j,node in enumerate(H):
      path_length.extend(nx.single_source_shortest_path_length(Gsub, node).values())
      pos[0:3]=voxels[str(node)]
      A[j,:]=numpy.dot(affine, numpy.transpose(pos))[0:3]
      Gsub.remove_node(node)
    D=spatial.distance.cdist(A,A)
    rmax.append(numpy.max(D))
    lmax.append(max(path_length))
    lmean.append(numpy.mean(path_length))
    del Gsub, D, A, pos
  json.dump({'N':N,'rmax':rmax,'lmax':lmax,'lmean':lmean},open(bn.cluster_properties_file,'w'))
  return {'N':N,'rmax':rmax,'lmax':lmax,'lmean':lmean} #(N,lmean,lmax,rmax)