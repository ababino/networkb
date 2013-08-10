# -*- coding: utf-8 -*-
"""
Created on Wed May 15 20:09:59 2013

@author: andres
"""
import networkb
import networkx
import itertools
import json

def weak_link_distribution(bn,N_clus=2):
  jumps=networkb.find_th_jumps(bn,N_clus)
  if len(jumps)==0:
    return []
  pc=max(jumps)
  G=bn.get_Graph(pc,correlation='positive')
  cluster_list=[x for x in networkx.connected_components(G) if x>0]
  if len(cluster_list)<2:
    return []
  thmin=(jumps[jumps.index(pc)-1]+pc)/2
  H=bn.get_Graph(thmin,th_up=pc,correlation='positive')
  if H.number_of_edges()<1:
    return []
  H=networkx.subgraph(H,itertools.chain.from_iterable(cluster_list))
  if H.number_of_edges()<1:
    return []
  d=[]  
  for e in H.edges_iter():
    d.append(bn.nodedistance(e))
  json.dump(d,open(bn.weak_link_distribution_file,'w'))
  return d
