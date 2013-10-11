# -*- coding: utf-8 -*-
"""
Created on Wed May 15 20:09:59 2013

@author: andres
"""

import networkb
import networkx
import itertools
import json

def weak_link_distribution(bn,N_clus=2,mcs=0,n_jumps=1):
  jumps=[j[0] for j in networkb.find_th_jumps(bn,N_clus)]
  jumps=sorted(list(set(jumps)),reverse=True)
  if len(jumps)==0:
    return []
  pcs=jumps[0:min(n_jumps,len(jumps)-1)]
  d=[]  
  for pc in pcs:
    G=bn.get_Graph(pc,correlation='positive')
    cluster_list=[x for x in networkx.connected_components(G) if x>mcs]
    if len(cluster_list)<2:
      continue
    thmin=(jumps[jumps.index(pc)+1]+pc)/2
    print thmin
    print pc
    H=bn.get_Graph(thmin,th_up=pc,correlation='positive')
    if H.number_of_edges()<1:
      continue
    H=networkx.subgraph(H,itertools.chain.from_iterable(cluster_list))
    if H.number_of_edges()<1:
      continue
    for e in H.edges_iter():
      d.append(bn.nodedistance(e))
  json.dump(d,open(bn.weak_link_distribution_file,'w'))
  return d
