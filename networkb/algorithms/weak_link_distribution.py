# -*- coding: utf-8 -*-
"""
Created on Wed May 15 20:09:59 2013

@author: andres
"""

import networkb
import networkx
import itertools
import json
import sys
sys.path.append('/home/andres/08-develop/nibabel')
sys.path.append('/home/andres/08-develop/snap-python/swig')
import snap

def weak_link_distribution(bn,N_clus=2,mcs=0,n_jumps=1):
    jumps = [j[0] for j in networkb.find_th_jumps(bn, N_clus)]
    jumps = sorted(list(set(jumps)), reverse=True)
    if len(jumps)==0:
        return []
    pcs = jumps[0:min(n_jumps, len(jumps)-1)]
    d = []
    for pc in pcs:
        G = bn.get_Graph(pc, correlation='positive')[0]
        CnComV  =  snap.TCnComV()
        MxWccGraph = snap.GetWccs(G, CnComV)
        cluster_list = [[n for n in clus] for clus in CnComV if clus.Len()>1]
        if len(cluster_list)<2:
            continue
        thmin = (jumps[jumps.index(pc)+1]+pc)/2
        print thmin
        print pc
        H = bn.get_Graph(thmin, th_up=pc, correlation='positive')[0]
        if H.GetEdges()<1:
            continue
        cluster_vector = snap.TIntV()
        for i in itertools.chain.from_iterable(cluster_list):
            cluster_vector.Add(i)
        H = snap.GetSubGraph(H, cluster_vector)
        if H.GetEdges()<1:
            continue
        for e in H.edges_iter():
            d.append(bn.nodedistance(e))
    json.dump(d, open(bn.weak_link_distribution_file, 'w'))
    return d
