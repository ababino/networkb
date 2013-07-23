# -*- coding: utf-8 -*-
"""
Created on Wed May 24 17:31:05 2013

@author: andres
"""
import numpy, networkb

def define_clusters(brain,jump_size, cluster_min):
  gc=numpy.array(brain.get_gc()[0])
  th=brain.get_th()
  cluster_dic=brain.get_cluster_dic()
  gc_diff=numpy.diff(gc)
  jump_thresholds=[]
  jump=min(gc_diff)

  while abs(jump)>=jump_size:
    min_ind=numpy.nanargmin(gc_diff)
    gc_diff[min_ind]=numpy.nan
    jump_thresholds.append(th[min_ind+1])
    jump=numpy.nanmin(gc_diff)
  """  
  clusters=[]
  if len(cluster_dic[str(max(jump_thresholds))][0])>=cluster_min:
    clusters.append(cluster_dic[str(max(jump_thresholds))][0])
  
  for jt in sorted(jump_thresholds):
    if len(cluster_dic[str(jt)])>1:
      if len(cluster_dic[str(jt)][1])>=cluster_min:
	clusters.append(cluster_dic[str(jt)][1])
  """
  clusters=[]
  if len(jump_thresholds)>0:
    for jt in sorted(jump_thresholds):
      c1=set(cluster_dic[str(jt)][0])
      th0=th[th.index(jt)-1]
      c2=set(cluster_dic[str(th0)][0])
      c3=list(c2.difference(c1))
      if len(c3)>cluster_min:
	clusters.append((jt,c3))
    last_jump=max(jump_thresholds)
    clusters.append((last_jump,cluster_dic[str(last_jump)][0]))
  return clusters
    


  

