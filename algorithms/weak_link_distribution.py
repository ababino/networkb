# -*- coding: utf-8 -*-
"""
Created on Wed May 15 20:09:59 2013

@author: andres
"""
from networkb.algorithms.utils import nodedistance, find_peaks, window_correlation
import json, scipy.stats, numpy, pprocess

def parallel_bimodal(X):
  N=X[1][3]
  temp=window_correlation(X[1][0],X[1][1],X[1][2]) 
  (count,bins)=numpy.histogram(temp,bins=20,normed=True)
  ku=scipy.stats.kurtosis(temp, axis=0, fisher=True, bias=True)
  sk=scipy.stats.skew(temp, axis=0, bias=True)
  bimodal_coeff=(sk**2+1) / ( ku + (3*(N-1)**2) / ((N-2)*(N-3)) )
  return (X[0],bimodal_coeff)

def weak_link_distribution(bn,ncores=4):
  voxels=bn.get_node2voxel()
  for key in voxels:
    voxels[key].append(1)
  
  D=bn.get_img_data()
  affine=bn.get_affine()
  gc=bn.get_gc()
  th=bn.get_th()

  peaks=find_peaks(th,gc)
  (th1,temp)=max(peaks,key=lambda x:x[1])
  i=peaks.index((th1,temp))
  if i>0:
    th2=(peaks[i-1][0]+th1)/2
  elif len(peaks)>i+1:
    th1=peaks[i+1][0]
    th2=(peaks[i][0]+th1)/2
  else:
    return (False,None,None,None)
  th2=0.6
  print '(th2,th1)=('+str(th2)+', '+str(th1)+')'
  G=bn.get_Graph(th2,th_up=th1)
  wl=[]
  wl=[(n1,n2) for (n1,n2) in G.edges_iter()]       
  pos=[nodedistance(affine,voxels,str(n1),str(n2)) for (n1,n2) in wl]
  n=20
  N=D.shape[3]-n
  bm=[0  for x in range(len(wl))]
  par_in=[]
  for (n1,n2) in wl:
    [i,j,k]=voxels[str(n1)][0:3]
    s1=D[i,j,k,:]
    [i,j,k]=voxels[str(n2)][0:3]
    s2=D[i,j,k,:]
    par_in.append((s1,s2,n,N))

  queue = pprocess.Queue(limit=ncores)
  calc = queue.manage(pprocess.MakeParallel(parallel_bimodal))
  for (i,x) in enumerate(par_in):
    calc((i,x))
  for (i,bmc) in queue:
    bm[i]=bmc
  """
  for (n1,n2) in wl:
    [i,j,k]=voxels[str(n1)][0:3]
    s1=D[i,j,k,:]
    [i,j,k]=voxels[str(n2)][0:3]
    s2=D[i,j,k,:]
    temp=window_correlation(s1,s2,n) 
    (count,bins)=numpy.histogram(temp,bins=20,normed=True)
    ku=scipy.stats.kurtosis(temp, axis=0, fisher=True, bias=True)
    sk=scipy.stats.skew(temp, axis=0, bias=True)
    bimodal_coeff=(sk**2+1) / ( ku + (3*(N-1)**2) / ((N-2)*(N-3)) )
    bm.append(bimodal_coeff)    
  del D
  """
  json.dump((True,wl,bm,pos),open(bn.weak_link_distribution_file,'w'))
  return (True,wl,bm,pos)
