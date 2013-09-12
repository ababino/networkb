# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 19:40:49 2013

@author: andres
"""

import nibabel as nib
import pylab 

def clusters_probability2img(bns,cluster_lists,filename):
  """
  Save clusters in img or nii file. filename is the name of the output file 
  without the extension. 
  """
  img=nib.load(bns[0].mask)
  D=img.get_data()
  sh=D.shape
  M=pylab.zeros((sh[0],sh[1],sh[2]))  
  for l,bn in enumerate(bns):
    cluster_list=cluster_lists[l]
    cluster_list=sorted(cluster_list,key=len,reverse=True)
    node2voxel=bn.get_node2voxel()
    for cluster in cluster_list:
      for node in cluster:
        (i,j,k)=node2voxel[str(node)]
        M[i,j,k]=M[i,j,k]+1
  img_out = nib.Nifti1Image(M, img.get_affine())
  nib.save(img_out,filename)
  