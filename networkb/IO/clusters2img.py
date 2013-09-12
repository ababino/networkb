# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:43:03 2013

@author: andres
"""

import nibabel as nib
import pylab 

def clusters2img(bn,cluster_list,filename):
  """
  Save clusters in img or nii file. filename is the name of the output file 
  without the extension. 
  """
  cluster_list=sorted(cluster_list,key=len,reverse=True)
  node2voxel=bn.get_node2voxel()
  img=nib.load(bn.mask)
  D=img.get_data()
  sh=D.shape
  M=pylab.zeros((sh[0],sh[1],sh[2]))  
  for ii,cluster in enumerate(cluster_list):
    for node in cluster:
      (i,j,k)=node2voxel[str(node)]
      M[i,j,k]=ii+1
    img_out = nib.Nifti1Image(M, img.get_affine())
    nib.save(img_out,filename)
  