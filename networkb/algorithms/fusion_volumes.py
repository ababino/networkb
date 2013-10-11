# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 16:56:53 2013

@author: andres
"""
import numpy

def fusion_volumes(Vol_list,weight=0.8):
  if weight>1 or weight<0:
    print 'weight must be a number in the (0,1) range'
    return
  sh=Vol_list[0].shape
  V1=Vol_list[0].astype(numpy.float32,copy=False)
  V1=V1/V1.max()
  V=numpy.zeros((sh[0],sh[1],sh[2],3))
  V[:,:,:,0]=V1
  V[:,:,:,1]=V1
  V[:,:,:,2]=V1
  for Vol in Vol_list[1:]:
    if not Vol.shape==sh:
      print 'Volume shapes not equal'
      return
    V2=Vol.astype(numpy.float32,copy=False)
    V2=V2/V2.max()
    V3=numpy.zeros((sh[0],sh[1],sh[2],3))
    V3[:,:,:,0]=V2
    V[V3.max(3)!=0,:]=(1-weight)*V[V3.max(3)!=0,:]+weight*V3[V3.max(3)!=0,:]
  return V