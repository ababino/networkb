# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 16:22:45 2013

@author: andres
"""
import pylab
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons
import nibabel as nib
import os
import numpy

def plot_clusters_in_brain(bn,background,cluster_list,actives=[True,False,False]):
  if len(cluster_list)<3:
    print 'just '+str(len(cluster_list))+ ' cluster'
  colors=[[1,0,0],[0,1,0],[0,0,1]]  
  node2voxel=bn.get_node2voxel()
  img=nib.load(bn.mask)
  imgB = nib.load(background)
  D=img.get_data()
  sh=D.shape
  M=pylab.zeros((sh[0],sh[1],sh[2],3))
  B=M
  B0=imgB.get_data().astype(numpy.float32, copy=False) 
  #B0img.get_data().astype(numpy.float32, copy=False)
  B[:,:,:,0]=(B0-B0.min())/(B0.max()-B0.min())
  B[:,:,:,1]=B[:,:,:,0]
  B[:,:,:,2]=B[:,:,:,0]
  M=B  
  #M[D>0,:]=pylab.ones_like(M[D>0,:])
  for ii in range(len(actives)):
    if actives[ii]:
      if len(cluster_list)>ii:
        for node in cluster_list[ii]:
          (i,j,k)=node2voxel[str(node)]
          M[i,j,k,:]=[B[i,j,k,l]*colors[ii][l] for l in range(3)]
  plt.figure()
  plt.subplot(221)
  j0 = round(M.shape[1]*0.5)
  lj = plt.imshow(M[:,j0,:,:].swapaxes(0,1),interpolation='None')
  plt.axis([0, sh[0], 0, sh[2]])
  plt.subplot(222)
  i0 = round(M.shape[0]*0.5)
  li = plt.imshow(M[i0,:,:,:].swapaxes(0,1),interpolation='None')
  plt.axis([0, sh[1], 0, sh[2]])
  plt.subplot(223)
  k0 = round(M.shape[2]*0.5)
  lk = plt.imshow(M[:,:,k0,:].swapaxes(0,1),interpolation='None')
  plt.axis([0, sh[0], 0, sh[1]])
  axcolor = 'lightgoldenrodyellow'
  axi = plt.axes([0.55, 0.3, 0.4, 0.03], axisbg=axcolor)
  axj = plt.axes([0.55, 0.2, 0.4, 0.03], axisbg=axcolor)
  axk = plt.axes([0.55, 0.1, 0.4, 0.03], axisbg=axcolor)
  si = Slider(axi, 'i', 0, sh[0]-1, valinit=i0,valfmt='%i')
  sj = Slider(axj, 'j', 0, sh[1]-1, valinit=j0,valfmt='%i')
  sk = Slider(axk, 'k', 0, sh[2]-1, valinit=k0,valfmt='%i')
  def update(val):
    i = int(si.val)
    j = int(sj.val)
    k = int(sk.val)
    lj.set_data(M[:,j,:,:].swapaxes(0,1))
    li.set_data(M[i,:,:,:].swapaxes(0,1))
    lk.set_data(M[:,:,k,:].swapaxes(0,1))
    plt.draw()
  si.on_changed(update)
  sj.on_changed(update)
  sk.on_changed(update)
  rax = plt.axes([0.025, 0.5, 0.1, 0.15], axisbg=axcolor)
  check = CheckButtons(rax, ('0', '1', '2'), actives=actives)

  def showcluster(label):
    i0=int(si.val)
    j0=int(sj.val)
    k0=int(sk.val)
    M=B
    for node in cluster_list[int(label)]:
      (i,j,k)=node2voxel[str(node)]
      if all(M[i,j,k,:]==colors[int(label)]):
        M[i,j,k,:]=[1,1,1]
      else:
        M[i,j,k,:]=colors[int(label)]
    lj.set_data(M[:,j0,:,:].swapaxes(0,1))
    li.set_data(M[i0,:,:,:].swapaxes(0,1))
    lk.set_data(M[:,:,k0,:].swapaxes(0,1))
    plt.draw()
  check.on_clicked(showcluster)
  plt.show()
  return