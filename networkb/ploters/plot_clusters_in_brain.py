# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 16:22:45 2013

@author: andres
"""
import pylab
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
import nibabel as nib

def plot_clusters_in_brain(bn,cluster_list):
  node2voxel=bn.get_node2voxel()
  img=nib.load(bn.mask)
  D=img.get_data()
  sh=D.shape
  M=pylab.zeros((sh[0],sh[1],sh[2],3))
  M[D>0,:]=pylab.ones_like(M[D>0,:])
  """
  for i in range(sh[0]):
    for j in range(sh[1]):
      for k in range(sh[2]):
        if D[i,j,k]>0:
          M[i,j,k,:]=[1,1,1]
  """
  for node in cluster_list[0]:
    (i,j,k)=node2voxel[str(node)]
    M[i,j,k,:]=[1,0,0]
  plt.subplot(221)
  j0 = 50
  lj = plt.imshow(M[:,j0,:,:].swapaxes(0,1))
  plt.axis([0, sh[0], 0, sh[2]])
  plt.subplot(222)
  i0 = 50
  li = plt.imshow(M[i0,:,:,:].swapaxes(0,1))
  plt.axis([0, sh[1], 0, sh[2]])
  plt.subplot(223)
  k0 = 50
  lk = plt.imshow(M[:,:,k0,:].swapaxes(0,1))
  plt.axis([0, sh[0], 0, sh[1]])
  axcolor = 'lightgoldenrodyellow'
  axi = plt.axes([0.55, 0.3, 0.4, 0.03], axisbg=axcolor)
  axj = plt.axes([0.55, 0.2, 0.4, 0.03], axisbg=axcolor)
  axk = plt.axes([0.55, 0.1, 0.4, 0.03], axisbg=axcolor)
  si = Slider(axi, 'i', 0, sh[0], valinit=i0,valfmt='%i')
  sj = Slider(axj, 'j', 0, sh[1], valinit=j0,valfmt='%i')
  sk = Slider(axk, 'k', 0, sh[2], valinit=k0,valfmt='%i')
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
  radio = RadioButtons(rax, ('0', '1', '2'), active=0)
  def colorfunc(label):
    i0=int(si.val)
    j0=int(si.val)
    k0=int(si.val)
    for cluster in cluster_list:
      for node in cluster:
        (i,j,k)=node2voxel[str(node)]
        if D[i,j,k]>0:
          M[i,j,k,:]=[1,1,1]
    for node in cluster_list[int(label)]:
      (i,j,k)=node2voxel[str(node)]
      M[i,j,k,:]=[1,0,0]
    lj.set_data(M[:,j0,:,:].swapaxes(0,1))
    li.set_data(M[i0,:,:,:].swapaxes(0,1))
    lk.set_data(M[:,:,k0,:].swapaxes(0,1))
    plt.draw()
  radio.on_clicked(colorfunc)
  plt.show()
  return