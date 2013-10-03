# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 12:56:02 2013

@author: andres
"""
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy
def plot_volume(V):
  sh=V.shape
  if len(sh)==3:
    V=V.astype(numpy.float32,copy=False)
    V=V/V.max()
    V2=numpy.zeros((sh[0],sh[1],sh[2],3))
    V2[:,:,:,0]=V
    V2[:,:,:,1]=V
    V2[:,:,:,2]=V
    V=V2
  elif len(sh)!=4:    
    print 'Must be 3D plus color dim'
    return
  plt.subplot(221)
  j0 = 50
  lj = plt.imshow(V[:,j0,:,:].swapaxes(0,1))
  plt.axis([0, sh[0], 0, sh[2]])
  plt.subplot(222)
  i0 = 50
  li = plt.imshow(V[i0,::-1,:,:].swapaxes(0,1))
  plt.axis([0, sh[1], 0, sh[2]])
  plt.subplot(223)
  k0 = 50
  lk = plt.imshow(V[:,:,k0,:].swapaxes(0,1))
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
      li.set_data(V[i,::-1,:,:].swapaxes(0,1))
      lj.set_data(V[:,j,:,:].swapaxes(0,1))
      lk.set_data(V[:,:,k,:].swapaxes(0,1))
      plt.draw()
  si.on_changed(update)
  sj.on_changed(update)
  sk.on_changed(update)
  return