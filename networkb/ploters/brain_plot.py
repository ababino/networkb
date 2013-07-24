# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 16:22:14 2013

@author: andres
"""
from pylab import get_cmap
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
def brain_plot(D):
  cm=get_cmap("Greys_r")
  sh=D.shape
  plt.subplot(221)
  j0 = 50
  lj = plt.imshow(D[:,j0,:].T,cmap=cm)
  plt.axis([0, sh[0], 0, sh[2]])
  plt.subplot(222)
  i0 = 50
  li = plt.imshow(D[i0,:,:].T,cmap=cm)
  plt.axis([0, sh[1], 0, sh[2]])
  plt.subplot(223)
  k0 = 50
  lk = plt.imshow(D[:,:,k0].T,cmap=cm)
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
      li.set_data(D[i,:,:].T)
      lj.set_data(D[:,j,:].T)
      lk.set_data(D[:,:,k].T)
      plt.draw()
  si.on_changed(update)
  sj.on_changed(update)
  sk.on_changed(update)
  return