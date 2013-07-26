# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 17:11:20 2013

@author: andres
"""

from pylab import get_cmap
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import nibabel as nib
import os.path

def plot_4d_brain(bn):
  img_dir=os.path.join(bn.func_dir,bn.name)
  img = nib.load(img_dir)
  D=img.get_data()
  cm=get_cmap("Greys_r")
  sh=D.shape
  t0=0
  plt.figure()
  plt.subplot(221)
  j0 = 50
  lj = plt.imshow(D[:,j0,:,t0].T,cmap=cm)
  plt.axis([0, sh[0], 0, sh[2]])
  plt.subplot(222)
  i0 = 50
  li = plt.imshow(D[i0,:,:,t0].T,cmap=cm)
  plt.axis([0, sh[1], 0, sh[2]])
  plt.subplot(223)
  k0 = 50
  lk = plt.imshow(D[:,:,k0,t0].T,cmap=cm)
  plt.axis([0, sh[0], 0, sh[1]])
  axcolor = 'lightgoldenrodyellow'
  axt = plt.axes([0.55, 0.4, 0.4, 0.03], axisbg=axcolor)
  axi = plt.axes([0.55, 0.3, 0.4, 0.03], axisbg=axcolor)
  axj = plt.axes([0.55, 0.2, 0.4, 0.03], axisbg=axcolor)
  axk = plt.axes([0.55, 0.1, 0.4, 0.03], axisbg=axcolor)
  si = Slider(axi, 'i', 0, sh[0]-1, valinit=i0,valfmt='%i')
  sj = Slider(axj, 'j', 0, sh[1]-1, valinit=j0,valfmt='%i')
  sk = Slider(axk, 'k', 0, sh[2]-1, valinit=k0,valfmt='%i')
  st = Slider(axt, 't', 0, sh[3]-1, valinit=t0,valfmt='%i')
  def update(val):
      i = int(si.val)
      j = int(sj.val)
      k = int(sk.val)
      t = int(st.val)
      li.set_data(D[i,:,:,t].T)
      lj.set_data(D[:,j,:,t].T)
      lk.set_data(D[:,:,k,t].T)
      plt.draw()
  si.on_changed(update)
  sj.on_changed(update)
  sk.on_changed(update)
  st.on_changed(update)
  plt.show()
  return