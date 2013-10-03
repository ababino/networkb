# -*- coding: utf-8 -*-
"""
Created on Thu May 24 18:25:28 2013

@author: andres
"""

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.pyplot as plt
import networkb
import numpy

def plot_percolation_curve(brain):
  fig, ax = plt.subplots(subplot_kw ={'axisbg':'k'})#subplot_kw ={'axisbg':'k'}
  gc=brain.get_gc()[0]
  th=brain.get_th()
  ymax=max(gc)
  plt.plot(th,gc,'y.-')
  plt.grid(True,color='w')#,color='w'
  
  colors=[[1,0,0],[0,1,0],[0,0,1],[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]]
  jumps=networkb.find_th_jumps(brain)
  while len(jumps)>len(colors):
    colors.extend(colors)

  half=round(brain.volume_shape[2]/2)
  B0=brain.get_mask_data().astype(numpy.float32, copy=False)[:,:,half]
  sh=B0.shape
  B=numpy.zeros((sh[0],sh[1],3))
  B[:,:,0]=(B0-B0.min())/(B0.max()-B0.min())
  B[:,:,1]=B[:,:,0]
  B[:,:,2]=B[:,:,0]
  y0=0.05#0.1
  x0=0.9#jumps[0][0]
  for i,jump in enumerate(jumps):
    M=networkb.nodelist2volumen(brain,jump[1],1).max(2)
    N=numpy.zeros_like(B)
    N[M==0,:]=B[M==0,:]
    N[M==1,:]=B[M==1,:]*colors[i]
    xy=(jump[0],gc[th.index(jump[0])])
    xy2=(x0,y0)
    y0=y0+0.1
    if y0>0.9:
      y0=0.9
      x0=x0-0.1
    imagebox = OffsetImage(N.swapaxes(0,1), zoom=0.5)
    ab = AnnotationBbox(imagebox, xy,
                        xybox=xy2,
                        xycoords='data',
                        boxcoords="axes fraction",#'data',#"axes fraction"
                        pad=0.0,
                        arrowprops=dict(arrowstyle="->",color='w')
                        )
    ax.add_artist(ab)
  return ax
