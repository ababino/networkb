# -*- coding: utf-8 -*-
"""
Created on Thu May 24 18:25:28 2013

@author: andres
"""
import pylab
def plot_percolation_curve(brain,i):
  pylab.figure(i)
  gc=brain.get_gc()
  th=brain.get_th()
  ax=pylab.plot(th,gc[0],'*-')
  pylab.grid(True)
  return ax
