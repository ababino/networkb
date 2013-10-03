# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:42:48 2013

@author: andres
"""
import pylab as pl

def custom_plot(x,y,op='linlin'):
  pl.figure()  
  ax=pl.subplot(111)
  pl.plot(x,y,'*')
  if op=='loglog':
    ax.set_xscale("log")
    ax.set_yscale("log")
  elif op=='linlog':
    ax.set_yscale("log")
  ax.grid(True)
  return ax