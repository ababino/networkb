# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:42:48 2013

@author: andres
"""
import pylab

def custom_plot(x,y,op='linlin'):
  pylab.figure()
  ax=pylab.subplot(111)
  pylab.plot(x,y,'*')
  if op=='loglog':
    ax.set_xscale("log")
    ax.set_yscale("log")
  elif op=='linlog':
    ax.set_yscale("log")
  ax.grid(True)
  return ax