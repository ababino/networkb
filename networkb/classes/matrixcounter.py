# -*- coding: utf-8 -*-
"""
Created on Wed May 15 18:52:15 2013

@author: andres
"""

class Matrix_Counter:
  """
  Iterator that returns a matrix without the first n columns in the n 
  itertation.
  """
  def __init__(self, M,ini,end):
    self.current = ini
    self.end = end
    self.high = M.shape[1]
    self.M=M
  def __iter__(self):
    return self
  def next(self):
    if self.current > self.high or self.current > self.end:
      raise StopIteration
    else:
      self.current += 1
      return self.M[:,self.current-1:]