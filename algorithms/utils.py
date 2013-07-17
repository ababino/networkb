# -*- coding: utf-8 -*-
"""
Created on Wed May 15 20:11:59 2013

@author: andres
"""
import numpy, networkx
from scipy import spatial
from scipy import stats

def find_peaks(th,gc):
  peaks=[]
  for i in range(1,len(th)-1):
    if gc[1][i-1]<gc[1][i] and gc[1][i]>gc[1][i+1]:
      peaks.append((th[i],gc[1][i]))
  return peaks

def nodedistance(affine,vP,n1,n2):
  """
  node distance in cm. (en general)
  """
  ind1=vP[n1]
  ind2=vP[n2]
  if len(ind1)==3:
    ind1.append(1)
  if len(ind2)==3:
    ind2.append(1)
  v1=numpy.dot(affine, numpy.transpose(ind1))[0:3]
  v2=numpy.dot(affine, numpy.transpose(ind2))[0:3]
  d=spatial.distance.euclidean(v1,v2)
  return d

def power_law_fit(x,y):
  pl = lambda A, d, x: A*x**d 
  a, b, r_value, p_value, std_err = stats.linregress(numpy.log(x),numpy.log(y))
  y_fit=pl(numpy.exp(b),a,x)
  return (a,y_fit)

def exp_fit(x,y):
  exp_fun = lambda A, x0, x: A*numpy.exp(x/x0) 
  a, b, r_value, p_value, std_err = stats.linregress(x,numpy.log(y))
  A=numpy.exp(b)  
  x0=1.0/a  
  y_fit=exp_fun(A,x0,x)
  return (A,x0,y_fit)

def gaussian_fit(x,y):
  pl = lambda A, x0, s, x: A*numpy.exp(((x-x0)**2)/s) 
  p = numpy.polyfit(x,numpy.log(y),2)
  s=1./p[0]
  x0=-p[1]/(2*p[0])
  A=numpy.exp(p[2]+(p[1]**2)/(4*p[0]))
  y_fit=pl(A,x0,s,x)
  return ((A,x0,s),y_fit)

def window_correlation(x,y,w):
  if len(x)!=len(y):
    print 'vector x and y must be of the same size'
    print 'len(x)='+str(len(x))
    print 'len(y)='+str(len(y))
    return
  if len(x)<w:
    print 'window mus be smaller than len(x)'
    print 'len(x)='+str(len(x))+' w='+str(w)
  N=len(x)-w
  return [stats.pearsonr(x[i:i+w],y[i:i+w])[0] for i in range(N)]

def find_th_jumps(bn,max_clus=3):
  NON=bn.get_non()
  nodes_dic={}
  ths=[]
  for (node,dat) in NON.nodes(data=True):
    ths.append(dat['th'])
    if dat['th'] in nodes_dic:
      nodes_dic[dat['th']].append((node,dat['cc']))
    else:
      nodes_dic[dat['th']]=[(node,dat['cc'])]
  
  ths=list(set(ths))
  subnodes=[]
  for th in ths:
    nodes_dic[th]=sorted(nodes_dic[th],key=lambda x: len(x[1]),reverse=True)[:min(len(nodes_dic[th]),max_clus)]
    subnodes.extend([x[0] for x in nodes_dic[th]])
  
  subNON=networkx.subgraph(NON,subnodes)
  cc=networkx.connected_components(subNON.to_undirected())
  subNON=networkx.subgraph(NON,cc[0])  
  jumps=[]
  for th in ths:
    node=nodes_dic[th][0][0]
    if node in subNON.nodes():
      if networkx.degree(subNON,node)==3:
        jumps.append(th)
      
  thresholds=bn.get_th()
  jumps=sorted([thresholds[thresholds.index(x)+1] for x in jumps])
  return jumps
    
    
    



  
