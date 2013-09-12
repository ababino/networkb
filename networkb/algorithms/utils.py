# -*- coding: utf-8 -*-
"""
Created on Wed May 15 20:11:59 2013

@author: andres
"""
import numpy
import networkx
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

def find_th_jumps(bn,max_clus=2):
  """
  Returns the thresholds where a jump occurs. A jump is defined as the 
  join of the biggest cluster with, up to, the max_clus cluster.
  """
  NON=bn.get_non()
  node_list=[node for node,dat in NON.nodes(data=True) if dat['order']==0]
  subNON=networkx.Graph()
  for n1,n2 in NON.edges_iter(nbunch=node_list):
    subNON.add_edge(n1,n2)
  node_list=networkx.connected_components(subNON)[0]
  subNON=NON.subgraph(node_list)
  max_th=max([dat['th'] for n,dat in subNON.nodes(data=True)])
  N=bn.number_of_nodes()
  jumps=[]
  first_cluster=(0,[])
  for node,data in NON.nodes(data=True):
    if NON.degree(node)>=3 and NON.node[node]['order']==0:
      for node2 in NON.neighbors(node):
        if 0<NON.node[node2]['order']<=max_clus:
          if 20*len(NON.node[node2]['cc'])>len(NON.node[node]['cc']) or 200*len(NON.node[node2]['cc'])>N:
            if NON.node[node2]['th']<max_th:
              jumps.append((NON.node[node2]['th'],NON.node[node2]['cc']))
              if NON.node[node2]['th']>first_cluster[0]:
                for node3 in NON.neighbors(node):
                  if NON.node[node3]['order']==0 and NON.node[node3]['th']==NON.node[node2]['th']:
                    first_cluster=((NON.node[node3]['th'],NON.node[node3]['cc']))
  jumps.append(first_cluster)
  jumps=sorted(jumps,key=lambda x: x[0],reverse=True)
  return jumps
    
def nodelist2volumen(bn,nodelist,element):
  node2voxel=bn.node2voxel
  B=numpy.zeros(bn.volume_shape)   
  for node in nodelist:
    (i,j,k)=node2voxel[str(node)]
    B[i,j,k]=element
  return B


  
