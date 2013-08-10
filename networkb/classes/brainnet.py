from networkb.classes.matrixcounter import Matrix_Counter
from networkb.algorithms.mass_path_distance import mass_path_distance
from networkb.algorithms.weak_link_distribution import weak_link_distribution
from glob import glob
import json 
import pprocess 
import itertools as itt
import os.path
import os
import networkx as nx
import nibabel as nib
import pylab as pl
import numpy
from scipy import spatial
from networkx.readwrite import json_graph
import logging
logger = logging.getLogger('networkb.classes.brainnet')

class BrainNet():
  def __init__(self,directory,name,mask=None,min_th=0.6,ncores=4,
               force_edgelist=False,force_percolation=False,
               correlation='both'):
    self.dir=directory
    self.subject=os.path.basename(directory)
    self.name=name
    self.mask=mask
    self.func_dir=os.path.join(directory,'func')
    self.network_dir=os.path.join(directory,'network')
    self.edgelist_file=os.path.join(self.network_dir,'edgelist.dat')
    self.node2voxel_file=os.path.join(self.network_dir,'node2voxel.json')
    self.affine=self.get_affine()
    
    if not os.path.exists(self.network_dir):
      os.makedirs(self.network_dir)
    fh=logging.FileHandler(self.network_dir+'/info.log')
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.info(self.subject)
    if not os.path.exists(self.edgelist_file) or force_edgelist:
      logger.info('Building edgelist')
      self.gen_edgelist(min_th,ncores)
    self.min_th=min_th
    self.node2voxel=self.get_node2voxel()
    if (
    not os.path.exists(os.path.join(self.network_dir,'cluster_dic.json')) or 
    force_percolation):
      logger.info('Percolating')   
      self.percolation_network(ncores,correlation=correlation)

    self.non_file=os.path.join(self.network_dir,'non.json')
    self.cluster_dic_file=os.path.join(self.network_dir,'cluster_dic.json')
    self.gc_file=os.path.join(self.network_dir,'gc.json')
    self.th_file=os.path.join(self.network_dir,'th.json')
    self.cluster_properties_file=os.path.join(self.network_dir,
                                              'cluster_properties.json')
    self.weak_link_distribution_file=os.path.join(self.network_dir,
                                                  'weak_link_distribution.json')
    return

  def check_dir_structure(self):
    nii_files=glob(self.dir+'*.nii.gz')
    nii_mask=''
    nii_data=''
    for nii_file in nii_files:
      if nii_file.find('mask')==-1:
        nii_mask=nii_file
      else:
        nii_data=nii_file
    if len(nii_mask)=='':
      print 'No mask nii.gz file'

    if len(nii_data)=='':
      print 'No mask nii.gz file'

    dat_files=glob(self.dir+'*.dat')
    if len(dat_files)!=1:
      print 'wrong number of .dat files: ' + str(len(dat_files))
    return

#===========================================================================
# Generate edge list
#===========================================================================
  def gen_edgelist(self,th,ncores):
    """
    Gennerates a list of correlations among voxel's time series of file 
    "filename". "mask" is the mask file for "filename". "ncores" is the 
    number of cores to use in the calculation. "th" is the minimun threshold 
    that will be save to the file.
    """
    img_dir=os.path.join(self.func_dir,self.name)
    img = nib.load(img_dir)
    D=img.get_data()
    network_dir=os.path.join(self.dir,'network')
    listname=os.path.join(network_dir,'edgelist.dat')
  
    sh=D.shape  
    if self.mask!=None:
      if os.path.isabs(self.mask):
        img2 = nib.load(self.mask)
      else:
        img2 = nib.load(os.path.join(self.func_dir,self.mask))
      M=img2.get_data()
    else:
      M=pl.ones((sh[0],sh[1],sh[2]))      
    
    Dr=pl.zeros((sh[0]*sh[1]*sh[2],sh[3]))  
    node2voxel={}
    count=0
    
    for (i,j,k) in itt.product(range(sh[0]),range(sh[1]),range(sh[2])):
      if M[i,j,k]!=0 and any(D[i,j,k,:]!=0):
        v=D[i,j,k,:]-pl.mean(D[i,j,k,:])
        node2voxel[count]=(i,j,k)
        Dr[count,:]=v/pl.linalg.norm(v)
        count=count+1
    
    logger.info('number of nodes in scan: %i', count)
    Dr=Dr[:count,:]
    Drt=Dr.copy()
    Drt=Drt.transpose()
    
    f=open(os.path.join(network_dir,'node2voxel.json'),'w')
    json.dump(node2voxel,f)
    f.close()
    
    T=itt.repeat(th)
    I=Matrix_Counter(Drt,1,count)
  
    queue = pprocess.Queue(limit=ncores)
    calc = queue.manage(pprocess.MakeParallel(self.correlate))
    for inp in itt.izip(Dr,I,xrange(count),T):
      calc(inp)
  
    f=open(listname,'w')
    for d in queue:
      for s in d:
        f.write(s)
    f.close()
    return

  def correlate(self,(v,M,i,th)):
    out=[]
    if not(M is None) and not(i is None):  
      C=pl.dot(v,M)
      l=pl.absolute(C)>th    
      Csub=C[l]
      ind=pl.find(l)
      for j,c in enumerate(Csub):
        out.append(str(i)+' '+str(ind[j]+i+1)+' '+str(c)+'\n')    
    return out  

  def percolation_network(self,ncores,correlation='both'):
    if correlation not in ['negative','positive','both']:
      raise Exception(
      'correlation must be one of: negative, positive or both')
    G=nx.Graph()
    with open(self.edgelist_file) as edgelist:
      for line in edgelist:
          s=line.split()
          w=float(s[2])
          if correlation=='both':
            G.add_edge(int(s[0]),int(s[1]),weight=abs(w))
          elif correlation=='positive':
            if w>0:
              G.add_edge(int(s[0]),int(s[1]),weight=w)
          elif correlation=='negative':
            if w<0:
              G.add_edge(int(s[0]),int(s[1]),weight=-w)
      logger.info('number of nodes in network: %i', G.number_of_nodes())
      mcs=max(int(0.001*G.number_of_nodes()),2)
      
      th=list(pl.arange(self.min_th,1,0.001))
      (gc,NON,cluster_dic)=self.percolation(G,th,mcs,ncores)
      self.save_data(NON,cluster_dic,gc,th)  
    return

  def save_data(self,NON,cluster_dic,gc,th):
    NONdata=json_graph.node_link_data(NON)
    
    f=open(os.path.join(self.network_dir,'non.json'),'w')
    f.write(json.dumps(NONdata))
    f.close()
    os.path.join(self.network_dir,'cluster_dic.json')
    f=open(os.path.join(self.network_dir,'cluster_dic.json'),'w')
    f.write(json.dumps(cluster_dic))
    f.close()
    
    f=open(os.path.join(self.network_dir,'gc.json'),'w')
    f.write(json.dumps(gc))
    f.close()
  
    f=open(os.path.join(self.network_dir,'th.json'),'w')
    f.write(json.dumps(th))
    f.close()
    return 

#===========================================================================
# Geters
#===========================================================================
  def get_node2voxel(self):
    node2voxel=json.load(open(self.node2voxel_file))
    return node2voxel

  def get_affine(self):
    img=nib.load(os.path.join(self.func_dir,self.name))
    affine=img.get_affine()
    return affine    

  def get_img(self):
    img=nib.load(os.path.join(self.func_dir,self.name))
    return img    

  def get_img_data(self):
    img=nib.load(os.path.join(self.func_dir,self.name))
    D=img.get_data()
    return D    
    
  def get_gc(self):
    gc=json.load(open(self.gc_file))
    return gc     
    
  def get_th(self):
    th=json.load(open(self.th_file))
    return th          

  def get_cluster_dic(self):
    cluster_dic=json.load(open(self.cluster_dic_file))
    return cluster_dic 

  def get_non(self):
    NON=json_graph.node_link_graph(json.load(open(self.non_file)))
    return NON 

  def get_cluster_properties(self,force=False,N_clus=3,op='first',
                             correlation='both'):
    if os.path.exists(self.cluster_properties_file) and not force:
      clusters_properties=json.load(open(self.cluster_properties_file))
    else:
      clusters_properties=mass_path_distance(self,N_clus,op,correlation)
    return clusters_properties

  def get_weak_link_distribution(self,force=False,N_clus=2,mcs=0):
    if os.path.exists(self.weak_link_distribution_file) and not force:
      d=json.load(open(self.weak_link_distribution_file))
    else:
      d=weak_link_distribution(self,Nclus=Nclus,mcs=mcs)
    return d
    
  def get_SubGraph(self,th,nodelist,edge_list=None,correlation='both'):
    if correlation not in ['negative','positive','both']:
      raise Exception(
      'correlation must be one of: negative, positive or both')
    G=nx.Graph()
    if edge_list is None:
      with open(self.edgelist_file) as edge_list:
        for line in edge_list:
          s=line.split()
          w=float(s[2])
          if correlation=='both':
            if (th<=abs(w) and int(s[0]) in nodelist and 
                int(s[1]) in nodelist):
              G.add_edge(int(s[0]),int(s[1])) 
          elif correlation=='positive':
            if th<=w and int(s[0]) in nodelist and int(s[1]) in nodelist:
              G.add_edge(int(s[0]),int(s[1])) 
          elif correlation=='negative':
            if w<=-th and int(s[0]) in nodelist and int(s[1]) in nodelist:
              G.add_edge(int(s[0]),int(s[1])) 
    else:
      for edge in edge_list:
        if edge[0] in nodelist and edge[1] in nodelist:
          G.add_edge(edge[0],edge[1])       
    return G

  def get_Graph(self,th,th_up=1.1,correlation='both'):
    if correlation not in ['negative','positive','both']:
      raise Exception(
      'correlation must be one of: negative, positive or both')
    G=nx.Graph()
    with open(self.edgelist_file) as edge_list:
      for line in edge_list:
        s=line.split()
        w=float(s[2])
        if correlation=='both':
          if th<abs(w)<th_up:
            G.add_edge(int(s[0]),int(s[1]))
        elif correlation=='positive':
          if th<w<th_up:
            G.add_edge(int(s[0]),int(s[1]))          
        elif correlation=='negative':
          if -th_up<w<-th:
            G.add_edge(int(s[0]),int(s[1]))          
    return G
  
  def number_of_nodes(self):
    img_dir=os.path.join(self.func_dir,self.name)
    img = nib.load(img_dir)
    D=img.get_data()
    sh=D.shape  
    if self.mask!=None:
      if os.path.isabs(self.mask):
        img2 = nib.load(self.mask)
      else:
        img2 = nib.load(os.path.join(self.func_dir,self.mask))
      M=img2.get_data()
    else:
      M=pl.ones((sh[0],sh[1],sh[2]))      
    count=0
    for i in range(sh[0]):
      for j in range(sh[1]):
        for k in range(sh[2]):
          if M[i,j,k]!=0 and any(D[i,j,k,:]!=0):
            count=count+1
    return count    

  def prune(self,G,th,mcs,cc_old,ncores):
    """
    Prunes G and returns a list of clusters biggers than mcs 
    """

    r_edges=[(n1,n2) for (n1,n2,w) in G.edges_iter(data=True) if
    w['weight']<th]
    G.remove_edges_from(r_edges)    

    cc=nx.connected_components(G)
    cc100=[cc[j] for j in range(len(cc)) if len(cc[j])>=mcs]

    #####remove small clusters#######
    nodes=set(G.nodes())
    nodes_in_c=set(sum(cc100,[]))
    r_nodes=list(nodes.difference(nodes_in_c))
    G.remove_nodes_from(r_nodes) 
    
    return [G,cc100]

  def percolation_data(self,gc,cc100,nn):
    cc_sizes=[float(len(cc100[j]))/nn for j in range(len(cc100))]
    for j in range(max(len(cc_sizes),len(gc))):
      if j>=len(gc):
        lengs=[len(x) for x in gc]
        gc.append([0 for k in range(max(lengs)-1)])
      if j>=len(cc_sizes):
        gc[j].append(0)
      else:
        gc[j].append(cc_sizes[j])
    return gc

  def update_NON(self,NON,cc100,th,th_old):
    check=True
    nodes=NON.nodes(data=True)
    if len(nodes)==0:
      n=0
      check=False
    else:
      n=max(NON.nodes())+1
    
    for j in range(len(cc100)):
      if len(cc100[j])>0:
        NON.add_node(n,th=th,cc=cc100[j])
        for (node,dat) in nodes:
          if set(dat['cc']).issuperset(set(cc100[j])) and dat['th']==th_old:
            NON.add_edge(node,n)
        if nx.degree(NON,n)==0 and check:
          logger.warn('node in NON without parent. (%i,%i,th=%f,th_old=%f)',
                      n,len(cc100[j]),th,th_old)
        n=n+1
    
    cond=True
    for (node,dat) in nodes:
      if dat['th']==th_old and nx.degree(NON,node)!=2 and node>0 :
        cond=False
  
    for (node,dat) in nodes:
      if dat['th']==th_old:
        if nx.degree(NON,node)==2 and node>0 and cond:
          neigs=NON.neighbors(node)
          hijo=max(neigs)
          NON.node[node]['th']=NON.node[hijo]['th']
          NON.remove_node(hijo)
        elif nx.degree(NON,node)==1 and node==0:
          neigs=NON.neighbors(node)
          NON.node[node]['th']=NON.node[neigs[0]]['th']
          NON.remove_node(neigs[0])
  
    return NON  

  def percolation(self,G,th,mcs,ncores):
    """
    Percolation process of G throu the thresholds th. Returns gc: a list of
    sizes of clusters for a given threshold; NON: Network of networks; 
    cluster_dic: a dictionary of clusters bigger than mcs for each 
    threshold. 
    """
    nn=self.number_of_nodes() #float(G.number_of_nodes())
    #giant components list
    gc=[[]]
    #network of networks          
    NON=nx.DiGraph() 
    cluster_dic={}
    cc100=[G.nodes()]
    for i in range(len(th)):
      logger.debug('threshold: %f', th[i])
      [G,cc100]=self.prune(G,th[i],mcs,cc100,ncores)
      cluster_dic[str(th[i])]=cc100
      gc=self.percolation_data(gc,cc100,nn)
      NON=self.update_NON(NON,cc100,th[i],th[i-1])#[:min(10,len(cc100))]
      
    return (gc,NON,cluster_dic)  

  def nodedistance(self,(n1,n2)):
    """
    node distance in cm. (en general)
    """
    ind1=self.node2voxel[str(n1)]
    ind2=self.node2voxel[str(n2)]
    if len(ind1)==3:
      ind1.append(1)
    if len(ind2)==3:
      ind2.append(1)
    v1=numpy.dot(self.affine, numpy.transpose(ind1))[0:3]
    v2=numpy.dot(self.affine, numpy.transpose(ind2))[0:3]
    d=spatial.distance.euclidean(v1,v2)
    return d

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return itt.izip_longest(fillvalue=fillvalue, *args)

