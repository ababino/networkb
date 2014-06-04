# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 17:06:12 2013

@author: andres
"""
import networkx as nx
import pylab as pl

def plot_NON(NON, colors=['b','r','g','c','m','y','k'], max_n_clus=7):
    """
    Plot network of networks.
    """
    pos = {}
    labels = {}
    nodes_dic = {}
    for (node, dat) in NON.nodes(data=True):
        labels[node] = node
        if dat['th'] in nodes_dic:
            nodes_dic[dat['th']].append((node, dat['cc']))
        else:
            nodes_dic[dat['th']] = [(node, dat['cc'])]
    pos = nx.graphviz_layout(NON, prog='dot')
    ymax = max([pos[a][0] for a in pos])
    ymin = min([pos[a][0] for a in pos])
    for key in pos:
        x = NON.node[key]['th']
        y = float(pos[key][0]-ymin)/ymax
        pos[key] = (x, y)
    clus_list = []
    for th in nodes_dic:
        sorted_nodes = sorted([x for x in nodes_dic[th]],
                              key=lambda x: len(x[1]), reverse=True)
        sorted_nodes = sorted_nodes[:min(len(sorted_nodes), max_n_clus)]
        for i in range(len(sorted_nodes)):
            if len(clus_list) >= i + 1:
                clus_list[i].append(sorted_nodes[i][0])
            else:
                clus_list.append([sorted_nodes[i][0]])
    subNON = nx.subgraph(NON, sum(clus_list, []))
    subNON = nx.subgraph(subNON,
                         nx.connected_components(subNON.to_undirected())[0])
    fig = pl.figure()
    i = 0
    for clus in clus_list:
        pos_clus = {}
        for node in clus:
            if node in subNON.nodes():
                pos_clus[node] = pos[node]
        nx.draw_networkx_nodes(subNON, pos_clus, nodelist=pos_clus.keys(),
                               node_color = colors[i], node_size = 100)
        i += 1
    nx.draw_networkx_edges(subNON, pos, arrows=False)
    pl.grid(True)
    pl.xlabel('p')
    ax = pl.gca()
    ax.axes.get_yaxis().set_visible(False)
    #pl.show()
    return fig
