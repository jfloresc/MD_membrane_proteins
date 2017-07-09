#!/usr/bin/env python
# Calculate allosteric networks in MD runs based on Dokholyan N.V. Chem. Rev. 2016, 116, 6463-6487.

from __future__ import print_function
import matplotlib.pyplot as plt
import networkx as nx
import os
import numpy as np
import mdtraj as md
from itertools import combinations
import pylab

def main():
  CUTOFF = 0.75 # 0.75 nm 
  OCCUPANCY = 0.5

  t = md.load('All_T310_restart_16_25.dcd',top='ADP65_ATP35_b0.2_ftrue_T310_16_restart/T0EG5.pdb')
  pdbref = md.load('ADP65_ATP35_b0.2_ftrue_T310_16_restart/T0EG5.pdb')
  topology = pdbref.topology
  atoms_contact = [a.index for a in topology.atoms if (a.name == 'CA' and a.residue.name == 'GLY') or (a.name == 'CB' and a.residue.name != 'GLY')]
  at_to_res_map = {a.index: a.residue.index for a in topology.atoms if (a.name == 'CA' and a.residue.name == 'GLY') or (a.name == 'CB' and a.residue.name != 'GLY')}
  #atoms_contact1 = topology.select('not resname GLY and name CB or resname GLY and name CA')
  pairs = np.array([(i,j) for (i,j) in combinations(atoms_contact, 2)])
  #print(atoms_contact)
  #print(at_to_res_map)
  dist = md.compute_distances(t, pairs, periodic=True)
  life_contacts = dict()

  for i in xrange(t.n_frames):
    contacts = pairs[dist[i]< CUTOFF] 
    for ipair in xrange(len(contacts)):
      key = tuple(contacts[ipair])
      if key in life_contacts:
        life_contacts[key] += 1
      elif key[::-1] in life_contacts:
        life_contacts[key] += 1
      else:
        life_contacts[key] = 1

  #print(life_contacts) 
   
  correl = np.loadtxt('covariance.dat')

  G = nx.Graph(name = 'Dynamic Network')
  
  for key in life_contacts:
    fraction = life_contacts[key]/float(t.n_frames)
    if fraction >= OCCUPANCY:
      i = at_to_res_map[key[0]]
      j = at_to_res_map[key[1]]
      w = 1.0 - np.fabs(correl[i,j])
      G.add_edge(i, j, weight= w)
      #print(i, j, w)
    #else:
    #  print(i, j, w, fraction)

  n_nodes = len(G)
  print("Number of nodes (after first filter):", n_nodes)

  connected_components = nx.connected_components(G)

  print('Number of connected components : ', \
        len(list(connected_components) ))

  print('Length of the connected components : ', \
       [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)])

  main_component = max(nx.connected_component_subgraphs(G), key=len)

  print("Average shortest path in G's main connected component : ", \
       nx.average_shortest_path_length(main_component,weight='weight'))


#  p=nx.shortest_path_length(G,weight='weight')
#  for key in p:
#    pp = p[key]
#    for k in pp:
#      print(key, k, pp[k])
#  print(p[53][52])
#  plt.figure()
#  pos = nx.nx_pydot.graphviz_layout(G, prog = 'dot')
#  nx.draw(G,pos, with_labels= True)
#  plt.show()
 
  steps = 40
  dc = 1.0/float(steps)
  cutoff_w = dc 
  for i in xrange(steps-1):
    removed_edges = [(gene1, gene2) for (gene1, gene2, w) in G.edges(data=True) if w['weight'] < cutoff_w]
    G.remove_edges_from(removed_edges)
    #connected_components = nx.connected_components(G)
    #print('Number of connected components : ', \
    #  len(list(connected_components)))

    #print('Length of the connected components : ', \
    #  [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)])

    main_component = max(nx.connected_component_subgraphs(G), key=len)
    print(len(main_component)*100.0/float(n_nodes), '% Correlation: ', 1.0 - cutoff_w )
   # print('i: ', i, 'cutoff: ', cutoff_w, 'Number of Edges: ', G.number_of_edges(), 'Number of nodes: ', len(G))
    cutoff_w += dc


  CUTOFF_W = 1.0 - 0.475
  G1 = nx.Graph(name = 'Dynamic Network')
  for key in life_contacts:
    fraction = life_contacts[key]/float(t.n_frames)
    if fraction >= OCCUPANCY:
      i = at_to_res_map[key[0]]
      j = at_to_res_map[key[1]]
      w = 1.0 - np.fabs(correl[i,j])
      if w < CUTOFF_W:
        continue
      else:
        G1.add_edge(i, j, weight= w)
  main_component = max(nx.connected_component_subgraphs(G1), key=len)

  print(len(main_component)*100.0/float(n_nodes), '% Correlation: ', 1.0 - CUTOFF_W )

  #plt.figure()
  #pos = nx.nx_pydot.graphviz_layout(G1, prog = 'dot')
  #nx.draw(G1,pos, with_labels= True)
  #plt.show()

  p=nx.shortest_path(G1,weight='weight')
  #p=nx.shortest_path_length(G1,weight='weight')
  for key in p:
    pp = p[key]
    for k in pp:
      print(key+17, k+17, np.array(pp[k])+17)
#  print(p[53][52])
  GMax = nx.Graph(data = main_component, name = 'Main Component')
  new_label = []
  for i in GMax.nodes():
    new_label.append(i+17)
  mapping = dict(zip(GMax.nodes(),new_label))
  #print(new_label)
  G1=nx.relabel_nodes(GMax, mapping, copy=True)
  #plt.figure()
  #pos = nx.nx_pydot.graphviz_layout(GMax, prog = 'fdp')
  #pos = nx.nx_agraph.graphviz_layout(GMax, prog = 'fdp', args='-Elen=0.5 -Gsize=2,2 -Gdpi=300')
  #nx.draw(GMax,pos, with_labels= True)
  #plt.savefig("Gmax_layout.png")

  ##graph_1 = nx.nx_agraph.to_agraph(G1)
  #graph_1.graph_attr.update(landscape='true')#,ranksep='1.0')
  #graph_1.graph_attr.update(size='12,10')#,ranksep='1.0')
  #graph_1.edge_attr.update(dpi=300)
  #graph_1.layout(prog='fdp')
  #graph_1.edge_attr.update(len=0.1)
  #graph_1.draw('test.png',format='png')
  #graph_1.edge_attr.update(labelfontsize=144)#,ranksep='1.0')
  ##graph_1.write('test.dot')
  #fdp  test.dot -Tpng -Elen=0.1 -Nfontsize=72 -o test.png
  #fdp  test.dot -Tpng -Elen=0.2 -Nfontsize=72 -Gdpi=300 -Gsize="6,4" -o test.png 
  #plt.show()
  #fdp  test.dot -Tpng -Elen=0.2  -Nfontsize=72 -Eweight=1.0  -o test.png

  bottleneck = {}
  for n,d in GMax.nodes_iter(data=True):
    Gtemp = nx.Graph(data = main_component, name = 'Main Component')
    Gtemp.remove_node(n)
    #print(len(Gtemp), Gtemp.edges(data=True))  
    n_comp =  len(list(nx.connected_components(Gtemp)))
    if n_comp >= 2:
      bottleneck[n]=[len(c) for c in sorted(nx.connected_components(Gtemp), key=len, reverse=True)]
      #print([len(c) for c in sorted(nx.connected_components(Gtemp), key=len, reverse=True)])
      #print("removed node: ", n)

  for key in sorted(bottleneck):
    print("EG5 Res. index: ", key+17, "index network: ", key, bottleneck[key])


  b = nx.betweenness_centrality(G1, k=None, normalized=True, weight='weight', endpoints=False, seed=None)
  x = [i for i in xrange(17,367)]
  between_array = []
  for i in x:
    try:
      between_array.append(b[i])
    except:
      between_array.append(np.nan)
  fig = plt.figure(figsize=(4,3))
  line1, = plt.plot(x,between_array,lw=1.5,color='k',label='Betweenness')
  leg = plt.legend(handles=[line1])
#leg = plt.legend(handles=[line3])
  for legobj in leg.legendHandles:
    legobj.set_linewidth(1.5)
  plt.legend(loc='best',frameon=True)
  plt.xlabel('Residue Index')
  plt.ylabel(r'Betweenness Centrality')
  plt.xlim((16,367))
  nameout='Betweenness_centrality'
  plt.subplots_adjust(left=0.2,bottom=0.15, top=0.95)
  plt.savefig(nameout+'.png',format='png',dpi=300)
  outf = open('betweenness.txt','w')
  for i,j in zip(x,between_array):
    outf.write('%s %s\n'%(i,j))
  outf.close()


if __name__ == '__main__':
    main()
