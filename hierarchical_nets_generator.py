"""
This module is to generate hierarchical networks 
using the framework of networkX.

You can import this model and call functions to generate graphs;
you can also run this script to check how it works
"""


import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def findij(n):
  n = int(n)
  if n<0:
    return NULL
  ij = [0, 0]
  ij[0] = 1
  ij[1] = 0
  while n%2==0:
    n = n/2
    ij[0] += 1
  ij[1] = (n-1)/2
  return ij


def find_hn3_edges(i, N, period):
  edges = []
  N_hn = 2**int(np.log2(N))
  if i<0 or i>N or (period==True and i==0):
    return 
  if i==0:
    return [(0,1), (0, N_hn/2), (0, N_hn)]
  if period:
    if i>1:
      edges.append((i, i-1))
    else:
      edges.append((i, N_hn))
    if i<N:
      edges.append((i, i+1))
    if i==N_hn/2:
      edges.append((i, N_hn))
    elif i==N_hn:
      edges.append((i, 1))
    else:
      ij = findij(i)
      if ij[1]%2 == 0:
        temp = 2**(ij[0]-1)*(2*(ij[1]+1)+1)
      else:
        temp = 2**(ij[0]-1)*(2*(ij[1]-1)+1)
      if temp<N:
        edges.append((i, temp))
  else:
    if i>0:
      edges.append((i, i-1))
    if i<N:
      edges.append((i, i+1))
    if i==N_hn/2:
      edges.append((i, N_hn))
    else:
      ij = findij(i)
      if ij[1]%2 == 0:
        temp = 2**(ij[0]-1)*(2*(ij[1]+1)+1)
      else:
        temp = 2**(ij[0]-1)*(2*(ij[1]-1)+1)
      if temp<N:
        edges.append((i, temp))
   
  return edges
  
def find_hn5_edges(i, N, period):
  edges = find_hn3_edges(i, N, period)
  ij = findij(i)
  N_hn = 2**int(np.log2(N))
  if period:
    temp = i-2**(ij[0]-1)
    if temp<1:
      edges.append((i, N_hn))
    else:
      edges.append((i, temp))
    temp = i+2**(ij[0]-1)
    if temp<(N+1):
      edges.append((i, temp))
  else:
    temp = i-2**(ij[0]-1)
    edges.append((i, N_hn))
    temp = i+2**(ij[0]-1)
    if temp<(N+1):
      edges.append((i, temp))
  return edges

def find_hnnp_edges(i, N, period):
  edges = []
  ij = findij(i)
  N_hn = 2**int(np.log2(N))
  if period:
    if i>1:
      edges.append((i, i-1))
    else:
      edges.append((i, N_hn))
    if i<N:
      edges.append((i, i+1))
    elif i==N_hn:
      edges.append((i, 1))
    if ij[1]%2 == 0:
      temp = i + 3*2**(ij[0]-1)
    else:
      temp = i - 3*2**(ij[0]-1)
    if temp==0:
      edges.append((i, N_hn))
    elif temp < N+1:
      edges.append((i, temp))
  else:
    if i>0:
      edges.append((i, i-1))
    if i<N:
      edges.append((i, i+1))
    if ij[1]%2 == 0:
      temp = i + 3*2**(ij[0]-1)
    else:
      temp = i - 3*2**(ij[0]-1)
    if temp<N+1:
        edges.append((i, temp))
  return edges

def find_hp6_edges(i, N, period):
  edges = find_hnnp_edges(i, N, period)
  ij = findij(i)
  N_hn = 2**int(np.log2(N))
  if period:
    temp = i-2**(ij[0]-1)
    if temp<1:
      edges.append((i, N_hn))
    else:
      edges.append((i, temp))
    temp = i+2**(ij[0]-1)
    if temp<(N+1):
      edges.append((i, temp))
  else:
    temp = i-2**(ij[0]-1)
    edges.append((i, temp))
    temp = i+2**(ij[0]-1)
    if temp<(N+1):
      edges.append((i, temp))
  return edges

def hn3(N, k = 0, period = True):
  N, k = int(N), int(k)
  if k != 0:
    N = 2**k
  G = nx.Graph()
  if period:
    G.add_nodes_from(np.arange(1, N+1))
  else:
    G.add_nodes_from(np.arange(0, N+1))
  for i in range(N):
    edges = find_hn3_edges(i+1, N, period)
    G.add_edges_from(edges)
  return G

def hn5(N, k = 0, period = True):
  N, k = int(N), int(k)
  if k != 0:
    N = 2**k
  G = nx.Graph()
  if period:
    G.add_nodes_from(np.arange(1, N+1))
  else:
    G.add_nodes_from(np.arange(0, N+1))
  for i in range(N):
    edges = find_hn5_edges(i+1, N, period)
    G.add_edges_from(edges)
  return G

def hnnp(N, k=0, period=True):
  N, k = int(N), int(k)
  if k!=0:
    N = 2**k
  G = nx.Graph()
  if period:
    G.add_nodes_from(np.arange(1, N+1))
  else:
    G.add_nodes_from(np.arange(0, N+1))
  for i in range(N):
    edges = find_hnnp_edges(i+1, N, period)
    G.add_edges_from(edges)
  return G

def hp6(N, k=0, period=True):
  N, k = int(N), int(k)
  if k!=0:
    N = 2**k
  G = nx.Graph()
  if period:
    G.add_nodes_from(np.arange(1, N+1))
  else:
    G.add_nodes_from(np.arange(0, N+1))
  for i in range(N):
    edges = find_hp6_edges(i+1, N, period)
    G.add_edges_from(edges)
  return G

##### test code ####
if __name__ =="__main__":
  hn_size = 16
  plt.figure(1, figsize=(10, 10))
  plt.subplot(221)
  G1 = hn3(hn_size)
  nx.draw_circular(G1, with_labels=True)
  plt.title("HN3")

  plt.subplot(222)
  G2 = hn5(hn_size)
  nx.draw_circular(G2, with_labels=True)
  plt.title("HN5")

  plt.subplot(223)
  G3 = hnnp(hn_size)
  nx.draw_circular(G3, with_labels=True)
  plt.title("HNNP")

  plt.subplot(224)
  G4 = hp6(hn_size)
  nx.draw_circular(G4, with_labels=True)
  plt.title("HP6")

  plt.show()
