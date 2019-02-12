#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 19:12:55 2017

Create ensembles of configuration-model scale-free networks under the given condition.

@author: KimHKyeong
"""

import networkx as nx
#import sys
import math
import numpy as np
#import matplotlib.pyplot as plt

def powerlaw_sequence(N,g):
    """
    extract degree sequence which follows power-law
    - N (int) : number of nodes
    - g (real) : degree exponent for power-law function
    - kmin (real) : minimum degree value
    - kmax (real) : maximum degree value
    """
    while(True):
        for x in range(N):
            while(True):
                r = np.random.random()
                k = (kmin**(-(g-1))+(kmax**(-(g-1))-kmin**(-(g-1)))*\
                     (1-r))**(-1./(g-1))  
                #extract a degree k under power-law form
                if k>kmin and k<kmax:
                    break
            s[x]=((int)(k))
        if sum(s)%2==0:
            break
        
def config_net(N,g,kmin,Ens,s):
    """
    Generate a network ensemble under the given condition
    - Ens (int) : number of networks in a ensemble
    - s (array) : degree sequence
    """
    fname2 = "/Users/KimHKyeong/Desktop/New/networkinfo/g{}_kmin{}.txt".format(g*10,kmin)
    fwr2 = open(fname2, 'w')
    print(kmin)
    
    for en in range(Ens):

        G= nx.configuration_model(s)  #generate a network using degree sequence
        G= nx.Graph(G)                #remove overlapped edge(link)
        G.remove_edges_from(G.selfloop_edges())  #remove self-loop
	
        fname = "/Users/KimHKyeong/edgelist/N%d_gamma%g_kmin%g_%d.txt" %(N, g, kmin, en)
        fwr = open(fname, 'wb')
        nx.write_edgelist(G, fwr, delimiter='\t') #write an edgelist of the network
        
        degree_values = list(G.degree().values())

        meank = 1.*sum(degree_values)/(2*N)
        maximum_k = max(degree_values)
        minimum_k = min(degree_values)
        Gc = max(nx.connected_component_subgraphs(G), key = len)
        giant = Gc.number_of_nodes()         
    
        wrt2 = "%g\t, %g\t, %g\t, %g\r\n" %(meank, maximum_k, minimum_k, giant)
        fwr2.write(wrt2)    #write the network's information
        
        fname3 = "/Users/KimHKyeong/degreedist/kmin{}_{}.txt".format(kmin, en)
        fwr3 = open(fname3, 'w')   
        for i in range(maximum_k+1):
            Pk = degree_values.count(i)/(N)
            wrt3= "%d\t, %g\r\n" %(i,Pk)
            fwr3.write(wrt3)         #write the degree distribution
        fwr3.close()
    fwr2.close()    
    G.clear()


N=10000
s = np.zeros(N, dtype=np.int)
g=3.0
Ens = 200
kmax=math.sqrt(N)
for i in range(1,201):
    kmin = i/100
    powerlaw_sequence()
    config_net(N,g,kmin,Ens,s)