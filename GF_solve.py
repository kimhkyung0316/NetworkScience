#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 12:00:41 2018

Generate a network from the given network file(edgelist file). Analyze the
network and solve the generating function numerically. Then, we can get the
size of the largest connected cluster numerically. We want compare it to the 
analytic solution.
@author: KimHKyeong
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import os
import math

def read_cal_network(fname,N):
    """
    get degree distribution, mean degree, number of links from the network file
    - fname (str) : path of the network file
    - N (int) : number of nodes
    """
    PkL=[]                 #Degree Distribution list
    mean_k=0               #degree dist of zero degree nodes
    NL=0                   #number of links
    if os.path.getsize(fname) == 0:
        PkL = [1]
        mean_k = 0
        NL = 0
    else:
        fread = open(fname, 'rb')
        G = nx.read_edgelist(fread, delimiter='\t')
        fread.close()               #read network file and make it as G
        
        NL = G.number_of_edges()    
        Nwl = G.number_of_nodes()    #number of nodes which have link        
        ZD = (N-Nwl)/N               #degree distribution of zero-degree nodes
        PkL.append(ZD)
   
        max_d = max(G.degree().values())      #maximum degree value
        DVlist = list(G.degree().values())    #list of degree values

        for i in range(1,max_d):       #get degree distribution and mean degree
            dd = DVlist.count(i)/N
            PkL.append(dd)
            mean_k += i*dd                    
        G.clear()                 
    return PkL,mean_k,NL


def Gen_GF_Fitting(Pk,mean_k):
    """
    Fitting the generating function to get a graph
    - Pk (array) : degree distribution
    - mena_k (real) : mean degree
    """
    X = np.arange(0,1.05,0.01)
    U = sum([j*Pk[j]/mean_k*pow(X,j-1) for j in range(1,len(Pk))])
    plt.plot(X,U)
    plt.plot(X,X)
    plt.xlim(0,1.01)
    plt.ylim(0,1.01)
    plt.grid(0.1)
    plt.show()
  
    
def Root_Find(Pk,mean_k):
    """
    solve the polynomial of generating function to get a solution
    """
    if len(Pk) <= 2:
        root = 1.0
    else:
        g1 = [j*Pk[j]/mean_k for j in range(1,len(Pk))]
        g1[1] = g1[1]-1
        g1.reverse()
        root = np.roots(g1)[-1]
    print(root.real)
    return root.real


def Giant_Com_Size(Pk,root):
    """
    get the largest connected cluster
    - root (real) : a solution from the polynomial
    """
    g0 = sum([Pk[k]*pow(root,k) for k in range(len(Pk))])
    gcs = 1 - g0
    print(gcs)
    return gcs

def Data_Extraction(N,g):
    """
    This function is made for extract the solution and the size of the LCC
    from network ensembles successively
    - N (int) : number of nodes in the network
    - g (real) : degree exponent of the network
    """
    fname = "/Users/KimHKyeong/Desktop/KHK/Python/GCSdata/GCS_Numerical_NL{}g{}.txt".format(math.log(N,10),g*10)
    f = open(fname, 'w')
    
    for j in range(150,151):    
        if j%100==0:
            j=int(j/100)
        else:
            j=j/100    
        print(j)
        Pk_list=[]
        K_MM = 0                      #mean of 200 mean degree values 
        NL = 0
        for i in range(0,200):
            file = "/Users/KimHKyeong/Desktop/KHK/Python/edgelist/NL4g3.0/N10000_gamma3_kmin{}_{}.txt".format(j,i)
            Pk,mean_k,Linknum = read_cal_network(file, N)
            K_MM += mean_k
            NL += Linknum
            Pk_list.append(len(Pk))
            
        leng = max(Pk_list)
        Pk_sum = np.zeros((1,leng))
        
        for k in range(0,200):
            file = "/Users/KimHKyeong/Desktop/KHK/Python/edgelist/NL4g3.0/N10000_gamma3_kmin{}_{}.txt".format(j,k)
            Pk = read_cal_network(file, N)[0]
        
            while len(Pk)<leng:
                Pk = Pk+[0]
            
            Pk_sum = Pk_sum + np.array(Pk)
        
        Pk_avg = Pk_sum/200
        Pk_avg = Pk_avg.tolist()[0]
        print(Pk_avg)
        K_MM = K_MM/200
        NL = NL/200
        
        root = Root_Find(Pk_avg,K_MM)
        gcs = Giant_Com_Size(Pk_avg,root)
        fwr = "{}\t, {}\r\n".format(NL/N, gcs) 
        f.write(fwr)
    f.close()
      
Data_Extraction(10000,3.0)