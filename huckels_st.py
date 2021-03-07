#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 14:46:58 2020

@author: benjamin
"""

import numpy as np


alpha = 0
beta = -1
        

# funtion with matrix input and output of sorted eigen values
def get_evals(Mh):
    evals, evecs = np.linalg.eig(Mh)                        ## get eigenvalues
    e = list(evals)
    for i in range(len(e)):
        e[i] = round(float(e[i]), 3)                        ## round eigenvalues to 3s.f
    
    e = sorted(e)  
    deg = 0
    tmp = 1000                                              ## set temporary value to unrealistic eigenvalue     
    for i in range(len(e)):
        if e[i] != tmp:                                     ## checking not a eigenvalue already printed
            for j in range(i,len(e)):
                
                                                            
                if j == len(e)-1 and e[i] == e[j]:          ## if reached end of list
                
                    deg += 1
                    print("degeneracy: " + str(deg) + ",   eigenvalue: " + str(e[i]))
                    tmp = e[i]
                    deg = 0
                
                    
                                                            ## if not at end of list
                else:
                    if e[i] == e[j]:
                        deg +=1                             ## finding the degeneracy
                                                            
                    else:
                        print("degeneracy: " + str(deg) + ",   eigenvalue: " + str(e[i]))
                        tmp = e[i]
                        deg = 0
                        break
    
            
#make huckel matrix for polyene case
def polyene(n):
    Mh = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            
            if (i-j)**2 == 1:                               #atoms adjacent only if difference between i and j is 1
                Mh[i,j] = beta
                
            elif i == j:
                Mh[i,j] = alpha
                
            else:
                Mh[i,j] = 0
                
    return Mh

#make huckel matrix for length n cyclic polyene
def c_polyene(n):
    Mh = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            
            if (i-j)**2 == 1 or (i-j)**2 == (n-1)**2:       #atoms adjacent if difference between i and j is 1 or n-1
                Mh[i,j] = beta
                
            elif i == j:
                Mh[i,j] = alpha
                
            else:
                Mh[i,j] = 0
                
    return Mh

## return Huckel matrix for any platonic solid, using hard coded adjacency matrices
def platonic_solid(n):
    if n==4:
        Mh = np.matrix([[0,1,0,0],
                        [1,0,1,0],
                        [0,1,0,1],
                        [0,0,1,0]])
    elif n==6:
        Mh = np.matrix([[0,1,1,1,1,0],
                        [1,0,1,1,0,1],
                        [1,1,0,0,1,1],
                        [1,1,0,0,1,1],
                        [1,0,1,1,0,1],
                        [0,1,1,1,1,0]])
    elif n==8:
        Mh = np.matrix([[0,0,0,0,0,1,1,1],
                        [0,0,0,0,1,0,1,1],
                        [0,0,0,0,1,1,0,1],
                        [0,0,0,0,1,1,1,0],
                        [0,1,1,1,0,0,0,0],
                        [1,0,1,1,0,0,0,0],
                        [1,1,0,1,0,0,0,0],
                        [1,1,1,0,0,0,0,0]])
    elif n==12:
        Mh = np.matrix([[0,1,0,0,1,0,1,1,0,0,1,0],
                        [1,0,0,0,1,1,0,1,0,1,0,0],
                        [0,0,0,1,0,0,1,0,1,0,1,1],
                        [0,0,1,0,0,0,1,1,0,1,0,1],
                        [1,1,0,0,0,1,0,0,1,0,1,0],
                        [0,1,0,0,1,0,0,0,1,1,0,1],
                        [1,0,1,1,0,0,0,1,0,0,1,0],
                        [1,1,0,1,0,0,1,0,0,1,0,0],
                        [0,0,1,0,1,1,0,0,0,0,1,1],
                        [0,1,0,1,0,1,0,1,0,0,0,1],
                        [1,0,1,0,1,0,1,0,1,0,0,0],
                        [0,0,1,1,0,1,0,0,1,1,0,0]])
    elif n==20:
        Mh = np.matrix([[0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                        [0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
                        [1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0],
                        [1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0],
                        [1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0],
                        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0],
                        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0],
                        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1],
                        [0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0],
                        [0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0],
                        [0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1],
                        [0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0],
                        [0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0],
                        [0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0],
                        [0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1],
                        [0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0],
                        [0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0],
                        [0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0],
                        [0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0],
                        [0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0]])
    
    Mh = beta*Mh
    return Mh
    
    
#get input from user
n = int(input("number of sp2 atoms: "))             #number of atoms
t = int(input("type of molecule 1 = polyene, 2 = cyclic polyene, 3 = platonic solid: "))         #input type of molecule

if t == 1:
    get_evals(polyene(n))

elif t == 2:
    get_evals(c_polyene(n))

elif t == 3:
    get_evals(platonic_solid(n))

else:
    print("invalid type")