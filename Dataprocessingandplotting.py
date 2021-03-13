#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 20:36:09 2020

@author: benjamin
"""

from numpy import *
import matplotlib.pyplot as plt
import re
import os
from scipy.optimize import curve_fit
import sys
import math

##define a function of the form expected to fit the data
def fit(X, eo, kr, ro, kt, to):
    x,y = X
    return eo + 0.5*kr*((x - ro)**2) + 0.5*kt*((y - to)**2)
                                       
##Make empty lists to be filled with data from files
r = []
theta = []
E = []

rp = []
thetap = []
Ep = []

dir = str(input("Directory name:   "))
for filename in os.listdir(dir):
    r_t_E = re.findall(r"[-+]?\d*\.\d+|\d+", filename)                     #extract floats from filename    
    r.append(float(r_t_E[1]))
    theta.append(float(r_t_E[2]))
    
    ##read Energy from each file into a separate list
    f = open(dir + filename, "r")
    for line in f:
        for line in f:
            if "SCF Done:" in line:
                l = line.split()
                E.append(float(l[4]))
    f.close()

## curvefitting data
X = (r, theta)
popt, pcov = curve_fit(fit, X, E)
for i in range(len(r)):                                                 ##extract data that is close to the minimum energy only
    if r[i] > popt[2]-0.15 and r[i] < popt[2]+0.15:
        if theta[i] > popt[4]-15 and theta[i] < popt[4]+15:
            rp.append(r[i])
            thetap.append(theta[i])
            Ep.append(E[i])

X = (rp, thetap)                                                        ##find new parameters using the extracted data only
popt, pcov = curve_fit(fit, X, Ep)

## plot the data
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlim([popt[2]-0.5, popt[2]+0.5])
ax.set_ylim([popt[4]-30, popt[4]+30])
#ax.plot_trisurf(rp, thetap, fit(X, *popt), cmap=plt.cm.CMRmap)
ax.plot_trisurf(r, theta, E, cmap="jet")
#ax.scatter(r, theta, E)
ax.set_xlabel("r/Angstroms")
ax.set_ylabel("Theta/degrees")
ax.set_zlabel("Energy/Hartrees")
ax.set_title("Energy surface")
plt.show()

##Find and print vibrational modes
##unit conversions
kr = popt[1]*4.35974*10**2
kt = popt[3]*((180/math.pi)**2)*4.35974*10**-18
r0 = popt[2]*10**-10

def sym_stretch(kr):
    return ((1/(2*math.pi)) * math.sqrt(kr/(2*1.66*10**-27)))* 3.3356*10**-11

def bend(kt, r0):
    return ((1/(2*math.pi))*math.sqrt(kt/(0.5*(r0**2)*1.66*10**-27)))*3.3356*10**-11

print ("Symmetric stretch:  " + str(sym_stretch(kr)) + " cm−1")
print ("Bend:  " + str(bend(kt, r0)) + " cm−1")


