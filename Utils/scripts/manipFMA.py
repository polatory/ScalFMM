# -*- coding: utf-8 -*-
"""
Created on Sun May 11 16:32:08 2014

@author: coulaud
"""

import numpy as np
from random import *
import matplotlib.pyplot as plt
#
#  Manipulation de fichier FMA
#
#Npart =0 

#def readFMA(inputFile,N) 
def readFMA(N, BoxSize) :
#
#part BoxSize,BoxCentre):
#
    print __doc__
    inputFile="/Users/coulaud/Dev/src/ScalFMM/scalfmmT/Data/prolate50.fma"
    fichier = open(inputFile, "r")
    A= fichier.readline().split()
    N = int(A[0])
    BoxSize = float(A[1])
    fichier.close()
    return N
    
def writeFMA(N, xyzq, outFile) :
#
#part BoxSize,BoxCentre):
#
    print __doc__
#    outFile="/Users/coulaud/Dev/src/ScalFMM/scalfmmT/Data/prolate50.fma"
    fichier = open(outFile, "r")
 #   A= fichier.readline().split()
 #   N = int(A[0])
 #   BoxSize = float(A[1])
    fichier.close()


def build1DDistribution(a, alpha,N):
    xyzq = np.zeros((N, 4))
    x = np.zeros(N)
    y = np.zeros(N)
    z = np.random.random_sample(N)
    
    sa = sin(alpha)
    ca = cos(alpha)
    print a,alpha,ca,sa
    for i in range(0, N):
        x[i] = a*sa*z[i]*cos(z[i])
        y[i] = a*sa*z[i]*sin(z[i])
        z[i] = a*ca*z[i]
        print i, x[i],y[i], z[i]
    return x,y,z

