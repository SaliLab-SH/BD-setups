"""
To calculate RDF from the ISG distance file.

"""
from __future__ import print_function, division
import IMP.atom
import IMP.algebra
import IMP.rmf
import IMP.core
import RMF
import IMP.container
import IMP.display
import sys
import random
from random import choice, uniform
import math
from math import trunc
import numpy as np
from numpy import zeros, sqrt, pi, mean, arange, histogram
import itertools as it

#----------NE RDF function----------
def NERDF(NA, NB, distance, R, R_NUCLEUS, R_GRANULES, dr, N):
    '''define the RDF function of B (Granules) with respect to A (Nucleus)'''
    # N is number of frames
    edges = arange(R_NUCLEUS + R_GRANULES, R - R_GRANULES + 0.1 * dr, dr)
    num_increments = len(edges) - 1
    radii = zeros(num_increments)
    numberDensity = NB / (4.0 / 3.0 * pi * ((R - R_GRANULES)**3-(R_NUCLEUS+R_GRANULES)**3))

    # Compute histogram over all pairs of distances
    (result, bins) = histogram(distance, bins=edges, normed=False)
    #print(len(result))
    #print(result)
    # Compute the average g(r) over all granules and all frames
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = result[i] / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3)) / numberDensity / N / NA
    print(np.mean(g_average))  # should be approx. 1 for random distributions
    return (radii, g_average)

#---------- Simulation parameters ----------

# I. Parts parameters
L = 50000 # Length of our bounding box, A
R = 20000 # PBC radius, A
R_NUCLEUS = 10000 # NE radius, A
N_GRANULES = 50 # Number of granules
R_GRANULES = 1500 # Radius of granules, A

#----------calculate and write RDF----------
dr = 100
k1 = np.loadtxt("random_ISG-NE_distance.xvg")
#print(len(k1))

distance=[]
for i in range(1, len(k1)):
    for j in range(0,N_GRANULES):
        distance.append(k1[i][j])
#print(distance)
    
f1=open("random_ISG-NE_RDF.xvg","w")

NERDF=NERDF(1, N_GRANULES, distance, R, R_NUCLEUS, R_GRANULES, dr, len(k1))

for x in zip(list(NERDF[0]), list(NERDF[1])):
    f1.write("{0}\t{1}\n".format(*x))

    


