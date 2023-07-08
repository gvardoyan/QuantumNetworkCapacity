#!/usr/bin/python
from gurobipy import *
from call_optimizer import optimize
import itertools
import numpy

# Node label data:
#  0 = Delft; 1 = Rotterdam; 2 = Leiden; 3 = Amsterdam; 4 = Hilversum; 5 = Utrecht; 6 = Amersfoort; 7 = Almere; 8 = Wageningen; 9 = Lelystad; 10 = Apeldoorn; 11 = Arnhem; 12 = Nijmegen; 13 = Zwolle; 14 = Deventer; 15 = Zutphen; 16 = Enschede

# number of nodes
numNodes = 17
# the number of "active" (present) edges in a snapshot;
# this is just one way to parallelize the overall computation
# (note that you will re-run this script for each k in [0,len(edges)])
k = int(sys.argv[1])
# path where snapshot capacity data will be stored
filename = "SURFNetExp/k"+str(k)+".txt"

# entanglement generation success probabilities
pvals = [0.4152, 0.2199, 0.0557, 0.0358, 0.2240, 0.1501, 0.1661, 0.1763, 0.1898, 0.1176, 0.0506, 0.0425, 0.2756, 0.0620, 0.1117, 0.2926, 0.1001, 0.1149, 0.0240, 0.0568]
nodes = range(0,numNodes)
# Create the pruned SURF network
# Note: VERY IMPORTANT that for any edge with s or t, s appears first and t appears second. Otherwise, the results will be different!
edges = [(0,1), (0,2), (2,3), (1,5), (3,4), (3,7), (4,5), (4,7), (5,6), (7,9), (6,8), (8,12), (11,12), (12,15), (10,11), (10,14), (9,13), (13,14), (13,16), (15,16)]
numEdges = len(edges)

# BSM measurement success probabilities
# for all nodes except s and t
qvals = [0.87, 0.74, 0.79, 0.62, 0.73, 0.98, 0.77, 0.76, 0.62, 0.74, 0.81, 0.84, 0.7, 0.68, 0.99]

edgeLabels = range(0,numEdges)
combo_iterator = itertools.combinations(edgeLabels,k)
n = int(sys.argv[2]);
cIdx = 1;
templ = '{0} {1:.9f}\n'
for combo in combo_iterator:
  if cIdx >= n:
    f = open(filename,"a+")
    snapEdges = [edges[i] for i in combo]
    snapProbsSuc = [pvals[i] for i in combo]
    snapProbsFail = [1 - pvals[i] for i in edgeLabels if i not in combo]
    snapProbs = snapProbsSuc + snapProbsFail
    snapP = numpy.prod(snapProbs)
    # Gurobi output
    obj = optimize(nodes,snapEdges,qvals)
    f.write(templ.format(obj, snapP))
    f.close()
  cIdx = cIdx + 1
