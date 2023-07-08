#!/usr/bin/python
from gurobipy import *
from call_optimizer import optimize
from computeSnapshotCapMult_revised import compute_snapshot_cap_mult
import itertools
import numpy
import pdb

# number of nodes
numNodes = 5
# entanglement generation success probabilities
pvals = [0.5679, 0.5179, 0.4723, 0.2479, 0.7839, 0.5423, 0.4308]
nodes = range(0,numNodes)
# create the five-node multiplexed network
# Note: VERY IMPORTANT that for any edge with s or t, s appears first and t appears second. Otherwise, the results will be different!
edges = [(0,1), (0,2), (0,3), (0,4), (1,2), (1,4), (3,4)]
# specify the edge capacities
edgeCapVals = [2, 4, 3, 1, 5, 3, 2]

# BSM measurement success probabilities for all nodes except s and t
qvals = [0.5, 0.27, 0.64]

# following are parameters used to parallelize the overall computation:
# file index to write to
k = int(sys.argv[1])
# start index of iterator
nStart = int(sys.argv[2])
# end index of iterator
nEnd = int(sys.argv[3])

# path where snapshot capacity data will be stored
filename = "FiveNodeExpMultiplexing/BruteForce/k"+str(k)+".txt"
templ = '{0} {1:.9e}\n'
# prepare the iterator
allCombos = []
for edgeCap in edgeCapVals:
  allCombos.append(range(edgeCap+1))
  
all_combos = itertools.product(*allCombos)
combo_iterator = itertools.islice(all_combos,nStart,nEnd)
for combo in combo_iterator:
  print(combo)
  # compute the probability of this snapshot
  snapProb = 1
  snapWeight = 1
  eIdx = 0
  for c in combo:
    pVal = pvals[eIdx]
    edgeCap = edgeCapVals[eIdx]
    numInactive = edgeCap-c

    snapWeight = snapWeight*math.comb(edgeCap,c)
    snapProb = snapProb*(pVal ** c)*((1-pVal) ** numInactive)
      
    eIdx = eIdx+1

  snapProb = snapWeight*snapProb
  # run the optimization and record the result
  f = open(filename,"a+")
  # Brute-force output
  obj = compute_snapshot_cap_mult(nodes,edges,qvals,combo)
  f.write(templ.format(obj, snapProb))
  f.close()
