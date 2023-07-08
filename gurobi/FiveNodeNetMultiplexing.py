#!/usr/bin/python
from gurobipy import *
from call_optimizer import optimize
from computeSnapshotCap import compute_snapshot_cap
import itertools
import numpy
from collections import defaultdict
import math
import sys

# number of nodes
numNodes = 5
# entanglement generation success probabilities
pvals = [0.5679, 0.5179, 0.4723, 0.2479, 0.7839, 0.5423, 0.4308]
nodes = range(0,numNodes)
# create the five-node (multiplexed) network
# Note: VERY IMPORTANT that for any edge with s or t, s appears first and t appears second. Otherwise, the results will be different!
edges = [(0,1), (0,2), (0,3), (0,4), (1,2), (1,4), (3,4)]
# specify the edge capacities
edgeCapVals = [2, 4, 3, 1, 5, 3, 2]

# BSM measurement success probabilities
# for all nodes except s and t
qvals = [0.5, 0.27, 0.64]

# file index to write to
k = str(sys.argv[1])
# index of iterator
nStart = int(sys.argv[2])
nEnd = int(sys.argv[3])
# path to file where snapshot capacity data will be stored
filename = "FiveNodeExpMultiplexing/k"+str(k)+".txt"
templ = '{0} {1:.9e}\n'

numNewNodes = 0;
edgeCaps = {}
idx = 0;
for edge in edges:
  numLinks = edgeCapVals[idx]
  edgeCaps[edge] = numLinks
  if numLinks > 1:
    numNewNodes = numNewNodes + numLinks
  idx = idx + 1
  
# node t needs to be relabeled in all data structures
t = numNodes-1+numNewNodes
nodes = [t if x==numNodes-1 else x for x in nodes]
edges = [(x,t) if (x,y)==(x,numNodes-1) else (x,y) for (x,y) in edges]
for edge in list(edgeCaps):
  # just need to check the "right" node in the edge, since node n-1
  # is last in the node list, so there are no edges of form (n-1,i)
  rnode = edge[1]
  if rnode == numNodes-1:
    newEdge = (edge[0],t)
    edgeCaps[newEdge] = edgeCaps.pop(edge)

# create pval dictionary, indexed by edges
pValDict = {}
for idx in range(len(edges)):
  pValDict[edges[idx]] = pvals[idx]

# make a copy of the original set of edges (need to keep track of singletons)
singleEdges = list(edges)
# make another copy so we can use pValDict
origEdges = list(edges)

# in the multiplexing variant, we transform the original graph
# as follows: for each edge (i,j) that has (strictly) more than one link
# (or equivalently, a higher-than-unit capacity), we create a new node
# for each such link and connect it to i and j.
# E.g.: if (i,j) has two links, create nodes k, l and edges (i,k), (k,j),
# (i,l), (l,j), and remove edge (i,j).
nextNodeNum = numNodes-1
newEdges = list(edges)
# keep track of the new nodes and which original edges they correspond to
newNodeDict = defaultdict(list)
i = 0;
for edge in edges:
  pVal = pvals[i]
  numLinks = edgeCaps[edge]
  if numLinks > 1:
    # first, remove this edge from the list of edges
    newEdges.remove(edge)
    # this isn't a singleton, so remove it from that list as well
    singleEdges.remove(edge)
    # next, create a new node for each link and add the new edges
    for link in range(numLinks):
      node = nextNodeNum + link
      newNodeDict[edge].append(node)
      # add this is a new node (but not in the 0th or last position,
      # since those are nodes s and t)
      nodes.insert(-1,node)
      qvals.append(1)
      newEdges.append((edge[0],node))
      newEdges.append((node,edge[1]))
    
    nextNodeNum = node+1
  i = i+1
    
edges = newEdges

#For the special (multiplexing-enabled) edges, we turn q's on/off to activate/deactivate these edges, and for the other (capacity = 1) edges, we add/remove them from the snapshot
# prepare the iterator
allCombos = []
for edgeCap in edgeCapVals:
  allCombos.append(range(edgeCap+1))
  
all_combos = itertools.product(*allCombos)
combo_iterator = itertools.islice(all_combos,nStart,nEnd)
for combo in combo_iterator:
  # create the snapshot based on the values in the combo
  snapEdges = list(edges)
  
  snapQs = list(qvals)
  eIdx = 0
  snapProb = 1
  snapWeight = 1
  for c in combo:
    # get the original edge
    origEdge = origEdges[eIdx]
    isSingleEdge = origEdge in singleEdges
    if isSingleEdge:
      if c > 0:
        snapProb = snapProb*pValDict[origEdge]
      else:
        snapEdges.remove(origEdge)
        snapProb = snapProb*(1-pValDict[origEdge])
    else:
      # this edge has a capacity higher than one;
      # get the list of nodes corresponding to this edge
      edgeCap = edgeCapVals[eIdx]
      qNodes = newNodeDict[origEdge]
      pVal = pValDict[origEdge]
      # activate/deactive the right number of parallel edges
      numInactive = edgeCap-c
      for i in range(numInactive):
        qNodeLabel = qNodes[i]
        snapQs[qNodeLabel-1] = 0
      
      snapWeight = snapWeight*math.comb(edgeCap,c)
      snapProb = snapProb*(pVal ** c)*((1-pVal) ** numInactive)

    eIdx = eIdx+1
  
  snapProb = snapWeight*snapProb
  # run the optimization and record the result
  f = open(filename,"a+")
  # Gurobi output
  obj = optimize(nodes,snapEdges,snapQs)
  f.write(templ.format(obj, snapProb))
  f.close()
