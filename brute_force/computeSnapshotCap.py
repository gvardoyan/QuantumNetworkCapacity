#!/usr/bin/python
import queue
from collections import Counter
import sys

# Compute the capacity of a single snapshot of an undirected graph.
# To do this, we first have to find all paths from s to t (we'll use BFS for this).
# Then, construct all sets of disjoint paths and compute the rate achieved by each, then output the max rate.

# assume that nodes[0] is src, nodes[end] is dst
# doing Breadth First Search
# return only the set of "good" paths
def searchPaths(nodes, edges):
  s = nodes[0]
  t = nodes[-1]
  
  # a list to store all valid paths
  allPaths = list()
  
  # a queue to keep track of all paths
  q = queue.Queue(maxsize=0)
  q.put([s])

  while not q.empty():
    cPath = q.get()
    lnode = cPath[-1]
    # get the neighbors of the last node in the path
    neighbors1 = [x for x, y in edges if y == lnode]
    neighbors2 = [x for y, x in edges if y == lnode]
    neighbors = neighbors1 + neighbors2
    # keep just the neighbors that don't already show up in the path
    # and create a new path for each of these neighbors
    for n in neighbors:
      if not n in cPath:
        newPath = list(cPath)
        newPath.append(n)
        if n == t:
          allPaths.append(newPath)
        else:
          q.put(newPath)
  
  return allPaths
    
# given the set of all paths in a graph from s to t,
# output maximal sets of edge-disjoint paths
def getDisjointPaths(paths):
  # first, determine whether each pair of paths is disjoint
  # and store the result in a dict
  disjointInfo = {}
  numPaths = len(paths)
  for ind1 in range(numPaths):
    for ind2 in range(ind1+1,numPaths):
      res = areTwoPathsDisjoint(paths[ind1],paths[ind2])
      disjointInfo[ind1,ind2] = res
  
  # now, construct the set of disjoint paths
  disjointPathSets = list()
  for ind in range(numPaths):
    # create a queue to keep track of disjoint sets with this path
    dsq = queue.Queue(maxsize=0)
    dsq.put([ind])
    while not dsq.empty():
      pathSet = dsq.get()
      # get the set of paths that are disjoint with this pathSet
      disjPaths = getPathsDisjointWithPathSet(pathSet,paths,disjointInfo)

      if not disjPaths:
        # if this paths set is not a duplicate, add it to the list
        if not duplicatePathsSet(pathSet,disjointPathSets):
          disjointPathSets.append(pathSet)
        continue
      for dpath in disjPaths:
        newDisjPathSet = list(pathSet)
        newDisjPathSet.append(dpath)
        dsq.put(newDisjPathSet)
    
  return disjointPathSets
 
# given a path set and the set of all paths, get the set of paths
# that are disjoint with the path set; disjointInfo is provided
def getPathsDisjointWithPathSet(pathSet,paths,disjointInfo):
  disjPaths = list()
  # don't need to look at paths that are already in pathSet
  allPaths = range(0,len(paths))
  relPaths = [p for p in allPaths if p not in pathSet]
  for pIdx in relPaths:
    pIsDisj = 1
    for dpIdx in pathSet:
      if pIdx != dpIdx:
        # the smaller idx always comes first in the dict
        idx1 = min(pIdx,dpIdx)
        idx2 = max(pIdx,dpIdx)
        if disjointInfo[idx1,idx2] == 0:
          pIsDisj = 0
          break
    if pIsDisj == 1:
      disjPaths.append(pIdx)
        
  return disjPaths

# check if a paths set is a duplicate of a list of paths sets
def duplicatePathsSet(paths, pathsSet):
  for pSet in pathsSet:
    if Counter(pSet) == Counter(paths):
      return 1
  return 0
    
# determine if two paths are disjoint
def areTwoPathsDisjoint(path1, path2):
  numNodes = len(path1)
  # don't need to check for the last node, since it's t;
  # sufficient to check just t's left neighbor
  for i in range(numNodes-1):
    node = path1[i]
    if node in path2:
      # get the neighbor to the right of node, in path1
      rightN = path1[i+1]
      # check if the neighbor is in path2 as well
      if rightN in path2:
        # get the indices of the node and the neighbor in path2
        j = path2.index(node)
        k = path2.index(rightN)
        if abs(j-k) == 1:
          # the paths are not disjoint
          return 0
  return 1

def compute_snapshot_cap(nodes, edges, qvals):
  s = nodes[0]
  t = nodes[-1]
  # get all paths in the graph
  allPaths = searchPaths(nodes,edges)
  
  # get sets of edge-disjoint paths in the graph
  disjointPathSets = getDisjointPaths(allPaths)
  # compute the objective function for each paths set
  objCosts = list()
  for pathSet in disjointPathSets:
    sumSoFar = 0
    for pathNum in pathSet:
      path = allPaths[pathNum]
      costSoFar = 1
      for node in path:
        if node != s and node != t:
          costSoFar = costSoFar*qvals[node-1]
      sumSoFar = sumSoFar + costSoFar
    objCosts.append(sumSoFar)

  # output the optimal objective value, and each paths set that yields it
  if not objCosts:
    return 0
  else:
    m = max(objCosts)
    return m
