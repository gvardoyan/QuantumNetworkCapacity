#!/usr/bin/python

from gurobipy import *

def optimize(nodes,edges,qvals):
  # Create the optimization model
  m = Model("netFlow")
  # comment out the line below if you want to see gurobi output
  #m.setParam( 'OutputFlag', False )
  # change gurobi solver parameters
  #m.setParam('MIPFocus', 3)
  #m.setParam('Cuts', 0)
  #m.setParam('Method',2)
  #m.setParam('PreMIQCPForm',0)

  s = nodes[0]
  t = nodes[-1]

  # Internal nodes only
  inodes = list(nodes)
  inodes.remove(s)
  inodes.remove(t)

  # Create flow variables, in both directions for each edge
  flowVars = {}
  for edge in edges:
    lnode = edge[0]
    rnode = edge[1]
    fNameDir1 = 'FN'+str(lnode)+'N'+str(rnode)
    fNameDir2 = 'FN'+str(rnode)+'N'+str(lnode)
    gurobiVar = m.addVar(ub=1.0, name=fNameDir1)
    flowVars[lnode,rnode] = gurobiVar
    gurobiVar = m.addVar(ub=1.0, name=fNameDir2)
    flowVars[rnode,lnode] = gurobiVar

  # Create binary variables that will be used to constrain the flow
  xVars = {}
  for node in inodes:
    # gather the list of nodes this node is adjacent to
    adjNodes1 = [x for x, y in edges if node == y]
    adjNodes2 = [x for y, x in edges if node == y]
    adjNodes = adjNodes1 + adjNodes2
    # create the variables, taking special care for s and t
    for adjnode1 in adjNodes:
      if adjnode1 == t:
        continue
      tempNodes = list(adjNodes)
      tempNodes.remove(adjnode1)
      for adjnode2 in tempNodes:
        if adjnode2 == s:
          continue
        varName = 'xN'+str(adjnode1)+'N'+str(node)+'N'+str(adjnode2)
        gurobiVar = m.addVar(vtype=GRB.BINARY, name=varName)
        xVars[adjnode1,node,adjnode2] = gurobiVar
        
  # add the constraints with x_{ijk}
  for j in inodes:
    qval = qvals[j-1]
    for i in nodes:
      relXvars = [] # sum over outgoing edges
      for k in nodes:
        # gather all information about node j
        if (i,j,k) in xVars:
          xvar = xVars[i,j,k]
          fvarin = flowVars[i,j]
          fvarout = flowVars[j,k]
          relXvars.append(xvar)
          m.addConstr(xvar*(fvarin*qval-fvarout) == 0)
        if (k,j,i) in xVars:
          xvar = xVars[k,j,i]
          fvarin = flowVars[k,j]
          fvarout = flowVars[j,i]
          relXvars.append(xvar)
          m.addConstr(xvar*(fvarin*qval-fvarout) == 0)
      # create the following constraints when necessary
      m.addConstr(quicksum(relXvars) <= 1)
        
  # if a flow is non-zero, the xvars should reflect that:
  for k, v in flowVars.items():
    i = k[0]
    j = k[1]
    if j == t:
      continue
    xvars = []
    for k in nodes:
      if (i,j,k) in xVars:
        xvars.append(xVars[i,j,k])
    m.addConstr(quicksum(xvars) >= v)

  # finally, add regular flow conservation constraints
  for i in inodes:
    flowin = []
    flowout = []
    qval = qvals[i-1]
    for j in nodes:
      if (i,j) in flowVars:
        # this is a flow out of i
        flowout.append(flowVars[i,j])
      if (j,i) in flowVars:
        # this is a flow into i
        flowin.append(flowVars[j,i])
    m.addConstr(quicksum(flowout) - qval*quicksum(flowin) == 0)

  # Set objective
  # gather edges with t as one of the nodes:
  flowTot = [flowVars[x,t] for x, y in edges if t == y]
  obj = quicksum(flowTot)
  m.setObjective(obj, GRB.MAXIMIZE)

  m.optimize()

  return obj.getValue()
