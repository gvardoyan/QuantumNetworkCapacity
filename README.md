# QuantumNetworkCapacity

- computeNetCap.m computes the overall network capacity. SURFnet is provided as an example; snapshot capacity data for this network is stored in SURFNetExp/. 

- Code files are divided into two sets: those that use gurobi optimization and brute-force code that was used to validate gurobi results. Brute-force code has two versions: one for multiplexed (computeSnapshotCapMult.py), and one for non-multiplexed networks (computeSnapshotCap.py).

- Data for the five-node multiplexed network is provided in FiveNodeExpMultiplexing/.

- Example network capacity computation using SURF:
  (1) Create a local SURFNetExp/ directory for storing data files,
  (2) Run SURFNetExp.py

- To run gurobi optimization software, a gurobi install is required.
