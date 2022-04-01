# Edge Removal

This subdirectory contains code to remove edges from a network in a semi-intelligent
way. These methods come from the paper "Finding missing edges and communities
in incomplete networks" by Yan and Gregory. They posit there are 4 reasons
edges are missing in data:

1. Random: Edges are randomly missing uniformly throughout the data
2. Snowball Effect: Due to a networks size, only a sample of data is collected
typically through a BFS or random walk style algorith. As a result, edges are
missing at the frontier of the sampled data.
3. Cold Ends: Vertices with low degrees are more likely to have missing edges
4. Right Censoring: Edges are missing because the degree of a vertex is bounded

This directory contains functions to remove edges according to each of the
above schemes. All functions assume that a complete network is passed to them
and they return 2 objects: a network with the known edges accrding to one of
the above schemes and a network with all edges that are unknown.

WARNING: The code in this directory is currently set for unweighted, undirected
networks only.

---
CREATED: Joshua Pickard, 4/1/2022
