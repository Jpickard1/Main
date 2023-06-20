# 3C with Kronecker Graphs

June 20, 2023

---

Hip-Hop Experiment:

Data: 200 sets of simulated coordinates with genes at HIGH, ON, OFF each (i.e. 600 total).

Process:
    1. for all coordinates construct a 3-way hypergraph representation (driverHipHop2HG)
    2. for each hypergraph representation learn the distribution theta
    3. for each theta compute:
        largest eigenvalues
        largest eigenvector
        corresponding laplacian
        corresponding fiedler vector

