### Structural Control of Polynomials

This directory contains code to determine the structurall controllability of polynomial systems with linear inputs through a hypergraph framework.
Homogeneous polynomial systems with linear inputs may be described in the form $\dot{\mathbf{x}}=\textscr{A}^{k-1}\mathbf{x}+\mathbf{B}\mathbf{u}$.
We define a hypergraph structure for the pair $(\textscr{A},\mathbf{B})$ upon which structural conditions are provided for the dynamical system to be controllable.
In particular, if the hypergraph of $(\textscr{A},\mathbf{B})$ contains no hyperedge dilations and no nonaccessible vertices, then the corresponding dynamics are controllable.

## Code Structure:

* readEdgeSet.m: reads a text file to a hypergraph struct
* HG2Star: constructs a star graph from a hypergraph
* HG2Clique: constructs a clique graph from a hypergraph
* strucCtrb: verifies if a polynomial is controllable


