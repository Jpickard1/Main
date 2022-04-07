# Random Networks
Auth: Joshua Pickard jpic@umich.edu

Date: April 6, 2022

---

This directory contains code to generate different types of random networks. There are a variety 

## Erdos-Renyi Network

An Erdos-Renyi graph is a graph with *V* vertices and *E* edges where the graph is chosen uniformly at random from the set of all graphs with *V* vertices and *E* edges. The vertices are ordered in the sense that graphs that can be obtained from one another via graph permutations are distinct.

## Erdos-Renyi-Gilbert Network

An Erdos-Renyi-Gilbert graph (ERG) is a graph with *V* vertices where every edge occurs with probability *p* independent of all other edges. The probability of an ERG graph with *V* vertices having *M* edges is 

<div align="center">
<img src="https://render.githubusercontent.com/render/math?math=p^M(1-p)^{\binom{n}{2}-M}.">
</div>

## Scale Free Network

A scale free network is a network where the degree distribution follows the Power Law. This means that the fraction of nodes with at least degree *k* is proportional to 
<img src="https://render.githubusercontent.com/render/math?math=k^{-\gamma}.">

### BA Algorithm
## Small World Network
