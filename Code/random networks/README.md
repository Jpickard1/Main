# Random Networks
Auth: Joshua Pickard jpic@umich.edu

Date: April 6, 2022

---

This directory contains code to generate different types of random networks. Currently we are able to generate 4 types of random networks. An example of each is seen below, and a more detailed explination of each type and the algorithms used to construct these networks is further down.

![Demo 4 random networks](https://github.com/Jpickard1/MissingData/blob/main/Code/random%20networks/random%20network%20demo.png?raw=true)

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


### Goals
The next steps for this directory include:
1. identify and build further types of random networks
2. fill in the details for each type in this README
3. learn about Ramanujin graphs
4. Write code to equate the size of each of these graphs by tuning their parameters so only the type, V, and maybe E need to be passed as parameters.
