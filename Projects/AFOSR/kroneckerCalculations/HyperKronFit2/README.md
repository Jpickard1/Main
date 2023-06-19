# HyperKronFit 2

June 18, 2023

---

This directory will contain all required code to execute the generalized KronFit algorithm on tensors rather than only graphs.

## Usage

Main Functions:

    1. HyperKronFit.m
    2. kronGen.m
    3. kHyperKronGen.m

## Functions

    1. hedgeLL.m: evaluates the log likelihood of an edge existing in a hypergraph
    2. hedgeDLL.m: evaluated the derivative of the log likeihood of an edge existing in a hypergraph
    3. nohedgeLL.m: evaluates the log likelihood of an edge not existing in a hypergraph
    4. nohedgeDLL.m: evaluates the derivative of the log likelihood of an edge not existing in a hypergraph
    5. genPermutation.m: generates the an alignment permutation of the Kronecker expansion and the data
    7. premutationProbabilityRatioTest.m: evaluates the ratio of likelihoods of 2 alignmnet permutations
    8. sampleGradient.m: samples the gradient of the log likelihood function with respect to theta
    9. emptyLL.m: evaluates the log likelihood of an empty hypergraph
    10. emptyDLL.m: evaluates the derivative of the log likelihood of an empty hypergraph with respect to theta

### Helper functions
    1. getCounts.m: returns the contributions of each element of theta to a particular hyperedge
    2. kronIndices.m:

### Generative Functions

    1. kronGen.m
    2. kHyperKronGen.m
