---
title: "Moment based estimation of stochastic Kronecker graph parameters"
layout: project
---

Moment based estimation of stochastic Kronecker graph parameters
=============

### David F. Gleich, Purdue University
### Art B. Owen, Stanford University

_These codes are research prototypes and may not work for you. No promises. But do email if you run into problems._

Download
--------

* [kgmoments.tar.gz](kgmoments.tar.gz) (updated 2012-01-16)
{: .nobullets}


Prereqs
-------

* [MatlabBGL](https://github.com/dgleich/matlab-bgl)
* [mcode](https://github.com/dgleich/mcode)
* A working C/C++ compiler
* A working Matlab mex compiler
* _Optional_ SNAP Codes (for liklihood calculations)

Instructions
------------

Will be prepared upon request.
 
Overview
--------

The package is organized by directory

`matlab`  
: All of the main matlab codes

`data`
: precomputed data for the experiments

`experiments`
: implementations of the experiments in the paper

`initial`
: initial versions of many of the codes

`skcoin`
: a python version of a coin-flipping random graph Kronecker graph generator

`snap`
: snap code to compute the log-likelihood of a kronecker graph

`web`
: this information and all the figures

Figures
-----------
    
|Experiment|Description|Figure|
|:------------------|:------------------------------------|:------------------|
|`experiments/fitting/objective_fits.m` | compute data for table on objective functions |  |
|`experiments/fitting/objectives_table.m` | output table | Tab. 2 |
|`experiments/fitting/kronecker_fits.m` | compute data for table on kronecker parameters |  |
|`experiments/fitting/kronecker_fits_table.m` | output table | Tab. 3 |
|`experiments/fitting/kronecker_partial_fits.m` | compute data for table on fitted kronecker parameters without certain features |  |
|`experiments/fitting/partial_fits_table.m` | output table | Tab. 4 |
|`experiments/identifiability/kron_identify_conflip.m` | compute data on variance of kronecker fits to kronecker model |  |
|`experiments/identifiability/parameter_variance.m` | output variance figures on a,b,c | Fig. 2  |
|`experiments/identifiability/feature_variance.m` | output variance figures | Fig. 3 |
|`experiments/identifiability/feature_misfit.m` | output misfit figures | Fig. 4 |

