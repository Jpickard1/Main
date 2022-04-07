# [Recovering Gener Interactions from Single-Cell Data Using Data Diffusion](https://www.sciencedirect.com/science/article/pii/S0092867418307244)

Auth: Joshua Pickard jpic@umich.edu

Created: April 7, 2022

---

This directory contains code to run the Markov Affinitiy based Graph Imputation of Cells (MAGIC) algorithm. The alorithm was developed by Smita Krishnaswamy and her lab. The majority of this code was downloaded from the Matlab directory at this git repository: https://github.com/KrishnaswamyLab/MAGIC, but I have written some code in `demo_magic.m` and modified the function `compute_optimal_t.m`.

## Selecting the Optimal *t*

The most interesting step of the MAGIC algorithm is multiplying the original data matrix *D* by the markov matrix *M* raised to some power *t*. We have been interested in determining the optimal value for *t,* and we think that it should be some function of the rank of the imputed data or the norm between the imputed data and the original data. That paper presents a different method for comparing the optimal value of *t* based on the following equation:

<div align="center">
<img src="https://render.githubusercontent.com/render/math?math=Error=R_{sq}(D_t, D_{t-1}) = 1 - SSE(D_t, D_{t-1})/SST(D_t, D_{t-1}),">
</div>

where SSE and SST stand for the sum of squared erros and totals respectively. Their function `compute_optimal_t.m` determines this error on line 30 using the MATLAB `procrustes.m` function, which I am not familiar with, but I will assume their code is correct. They say the optimal *t* occurs when this error is below a prespecified threshold value.

The figure below plots their error as a function of *t* as well as the rank of the imputed data and the norm between the imputed data and the original measurements. Note that the norm considers all points in the matrix, not only those points that are "known"; however the general trend is similar to if only "known" edges were considered. This plot was generated using a random matrix with 200 cells and 500 genes.

I find it interesting that the optimal *t* value according to the paper would occur after the norm is maximized and after the rank is minimized. It would seem better to select a good *t* value at least where the rank is greater than 1, especially if it won't increase the distance from the original data any further.

![Optimal t by MAGIC, Norm, and Rank](https://github.com/Jpickard1/MissingData/blob/main/Code/reproductions/Recovering%20Gener%20Interactions%20from%20Single-Cell%20Data%20Using%20Data%20Diffusion/Optimal%20t%20with%20MAGIC%20Rank%20and%20Norms.png?raw=true)
