# Structure Regularized Learning of Higher Order Dynamics

---

A central problem in nonlinear system identification of high dimensional systems is the super-exponential number of nonlinear terms that may exist in the governing, dynamic equations.
Orthognal to dynamic system idnetification, many experimental techniques produces static -- rather than dynamic -- information about the structure of a system.
Here, we propose a method for nonlinear system identification where the structural information is utilized to prune the set of nonlinear interactions into a tractible problem that may be solved with classical approaches that work well on low dimensional systems.

## Applications
Cell dynamics -- idealized as occuring on Waddington's epigenetic landscape -- are high-dimenaional, nonlinear, Hamiltonian dynamical systems.
Identifying the global governing dynamics of cellular identity remains the Holy Grail of cell reprogramming.
Reprogramming algorithm have taken two approaches in recent years:
1. Statistical Framework: Deep learning models are trained on a collection of data about the source and target cell types and predict control inputs to steer cells toward the target.
2. Dynamic Framework: Cell dynamics are approximated as a system identification task (DMD, DGC, etc.) and control inputs are selected correspondingly.
The dynamic framework is advantageous because it provides insight to the underlying system, in addition to a predicted output; however, deep models are advantageous due to their ability to be trained on a much larger set of experimental data.
In addition to time series gene expression, we can utilize this framework the improve the type of dynamic models we learn for gene expression with the aim of improving our performance in cell reprogramming.
