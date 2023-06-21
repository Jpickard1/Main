# Big Fit

June 21, 2023

---

In this experiment, we demonstrate our ability to fit Stochastic Kronecker Hypergraph parameters for a variety of hypergraphs that varry in terms of their size, domain, and other features.

## Data sets of interest:
    1. Enron Email (143V, 1459E)
    2. DAWN: Drug Abuse Warning Network (2558V, 143523E)
    3. Math Stack Exchange (176,445V, 595,778E)
    4. coauth-DBPL (2599087V, 2599087E)

## Experiment Process:
    1. Construct 3, 4, and 5 uniform hypergraphs for each data set (dawn this is in progress)
    2. Learn SKH distribution theta for each hypergraph
    3. Generate synthetic/null hypergraphs using HyperKron model
    4. Compare descriptive statistics of synthetic versus original hypergraphs
