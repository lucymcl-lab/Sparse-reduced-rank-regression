# Sparse-reduced-rank-regression
Simulation study in R comparing OLS, RRR, and Sparse Reduced Rank Regression (SRRR) in multivariate linear regression settings. Simulations evaluated prediction error, estimation accuracy, and variable selection under varying sparsity and correlation structures.

## Repository Contents

- `r_sim_stuff_F.R` - R code used to generate the simulation study
- `Sparse_Reduced_Rank_Regression_Simulation_Study.pdf` - LaTeX report summarising the simulation design, evaluation criteria and results.

## Methods Used
- Statistical computing in R
- Simulation study
- Multivariate linear regression
- Ordinary Least Squares
- Reduced Rank Regression
- Sparse Reduced Rank Regression

## Key Findings 
Across the scenarios considered, SRRR consistently achieved the lowest or near-lowest prediction and estimation errors, while also identifying the majority of relevant predictors. The benefits of both SRRR and RRR over OLS were clear in all scenarios, while the advantages of SRRR over RRR became more evident as the level of sparsity in the true coefficient matrix
increased.
