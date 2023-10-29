# TT-SIR

Matlab codes implementing the algorithms and numerical examples from [Dolgov, Savostyanov, [Tensor product approach to modelling epidemics on networks](https://doi.org/10.1016/j.amc.2023.128290), _Applied Mathematics and Computation_ 460, 128290, 2024].

## Main files

The two main script files are

 - `test_cme_sir.m` for setting up and running low-dimensional experiments, for example with the chain and Austria networks
 - `test_cme_smallworld.m` for setting up and running Small World experiments.

Each script will check for [TT-Toolbox](https://github.com/oseledets/TT-Toolbox) and [tAMen](https://github.com/dolgov/tamen) packages in Matlab paths (downloading them if necessary), ask interactively for model and approximation parameters, and run the corresponding methods. Currently available are TT and SSA methods. Each method is run when it's corresponding approximation parameter (TT truncation `tol`erance for TT, Number of samples for SSA) is nonzero.

`test_cme_sir.m` uses `tamen` with adaptive discretization in time. This requires compressing several snapshots of the probability distribution function in the same TT decomposition, as explained in Section 3.3 of the paper. However, this may increase the error in the output statistics for very high dimensions, as explained in Section 5. To alleviate this issue, `test_cme_smallworld.m` uses an implicit Euler discretization with a reasonably fine time step (0.01), which requires storing only one snapshot of the solution at a time. This reduces the TT ranks and the TT approximation error. The linear system of the implicit Euler scheme is solved via `amen_solve` from the tAMEn package.

### Outputs

Both scripts can compute and plot the following statistics of interest:
 - mean and standard deviation of the total number of infected individuals. Note that both TT and SSA methods compute initially the mean and the second moment of the number of infected individuals. However, in the end of the `test_cme_sir.m` and `test_cme_smallworld.m` scripts the second moment is replaced by the standard deviation.
 - Exceedance probability. For convenience of reproducing the experiments in the paper, both methods compute P(I > Icritical-1) and P(I > Icritical+1).
 - Mean total number of susceptible individuals.


## Method files

Each file starts with a description of its purpose, inputs and outputs. Type e.g. `help tt_cme_sir` or open the file in the editor.

 - `tt_cme_sir.m` Constructs and solves the SIR CME using tamen.
 - `ssa.m` Stochastic Simulation Algorithm sampling SIR trajectories.

## Auxiliary files

 - `check_tt.m` Check/download/add-to-path for TT-Toolbox and tAMEn
 - `expectations_ssa.m` A function to average quantities of interest over pre-sampled SSA trajectories
 
