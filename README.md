# TT-SIR

Matlab codes implementing the algorithms and numerical examples from [Dolgov, Savostyanov, [Tensor product approach to modelling epidemics on networks](https://doi.org/10.1016/j.amc.2023.128290), _Applied Mathematics and Computation_ 460, 128290, 2024] and [Dolgov, Savostyanov, [Tensor product algorithms for inference of contact network from epidemiological data](https://arxiv.org/abs/2401.15031), _arXiv_ 2401.15031, 2024]

## Main files

The main script files are

### Chemical Master Equation (CME) and rare event simulation

- `test_cme_sir.m` for setting up and running low-dimensional CME experiments, for example with the chain and Austria networks
- `test_cme_smallworld.m` for setting up and running Small World CME experiments

### Network inference

- `test_infer_chain.m` Set up the inference problem on the linear chain, (optionally) simulate new data, and run MCMC algorithm for sampling and maximizing the likelihood
- `test_infer_austria.m` Set up and run the inference of the Austria network
- `test_infer_smallworld.m` Set up and run the inference of one link in the Small World network

 ```markdown
 The full inference experiments above take A LOT of time, e.g. 3 days for the 15-dimensional Small World network. For a quick reproduction of the plots in the paper, the following scripts parse precomputed data.
 ```

- `plot_infer_network_stats.m` Plot mean $\pm$ standard deviation of the network configuration index error for chain or Austria data saved in `inferdata/`
- `plot_infer_loglike_stats.m` Plot mean $\pm$ standard deviation of the log-likelihood for chain or Austria data saved in `inferdata/`
- `plot_infer_smallworld.m` Plot log-likelihoods for Small World data saved in `inferdata/`

Each script will check for [TT-Toolbox](https://github.com/oseledets/TT-Toolbox) and [tAMen](https://github.com/dolgov/tamen) packages in Matlab paths (downloading them if necessary), ask interactively for model and approximation parameters, and run the corresponding methods. Currently available are TT (Tensor-Train) and SSA (Gillespie's Stochastic Simulation Algorithm) methods. Each method is run when it's corresponding approximation parameter (TT truncation `tol`erance for TT, Number of samples for SSA) is nonzero.

TT methods use `tamen` with adaptive discretization in time. This requires compressing several snapshots of the probability distribution function in the same TT decomposition, as explained in Section 3.3 of the paper. However, this may increase the error in the output statistics for very high dimensions, as explained in Section 5. To alleviate this issue, `test_cme_smallworld.m` uses an implicit Euler discretization with a reasonably fine time step (0.01), which requires storing only one snapshot of the solution at a time. This reduces the TT ranks and the TT approximation error. The linear system of the implicit Euler scheme is solved via `amen_solve` from the tAMEn package.

The forward CME solver for the inference problem can be chosen in the `cme_likelihood_si_bath.m` file:

- The `method` variable selects TT (tamen) or SSA method. _Note: SSA may produce zero likelihoods with insufficient number of trajectories._
- The `nt` variable sets the number of Chebyshev time grid points in tamen.
- The `Nsim` variable sets the number of SSA trajectories.

### Outputs

The `test_` scripts can compute and plot the following statistics of interest:

- mean and standard deviation of the total number of infected individuals. Note that both TT and SSA methods compute initially the mean and the second moment of the number of infected individuals. However, in the end of the `test_cme_sir.m` and `test_cme_smallworld.m` scripts the second moment is replaced by the standard deviation.
- Exceedance probability. For convenience of reproducing the experiments in the paper, both methods compute P(I > Icritical-1) and P(I > Icritical+1).
- Mean total number of susceptible individuals.

Inference scripts (`test_infer_`) save the entire current workspace to a `.mat` file (which is named with certain parameters such as the dimension and network type) in each MCMC step where the candidate maximum likelihood increased.

## Method files

Each file starts with a description of its purpose, inputs and outputs. Type e.g. `help tt_cme_sir` or open the file in the editor.

- `tt_cme_sir.m` Constructs and solves the SIR CME using tamen.
- `ssa.m` Stochastic Simulation Algorithm sampling SIR trajectories.
- `cme_likelihood_si_bath.m` Forward CME solvers for computing the log-likelihood of the given data in the inference problem.
- `infer_mcmc.m` MCMC methods for the inference problem (simple and no-Replacement)
- `scoring_initial_net.m` Score-based initial guess for network inference
- `simulate_load_infer_data.m` Simulate new observation data for the inference problem, or load existing data from a file.

## Auxiliary files

- `check_tt.m` Check/download/add-to-path for TT-Toolbox and tAMEn.
- `parse_parameter.m` Input a parameter interactively from a keyboard, offering a default value.
- `expectations_ssa.m` A function to average quantities of interest over pre-sampled SSA trajectories.
- `WattsStrogatz.m` Creates the Watts-Strogatz Small World network with random rewiring of some links.
- `WattsStrogatzDet.m` Creates the Watts-Strogatz Small World network with deterministic rewiring of the specified link.
- `adj_to_ind.m` Converts the network adjacency matrix into the network configuration index ($d(d-1)/2$ entries which take values 2 when the link is present, or 1 when the link is absent).
- `ind_to_adj.m` Converts the network configuration index into the adjacency matrix.

## Data

The folder `inferdata/` contain precomputed data for the inference experiments. The file names have the following format:

- `Xobs-Net*-d*-irun*.mat` Observed state trajectory data. Each wildcard (`*`) corresponds to the value of the corresponding parameter:
  - `Net`: network type: `11` for the linear chain, `14` for Austria, `17` for Small World.
  - `d`: dimension, i.e. the number of network nodes.
  - `irun`: the index of the experiment (realisation of the dataset and initial network configuration).
- `Net*-d*-irun*-methods*.mat` Matlab workspace from inference experiments. The `methods` wildcard can take values `2` for the simple MCMC method, or `4` for the MCMC no-Replacement method. Some of the useful variables are:

  - `Lmax_mcmc`, `Lmax_mcmc2` the history of candidate maximal log-likelihoods for MCMC and MCMC no-Replacement, respectively
  - `eval_mcmc`, `eval_mcmc2` the cumulative history of the number of forward CME solves where the log-likelihood has increased.
  - `ttimes_mcmc`, `ttimes_mcmc2` the cumulative CPU time where the log-likelihood has increased.
  - `W_mcmc`, `W_mcmc2` final maximum-likelihood adjacency matrices.
  - `W_ex` ground truth adjacency matrix.

```markdown
These data were produced with the default model and approximation parameters (e.g. T, the observation time) set up in the test_ scripts. If you change these default values, you may need to rerun the entire inference experiment with new simulated data.
```
