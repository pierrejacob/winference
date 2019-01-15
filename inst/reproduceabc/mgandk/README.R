Reproduce results for the multivariate g-and-k model in Section 5.1.2 of the WABC article.
===
  
mgandk_generate_data.R: generate and save a data set from the model

mgandk_wsmc_wasserstein.R: loads a data set of size 500 and runs SMC to approximate the WABC
posterior using the exact Wasserstein distance with a budget of 2*10^6 model simulations.

mgandk_wsmc_swap.R: loads the data set and runs SMC to approximate the WABC
posterior using the swapping distance with a budget of 2*10^6 model simulations.

mgandk_wsmc_hilbert.R: loads the data set and runs SMC to approximate the WABC
posterior using the Hilbert distance with a budget of 2*10^6 model simulations.

mgandk_wsmc_mmd.R: loads the data set and runs SMC to approximate the WABC
posterior using an approximation of  the MMD distance with a budget of 2*10^6
model simulations.

mgandk_mcmc.R: loads the data the WABC using the Hilbert approximation to tune
the MCMC proposal used to approximate the posterior distribution. Uses numerical
approximations of the likelihood. Performs 150,000 iterations.

mgandk_plots_ncomputed.R: loads the data and output from the scripts above.
Plot corresponds to Fig 6 a). number of sims vs number of smc steps

mgandk_plots_threshold.R: loads the data and output from the scripts above.
Plot corresponds to Fig 6 b). threshold vs number of sims.

mgandk_plots_w_to_posterior.R: loads the data and output from the scripts above.
Approximates the W1 distance between the posterior and the different WABC approximations,
plots W1 distance to posterior as function of model simulations. Corresponds to Fig 6 c).

mgandk_timings.R: estimate the average time it takes to compute the different distances

mgandk_plots.R: produces plots of WABC marginals.