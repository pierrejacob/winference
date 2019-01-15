Reproduce results for the multivariate model in Section 2.3 of the WABC article.
===

mvnormal_generate_data.R: generate and save a data set from the model

mvnormal_wsmc_wasserstein.R: loads a data set of size 100 and runs SMC to approximate the WABC
posterior using the exact Wasserstein distance with a budget of 10^6 model simulations.

mvnormal_wsmc_summary.R: loads the data set and runs SMC to approximate the ABC
posterior based on the sample mean with a budget of 10^6 model simulations.

mvnormal_wsmc_euclidean.R: loads the data set and runs SMC to approximate the ABC
posterior based on the Euclidean distance with a budget of 10^6 model simulations.

mvnormal_rejection_summary.R: loads the data set and runs a rejection sampler to approximate the ABC
posterior based on the sample mean with a budget of 10^6 model simulations.

mvnormal_rejection_wasserstein.R: loads the data set and runs a rejection sampler to approximate the WABC
posterior based on the exact Wasserstein distance with a budget of 10^6 model simulations.

mvnormal_plots.R: loads the data and the output from the scripts above to plot the marginal ABC
posteriors as well as their W1 distances to the posterior. Corresponds to Fig 1 of the paper.

mvnormal_timings.R: estimates the average time it takes to compute the different distances
and simulate data
