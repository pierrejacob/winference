Reproduce results for the cosine model in Section 4.1 of the article.
===

cosine_generate_data.R: generate a time series from a cosine trend + noise model

cosine_mcmc.R: loads data, take first 100 observations,
and runs a Metropolis--Hastings algorithm to approximate the posterior distribution.

cosine_wsmc_euclidean.R: loads data, take first 100 observations,
and compute ABC posterior based on Euclidean distances between
synthetic and observed data sets.

cosine_wsmc_curvematching_wasserstein.R: loads data, take first 100 observations,
and compute WABC posterior based on curve matching, with lambda = 1,
and Wasserstein distance computed exactly using the "transport" package
(which needs to be installed, i.e. install.package("transport"))

cosine_plots.R: loads results from above scripts and creates the four plots
of Figure 2 (a,b,c,d).
