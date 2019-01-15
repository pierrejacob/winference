Reproduce results for the AR(1) model in Section 4.2 of the article.
===

ar1_generate_data.R: generate a time series from an AR(1) model

ar1_wsmc_marginal.R: loads data, take first 1,000 observations,
and compute WABC posterior based on Wasserstein distance between
the marginal distributions of synthetic and observed data sets.

ar1_wsmc_delay1.R: loads data, take first 1,000 observations,
and compute WABC posterior based on Wasserstein distance between
the delay reconstructions of synthetic and observed data sets, with a lag of one.

ar1_plots.R: loads results from above scripts and creates the two plots
of Figure 3 (a) and (b).
