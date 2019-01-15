Reproduce results for the Levy-driven stochastic volatility model in Section 5.4 of the article.
===

levydriven_generate_data.R: generate a time series from the model.

levydriven_wsmc_hilbert.R: loads data, takes first 10,000 observations,
approximates WABC with delay reconstruction (lag of one), using the Hilbert
distance.

levydriven_wsmc_with_summary.R: loads results from above scripts,
and defines new distance that uses both the Hilbert-delay reconstruction one,
and a distance on summary statistics, and then approximates the corresponding
WABC posterior.

levydriven_plots.R: loads results from above scripts and creates the three plots
of Figure 9 (a,b,c) and the three plots of Figure 10 (a,b,c).

===
Additionally, the folder contains files that could be useful but not necessary
to reproduce the figures:

levydriven_timings.R: runs particle filters and compute variance of likelihood
estimator as well as record the timings; this is to convince oneself that a
classic PMMH approach to this problem would be time consuming, due to the length
of the time series.

levydriven_mh.R: implements PMMH on a subset of the data.

levydriven_is_correction.R: implements an IS correction step to
go from the WABC posterior to the actual posterior.



