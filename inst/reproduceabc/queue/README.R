Reproduce results for the queueing model in Section 5.3 of the article.
===

queue_wsmc_marginal_intermediate.R: load data set (RData file in the folder),
and approximates WABC posterior based on Wasserstein distance between marginal
distributions of synthetic and observed data sets.

queue_pmmh_intermediate.R: loads results from above script, uses them to tune the parameters
of a PMMH algorithm to target the exact posterior distribution, and runs PMMH.

queue_abctools.R: runs semi-automatic ABC using abctools package
(install.packages("abctools"))

queue_plots_compare.R: loads results from above scripts,
and creates the three plots in Figure 8 (a,b,c).
