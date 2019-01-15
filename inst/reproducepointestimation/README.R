These scripts produce the figure for the manuscript
"On parameter estimation with the Wasserstein distance"
by
Espen Bernton, Pierre E. Jacob, Mathieu Gerber, Christian P. Robert

The instructions below describe the different files.

Some of the scripts require having loaded packages such as ggplot2 for plotting,
doParallel, foreach, doRNG for parallel computation and some 'dplyr' for data.frame manipulations.

#######

Section 4.1 Quantile “g-and-k” distribution

gandk_functions.R: contains function definitions, in particular the inverse CDF of the g-and-k distribution in C++

gandk.cluster.script.R: computes the MEWE on B=1000 bootstrap versions of the data, and saves the output; we ran this 400 times in parallel on a cluster.
(This function will have to be adapted to different machines, and file paths changed appropriately.)

gandk_coverage.R: postprocesses the outputs of the above script to compute the coverage.

gandk_plots.R: computes MEWE for different values of data size n, and a number of independent runs,
in order to produce the scatter plots and rescaled plots of the section, as well as a histogram of a data set.
These plots are in Figure 1 of the manuscript.

gandk_correlated.R: same but for a data-generating process that generates dependent
observations. The resulting plots are in Figure 2 of the manuscript.


#######

Section 4.2 Sum of log-Normal random variables

same as for g-and-k, with files starting with "lognormal_"

The resulting plots are in Figure 3 of the manuscript.

#######

Section 4.3 Gamma data fitted with a Normal model

gamma_normal_functions.R: contains function definitions, in particular those that generate data from the model

gamma_normal_fixedmk_diffn.R: produces plots to show the behavior of the MEWE as the number of observations increase
(plots of Figure 4 in the manuscript)

gamma_normal_fixedn_diffmk.R: produces plots to show the effect of different m and k on the quality of the MEWE approximation of the MWE, for a fixed data set.
(plots of Figure 5 in the manuscript)

gamma_normal_bootstrap.R: runs 400 independent experiments in which the MEWE is computed on B=1000 bootstrap versions of the data, and saves the the corresponding bootstrap intervals. It also calculates the coverage rate of the intervals.

#######

Section 4.4 Cauchy data fitted with a Normal model

cauchydata_normalfit_fixedmk_diffn.R: computes MEWE repeatedly, on independent data sets of various sizes,
and produces scatter plots and rescaled plots.

cauchydata_normalfit_fixedmk_diffn_KL.R: same for a minimum KL estimator, using the
FNN package.

These two scripts produce the plots presented in Figure 6 of the manuscript.
