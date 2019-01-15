Reproduce results for the univariate g-and-k model in Section 5.1.1 of the WABC article.
===
  
gandk_generate_data.R: generate and save a data set from the model

gandk_wsmc.R: loads a data set of size 250 and runs SMC to approximate the WABC
posterior using the exact Wasserstein distance with a budget of 2.4*10^6 model simulations.
Can be continued for a long run (e.g. 10^8 simulations in as Fig 5 of the paper),
using wsmc_continue.

gandk_abctools.R: loads the abctools package, loads the data. Runs the semi-automatic
ABC procedure of Fearnhead and Prangle (2012) using 2.4*10^6 model simulations. Might 
require a lot of memory.

gandk_mcmc.R: load the data and run MCMC to approximate the posterior distribution,
using numerical approximations of the likelihood. Proposal adapted to target using a
previous run of the algorithm. Performs 75,000 iterations.

gandk_plots_compare.R: loads the data and the output of the scripts above. Plots 
the marginal distributions corresponding to Fig 4 of the paper. 

gandk_plots_convergence.R: loads the data and output from gandk_wsmc.R and gandk_mcmc.R.
Plots correspond WABC marginal posteriors in Fig 5 a)-d).

gandk_plots_ncomputed.R: loads the data and output from gandk_wsmc.R and gandk_mcmc.R.
Plot corresponds to Fig 5 e). number of sims vs number of smc steps

gandk_plots_threshold.R: loads the data and output from gandk_wsmc.R and gandk_mcmc.R.
Plot corresponds to Fig 5 f). threshold vs number of model sims.

gandk_plots_w_to_posterior.R: loads the data and output from gandk_wsmc.R and gandk_mcmc.R.
Approximates the W1 distance between the posterior and the WABC approximation, plots W1
distance as function of model simulations. Plot corresponds to Fig 5 g).