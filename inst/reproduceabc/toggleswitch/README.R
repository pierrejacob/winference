Reproduce results for the toggle switch model in Section 5.2 of the WABC article.
===
  
toggle_switch_generate.R: generate and save a data set from the model

toggle_switch_wsmc.R: loads a data set of size 2000 and runs SMC to approximate the WABC
posterior using the exact Wasserstein distance with a budget of 10^6 model simulations.

toggle_switch_summary.R: construct the summary statistic from Bonassi et al.

toggle_switch_load_summary.R: after toggle_switch_summary.R has been run, use this
script to load the summary statistic.

toggle_switch_summary_smc.R: loads the data set and the summary statistic and runs SMC to
approximate the ABC posterior based on the summary statistic with a budget of 10^6 model
simulations.

toggle_switch_plots.R: loads the data and the output from the SMC runs and plots the marginal 
distributions of the resulting ABC posteriors. Also plots a histogram of the data. Corresponds
to Fig 7 of the article.