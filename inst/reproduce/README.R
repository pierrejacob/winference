Reproduce results of the article
===

This folder contains files to reproduce the figures in the article:
"Inference in generative models using the Wasserstein distance"
by Espen Bernton, Pierre E. Jacob, Mathieu Gerber, Christian P. Robert

They are part of the 'winference' package:
https://github.com/pierrejacob/winference

The files are organized by name: e.g. all files relating to the queueing model
start with 'queue_'. They have to be executed in a certain order, starting with
'_generate_data', then '_wsmc', and finally '_plots'. For instance, for the
bivariate normal example, execute:

mvnormal_generate_data.R
mvnormal_wsmc_euclidean.R
mvnormal_wsmc_summary.R
mvnormal_wsmc_wasserstein.R
mvnormal_plots.R

The files will save the results in the current working directory, which you might want to change
with the 'setwd' function. The scripts might take hours to run, depending on the machine; if this is a problem,
make sure you take a look at what you run before you run it!
