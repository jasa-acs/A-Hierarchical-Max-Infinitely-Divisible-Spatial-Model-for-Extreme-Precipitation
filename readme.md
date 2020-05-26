# A Hierarchical Max-Infinitely Divisible Spatial Model for Extreme Precipitation

# Author Contributions Checklist Form

## Data

### Abstract

Annual maximum precipitation data at gauge stations in the northeastern United States and Canada.

### Availability 

The data, which consist of annual maximum daily precipitation accumulations observed at rain gauge stations throughout the northeastern US and Canada are publically available through the National Oceanic Atmospheric Administration (NOAA). 

### Description 

Permissions: the data are publically available
Licensing information: NA
Link to data: https://hdsc.nws.noaa.gov/hdsc/pfds/pfds_series.html
Data provenance, including identifier or link to original data if different than above
File format: text file
Metadata: https://hdsc.nws.noaa.gov/hdsc/pfds/pfds_series.html
Version information: none available


## Code

### Abstract

The code provided contains a package for simulating from and fitting, via a Metropolis-Hastings algorithm, the proposed hierarchical max-id model for extreme precipitation. The R package is called stablemix. In addition to the R package, additional scripts are provided that perform the analyses of annual maximum precipitation observed over the Northeaster US and Canada.


### Description

How delivered: R package
Licensing information: MIT License
Link to code/repository: Submitted with manuscript, but can put on github
Version information 

### Description



### Optional Information 

Hardware requirements: access to a computing cluster
Supporting software requirements:
	R version: 3.6.0
	R packages: 
	extRemes (2.0.10),
	RandomFields (3.3.6),
	RcppGSL (0.3.6),
	Rcpp (1.0.1),
	abind (1.4.5),
	dplyr (0.8.1),
	tidyr (0.8.3), 
    	testthat (2.2.1)
	stabledist (0.7.1),
	texmex (2.4.2),
	foreach (1.4.4),
	parallel (3.6.0),
	doMC (1.3.5),
	SALTSampler (1.1.0),
	ggplot2 (3.1.1),
	ggthemes (4.2.0),
	cowplot (0.9.4),
	fields (9.8.3),
	RColorBrewer (1.1.2),
	mvtnorm (1.0.10),
	ggmap (3.0.0),
	utils (3.6.0),
	gdata (2.18.0),
	gridExtra (2.3.0),
	rgdal (1.4.6),
	raster (3.0.7)
	Other software:
	GSL (2.5)

## Instructions for Use

### Reproducibility 

What is to be reproduced 
-	Figures 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 
-	Table 1
How to reproduce analyses (e.g., workflow information, makefile, wrapper scripts)
Setup/Install
Before performing any subsequent analyses, install ‘stablemix’ R package 
R CMD install stablemix_0.1.0.tar.gz
MCMC Analysis
Within the /analysis/mcmc/ directory, the subdirectories are further divided (and labeled) by max-id (maxid), max-stable (maxstable), Gaussian density basis (fixed), and log-Gaussian process basis (lnorm). 

For each of these four models, the mcmc sampler is initialized with the initialize_mcmc.R script. The sampler is then run using mcmc.R. The analysis in the paper took 2 months to complete, and was run on a cluster in sequences of chunks of 48 hour wall times. The resulting samples were combined. To run for shorter periods reduce the number of mcmc samples in the mcmc.R scripts.

Posterior Predictive Sampling (<24 hours on computing cluster)
First run MCMC Analysis

To make posterior predictive samples, for each respective model, run simulate_conditional_psurfaces.R and simulate_conditional_psurfaces_low_level.R in the /precipitation/summarize/<basis_type/ model_name> directories, where basis_type is one of lnorm or fixed, and model_name is one of maxid or maxstable.

Figure 1 (<5 min)
Run /analysis/model_properties/plot_example_surface.R
Figure 2 (<10 min)
Run 
1.	/analysis/model_properties/chi/calc_chi.R
2.	/analysis/model_properties/chi/plot_chi.R
Figure 3 (10 hour walltime on a cluster)
Run
1.	/analysis/model_properties/chi/calc_gauss_density_spatial_chi.R
2.	/analysis/model_properties/chi/calc_gauss_proc_spatial_chir.R
3.	/analysis/model_properties/chi/plot_spatial_chi.R
Figure 4 (<5 min)
Run
1.	/analysis/precipitation/summarize/plot_results/plot_gauge_knot_locations_together.R
Figure 5 (<5 min)
First run MCMC Analysis and Posterior Predictive Sampling sections
Run
1.	/analysis/precipitation/summarize/chi_u/plot_chi.R
Figure 6 (<10 min)
First run MCMC Analysis and Posterior Predictive Sampling sections
Run
1.	/analysis/precipitation/summarize/plot_results/plot_lnmid_qq_transform_gumbel_holdout.R
Figure 7 (<5 min)
First run MCMC Analysis and Posterior Predictive Sampling sections
Run
/analysis/precipitation/summarize/plot_results/plot_lnmid_psurface_paper.R
Figure 8 (<5 min)
First run MCMC Analysis and Posterior Predictive Sampling sections
Run
1.	/analysis/precipitation/summarize/plot_results/plot_lnmid_postpred_mean_draw.R
Figure 9 (<10 min)
Run
1.	/analysis/precipitation/summarize/plot_results/factor_analysis.R
Figure 10 (<2 hours) 
Run
1.	/analysis/precipitation/summarize/watershed/lnorm_merrimack_predict_quantiles.R

Table 1
First run the Posterior Predictive Sampling scripts then
Run
1.	/analysis/precipitation/<basis_type/ model_name>/log_scores 

where basis_type is one of lnorm or fixed, and model_name is one of maxid or maxstable


Detailed description of which scripts produce figures and tables in the paper

General descriptions of what each script produces and their connections to the paper are given below. It is assumed that each script is called from the directory in which it resides.

setup.R
	# This sets the paths to data and figure directories for loading and storing data and figures. This script should always be called before any other analysis to set paths. Throughout it is assumed that all scripts are run from the directory they reside in.

./model_properties:
	plot_example_surface.R
		# Description: Simulate a realization of the max-id and max-stable models
		#              with Gaussian density and log-Gaussian process basis functions
		#
		# Output: 
		#     - Figure in ./fig containing sample realizations of the four models
		#       described in the paper. See, e.g., Figure 1 of the paper.

./model_properties/chi:
	calc_chi.Rdata
		# Description: Generate Monte Carlo estimates of chi_u  for all four models
		#              discussed in the paper at two fixed spatial locations.
		#              The output is saved and used in plot chi.
		#
		# Output: Two datasets ./data/fixed_chi.Rdata and ./data/lnorm_chi.Rdata
		#         which are called by plot_chi.R to create Figure 2 of the paper.
	calc_gauss_density_spatial_chi.R
		# Description: Generate Monte Carlo estimates of chi_u as a function of distance
		#              for the gaussian density basis models discussed in the paper.
		#
		# Output: Produces a dataset of monte carlo samples of chi_u and stores it in
		#         ./data/fixed_chi_spatial.Rdata which is called by 
		#          plot_spatial_chi.R.
	calc_gauss_proc_spatial_chi.R
		# Description: Generate Monte Carlo estimates of chi_u as a function of distance
		#              for the log-gaussian process basis models discussed in the paper.
		#
		# Output: Produces a dataset of monte carlo samples of chi_u and stores it in
		#         ./data/lnorm_chi_spatial.Rdata which is called by 
		#          plot_spatial_chi.R.
	plot_chi.R
		# Description: plot the Monte Carlo estimates of chi_u generated in the script
		#              calc_chi.R as a function of u for fixed locations. calc_chi.R
		#              must be run prior to this script.
		#
		# Output: Produces Figure 2 of the paper, plots of of the functions chi_u
		#         for all four models described in the paper. Figure is stored in
		#         ./fig/chi_plot.pdf
	plot_spatial_chi.R
		# Description: plot the Monte Carlo estimates of chi_u generated in the script
		#              calc_spatial_chi.R as a function spatial distance for select
		#              quantiles. calc_gauss_density_spatial_chi.R and
		#              calc_gauss_proc_spatial_chi.R must be run before calling this
		#              script
		# 
		# Output: Produces Figure 3 of the paper. Plots of chi_u as a function of distance
		#         for the four models described in the paper. Figure is stored in 
		#         ./fig/spatial_chi_plots.pdf

./precipitation/data:
	map_polygon.Rdata
		# A polygon dataset for plotting maps
	mprecip.Rdata
		# data set containing spatial locations and annual precipitation maxima used in analysis
	ne_states.Rdata
		# definitions of north eastern US polygon dataset for plotting
	predcoord.Rdata
		# prediction spatial coordinate matrix
	subregion_groups.Rdata
		# a dataset specifying which subregion each gauge station in mprecip belongs to


./precipitation/mcmc:
	combine_mcmc_output.R
		# Description: If MCMC sampler is run in chunks (to deal with e.g. max wall times)
		#			   this script will combine the output chunks and store the combined data.
		#
		# Output: 
		#     - Combined datasets for each of the four models are stored in 
		#		lnorm_maxid.Rdata
		#		lnorm_maxstable.Rdata
		#		fixed_maxid.Rdata
		#		fixed_maxstable.Rdata

./precipitation/mcmc/fixed/maxid:
	initialize_mcmc.R
		# Description: This code initializes the parameters sampled in the mcmc
		#              algorithm for the Gaussian density basis, theta > 0 model
		#              described in the paper. Run this before mcmc.R
		# 
		# Output: A dataset called initial_values1.Rdata containing initial values 
		#         for the mcmc for the theta > 0 Gaussian density basis model.
	mcmc.R
		# Description: This code runs the mcmc sampler forthe Gaussian density basis,
		#              theta > 0 model described in the paper. This should be run 
		#              this should be run in chunks for long jobs.
		# Output: A dataset called fixed/maxid<chunk>.RData containing the posterior 
		#         samples for the Gaussian density theta > 0 model, where chunk is
		#         the number of times this script has been run sequentially if broken
		#         up into smaller chunks. Initial values for the next chunk are also
		#         stored in a dataset called initial_values<chunk>.Rdata

./precipitation/mcmc/fixed/maxstable:
	initialize_mcmc.R
		# Description: This code initializes the parameters sampled in the mcmc
		#              algorithm for the Gaussian density basis, theta = 0 model
		#              described in the paper. Run this before mcmc.R
		# Output: A dataset called initial_values1.Rdata containing initial values 
		#         for the mcmc for the theta = 0 Gaussian density basis model.
	mcmc.R
		# Description: This code runs the mcmc sampler forthe Gaussian density basis,
		#              theta = 0 model described in the paper. This should be run 
		#              this should be run in chunks for long jobs.
		# Output: A dataset called fixed/maxid<chunk>.RData containing the posterior 
		#         samples for the Gaussian density theta = 0 model,where chunk is
		#         the number of times this script has been run sequentially if broken
		#         up into smaller chunks. Initial values for the next chunk are also
		#         stored in a dataset called initial_values<chunk>.Rdata


./precipitation/mcmc/lnorm/maxid:
	initialize_mcmc.R
		# Description: This code initializes the parameters sampled in the mcmc
		#              algorithm for the Gaussian process basis, theta > 0 model
		#              described in the paper. Run this before mcmc.R
		# Output: A dataset called initial_values1.Rdata containing initial values 
		#         for the mcmc for the theta > 0 Gaussian process basis model.
	mcmc.R
		# Description: This code runs the mcmc sampler forthe Gaussian process,
		#              theta > 0 model described in the paper. This should be run 
		#              this should be run in chunks for long jobs.
		# Output: A dataset called fixed/maxid<chunk>.RData containing the posterior 
		#         samples for the Gaussian density theta = 0 model,where chunk is
		#         the number of times this script has been run sequentially if broken
		#         up into smaller chunks. Initial values for the next chunk are also
		#         stored in a dataset called initial_values<chunk>.Rdata

./precipitation/mcmc/lnorm/maxstable:
	initialize_mcmc.R
		# Description: This code initializes the parameters sampled in the mcmc
		#              algorithm for the Gaussian process basis, theta = 0 model
		#              described in the paper. Run this before mcmc.R
		# Output: A dataset called initial_values1.Rdata containing initial values 
		#         for the mcmc for the theta = 0 Gaussian process basis model.
	mcmc.R
		# Description: This code runs the mcmc sampler forthe Gaussian process,
		#              theta = 0 model described in the paper. This should be run 
		#              this should be run in chunks for long jobs.
		# Output: A dataset called fixed/maxid<chunk>.RData containing the posterior 
		#         samples for the Gaussian process theta = 0 model, where chunk is
		#         the number of times this script has been run sequentially if broken
		#         up into smaller chunks. Initial values for the next chunk are also
		#         stored in a dataset called initial_values<chunk>.Rdata

./precipitation/summarize/chi_u:
	calc_chi_emp_holdout.R
		# Description: Calculate an empirical estimate of chi for annual maximum 
		# 			   precipitation at the holdout gauge locations. Run this
		# 			   before running plot_chi.R
		# Output: Creates a dataset called emp_chi.Rdata in the precipitation/data
		#		  directory with emprical estimates of chi_u for different u for the
		# 		  annual max precip data. This dataset will be used by plot_chi.
	calc_chi_maxid.R
		# Description: Calculate max-id model estimate of chi for annual maximum 
		#              precipitation. This should be run after all scripts in
		#              /precipitation/mcmc/ and precipitation/mcmc_post_proc/
		#              that generate posterior predictive samples
		# Output: A dataset in /precipitation/data/ called lnorm_maxid_chi_low.Rdata that 
		#         contains estimates of chi_u for the theta > 0, log-Gaussian
		#         process basis model for varying spatial lags.
	functions.R
		# Functions for calculating chi_u used by calc_chi_emp_holdout.R and
		# calc_chi_maxid.R. Must be in the same directory as those scripts.
	plot_chi.R
		# Description: Plot the empirical and model estimates of chi_u as a function
		#              of distance.
		# Output: A figure of chi_u as a function of spatial distance for the log-
		#         Gaussian process basis theta >0 model. Figure 5 of the paper.
		#         The output figure is called chi_10mi.pdf

./precipitation/summarize/mcmc_post_proc/fixed/maxid:
	log_scores.R
		# Description: Calculate the log-scores for the Gaussian density basis
		#              theta > 0 model. This should be run after all scripts in
		#              /precipitation/mcmc/ and precipitation/mcmc_post_proc/
		#              that generate posterior predictive samples
		# Output: A dataset containing the log-scores for the Gaussian density basis
		#         theta >0 model. These log-scores are shown in Table 1 of the paper.
		#         The dataset is called fixed_maxid_lscores.Rdata
	simulate_conditional_psurfaces.R
		# Description: Simulate posterior predictive surfaces for the Gaussian
		#			   density basis, theta > 0 model. This should be run after all scripts in
		#              /precipitation/mcmc/
		# Output: A dataset containing posterior predictive samples at prediction 
		# 		  locations. This is not used in the paper, but may be useful 
		#         for visualizing predictions.The dataset is called fixed_maxid_csims.Rdata

./precipitation/summarize/mcmc_post_proc/fixed/maxstable:
	log_scores.R
		# Description: Calculate the log-scores for the Gaussian density basis
		#              theta = 0 model. This should be run after all scripts in
		#              /precipitation/mcmc/ and precipitation/mcmc_post_proc/
		#              that generate posterior predictive samples
		# Output: A dataset containing the log-scores for the Gaussian density basis
		#         theta = 0 model. These log-scores are shown in Table 1 of the paper.
		#         The dataset is called fixed_maxstable_lscores.Rdata
	simulate_conditional_psurfaces.R
		# Description: Simulate posterior predictive surfaces for the Gaussian
		#			   density basis, theta = 0 model 
		# Output: A dataset containing posterior predictive samples at prediction 
		# 		  locations. This is not used in the paper, but may be useful 
		#         for visualizing predictions.The dataset is called fixed_maxstable_csims.Rdata

./precipitation/summarize/mcmc_post_proc/lnorm/maxid:
	cond_sim_holdout.R
		# Description: Make posterior predictie draws from the model at holdout 
		#  			   locations. Model: log-Gaussian process basis, theta> 0
		# Output: A dataset called lnorm_maxid_holdout_yhat.Rdata containing 
		# 		  posterior predictive draws. Used in QQ plots in 
		#		 /precipitation/summarize/plot_results/plot_lnmid_postpred_qq_transform_gumbel.R
	cond_sim_holdout_gev_gps.R
		# Description: Simulate the GEV marginal parameter Gaussian
		#              processes at holdout locations, conditionally on the posterior
		#              mcmc samples. Model: log-Gaussian process basis, theta> 0
		# Output: A dataset containing conditional samples of GEV GPs at holdout locations
		#         given posterior samples of GEV GPs at observation locations. Dataset
		#         is called lnorm_maxid_holdout_gev_par.Rdata. Used in QQ plots in 
		#    /precipitation/summarize/plot_results/plot_lnmid_postpred_qq_transform_gumbel.R
	log_scores.R
		# Description: Calculate the log-scores for the log-Gaussian process basis
		#              theta > 0 model. This should be run after all scripts in
		#              /precipitation/mcmc/ and precipitation/mcmc_post_proc/
		#              that generate posterior predictive samples
		# Output: A dataset containing the log-scores for the log Gaussian process basis
		#         theta > 0 model. These log-scores are shown in Table 1 of the paper.
		#         The dataset is called fixed_maxid_lscores.Rdata
	simulate_conditional_psurfaces.R
		# Description: Simulate posterior predictive surfaces for the log-Gaussian
		#			   process basis, theta > 0 model 
		# Output: A dataset containing posterior predictive samples at prediction 
		# 		  locations. This is not used in the paper, but may be useful 
		#         for visualizing predictions.The dataset is called lnorm_maxid.Rdata
	simulate_conditional_psurfaces_low_level.R
		# Description: Simulate posterior predictive surfaces for the log-Gaussian
		#			   process basis, theta > 0 model. Condition on low-level parameters
		#             (e.g.
		# 			   not basis function, but yes basis function hyper prior parameters)
		# Output: A dataset called cond_sim_low_lnorm_maxid_csims.Rdata. This is used
		#		  by the /summarize/chi_u/ scripts to make Figure 5.

./precipitation/summarize/mcmc_post_proc/lnorm/maxstable:
	log_scores.R
		# Description: Calculate the log-scores for the log-Gaussian process basis
		#              theta = 0 model. This should be run after all scripts in
		#              /precipitation/mcmc/ and precipitation/mcmc_post_proc/
		#              that generate posterior predictive samples
		# Output: A dataset containing the log-scores for the log Gaussian process basis
		#         theta > 0 model. These log-scores are shown in Table 1 of the paper.
		#         The dataset is called fixed_maxid_lscores.Rdata
	simulate_conditional_psurfaces.R
		# Description: Simulate posterior predictive surfaces for the log-Gaussian
		#			   process basis, theta = 0 model 
		# Output: A dataset containing posterior predictive samples at prediction 
		# 		  locations. This is not used in the paper, but may be useful 
		#         for visualizing predictions.The dataset is called lnorm_maxstable.Rdata

./precipitation/summarize/plot_results:
	factor_analysis.R
		# Description: Plot model posterior spatial factor means, ranked by the 
		#              year-to-year variance of the basis scaling factors. All
		#              scripts in /mcmc/lnorm/maxid and /mcmc_post_proc/lnorm/maxid
		#              must be run before this script.
		# Output: Figure 9 of the paper. The top basis factor means are plotted in a
		#         figure called top_K.jpeg
	plot_gauge_knot_locations.R
		# Description: Make a plot of the spatial gauge locations and knot coordinates.
		#
		# Output: Figure 4 of the paper. The figure is stored in a file called gauge_knot_locations_together.pdf
	plot_lnmid_postpred_mean.R
		# Description: Plot the posterior predictive mean, sd, and single draw 
		#              surfaces for the log-Gaussian process basis theta > 0 model.
		#              Scripts in /mcmc/lnorm/maxid and /mcmc_post_proc/lnorm/maxid
		#              must be run before this script.
		# Output: Figure 8 of the paper. The figure is stored in ppred_mean.jpeg
	plot_lnmid_postpred_qq_transform_gumbel_holdout.R
		# Description: Calculate group-wise min, mean, and maxima for the annual 
		#              maxima transformed to the Gumbel scale for the holdout 
		#              locations. Compare model and empirical quantiles of each.
		#              Scripts in /mcmc/lnorm/maxid and /mcmc_post_proc/lnorm/maxid
		#              must be run before this script.
		# Output: Figure 6 of the paper. The figure is saved in lnmid_qqplots_gumbel_hold.jpeg
	plot_lnmid_psurface.R
		# Description: Plot posterior predictive estimate and uncertainty 
		#              of 99th percentile using the log-Gaussian process basis,
		#              theta > 0 model. 
		#              Scripts in /mcmc/lnorm/maxid and /mcmc_post_proc/lnorm/maxid
		#              must be run before this script.
		# Output: Figure 7 of the paper. The figure is stored in lnmidQ99.jpeg

## Replication (Optional)
How to use software in other settings (or links to such information, e.g., R package vignettes, demos or other examples)

## Notes
The code and data have been attached with the submission. 

