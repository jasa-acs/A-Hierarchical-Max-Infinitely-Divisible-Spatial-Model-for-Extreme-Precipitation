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




## Replication (Optional)
How to use software in other settings (or links to such information, e.g., R package vignettes, demos or other examples)

## Notes
The code and data have been attached with the submission. 

