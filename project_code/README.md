# Companion code to 'A Hierarchical Max-infinitely Divisible Process for Extreme Precipitation'
To install the R package that can be used to simulate data and fit the model described in the paper, run

```R
R CMD INSTALL stablemix_0.1.0.tar.gz
```

For detailed information on which scripts to run to produce plots and tables from the paper, see the README file in the /analysis folder

# R Package
An R package with code for simulating from, fitting, and analyzing the models described in the paper. This code can be found in the `/stablemix` directory.

#Setup/Install
Before performing any subsequent analyses, install 'stablemix' R package 
R CMD install stablemix_0.1.0.tar.gz
MCMC Analysis
Within the /analysis/mcmc/ directory, the subdirectories are further divided (and labeled) by max-id (maxid), max-stable (maxstable), Gaussian density basis (fixed), and log-Gaussian process basis (lnorm). 

For each of these four models, the mcmc sampler is initialized with the initialize_mcmc.R script. The sampler is then run using mcmc.R. The analysis in the paper took 2 months to complete, and was run on a cluster in sequences of chunks of 48 hour wall times. The resulting samples were combined. To run for shorter periods reduce the number of mcmc samples in the mcmc.R scripts.

# Posterior Predictive Sampling (<24 hours on computing cluster)
First run MCMC Analysis

To make posterior predictive samples, for each respective model, run simulate_conditional_psurfaces.R and simulate_conditional_psurfaces_low_level.R in the /precipitation/summarize/<basis_type/ model_name> directories, where basis_type is one of lnorm or fixed, and model_name is one of maxid or maxstable.

# Figure 1 (<5 min)
Run /analysis/model_properties/plot_example_surface.R
# Figure 2 (<10 min)
Run 
1.	/analysis/model_properties/chi/calc_chi.R
2.	/analysis/model_properties/chi/plot_chi.R
# Figure 3 (10 hour walltime on a cluster)
Run
1.	/analysis/model_properties/chi/calc_gauss_density_spatial_chi.R
2.	/analysis/model_properties/chi/calc_gauss_proc_spatial_chir.R
3.	/analysis/model_properties/chi/plot_spatial_chi.R
# Figure 4 (<5 min)
Run
1.	/analysis/precipitation/summarize/plot_results/plot_gauge_knot_locations_together.R

# Figure 5 (<5 min)
First run MCMC Analysis and Posterior Predictive Sampling sections
Run
1.	/analysis/precipitation/summarize/chi_u/plot_chi.R

# Figure 6 (<10 min)
First run MCMC Analysis and Posterior Predictive Sampling sections
Run
1.	/analysis/precipitation/summarize/plot_results/plot_lnmid_qq_transform_gumbel_holdout.R

# Figure 7 (<5 min)
First run MCMC Analysis and Posterior Predictive Sampling sections
Run
1. /analysis/precipitation/summarize/plot_results/plot_lnmid_psurface_paper.R

# Figure 8 (<5 min)
First run MCMC Analysis and Posterior Predictive Sampling sections
Run
1.	/analysis/precipitation/summarize/plot_results/plot_lnmid_postpred_mean_draw.R

# Figure 9 (<10 min)
Run
1.	/analysis/precipitation/summarize/plot_results/factor_analysis.R


# Table 1
First run the Posterior Predictive Sampling scripts then
Run
1.	/analysis/precipitation/<basis_type/model_name>/log_scores 

where basis_type is one of lnorm or fixed, and model_name is one of maxid or maxstable

