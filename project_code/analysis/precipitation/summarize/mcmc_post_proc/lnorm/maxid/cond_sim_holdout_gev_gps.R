library(stablemix)

#   -----------------------------------------------------------------------
# Description: Simulate the GEV marginal parameter Gaussian
#              processes at holdout locations, conditionally on the posterior
#              mcmc samples. Model: log-Gaussian process basis, theta> 0
# Output: A dataset containing conditional samples of GEV GPs at holdout locations
#         given posterior samples of GEV GPs at observation locations. Dataset
#         is called lnorm_maxid_holdout_gev_par.Rdata. Used in QQ plots in 
#    /precipitation/summarize/plot_results/plot_lnmid_postpred_qq_transform_gumbel.R
#   -----------------------------------------------------------------------

# Parameters --------------------------------------------------------------
stub <- "lnorm_maxid"

# Load data ---------------------------------------------------------------
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(samples.dir, paste0(stub, ".Rdata")))

# Sample gev GPs at holdout locations --------------------------------------
csim_gev <-
  lnorm_cond_sim_all(out,
                     obs_coords,
                     holdout_coords,
                     thin_int = 1,
                     sim_y = FALSE)
dir.create(file.path(samples.dir, "holdout"),
           showWarnings = FALSE,
           recursive = TRUE)
save(csim_gev, file = file.path(
  samples.dir,
  "holdout",
  paste0(stub, "_holdout_gev_par.Rdata")
))
