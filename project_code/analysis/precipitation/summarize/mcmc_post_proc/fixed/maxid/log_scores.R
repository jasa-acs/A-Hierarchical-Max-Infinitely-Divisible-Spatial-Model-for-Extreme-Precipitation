library(stablemix)


#   -----------------------------------------------------------------------
# Description: Calculate the log-scores for the Gaussian density basis
#              theta > 0 model. This should be run after all scripts in
#              /precipitation/mcmc/ and precipitation/mcmc_post_proc/
#              that generate posterior predictive samples
# Output: A dataset containing the log-scores for the Gaussian density basis
#         theta >0 model. These log-scores are shown in Table 1 of the paper.
#         The dataset is called fixed_maxid_lscores.Rdata
#   -----------------------------------------------------------------------

# Parameters --------------------------------------------------------------
stub <- "fixed_maxid"             # Model name 

# Load data ---------------------------------------------------------------
load(file.path(data.dir, "mprecip.Rdata"))

# Sample latent path at holdout locations ---------------------------------
load(file.path(samples.dir, paste0(stub, ".Rdata")))
csims <-
  fixed_cond_sim_all(
    out,                                  # mcmc posterior predictive samples
    obs_coords,                           # observation coordinates
    holdout_coords,                       # holdout coordinates
    knot_coords = knot_coords,            # knot coordinates
    thin_int = 1,                  # thinning interval
    sim_y = FALSE                         # should response be simulated
  )
dir.create(file.path(samples.dir, "holdout"),
           showWarnings = FALSE,
           recursive = TRUE)
save(csims, file = file.path(samples.dir, "holdout", paste0(stub, "_holdout_latent.Rdata")))

# Calculate log-scores ----------------------------------------------------
# Extract latent variables
load(file.path(samples.dir, paste0(stub, ".Rdata")))
load(file.path(samples.dir, "holdout", 
  paste0(stub, "_holdout_latent.Rdata")))
nmcmc <- nrow(out$smp_par)
thin_id <- which((1:nmcmc))
lB_smps <- out$smp_lB[thin_id,]
theta_smps <- out$smp_par[thin_id, "theta"]
alpha_smps <- out$smp_par[thin_id, "alpha"]

# Calculate log-scores
lscores <-
  calc_log_scores(yhold,
                  alpha_smps,
                  theta_smps,
                  lB_smps,
                  csims[["lK"]],
                  csims[["mu"]],
                  csims[["sigma"]],
                  csims[["xi"]])

# Save result -------------------------------------------------------------
save(lscores, file = file.path(samples.dir, "holdout", paste0(stub, "_lscores.Rdata")))
