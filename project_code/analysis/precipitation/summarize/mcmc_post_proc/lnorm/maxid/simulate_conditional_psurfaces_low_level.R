library(stablemix)

#   -----------------------------------------------------------------------
# Description: Simulate posterior predictive surfaces for the log-Gaussian
#			   process basis, theta > 0 model. Condition on low-level parameters
#             (e.g.
# 			   not basis function, but yes basis function hyper prior parameters)
# Output: A dataset called cond_sim_low_lnorm_maxid_csims.Rdata. This is used
#		  by the /summarize/chi_u/ scripts to make Figure 5.
#   -----------------------------------------------------------------------

# Parameters
stub <- "lnorm_maxid"  						# Model name

# Load data ---------------------------------------------------------------
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(data.dir, "predcoord.Rdata"))
load(file.path(samples.dir, paste0(stub, ".Rdata")))

# Make conditional draws
csims <-
  lnorm_cond_sim_low_all(out, pred_coords, thin_int = 1)  # Conditional simulation
dir.create(
  file.path(samples.dir, "cond_sim", paste0(stub, "initial_values")),
  showWarnings = FALSE,
  recursive = TRUE
)

# Save conditional simulations --------------------------------------------
save(csims, file = file.path(samples.dir, "cond_sim_low", paste0(stub, "_csims.Rdata")))
