library(stablemix)

#   -----------------------------------------------------------------------
# Description: Simulate posterior predictive surfaces for the log-Gaussian
#			   process basis, theta = 0 model 
# Output: A dataset containing posterior predictive samples at prediction 
# 		  locations. This is not used in the paper, but may be useful 
#         for visualizing predictions.The dataset is called lnorm_maxstable.Rdata
#   -----------------------------------------------------------------------

# Parameters
stub <- "lnorm_maxstable" 	# model name

# Load data --------------------------------------------------------------
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(data.dir, "predcoord.Rdata"))
load(file.path(samples.dir, paste0(stub, ".Rdata")))

# Conditional simulations ------------------------------------------------
csims <-
  lnorm_cond_sim_all(out, obs_coords, pred_coords, thin_int = 1)  # Conditional simulation
dir.create(
  file.path(samples.dir, "cond_sim", paste0(stub, "initial_values")),
  showWarnings = FALSE,
  recursive = TRUE
)

# Save conditional simulations --------------------------------------------
save(csims, file = file.path(samples.dir, "cond_sim", paste0(stub, "_csims.Rdata")))
