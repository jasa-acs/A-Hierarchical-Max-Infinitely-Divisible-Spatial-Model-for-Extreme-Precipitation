library(stablemix)

#   -----------------------------------------------------------------------
# Description: Simulate posterior predictive surfaces for the Gaussian
#			   density basis, theta = 0 model 
# Output: A dataset containing posterior predictive samples at prediction 
# 		  locations. This is not used in the paper, but may be useful 
#         for visualizing predictions.The dataset is called fixed_maxstable_csims.Rdata
#   -----------------------------------------------------------------------

# Parameters

stub <- "fixed_maxstable"
# Load data ---------------------------------------------------------------
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(data.dir, "predcoord.Rdata"))
load(file.path(samples.dir, paste0(stub, ".Rdata")))

# Make conditional simulations --------------------------------------------
csims <-
  fixed_cond_sim_all(out,										# MCMC posterior samples		
                     obs_coords,								# Observation coordinates				
                     pred_coords,								# Prediction coordinates				
                     knot_coords = knot_coords,					# Basis knot coordiantes							
                     thin_int = 1)  							# Thinning interval						

# Save output -------------------------------------------------------------
dir.create(
  file.path(samples.dir, "cond_sim", paste0(stub, "initial_values")),
  showWarnings = FALSE,
  recursive = TRUE
)
save(csims, file = file.path(samples.dir, "cond_sim", paste0(stub, "_csims.Rdata")))
