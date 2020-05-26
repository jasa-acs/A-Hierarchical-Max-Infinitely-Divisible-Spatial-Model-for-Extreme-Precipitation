library(stablemix)

#   -----------------------------------------------------------------------
# Description: Make posterior predictie draws from the model at holdout 
#  			   locations. Model: log-Gaussian process basis, theta> 0
# Output: A dataset called lnorm_maxid_holdout_yhat.Rdata containing 
# 		  posterior predictive draws. Used in QQ plots in 
#		 /precipitation/summarize/plot_results/plot_lnmid_postpred_qq_transform_gumbel.R
#   -----------------------------------------------------------------------


# Parameters --------------------------------------------------------------
stub <- "lnorm_maxid"

# Load data ---------------------------------------------------------------
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(samples.dir, paste0(stub, ".Rdata")))
# Make posterior predictive draws-------------------------------------------
csims <-
  lnorm_cond_sim_all(out,
                     obs_coords,
                     holdout_coords,
                     thin_int = 1,
                     sim_y = TRUE)
dir.create(file.path(samples.dir, "holdout"),
           showWarnings = FALSE,
           recursive = TRUE)
save(csims, file = file.path(samples.dir, "holdout", paste0(stub, "_holdout_yhat.Rdata")))
