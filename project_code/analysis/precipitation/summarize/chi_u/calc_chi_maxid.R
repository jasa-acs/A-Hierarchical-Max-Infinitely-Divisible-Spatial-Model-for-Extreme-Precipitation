library(stablemix)
library(dplyr)
library(foreach)
library(doMC)
library(parallel)
library(extRemes)
source("./functions.R")
#   -----------------------------------------------------------------------
# Description: Calculate max-id model estimate of chi for annual maximum 
#              precipitation. This should be run after all scripts in
#              /precipitation/mcmc/ and precipitation/mcmc_post_proc/
#              that generate posterior predictive samples
# Output: A dataset in /precipitation/data/ called lnorm_maxid that 
#         contains estimates of chi_u for the theta > 0, log-Gaussian
#         process basis model for varying spatial lags.
#   -----------------------------------------------------------------------

# Load precip data ---------------------------------------------------------------
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(data.dir, "predcoord.Rdata"))
stub <- "lnorm_maxid"

# Define bins -------------------------------------------------------------
maxrng <- 320                                        # Maximum distance to consider (miles)                               
bin_centers <- seq(20, maxrng, by = 10)              # Bin centers to calculate chi for                                                         
rad <- (diff(bin_centers) / 2)[1]                    # half-width of bins                                                   
bins <-  c(bin_centers[-length(bin_centers)] - rad,  #                                                                      
            bin_centers[length(bin_centers)] + rad)  #                                                                     
bins <- bins * 1609.34                               # Convert from miles to meters                                         
bins[1] <- 1                                         # lower bound is some small number greate                               
ps <- seq(0.01, 0.98, by = 0.02)                     # upper quantiles to calculate chi for                                                   
ps <- c(ps, 0.25, 0.5, 0.75, 0.9, 0.98)                                                                       
ps <- sort(unique(ps))                                                                                                                  

# Calc Chi ----------------------------------------------------------------
load(file.path(samples.dir, "cond_sim_low", 
               paste0(stub, "_csims.Rdata")))
ysim <- csims[["y"]]                                # Conditional samples of response
pred_coord_list <- list()       
pred_id_list <- list()
nmc <- dim(ysim)[1]                                
for (k in 1:10) {
  pred_coord_list[[k]] <- list()
  pred_id_list[[k]] <- list()
  for (i in 1:nmc) {
    pred_id_list[[k]][[i]] <-
      sample(1:nrow(pred_coords), 100, replace = F)
    pred_coord_list[[k]][[i]] <-
      pred_coords[pred_id_list[[k]][[i]],]
  }
}
n.core = parallel::detectCores()
doMC::registerDoMC(cores = n.core)
moddfs <- foreach(i = 1:dim(ysim)[1]) %dopar% {
  chidf_list <- list()
  for (k in 1:10) {
    chidf_one <- calc_chi_u_spatial(ysim[i, , pred_id_list[[k]][[i]]],
                                    pred_coord_list[[k]][[i]],
                                    bins, ps, bin_centers = bin_centers *
                                      1609.34)
    chidf_list[[k]] <- data.frame(chidf_one, mcit = i, wi = k)
  }
  chidf <- do.call(rbind.data.frame, chidf_list)
  chidf
}
chimoddf <- do.call(rbind.data.frame, moddfs)
dir.create(file.path(samples.dir,  "chi_u"),
           showWarnings = FALSE,
           recursive = TRUE)

# Save output -------------------------------------------------------------
save(chimoddf,
     pred_id_list,
     pred_coord_list,
     file = file.path(samples.dir, "chi_u",
                      paste0(stub, "_chi_low.Rdata")))
