library(stablemix)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(RColorBrewer)
library(extRemes)
library(mvtnorm)
source("./functions.R")

#   -----------------------------------------------------------------------
# Description: Calculate an empirical estimate of chi for annual maximum 
# 			   precipitation at the holdout gauge locations. Run this
# 			   before running plot_chi.R
# Output: Creates a dataset called emp_chi.Rdata in the precipitation/data
#		  directory with emprical estimates of chi_u for different u for the
# 		  annual max precip data. This dataset will be used by plot_chi.
#   -----------------------------------------------------------------------

# Load data ---------------------------------------------------------------
load(file.path(data.dir, "mprecip.Rdata"))

# Define bins -------------------------------------------------------------
maxrng <- 320												# Maximum distance to consider (miles)
bin_centers <- seq(20, maxrng, by = 20)						# Bin centers to calculate chi for
rad <- (diff(bin_centers)/2)[1]								# half-width of bins
bins <- c(bin_centers[-length(bin_centers)] -rad, 			#
	bin_centers[length(bin_centers)] + rad)					#
bins <- bins * 1609.34										# Convert from miles to meters
bins[1] <- 1										        # lower bound is some small number greater than 0
ps <- seq(0.01, 0.98, by = 0.02)							# upper quantiles to calculate chi for
ps <- c(ps, 0.25, 0.5, 0.75, 0.9, 0.98)
ps <- sort(unique(ps))

# Calculate chi -----------------------------------------------------------
emp_chi <- calc_chi_u_spatial(yhold, holdout_coords, bins, ps, bin_centers = bin_centers*1609.34)

# Save output -------------------------------------------------------------
save(emp_chi, file = file.path(samples.dir, "chi_u",
                                paste0("emp_chi.Rdata")))

