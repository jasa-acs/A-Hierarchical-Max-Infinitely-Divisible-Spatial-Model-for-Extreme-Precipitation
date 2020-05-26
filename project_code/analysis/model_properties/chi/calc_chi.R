library(stablemix)


#   -----------------------------------------------------------------------
# Description: Generate Monte Carlo estimates of chi_u  for all four models
#              discussed in the paper at two fixed spatial locations.
#              The output is saved and used in plot chi.
#
# Output: Two datasets ./data/fixed_chi.Rdata and ./data/lnorm_chi.Rdata
#         which are called by plot_chi.R to create Figure 2 of the paper.
#   -----------------------------------------------------------------------

# Parameters --------------------------------------------------------------
nq <- 1000                            # number of quantiles to evaluate chi and chibar at
qlim <- c(0.0001, 0.9999)             # limits of quantiles
obs_coord <- c(0, 1/4)                # observation coordinates
knot_coord <- seq(0, 1,  l= 25)       # knot coordinates
kern_bws  <- c(1/6, 1/4, 1/3)         # kernel bandwidths
alphas <- c(0.1, 0.25)                # Positive Stable index
thetas <- c(0, 1e-10, 1e-8, 1e-4, 1e-2, 1e-1)    # Exponential tilting parameter

# Make chis ---------------------------------------------------------------
chilist_all <- list()
for(k in 1:length(kern_bws)) {
  kern_bw <- kern_bws[k]
  chilist <- prepare_fixed_chis(nq, qlim, obs_coord, knot_coord, kern_bw, alphas, thetas)
  chilist_all[[k]] <- chilist
}
# Save output
save(nq, qlim, obs_coord, knot_coord, kern_bws, alphas, thetas, chilist_all,
     file = "./data/fixed_chi.Rdata")


# Log-GP basis
# Parameters --------------------------------------------------------------
nq <- 1000                             # number of quantiles to evaluate chi and chibar at
qlim <- c(0.0001, 0.9999)              # limits of quantiles
obs_coord <- c(0, 1/4)                 # observation coordinates
gvars <- 25                            # Gaussian Process variance
gscls <- 0.75                          # GP Scale
alphas <- c(0.1, 0.25)                 # Positive stable index
thetas <- c(0, 1e-10, 1e-8, 1e-4, 1e-2, 1e-1) # Exponential tilting parameter
L <- 15                                # Number of basis functions
gppars <- expand.grid(gvars, gscls)
colnames(gppars) <- c("gvar","gscl")

# Make chis ---------------------------------------------------------------
chilist_all <- list()
for(k in 1:nrow(gppars)) {
  gvar <- gppars[k,1]
  gscl <- gppars[k,2]
  chilist <- prepare_lnorm_chis(nq, qlim, obs_coord,
                                L, gvar, gscl, alphas, thetas, nmc = 1000)
  chilist_all[[k]] <- chilist
}
# Save output
save(nq, qlim, obs_coord, L, gvars, gscls, alphas, thetas, chilist_all,
     gppars, file = "./data/lnorm_chi.Rdata")
