library(stablemix)

#   -----------------------------------------------------------------------
# Description: Generate Monte Carlo estimates of chi_u as a function of distance
#              for the gaussian density basis models discussed in the paper.
#
# Output: Produces a dataset of monte carlo samples of chi_u and stores it in
#         ./data/fixed_chi_spatial.Rdata which is called by 
#          plot_spatial_chi.R.
#   -----------------------------------------------------------------------


# Function to generate the data for spatial plot of chi_u for gauss density basis
# Parameters:
#   us: vector of probabilities
#   hs: vecor of distances from origin
#   knot_coord: d x 2 matrix of basis knot locations
#   kern_bw: Gaussian density kernel bandwidth
#   alpha: model dependence parameter
#   theta: model dependence parameter
fixed_chi_spatial <-
  function(us,
           hs,
           knot_coord,
           kern_bw,
           alpha,
           theta) {
    delta <- alpha
    n <- length(hs)
    chibarm <- chim <- matrix(NA, nrow = length(us), ncol = n)
    for (i in 1:length(us)) {
      for (j in 1:n) {
        qlim <- c(us[i], us[i])
        obs_coord <- c(0, hs[j])
        chis <- fixed_chi_exact(1, qlim, alpha, theta, obs_coord,
                                knot_coord, kern_bw)
        chim[i, j] <- chis$chi
        chibarm[i, j] <- chis$chibar
      }
    }
    dh <- data.frame(dummy = paste0("V", 1:n), h = hs)
    chiw <- cbind.data.frame(chim, us)
    colnames(chiw)[1:n] <- paste0("V", 1:n)
    chilong <- gather(chiw, dummy, chi, V1:V100, factor_key = FALSE)
    chid <- inner_join(chilong, dh, by = "dummy") %>% select(-dummy)
    chibarw <- cbind.data.frame(chibarm, us)
    colnames(chibarw)[1:n] <- paste0("V", 1:n)
    chibarlong <-
      gather(chibarw, dummy, chibar, V1:V100, factor_key = FALSE)
    chibard <-
      inner_join(chibarlong, dh, by = "dummy") %>% select(-dummy)
    return(list(chid = chid, chibard = chibard))
  }


# Parameters --------------------------------------------------------------
us <- c(0.25, 0.5, 0.75, 0.9, 0.98)                # quantiles
hs <- seq(0, 1, l = 100)                           # distances between spatial locs
kern_bws  <-  1 / 6                                # gaussian density sds
knot_coord <- seq(0, 1, l = 36)                    # kernel knot locations
alphas <- c(0.1, 0.25)                             # PS index
thetas <- c(0, 1 / 10, 1e-4)                       # exponential tilting parameter
pars <- expand.grid(alphas, thetas, kern_bws)
colnames(pars) <- c("alpha", "theta", "kern_bw")

# Make Chi data -----------------------------------------------------------
chilist_all <- list()
for (k in 1:nrow(pars)) {
  alpha <- pars[k, "alpha"]
  theta <- pars[k, "theta"]
  kern_bw <- pars[k, "kern_bw"]
  chilist <-
    fixed_chi_spatial(us, hs,  knot_coord, kern_bw, alpha, theta)
  chilist_all[[k]] <- chilist
}

# Save chi ----------------------------------------------------------------
save(us,
     hs,
     knot_coord,
     kern_bws,
     alphas,
     thetas,
     chilist_all,
     pars,
     file = "./data/fixed_chi_spatial.Rdata")
