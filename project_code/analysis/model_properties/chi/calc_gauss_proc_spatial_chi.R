library(stablemix)
library(parallel)
library(doMC)
library(foreach)

#   -----------------------------------------------------------------------
# Description: Generate Monte Carlo estimates of chi_u as a function of distance
#              for the log-gaussian process basis models discussed in the paper.
#
# Output: Produces a dataset of monte carlo samples of chi_u and stores it in
#         ./data/lnorm_chi_spatial.Rdata which is called by 
#          plot_spatial_chi.R.
#   -----------------------------------------------------------------------

# Function to generate the data for spatial plot of chi_u for log-GP basis
# Parameters:
#   us: vector of probabilities
#   hs: vecor of distances from origin
#   L: number of basis functions to use
#   gvar: GP variance parameter
#   gscl: GP range parameter
#   alpha: model dependence parameter
#   theta: model dependence parameter
#   nmc: Number of monte carlo replicates to simulate
lnorm_chi_spatial <-
  function(us,
           hs,
           L,
           gvar,
           gscl,
           alpha,
           theta,
           nmc = 50000) {
    delta <- alpha
    n <- length(hs)
    nq <- length(us)
    chibarm <- chim <- matrix(NA, nrow = length(us), ncol = n)
    chis <- chibars <- array(NA, dim = c(nmc,  n, nq))
    for (j in 1:n) {
      qlim <- range(us)
      obs_coord <- c(0, hs[j])
      for (k in 1:nmc) {
        chicomb <- lnorm_cond_chi_exact(length(us),
                                        qlim,
                                        alpha,
                                        theta,
                                        obs_coord,
                                        L,
                                        gvar,
                                        gscl,
                                        us = us)
        chis[k, j, ] <- chicomb$chi
        chibars[k, j, ] <- chicomb$chibar
      }
      chim[, j] <- apply(chis[, j,], 2, median)
      chibarm[, j] <- apply(chibars[, j,], 2, median)
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
us <- c(0.25, 0.5, 0.75, 0.9, 0.98)                    # quantiles
hs <- seq(0.01, 1, l = 100)                            # distances between spatial locs
gvar <- 25                                             # GP var
gscl <- 3 / 4                                          # GP scale
alphas <- c(0.1, 0.25)                                 # PS index
thetas <- c(0, 1e-4)                                   # Exponential tilting par
L <- 15                                                # Number of basis functions
pars <- expand.grid(alphas, thetas, gvar, gscl)
colnames(pars) <- c("alpha", "theta", "gvar", "gscl")

# Make Chi data -----------------------------------------------------------
chilist_all <- list()
registerDoMC(cores = detectCores())
chilist_all <- foreach(k = 1:nrow(pars)) %dopar% {
  alpha <- pars[k, "alpha"]
  theta <- pars[k, "theta"]
  gvar <- pars[k, "gvar"]
  gscl <- pars[k, "gscl"]
  lnorm_chi_spatial(us, hs,  L, gvar, gscl, alpha, theta, nmc = 50000)
}

# Save chi ----------------------------------------------------------------
dir.create("./data/",
           recursive = T,
           showWarnings = F)
save(us, hs, L, gvar, gscl, alphas, thetas, chilist_all, pars,
     file = "./data/lnorm_chi_spatial.Rdata")
