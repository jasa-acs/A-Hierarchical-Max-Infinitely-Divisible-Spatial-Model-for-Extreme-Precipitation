library(stablemix)

#   -----------------------------------------------------------------------
# Description: This code runs the mcmc sampler forthe Gaussian density basis,
#              theta = 0 model described in the paper. This should be run 
#              this should be run in chunks for long jobs.
# Output: A dataset called fixed/maxid<chunk>.RData containing the posterior 
#         samples for the Gaussian density theta = 0 model,where chunk is
#         the number of times this script has been run sequentially if broken
#         up into smaller chunks. Initial values for the next chunk are also
#         stored in a dataset called initial_values<chunk>.Rdata
#   -----------------------------------------------------------------------

# File paths --------------------------------------------------------------
stub <- "fixed/maxstable"

# Load data ---------------------------------------------------------------
load(file.path(data.dir,"mprecip.Rdata"))             # Precipitation data                        
load(file.path(data.dir,stub,"job_list.Rdata"))       # Manage job list                              
last_it <- max(job_list)                              # Increment the job number        
load(file.path(data.dir,stub,                         # Load the mcmc initial values            
  "initial_values", 
  paste0("initial_values",last_it,".Rdata")))

it <- last_it + 1                                     # Increment the job number
colnames(obs_coords) <- NULL
colnames(knot_coords) <- NULL

# Run MCMC ----------------------------------------------------------------
runtime <- system.time(
  out <-
    mcmc_fixed_basis(
      1000,                             # MCMC iterations per chunk                                            
      obs_coords,                       # Observation coordinates                    
      knot_coords,                      # Knot coordinates                      
      y,                                # Precip observations            
      lB,                               # Log scaling factors            
      init,                             # Initial values for mcmc              
      tune_var,                         # Tuning variances                  
      df,                               # Data frame with spatial covariates (if            
      mar_gp_par,                       # Marginal GEV GP cov parameter list                    
      mar_gp_init = mar_gp_init,        # GEV Initial values                      
      mar_gp_which = mar_gp_which,      # Which GEV parameters are modeled using a GP                        
      loc = loc_form,                   # R formula for GEV location          
      scale = scale_form,               # R formula for GEV log scale              
      shape = shape_form,               # R formula for GEV shape              
      tau_tol = tau_min,                # Minimum value for kernel bw              
      thin_int = 10,                    # Thinning interval          
      clusters = clusters,              # Number of clusters to use when block sampling in space                
      infer_general_theta = FALSE,      # Not used                        
      parallelize = TRUE                # Should parallelization be used where possible              
    )
)

# Save output as initial values for next chunk -------------------------------------
clusters <- out$clusters
nmc <- nrow(out$smp_par)
pars <- out$smp_par[nmc, ]
pnames <- colnames(out$smp_par)
for (i in 1:length(pnames))
  assign(pnames[i], pars[i])
mar_gp_par <-
  list(
    mu = list(gvar = out$smp_mu_gvar[nmc], gscl = out$smp_mu_gscl[nmc]),
    sigma = list(
      gvar = out$smp_sigma_gvar[nmc],
      gscl = out$smp_sigma_gscl[nmc]
    )
  )
init <- list(
  alpha = alpha,
  theta = theta,
  tau = tau,
  loc = mu,
  scale = sigma,
  shape = xi
)
lB <- matrix(out$smp_lB[nmc, ], nrow = L, ncol = n)
npv <- nrow(out$prop_var)
tune_var <- list(
  alpha = out$prop_var[npv, "alpha"],
  theta = out$prop_var[npv, "theta"],
  tau = out$prop_var[npv, "tau"],
  loc = out$loc_prop_var[npv,],
  scale = out$scale_prop_var[npv,],
  shape = out$shape_prop_var[npv,],
  lB = out$final_lB_prop_var,
  lK = out$final_lK_prop_var,
  gps = list(gvar = 0.1, gscl = 0.1, gp = 0.1)
)
mar_gp_init <- list(
  mu = out$smp_mu_gp[nmc,],
  sigma = out$smp_sigma_gp[nmc,],
  xi = rep(0, length(out$smp_sigma_gp[nmc,]))
)

# Save values -------------------------------------------------------------
dir.create(file.path(data.dir,stub,"initial_values"), showWarnings = FALSE, recursive = TRUE)
save(it, n, nloc, L, alpha, theta, tau, mar_gp_par, df, mar_gp_which,
     loc_form, scale_form, shape_form, init, lB, clusters,
     tune_var, mar_gp_init, clusters, tau_min,
     file = file.path(data.dir,stub,"initial_values",paste0("initial_values",it,".Rdata"))
)
# Manage job list
job_list <- c(job_list,it)
save(job_list, file = file.path(data.dir,stub,"job_list.Rdata"))
# MCMC samples
dir.create(file.path(samples.dir,stub), showWarnings = FALSE, recursive = TRUE)
save(out, runtime, file = file.path(samples.dir, stub, paste0(it,".Rdata")))
