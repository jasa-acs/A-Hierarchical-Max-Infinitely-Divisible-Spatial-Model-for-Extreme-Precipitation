library(stablemix)

#   -----------------------------------------------------------------------
# Description: This code runs the mcmc sampler forthe Gaussian process,
#              theta = 0 model described in the paper. This should be run 
#              this should be run in chunks for long jobs.
# Output: A dataset called fixed/maxid<chunk>.RData containing the posterior 
#         samples for the Gaussian process theta = 0 model, where chunk is
#         the number of times this script has been run sequentially if broken
#         up into smaller chunks. Initial values for the next chunk are also
#         stored in a dataset called initial_values<chunk>.Rdata
#   -----------------------------------------------------------------------

# File paths --------------------------------------------------------------
stub <- "lnorm/maxstable"

# Load data ---------------------------------------------------------------
load(file.path(data.dir,"mprecip.Rdata"))                # Precipitation maxima                                         
load(file.path(data.dir,stub,"job_list.Rdata"))          # Job/ chunk iteration info                                         
last_it <- max(job_list)                                 # Get last iteration                   
load(file.path(data.dir,stub,"initial_values",           # Get initial values                                          
  paste0("initial_values",last_it,".Rdata")))            #                                        
it <- last_it + 1                                        # Update iteration           

# Run MCMC ----------------------------------------------------------------
runtime <- system.time(
  out <-
    mcmc_lnorm_basis(                                     
      1000,                                               # MCMC iterations per chunk                                
      obs_coords,                                         # Observation coordinates                          
      y,                                                  # Precipitation maxima                  
      lB,                                                 # log basis scaling factors                  
      lK,                                                 # log basis function matrix                  
      init,                                               # Initial values                    
      tune_var,                                           # Initial tuning variance                        
      df,                                                 # Spatial covariates dataframe (optional)                  
      mar_gp_par,                                         # Marginal GEV parameter GP initial values                          
      clusters = clusters,                                # Spatial block definitions to use for block sampling                                    
      mar_gp_init = mar_gp_init,                          # Marginal GEV GP initial values                                          
      mar_gp_which = mar_gp_which,                        # Which GEV parameters to sample GPs for                                            
      loc = loc_form,                                     # GEV location linear function form (R formula)                              
      scale = scale_form,                                 # GEV log-scale linear function form (R formula)                                  
      shape = shape_form,                                 # GEV shape linear function form (R formula)                                  
      thin_int = 10,                                      # Thinning interval                              
      infer_general_theta = FALSE                         # Not used                                          
    )
)

# Save output as initial values -------------------------------------------
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
# Initial values for next chunk
init <- list(
  alpha = alpha,
  theta = theta,
  gvar = gvar,
  gscl = gscl,
  loc = mu,
  scale = sigma,
  shape = xi
)
lK <- out$smp_lK[nmc, , ]
lB <- matrix(out$smp_lB[nmc, ], nrow = L, ncol = n)
npv <- nrow(out$prop_var)
tune_var <- list(
  alpha = out$prop_var[npv, "alpha"],
  theta = out$prop_var[npv, "theta"],
  gvar = out$prop_var[npv, "gvar"],
  gscl = out$prop_var[npv, "gscl"],
  loc = out$loc_prop_var[npv, ],
  scale = out$scale_prop_var[npv, ],
  shape = out$shape_prop_var[npv, ],
  lB = out$final_lB_prop_var,
  lK = out$final_lK_prop_var,
  gps = list(gvar = 0.1, gscl = 0.1, gp = 0.1)
)
mar_gp_init <- list(
  mu = out$smp_mu_gp[nmc, ],
  sigma = out$smp_sigma_gp[nmc, ],
  xi = rep(0, length(out$smp_sigma_gp[nmc, ]))
)

# Save values -------------------------------------------------------------
# Initial values for next chunk
dir.create(file.path(data.dir,stub,"initial_values"),
           showWarnings = FALSE, recursive = TRUE)
save(it, n, nloc, L, alpha, theta, gvar, gscl, mar_gp_par, df, mar_gp_which,
     loc_form, scale_form, shape_form, init, lK, lB, clusters,
     tune_var, mar_gp_init, clusters,
     file = file.path(data.dir,stub,"initial_values",paste0("initial_values",it,".Rdata"))
)
# Manage job list
job_list <- c(job_list, it)
save(job_list, file = file.path(data.dir,stub,"job_list.Rdata"))
# MCMC samples
dir.create(file.path(samples.dir,stub), showWarnings = FALSE, recursive = TRUE)
save(out, runtime, file = file.path(samples.dir, stub, paste0(it,".Rdata")))
