library(stablemix)

#   -----------------------------------------------------------------------
# Description: This code initializes the parameters sampled in the mcmc
#              algorithm for the Gaussian process basis, theta = 0 model
#              described in the paper.
# Output: A dataset called initial_values1.Rdata containing initial values 
#         for the mcmc for the theta = 0 Gaussian process basis model.
#   -----------------------------------------------------------------------


# File paths --------------------------------------------------------------
stub <- "lnorm/maxstable"

# Load data ---------------------------------------------------------------
load(file.path(data.dir,"mprecip.Rdata"))

# Set initial values ------------------------------------------------------
it <- 1                                                        # Chunk number                               
colnames(obs_coords) <- NULL                                   # Remove any names                                                                               
alpha <- 0.4                                                   # Positive stable parameter                                    
theta <- 0                                                     # Exponential tilting parameter                                  
D <- dist(obs_coords)                                          # Distance of matrix between obs coordinates                                             
gvar <- 20                                                     # GP variance                                            
gscl <- max(D) / 2                                             # GP scale                                               
mar_gp_par <- list(mu = list(gvar = gvar, gscl = gscl),        # GEV Marginal params' GP hyper params                                                                                 
                   sigma = list(gvar = gvar, gscl = gscl))                                                                                        
mar_betas <- list(loc = 1, scale = 0.7, shape = 0)             # Linear coef for GEV GP means                                                                              
                                                                                        
n <- nrow(y)                                                   # number of spatial locations                                     
nloc <- ncol(y)                                                # number of random effects per year                                          
L <- 15                                                        # Sample a process for init of lK and lB terms                                              
zall <- rstabmix(n, obs_coords, nbasis = L, alpha = alpha,                                                                                        
                 delta = alpha, theta = 0.5,                                                                                        
                 gvar = gvar, gscl = gscl, return_all = T,                                                                                        
                 type = "br")                                                                                                                                            
lK <- zall[["lK"]]                                            # Log spatial factors                                                                                                 
lB <- zall[["lB"]]                                            # Log spatial factor coefficients                                                                                                                                                                                     
loc_form = ~ 1                                                # R formula for GEV location GP mean      
scale_form = ~1                                               # R formula for GEV log scale GP mea      
shape_form = ~1                                               # R formula for GEV scale form      
mar_gp_which = c("mu","sigma")                                                      
lB[] <- 0   
df <- data.frame(x = rep(1,n*nloc))
# GEV location
mu <- sample_marpar_gp(loc_form,
                       df,
                       mar_betas[["loc"]],
                       n,
                       obs_coords,
                       list(gvar = mar_gp_par[["mu"]][["gvar"]], 
                            gscl = mar_gp_par[["mu"]][["gscl"]]))
# GEV scale
sigma <- sample_marpar_gp(
  scale_form,
  df,
  mar_betas[["scale"]],
  n,
  obs_coords,
  list(gvar = mar_gp_par[["sigma"]][["gvar"]], 
       gscl = mar_gp_par[["sigma"]][["gscl"]]),
  log.link = TRUE
)
# MCMC initial values
init <- list(
  alpha = alpha,
  theta = theta,
  gvar = gvar,
  gscl = gscl,
  loc = mar_betas[["loc"]],
  scale = mar_betas[["scale"]],
  shape = mar_betas[["shape"]]
)
# Initial tuning variances
tune_var <-
  list(
    alpha = 0.1,
    theta = 0.1,
    gvar = 0.1,
    gscl = 0.1,
    loc = 0.1,
    scale = 0.1,
    shape = 0.001,
    lB = 1e-10,
    lK = 0.1,
    gps = list(gvar = 0.1, gscl = 0.1, gp = 0.1)
  )
# Initialize MCMC ---------------------------------------------------------
out <-
  mcmc_lnorm_basis(
    10,
    obs_coords,
    y,
    lB,
    lK,
    init,
    tune_var,
    df,
    mar_gp_par,
    mar_gp_which = mar_gp_which,
    loc = loc_form,
    scale = scale_form,
    shape = shape_form,
    thin_int = 1
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
dir.create(file.path(data.dir,stub, "initial_values"), showWarnings = FALSE, recursive = TRUE)
save(it, n, nloc, L, alpha, theta, gvar, gscl, mar_gp_par, df, mar_gp_which,
     loc_form, scale_form, shape_form, init, lK, lB, clusters,
     tune_var, mar_gp_init, clusters,
     file = file.path(data.dir,stub,"initial_values",paste0("initial_values",it,".Rdata"))
)

# Save job list for managing jobs ------------------------------------------
job_list <- it
save(job_list, file = file.path(data.dir,stub,"job_list.Rdata"))