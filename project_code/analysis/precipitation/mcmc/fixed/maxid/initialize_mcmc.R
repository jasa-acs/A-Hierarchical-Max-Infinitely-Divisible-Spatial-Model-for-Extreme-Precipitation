library(stablemix)


#   -----------------------------------------------------------------------
# Description: This code initializes the parameters sampled in the mcmc
#              algorithm for the Gaussian density basis, theta > 0 model
#              described in the paper.
# 
# Output: A dataset called initial_values1.Rdata containing initial values 
#         for the mcmc for the theta > 0 Gaussian density basis model.
#   -----------------------------------------------------------------------

# File paths --------------------------------------------------------------
stub <- "fixed/maxid"

# Load data ---------------------------------------------------------------
load(file.path(data.dir, "mprecip.Rdata"))

# Set initial values ------------------------------------------------------
it <- 1                                 # Chunk number
colnames(obs_coords) <- NULL            
colnames(knot_coords) <- NULL           
alpha <- 0.1                            # Positive stable index
theta <- 0.1                            # Exponential tilting parameter
D <- calc_dist(obs_coords, knot_coords) # Distance between observation and knot coordinates
tau_min <- max(apply(D, 1, min)) / 2    # Kernel bandwidth min
tau <- tau_min * 2                      # Kernel bandwidth
gvar <- 1                               # GP var for GEV parameters
gscl <- tau_min * 2                     # GP scale for GEV parameters
mar_gp_par <- list(mu = list(gvar = gvar, gscl = gscl), # GEV hyperparameters
                 sigma = list(gvar = gvar, gscl = gscl))
mar_betas <- list(loc = 0, scale = 0, shape = 0) # Linear coefs for GEV parameters
n <- nrow(y)                              # Number of time replicates
nloc <-  ncol(y)                          # number of spatial locations
L <- nrow(knot_coords)                    # number of random effects per year
zall <-                                   # Used to get initial values
  rstabmix(
    n,
    obs_coords,
    knot_coord = knot_coords,
    nbasis = L,
    alpha = alpha,
    delta = alpha,
    theta = theta,
    kern_bw = tau,
    return_all = T,
    type = "smith"
  )
lK <- zall[["lK"]]                       # Log basis matrix
lB <- zall[["lB"]]                       # Log basis scaling factors
loc_form = ~ 1
scale_form = ~ 1
shape_form = ~ 1
mar_gp_which = c("mu", "sigma")
df <- data.frame(x = rep(1, n * nloc))
mu <- sample_marpar_gp(loc_form,
                       df,
                       mar_betas[["loc"]],
                       n,
                       obs_coords,
                       list(gvar = mar_gp_par[["mu"]][["gvar"]],
                            gscl = mar_gp_par[["mu"]][["gscl"]]))
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
# MCMC initial values list
init <- list(
  alpha = alpha,
  theta = theta,
  tau = tau,
  loc = mar_betas[["loc"]],
  scale = mar_betas[["scale"]],
  shape = mar_betas[["shape"]]
)
# Initial values for tuning variance in mcmc
tune_var <- list(
  alpha = 0.1,
  theta = 0.1,
  tau = 0.1,
  loc = 0.1,
  scale = 0.5,
  shape = 0.1,
  lB = 0.1,
  lK = 0.1,
  gps = list(gvar = 0.1, gscl = 0.1, gp = 0.1)
)
# Initialize MCMC ---------------------------------------------------------
out <-
  mcmc_fixed_basis(
    2,
    obs_coords,
    knot_coords,
    y,
    lB,
    init,
    tune_var,
    df,
    mar_gp_par,
    mar_gp_which = mar_gp_which,
    loc = loc_form,
    scale = scale_form,
    shape = shape_form,
    tau_tol = tau_min,
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
    mu = list(gvar = out$smp_mu_gvar[nmc],
              gscl = out$smp_mu_gscl[nmc]),
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
dir.create(file.path(data.dir, stub, "initial_values"), showWarnings = FALSE, recursive = TRUE)
save(it, n, nloc, L, alpha, theta, tau, mar_gp_par, df, mar_gp_which,
     loc_form, scale_form, shape_form, init, lB, clusters,
     tune_var, mar_gp_init, clusters, tau_min,
     file = file.path(data.dir,stub,"initial_values", paste0("initial_values",it,".Rdata"))
)
# Save job list for managing jobs ------------------------------------------
job_list <- it
save(job_list, file = file.path(data.dir, stub,"job_list.Rdata"))
