library(stablemix)
library(fields)
library(cowplot)

#   -----------------------------------------------------------------------
# Description: Simulate a realization of the max-id and max-stable models
#              with Gaussian density and log-Gaussian process basis functions
#
# Output: 
#     - Figure in ./fig containing sample realizations of the four models
#       described in the paper. See, e.g., Figure 1 of the paper.
#   -----------------------------------------------------------------------

# Simulate Realizations ---------------------------------------------------
seed.start <- 2019
# Gaussian density basis --------------------------------------------------
n <- 1                                    # Number of processes to simulate
nloc <- 100^2                             # Number of spatial locations
alphas <- c(0.1, 0.25)                    # Positive Stable index
thetas <- c(0, 1e-4)                      # Exponential tilting parameter
pars <- expand.grid(alphas, thetas)
colnames(pars) <- c("alpha", "theta")
L <- 25                                   # Number of basis functions to use
tau <-  1/6                               # Gaussian density sd
s <- seq(0, 1, l = sqrt(nloc))            # Spatial locations
s <- as.matrix(expand.grid(s,s))
mar_betas <-list(loc = 0, scale = 0, shape = 0)  # linear coef. for marginal GEV param.
v <- seq(0, 1, l = sqrt(L))               # knot locations
v <- as.matrix(expand.grid(v,v))
ylist_sm <- list()
for(k in 1:nrow(pars)){
  seed <- seed.start
  alpha <- pars[k,"alpha"]
  theta <- pars[k,"theta"]
  # Simulate process
  zall <- rstabmix_fix_seed(n, s, v, nbasis = L, alpha = alpha,
                            delta = alpha, theta = theta,
                            kern_bw = tau, return_all = T,
                            type = "smith", seed = seed)
  lZ <- zall[["lZ"]]
  lK <- zall[["lK"]]
  mu <- matrix(mar_betas[["loc"]], nrow = 1, ncol = nloc)
  sigma <- exp(matrix(mar_betas[["scale"]], nrow = 1, ncol = nloc))
  xi <- matrix(mar_betas[["shape"]], nrow = 1, ncol = nloc)
  # Transform to desired GEV margins
  ylist_sm[[k]] <- lsm2gev(lZ, alpha, theta, lK, mu, sigma, xi)
}

# Log GP Basis -----------------------------------------------------------
n <- 1                                      # Number of processes to simulate
nloc <- 100^2                               # Number of spatial locations
L <- 15                                     # Positive Stable index
s <- seq(0, 1, l = sqrt(nloc))              # Exponential tilting parameter
s <- as.matrix(expand.grid(s,s))
mar_betas <-list(loc = 0, scale = 0, shape = 0)  # linear coef. for marginal GEV param.
gvar <- 25                                  # GP variance
gscl <- 0.75                                # GP scale
alphas <- c(0.1, 0.25)                      # Positive stable index
thetas <- c(0, 1e-4)                        # Exponential tilting parameter
L <- 15
pars <- expand.grid(alphas, thetas)
colnames(pars) <- c("alpha", "theta")

ylist_br <- list()
for(k in 1:nrow(pars)){
  seed <- seed.start
  alpha <- pars[k,"alpha"]
  theta <- pars[k,"theta"]

  zall <- rstabmix_fix_seed(n, s, nbasis = L, alpha = alpha,
                   delta = alpha, theta = theta,
                   gvar = gvar, gscl = gscl, return_all = T,
                   type = "br", seed = seed)
  lZ <- zall[["lZ"]]
  lK <- zall[["lK"]]
  mu <- matrix(mar_betas[["loc"]], nrow = 1, ncol = nloc)
  sigma <- exp(matrix(mar_betas[["scale"]], nrow = 1, ncol = nloc))
  xi <- matrix(mar_betas[["shape"]], nrow = 1, ncol = nloc)
  ylist_br[[k]] <- lsm2gev(lZ, alpha, theta, lK, mu, sigma, xi)
}

# Plotting ----------------------------------------------------------------
# Gaussian density basis -------------------------------------------------------------------
pl_sm <- list()
text.size <- 5
for(i in seq_along(ylist_sm)){
  df <- data.frame(s, y = t(ylist_sm[[i]]))
  colnames(df) <- c("lon","lat","y")
  p <- df %>%
    ggplot(aes(
      x = lon,
      y = lat
    )) +
    geom_raster(aes(fill = y), interpolate = FALSE) +
    scale_fill_gradientn(colours = tim.colors(n = 1000))+
    theme_minimal() +
    theme(
      axis.text = element_text(size = text.size),
      axis.title = element_text(size = text.size),
      text = element_text(size = text.size),
      legend.key.size = unit(3/4, "line"),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_blank()
    ) +

    guides(colour = guide_legend(override.aes = list(size = 1))) +
    ylim(c(0,1)) +
    scale_colour_gradientn(colours = tim.colors(101), lim = c(-2,8)) +
    labs(x = "", y = "") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme_void()
  legend <- get_legend(p)
  p <- p +  theme(legend.position = "none")
  pl_sm[[i]] <- p
}

# Log GP Basis -----------------------------------------------------------------
pl_br <- list()
text.size <- 5
for(i in seq_along(ylist_br)){
  df <- data.frame(s, y = t(ylist_br[[i]]))
  colnames(df) <- c("lon","lat","y")
  p <- df %>%
    ggplot(aes(
      x = lon,
      y = lat
    )) +
    geom_raster(aes(fill = y), interpolate = FALSE) +
    scale_fill_gradientn(colours = tim.colors(n = 1000))+
    theme_minimal() +
    theme(
      axis.text = element_text(size = text.size),
      axis.title = element_text(size = text.size),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_blank(),
      axis.ticks = element_blank(),
      text = element_text(size = text.size),
      legend.key.size = unit(1, "line"),
      legend.text = element_text(size = 30),
      legend.title = element_text(size = 20)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 1))) +
    ylim(c(0,1)) +
    scale_colour_gradientn(colours = tim.colors(101), lim = c(-2,8)) +
    labs(x = "", y = "") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.ticks = element_blank()) +
    theme_void()
  p$labels$fill <- expression(tilde(Z))
  legend <- get_legend(p)
  p <- p +  theme(legend.position = "none")
  pl_br[[i]] <- p
}
pgrd <- plot_grid(pl_sm[[1]],pl_sm[[2]], pl_sm[[3]],pl_sm[[4]], legend,
                  pl_br[[1]],pl_br[[2]], pl_br[[3]],pl_br[[4]], NULL,
                  ncol = 5, nrow = 2, rel_widths = c(1,1,1,1, 0.2))#, labels = c("A", "B"))


# Save plots --------------------------------------------------------------
dir.create("./fig/",recursive = T, showWarnings = F)
save_plot(pgrd, base_aspect_ratio = 4.3/2,base_height = 5,
          filename = "./fig/realizations.pdf")
