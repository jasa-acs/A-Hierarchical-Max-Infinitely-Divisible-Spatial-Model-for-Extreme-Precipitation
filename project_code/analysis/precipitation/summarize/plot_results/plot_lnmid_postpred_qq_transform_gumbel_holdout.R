library(stablemix)
library(ggplot2)
library(scales)
library(ggmap)
library(dplyr)
library(utils)
library(gdata)
library(RColorBrewer)
require(gridExtra)
library(abind)
library(fields)
library(cowplot)
library(extRemes)

#   -----------------------------------------------------------------------
# Description: Calculate group-wise min, mean, and maxima for the annual 
#              maxima transformed to the Gumbel scale for the holdout 
#              locations. Compare model and empirical quantiles of each.
#              Scripts in /mcmc/lnorm/maxid and /mcmc_post_proc/lnorm/maxid
#              must be run before this script.
# Output: Figure 6 of the paper. The figure is saved in lnmid_qqplots_gumbel_hold.jpeg
#   -----------------------------------------------------------------------

# Functions ---------------------------------------------------------------
# Group-wise minima
calc_multi_min <- function(y, station_ids) {
  sort(apply(y[, station_ids], 1,
             function(x) {
               if (all(is.na(x))) {
                 return(NA)
               }
               else{
                 min(x, na.rm = T)
               }
             }), na.last = T)
}
# Groupwise means
calc_multi_mean <- function(y, station_ids) {
  sort(apply(y[, station_ids], 1,
             function(x) {
               if (all(is.na(x))) {
                 return(NA)
               }
               else{
                 mean(x, na.rm = T)
               }
             }), na.last = T)
}
# Group-wise maxima
calc_multi_max <- function(y, station_ids) {
  sort(apply(y[, station_ids], 1,
             function(x) {
               if (all(is.na(x))) {
                 return(NA)
               }
               else{
                 max(x, na.rm = T)
               }
             }), na.last = T)
}

# Load data ---------------------------------------------------------------
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(data.dir, "predcoord.Rdata"))
rand_station <- 1:ncol(yhold)

# Load samples ------------------------------------------------------------
load(file.path(samples.dir, "lnorm_maxid_holdout_yhat.Rdata"))
ymc <- csims[["y"]]
nmc <- dim(ymc)[1]
nyr <- dim(ymc)[2]
nloc <- dim(ymc)[3]

# Load GEV samples --------------------------------------------------------
load(file.path(samples.dir, "lnorm_maxid_holdout_gev_par.Rdata"))
mus <- csim_gev$mu
sigmas <- csim_gev$sigma
xis <- csim_gev$xi
nmcmc <- nrow(xis)

# Transform to Gumbel scale using fitted GEVs
gyhold <- array(NA_real_, dim = c(nmcmc, nyr, nloc))
for (i in 1:nmcmc) {
  gyhold[i, ,] <- qevd(
    stablemix::pevdM(
      yhold,
      loc = matrix(
        mus[i,],
        nrow = nyr,
        ncol = nloc,
        byrow = T
      ),
      scale = matrix(
        sigmas[i,],
        nrow = nyr,
        ncol = nloc,
        byrow = T
      ),
      shape = matrix(
        xis[i,],
        nrow = nyr,
        ncol = nloc,
        byrow = T
      )
    ),
    loc = 0,
    scale = 1,
    shape = 0
  )
}

# Transform mcmc postpred samples to Gumbel ----------------------------------------
gymc <- array(NA_real_, dim = c(nmcmc, nyr, nloc))
for (i in 1:nmcmc) {
  gymc[i, ,] <- qevd(
    stablemix::pevdM(
      ymc[i, ,],
      loc = matrix(
        mus[i,],
        nrow = nyr,
        ncol = nloc,
        byrow = T
      ),
      scale = matrix(
        sigmas[i,],
        nrow = nyr,
        ncol = nloc,
        byrow = T
      ),
      shape = matrix(
        xis[i,],
        nrow = nyr,
        ncol = nloc,
        byrow = T
      )
    ),
    loc = 0,
    scale = 1,
    shape = 0
  )
}

# Minima
pred_pmin <- obs_pmin <- matrix(NA_real_, nmcmc, nyr)
for (i in 1:nmcmc) {
  obs_pmin[i,] <- calc_multi_min(gyhold[i, ,], rand_station)
  pred_pmin[i,] <- calc_multi_min(gymc[i, ,], rand_station)
}
obs <- apply(obs_pmin, 2, mean)
est <- apply(pred_pmin, 2, mean)
lb <- apply(pred_pmin,2, quantile, probs = 0.025)
ub <- apply(pred_pmin,2, quantile, probs = 0.975)
mindf <- data.frame(obs = obs, est = est, lb = lb, ub = ub, type = "min")


# Means
pred_pmean <- obs_pmean <- matrix(NA_real_, nmcmc, nyr)
for (i in 1:nmcmc) {
  obs_pmean[i,] <- calc_multi_mean(gyhold[i, ,], rand_station)
  pred_pmean[i,] <- calc_multi_mean(gymc[i, ,], rand_station)
}
obs <- apply(obs_pmean, 2, mean)
est <- apply(pred_pmean, 2, mean)
lb <- apply(pred_pmean, 2, quantile, probs = 0.025)
ub <- apply(pred_pmean, 2, quantile, probs = 0.975)
meandf <- data.frame(obs = obs, est = est, lb = lb, ub = ub, type = "mean")

# Maxima
pred_pmax <- obs_pmax <- matrix(NA_real_, nmcmc, nyr)
for (i in 1:nmcmc) {
  obs_pmax[i,] <- calc_multi_max(gyhold[i, ,], rand_station)
  pred_pmax[i,] <- calc_multi_max(gymc[i, ,], rand_station)
}
obs <- apply(obs_pmax, 2, mean)
est <- apply(pred_pmax, 2, mean)
lb <- apply(pred_pmax, 2, quantile, probs = 0.025)
ub <- apply(pred_pmax, 2, quantile, probs = 0.975)
maxdf <- data.frame(obs = obs, est = est, lb = lb, ub = ub, type = "max")


# Plotting ----------------------------------------------------------------
p1 <- ggplot(mindf, aes(x = est, y = obs)) +
  geom_point(size = 0.8) +
  geom_line(aes(x = est, y = lb), linetype = 2) +
  geom_line(aes(x = est, y = ub), linetype = 2) +
  geom_abline(intercept = 0,
              slope = 1,
              size = 0.2) +
  xlab("Model") +
  ylab("Observed") +
  ylim(c(-3.5, 0)) +
  xlim(c(-3.5, 0)) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    text = element_text(size = 18)
  )

p2 <- ggplot(meandf, aes(x = est, y = obs)) +
  geom_point(size = 0.8) +
  geom_line(aes(x = est, y = lb), linetype = 2) +
  geom_line(aes(x = est, y = ub), linetype = 2) +
  geom_abline(intercept = 0,
              slope = 1,
              size = 0.2) +
  xlab("Model") +
  ylab("Observed") +
  ylim(c(-0.8, 2.8)) +
  xlim(c(-0.8, 2.8)) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    text = element_text(size = 18)
  )

p3 <- ggplot(maxdf, aes(x = est, y = obs)) +
  geom_point(size = 0.8) +
  geom_line(aes(x = est, y = lb), linetype = 2) +
  geom_line(aes(x = est, y = ub), linetype = 2) +
  geom_abline(intercept = 0,
              slope = 1,
              size = 0.2) +
  xlab("Model") +
  ylab("Observed") +
  ylim(c(1, 14)) +
  xlim(c(1, 14)) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    text = element_text(size = 18)
  )
# Combine plots
pgrd <-
  plot_grid(
    p1,
    p2,
    p3,
    nrow = 1,
    ncol = 3,
    rel_widths = c(1, 1, 1, 0.01)
  )


# Save output -------------------------------------------------------------
save_plot(
  pgrd,
  base_aspect_ratio = 3.1,
  base_height = 5,
  filename = file.path(fig.dir, "lnmid_qqplots_gumbel_hold.jpeg")
)
