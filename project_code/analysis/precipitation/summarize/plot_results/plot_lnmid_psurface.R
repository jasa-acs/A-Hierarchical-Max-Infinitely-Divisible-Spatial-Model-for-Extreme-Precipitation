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
library(raster)
library(rgdal)

#   -----------------------------------------------------------------------
# Description: Plot posterior predictive estimate and uncertainty 
#              of 99th percentile using the log-Gaussian process basis,
#              theta > 0 model. 
#              Scripts in /mcmc/lnorm/maxid and /mcmc_post_proc/lnorm/maxid
#              must be run before this script.
# Output: Figure 7 of the paper. The figure is stored in lnmidQ99.jpeg
#   -----------------------------------------------------------------------


# Load data ---------------------------------------------------------------
load(file.path(data.dir, "map_polygon.Rdata"))
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(data.dir, "predcoord.Rdata"))
load(file.path(samples.dir, "lnorm_maxid_csims.Rdata"))


# Q99 ---------------------------------------------------------------------
u <- 0.99                                     # Quantile to plot
nmc <- nrow(csims[["mu"]])                    # Number of monte carlo samples
nloc <- ncol(csims[["mu"]])                   # Number of spatial locations
qs_all <- matrix(NA, nrow = nmc, ncol = nloc) # Initialize quantile matrix
for (j in 1:nmc) {      
  qs_all[j, ] <- qevdM(                       # Estimate quantile
    matrix(u, nrow = 1,
           ncol = nloc),
    loc = matrix(csims[["mu"]][j, ], nrow = 1),
    scale = matrix(csims[["sigma"]][j, ], nrow = 1),
    shape = matrix(csims[["xi"]][j, ], nrow = 1)
  )
}
qs <- apply(qs_all, 2, mean)
qs_sd <- apply(qs_all, 2, sd)
df <-
  data.frame(
    est = c(qs),
    sd = c(qs_sd),
    lon = pred_latlon$x,
    lat = pred_latlon$y
  )
# Plot the result
pest <- ggplot(data = df, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = est), interpolate = FALSE) +
  geom_polygon(
    data = state,
    aes(x = long, y = lat, group = group),
    fill = NA,
    color = "black",
    size = 0.2
  ) +
  scale_fill_gradientn(colours = tim.colors(n = 1000)) +
  ggtitle(expression(Q[99] ~ Posterior ~ Mean)) +
  theme(panel.background = element_blank())  +
  xlab("") +
  ylab("") +
  theme_void() +
  theme(
    plot.title = element_text(size = 25),
    text = element_text(size = 18),
    legend.key.size = unit(2 / 3, "cm"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.position = c(0.8, 0.2)
  ) +
  coord_quickmap(xlim = c(-80.8,-65.6333),
                 ylim = c(39.7358, 48.25))
pest$labels$fill <- "in"

psd <- ggplot(data = df, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = sd), interpolate = FALSE) +
  geom_polygon(
    data = state,
    aes(x = long, y = lat, group = group),
    fill = NA,
    color = "black",
    size = 0.2
  ) +
  scale_fill_gradientn(colours = tim.colors(n = 1000)) +
  ggtitle(expression(Q[99] ~ Posterior ~ SD)) +
  theme(panel.background = element_blank())  +
  xlab("") +
  ylab("") +
  theme_void() +
  theme(
    plot.title = element_text(size = 25),
    text = element_text(size = 18),
    legend.key.size = unit(2 / 3, "cm"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.position = c(0.8, 0.2)
  ) +
  coord_quickmap(xlim = c(-80.8,-65.6333),
                 ylim = c(39.7358, 48.25))
psd$labels$fill <- "in"

# Combine plots -----------------------------------------------------------
pgrd <-
  plot_grid(
    pest,
    psd,
    nrow = 1,
    ncol = 3,
    rel_widths = c(1, 1, 0.01)
  )


# Save output -------------------------------------------------------------
save_plot(
  pgrd,
  base_aspect_ratio = 2.1,
  base_height = 9,
  filename = file.path(fig.dir, "lnmidQ99.jpeg")
)
