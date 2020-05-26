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
# Description: Plot the posterior predictive mean, sd, and single draw 
#              surfaces for the log-Gaussian process basis theta > 0 model.
#              Scripts in /mcmc/lnorm/maxid and /mcmc_post_proc/lnorm/maxid
#              must be run before this script.
# Output: Figure 8 of the paper. The figure is stored in ppred_mean.jpeg
#   -----------------------------------------------------------------------

# Load data ---------------------------------------------------------------
load(file.path(data.dir, "map_polygon.Rdata"))
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(data.dir, "predcoord.Rdata"))


load(file.path(samples.dir, "lnorm_maxid_csims.Rdata"))
mcit <- 40000                   # Iteration number to plot draw for a single iteration for sanity check
j <- 52                         # Choose a year to plot (index from earliest year)

# Response ----------------------------------------------------------------
yimp <- apply(csims[["y"]], c(2, 3), mean, na.rm = T)
col <-
  c(min(c(csims[["y"]][mcit, j, ], y[j, ]), na.rm = T), 
    quantile(c(csims[["y"]][mcit, j, ], y[j, ]), 0.9992, na.rm = T))
lon = obs_latlon@data$LON
lat = obs_latlon@data$LAT
df <- data.frame(y = y[j, ], lon = lon, lat = lat)

# Plot observation --------------------------------------------------------
p1 <- ggplot(data = df, aes(x = lon, y = lat)) +
  geom_point(aes(colour = y), size = 4) +
  geom_polygon(
    data = state,
    aes(x = long, y = lat, group = group),
    fill = NA,
    color = "black",
    size = 0.2
  ) +
  scale_color_gradientn(colours = tim.colors(n = 1000), limits = col) +
  ggtitle(paste0("Observed Maxima (", 1960 + j, ")")) +
  theme(panel.background = element_blank())  +
  xlab("Lon") +
  ylab("Lat") +
  theme_void() +
  theme(
    plot.title = element_text(size = 35,
                              face = "plain"),
    text = element_text(size = 18),
    legend.key.size = unit(2 / 3, "cm"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.position = c(0.8, 0.2)
  ) +
  coord_quickmap(xlim = c(-80.8,-65.6333),
                 ylim = c(39.7358, 48.25))
p1$labels$fill <- "in"


# Plot mean surface -------------------------------------------------------
df <-
  data.frame(gp = yimp[j, ], lon = pred_latlon$x, lat = pred_latlon$y)
p2 <- ggplot(data = df, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = gp), interpolate = FALSE) +
  geom_polygon(
    data = state,
    aes(x = long, y = lat, group = group),
    fill = NA,
    color = "black",
    size = 0.2
  ) +
  scale_fill_gradientn(colours = tim.colors(n = 1000), limits = col) +
  ggtitle(paste0("Posterior Predictive Mean (", 1960 + j, ")")) +
  theme(panel.background = element_blank())  +
  xlab("") +
  ylab("") +
  theme_void() +
  theme(
    plot.title = element_text(size = 35,
                              face = "plain"),
    text = element_text(size = 18),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.position = c(0.8, 0.2)
  ) +
  coord_quickmap(xlim = c(-80.8,-65.6333),
                 ylim = c(39.7358, 48.25))
p2$labels$fill <- "in"


# Plot single draw --------------------------------------------------------
df <-
  data.frame(gp = csims[["y"]][mcit, j,],
             lon = pred_latlon$x,
             lat = pred_latlon$y)
p3 <- ggplot(data = df, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = gp), interpolate = FALSE) +
  geom_polygon(
    data = state,
    aes(x = long, y = lat, group = group),
    fill = NA,
    color = "black",
    size = 0.2
  ) +
  scale_fill_gradientn(colours = tim.colors(n = 1000), limits = col) +
  ggtitle(paste0("Posterior Predictive Draw (", 1960 + j, ")")) +
  theme(panel.background = element_blank())  +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_void() +
  theme(
    plot.title = element_text(size = 35,
                              face = "plain"),
    text = element_text(size = 18),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.position = c(0.8, 0.2)
  ) +
  coord_quickmap(xlim = c(-80.8,-65.6333),
                 ylim = c(39.7358, 48.25))
p3$labels$fill <- "in"

p1 <- p1 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

pgrd <-
  plot_grid(
    p1,
    p3,
    p2,
    nrow = 1,
    ncol = 3,
    rel_widths = c(1, 1, 1, 0.01)
  )

# Save output
save_plot(
  pgrd,
  base_aspect_ratio = 3.1,
  base_height = 9,
  filename = file.path("./fig", paste0("ppred_mean.jpeg"))
)
