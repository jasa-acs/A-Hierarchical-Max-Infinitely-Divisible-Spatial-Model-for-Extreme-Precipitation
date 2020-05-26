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
# Description: Make a plot of the spatial gauge locations and knot coordinates.
#
# Output: Figure 4 of the paper. The figure is stored in a file called gauge_knot_locations_together.pdf
#   -----------------------------------------------------------------------


# Load data ---------------------------------------------------------------
load(file.path(data.dir, "map_polygon.Rdata"))
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(data.dir, "predcoord.Rdata"))

# Combine knot and observation locations-----------------------------------
df1 <-
  data.frame(
    lon = all_latlon@data$LON,
    lat = all_latlon@data$LAT,
    col = "black",
    sze = 1,
    id = "a",
    stk = 1
  )
df2 <-
  data.frame(
    lon = knot_latlon$x,
    lat = knot_latlon$y,
    col = "red",
    sze = 3,
    id = "b",
    stk = 1.5
  )
df <- rbind.data.frame(df1, df2)

# Plot knot and observation locations -------------------------------------
p <- ggplot(df) +
  geom_point(aes(
    x = lon,
    y = lat,
    group = id,
    shape = id,
    colour = col,
    size = sze,
    stroke = stk
  )) +
  geom_polygon(
    data = state,
    aes(x = long, y = lat, group = group),
    fill = NA,
    color = "black",
    size = 0.1
  ) +
  theme(panel.background = element_blank())  +
  scale_color_manual(
    values = c("black", "red"),
    labels = c("Gauge", "Knot"),
    guide = FALSE
  ) +
  scale_shape_manual(values = c(20, 4), labels = c("Gauge", "Knot")) +
  scale_size(guide = FALSE) +
  guides(shape = guide_legend(override.aes = list(
    shape = c(20, 4),
    colour = c("black", "red"),
    size = c(3, 3)
  ))) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(
    plot.title = element_text(size = 25,
                              face = "plain"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    text = element_text(size = 18),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.position = c(0.8, 0.2)
  ) +
  coord_quickmap(xlim = c(-80.8,-65.6333),
                 ylim = c(39.7358, 48.25))
p$labels$shape <- "Type"

pgrd <- plot_grid(p,
                  nrow = 1,
                  ncol = 2,
                  rel_widths = c(1, 0.01))

# Save plot ---------------------------------------------------------------
save_plot(
  pgrd,
  base_aspect_ratio = 1.2,
  base_height = 9,
  filename = file.path(fig.dir, "gauge_knot_locations_together.pdf")
)
