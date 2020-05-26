library(stablemix)
library(ggplot2)
library(scales)
library(ggmap)
library(dplyr)
library(utils)
library(gdata)
library(RColorBrewer)
library(fields)
require(gridExtra)
library(cowplot)
library(raster)
library(rgdal)

#   -----------------------------------------------------------------------
# Description: Plot model posterior spatial factor means, ranked by the 
#              year-to-year variance of the basis scaling factors. All
#              scripts in /mcmc/lnorm/maxid and /mcmc_post_proc/lnorm/maxid
#              must be run before this script.
# Output: Figure 9 of the paper. The top basis factor means are plotted in a
#         figure called top_K.jpeg
#   -----------------------------------------------------------------------

# Parameters --------------------------------------------------------------
ntop_var <- 6     # number of top spatial factors to plot      

# Load data ---------------------------------------------------------------
load(file.path(data.dir, "map_polygon.Rdata"))
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(data.dir, "predcoord.Rdata"))
load(file.path(samples.dir, "lnorm_maxid.Rdata"))
load(file.path(samples.dir, "cond_sim", "lnorm_maxid_csims.Rdata"))

# Rank the factors --------------------------------------------------------
nmcmc <- nrow(out$smp_par)                                  # Number of mcmc samples
L <- dim(out$smp_lK)[3]                                     # Number of spatial factors
nrep <- dim(y)[1]                                           # Number of time replicates
vars <- matrix(NA, nrow = nmcmc, ncol = L)                  # Initialize variance matrix
for (i in 1:nmcmc) {                                        # Calculate sample variance of poserior basis scaling facotrs
  B <- matrix(exp(out$smp_lB[i, ]), nrow = L, ncol = nrep)
  vars[i, ] <- apply(B, 1, var)
}
vmn <- apply(vars, 2, mean)                                 # Estimate variances
svmn <- sort(vmn, decreasing = T)                           # Sum of variances  
print(round(svmn / sum(svmn), 2))                           # Proportion of variance sum in each factor


# Plot factor variance order ----------------------------------------------
jpeg(
  file.path(fig.dir, "factor_var_plots.jpeg"),
  height = 500,
  width = 2 * 500
)
par(mfrow = c(1, 2))
plot(vmn, ylab = "B Variance", xlab = "K Index")
plot(rev(sort(vmn)), ylab = "B Variance", xlab = "Ordered Factor by Variance")
dev.off()
which.top <- order(vmn, decreasing = T)[1:ntop_var]

df <- data.frame(rank = 1:length(vmn), Bvar = rev(sort(vmn)))
varplt <- df %>%
  ggplot() +
  geom_point(aes(x = rank,
                 y = Bvar)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  labs(list(x = "Ordered Factor Indices", y = "B Variance"))


# Plot factor means -------------------------------------------------------
Kimp <- apply(exp(csims[["lK"]]), c(2, 3), mean)
p <- list()
col <- range(Kimp[, which.top])
for (i in 1:ntop_var) {
  df <-
    data.frame(gp = Kimp[, which.top[i]],
               lon = pred_latlon$x,
               lat = pred_latlon$y)
  p[[i]] <- ggplot(data = df, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = gp), interpolate = FALSE) +
    geom_polygon(
      data = state,
      aes(x = long, y = lat, group = group),
      fill = NA,
      color = "black",
      size = 0.2
    ) +
    scale_fill_gradientn(colours = tim.colors(n = 1000), lim = col) +
    ggtitle("") +
    theme(panel.background = element_blank()) +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      text = element_text(size = 20),
      
      legend.key.size = unit(2 / 3, "cm"),
      legend.text = element_text(size = 10),
      legend.position = "none"
    ) +
    xlab("") +
    ylab("") +
    theme_void() +
    theme(legend.position = "none") +
    coord_quickmap(xlim = c(-80.8,-65.6333),
                   ylim = c(39.7358, 48.25))
  p[[i]]$labels$fill <- "K(s)"
}

# Grab legend -------------------------------------------------------------
df <-
  data.frame(gp = Kimp[, which.top[i]], lon = pred_latlon$x, lat = pred_latlon$y)
pnull <- ggplot(data = df, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = gp), interpolate = FALSE) +
  geom_polygon(
    data = state,
    aes(x = long, y = lat, group = group),
    fill = NA,
    color = "black",
    size = 0.2
  ) +
  scale_fill_gradientn(colours = tim.colors(n = 1000), lim = col) +
  ggtitle("") +
  theme(panel.background = element_blank()) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    text = element_text(size = 20),
    legend.key.size = unit(2 / 3, "cm"),
    legend.text = element_text(size = 20),
    legend.position = c(0, 0.5)
  ) +
  xlab("Lon") +
  ylab("Lat") +
  coord_quickmap()
pnull$labels$fill <- "K(s)"

# Combine plots and save --------------------------------------------------
legend <- get_legend(pnull)
legend$vp$justification <- c(0, 0.5)
pgrd <- plot_grid(
  p[[1]],
  p[[2]],
  p[[3]],
  legend,
  p[[4]],
  p[[5]],
  p[[6]],
  NULL,
  ncol = 4,
  nrow = 2,
  rel_widths  = c(1, 1, 1, .23),
  labels = c('1', '2', '3', '', '4', '5', '6', '')
)
save_plot(
  pgrd,
  base_aspect_ratio = 3.3 / 2,
  base_height = 9,
  filename = file.path(fig.dir, "top_K.jpeg")
)
