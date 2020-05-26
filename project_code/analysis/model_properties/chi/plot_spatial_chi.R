library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(stablemix)
library(fields)
library(RColorBrewer)

#   -----------------------------------------------------------------------
# Description: plot the Monte Carlo estimates of chi_u generated in the script
#              calc_spatial_chi.R as a function spatial distance for select
#              quantiles. calc_gauss_density_spatial_chi.R and
#              calc_gauss_proc_spatial_chi.R must be run before calling this
#              script
# 
# Output: Produces Figure 3 of the paper. Plots of chi_u as a function of distance
#         for the four models described in the paper. Figure is stored in 
#         ./fig/spatial_chi_plots.pdf
#   -----------------------------------------------------------------------

# log-GP Basis ------------------------------------------------------------
# Load data
load("./data/lnorm_chi_spatial.Rdata")

# Combine data ------------------------------------------------------------
chidfs <- list()
chibardfs <- list()
for (k in 1:length(chilist_all)) {
  chilist <- chilist_all[[k]]
  chid <- chilist[["chid"]]
  chibard <- chilist[["chibard"]]
  theta <- pars[k, "theta"]
  gvar <- pars[k, "gvar"]
  gscl <- pars[k, "gscl"]
  alpha <- pars[k, "alpha"]
  chidfs[[k]] <-
    data.frame(
      chid,
      theta = theta,
      alpha = alpha,
      gvar = gvar,
      gscl = gscl
    )
  chibardfs[[k]] <-
    data.frame(
      chibard,
      theta = theta,
      alpha = alpha,
      gvar = gvar,
      gscl = gscl
    )
}
chidf <- do.call(rbind.data.frame, chidfs)
chibardf <- do.call(rbind.data.frame, chibardfs)


# Make plots --------------------------------------------------------------
pal <- brewer.pal(n = 5, name = "Set1")
p1 <- chidf %>% filter(theta == 0) %>%
  ggplot(aes(
    x = h,
    y = chi,
    colour = factor(us),
    linetype = factor(alpha)
  )) +
  geom_line(size = 0.7) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  ylim(c(0, 1)) +
  scale_color_manual(values = pal) +
  labs(list(x = "h", y = expression(chi[u]))) +
  labs(colour = "u") +
  geom_hline(yintercept = 1 - us,
             linetype = "twodash",
             colour = "grey") +
  theme(legend.position = "none")

p2 <- chidf %>% filter(theta == 1e-4) %>%
  ggplot(aes(
    x = h,
    y = chi,
    colour = factor(us),
    linetype = factor(alpha)
  )) +
  geom_line(size = 0.7) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  ylim(c(0, 1)) +
  scale_color_manual(values = pal) +
  labs(list(x = "h", y = expression(chi[u]))) +
  labs(colour = "u") +
  geom_hline(yintercept = 1 - us,
             linetype = "twodash",
             colour = "grey")
p2$labels$linetype <- expression(alpha)
legend <- get_legend(p2)
p2 <- p2 +  theme(legend.position = "none")


# Gaussian Density Basis --------------------------------------------------
load("./data/fixed_chi_spatial.Rdata")
# Combine data ------------------------------------------------------------
chidfs <- list()
chibardfs <- list()
for (k in 1:length(chilist_all)) {
  chilist <- chilist_all[[k]]
  chid <- chilist[["chid"]]
  chibard <- chilist[["chibard"]]
  theta <- pars[k, "theta"]
  kern_bw <- pars[k, "kern_bw"]
  alpha <- pars[k, "alpha"]
  chidfs[[k]] <-
    data.frame(chid,
               theta = theta,
               alpha = alpha,
               bw = kern_bw)
  chibardfs[[k]] <-
    data.frame(chibard,
               theta = theta,
               alpha = alpha,
               bw = kern_bw)
}
chidf <- do.call(rbind.data.frame, chidfs)
chibardf <- do.call(rbind.data.frame, chibardfs)

# Make plots --------------------------------------------------------------
p3 <- chidf %>% filter(theta == 0) %>%
  ggplot(aes(
    x = h,
    y = chi,
    colour = factor(us),
    linetype = factor(alpha)
  )) +
  geom_line(size = 0.7) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  ylim(c(0, 1)) +
  scale_color_manual(values = pal) +
  labs(list(x = "h", y = expression(chi[u]))) +
  labs(colour = "u") +
  geom_hline(yintercept = 1 - us,
             linetype = "twodash",
             colour = "grey") +
  theme(legend.position = "none")


p4 <- chidf %>% filter(theta == 1e-4) %>%
  ggplot(aes(
    x = h,
    y = chi,
    colour = factor(us),
    linetype = factor(alpha)
  )) +
  geom_line(size = 0.7) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  ylim(c(0, 1)) +
  scale_color_manual(values = pal) +
  labs(list(x = "h", y = expression(chi[u]))) +
  labs(colour = "u") +
  geom_hline(yintercept = 1 - us,
             linetype = "twodash",
             colour = "grey") +
  theme(legend.position = "none")

# Combine plots -----------------------------------------------------------
pgrd <- plot_grid(p3, p4, legend, p1, p2, NULL, nrow = 2, ncol = 3, rel_widths = c(1,1,0.2))
# Save plots
dir.create("./fig", showWarnings = F, recursive = T)
save_plot(pgrd, base_aspect_ratio = 1.1,base_height = 9,
          filename = "./fig/spatial_chiplots.pdf")
