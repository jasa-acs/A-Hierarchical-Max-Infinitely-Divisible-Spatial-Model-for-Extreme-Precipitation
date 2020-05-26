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
#              calc_chi.R as a function of u for fixed locations. calc_chi.R
#              must be run prior to this script.
#
# Output: Produces Figure 2 of the paper, plots of of the functions chi_u
#         for all four models described in the paper. Figure is stored in
#         ./fig/chi_plot.pdf
#   -----------------------------------------------------------------------

# Parameters -------------------------------------------------------------------
pal <- brewer.pal(n = 7, name = "Set1")[-6]        # Plotting color pal
alpha_plt <- c(0.1, 0.25)                          # Values of PS tail dep to plot
kern_bw_plt <- 1 / 6                               # Gaussian density sd parameter to plot
gvar_plt <- 25                                     # Gaussian
gscl_plt <- 0.75                                   # GP Scale


# Gaussian density basis -------------------------------------------------------
# Load the Monte Carlo samples
load("./data/fixed_chi.Rdata")
# Reformat data for plotting
chidfs <- list()
chibardfs <- list()
for (k in 1:length(chilist_all)) {
  chilist <- chilist_all[[k]]
  chid <- chilist[["chid"]]
  chibard <- chilist[["chibard"]]
  kern_bw <- kern_bws[[k]]
  chidfs[[k]] <- data.frame(chid, bw = kern_bw)
  chibardfs[[k]] <- data.frame(chibard, bw = kern_bw)
}
chidf <- do.call(rbind.data.frame, chidfs)
chibardf <- do.call(rbind.data.frame, chibardfs)

# Max-stable case (Gaussian density) -------------------------------------------
# Make plots
p1 <- chibardf %>%
  filter((alpha == alpha_plt[1]) & (bw == kern_bw_plt)) %>%
  ggplot(aes(
    x = u,
    y = chibar,
    group = factor(theta),
    colour = factor(theta)
  )) +
  ylim(0, 1) +
  geom_line() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 20),
    text = element_text(size = 20)
    ,
    legend.position = "none"
  ) +
  scale_color_manual(values = pal) +
  labs(list(x = "u", y = expression(bar(chi)[u]))) +
  guides(colour = FALSE)

p2 <- chidf %>%
  filter((alpha == alpha_plt[1]) & (bw == kern_bw_plt)) %>%
  ggplot(aes(
    x = u,
    y = chi,
    group = factor(theta),
    colour = factor(theta)
  )) +
  ylim(0, 1) +
  geom_line() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 20),
    text = element_text(size = 20)
    ,
    legend.position = "none"
  ) +
  scale_color_manual(values = pal) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  labs(list(x = "u", y = expression(chi[u]))) +
  labs(colour = expression(theta))


# Max-id case (Gaussian density) -----------------------------------------------
# Make plots
p3 <- chibardf %>%
  filter((alpha == alpha_plt[2]) & (bw == kern_bw_plt)) %>%
  ggplot(aes(
    x = u,
    y = chibar,
    group = factor(theta),
    colour = factor(theta)
  )) +
  ylim(0, 1) +
  geom_line() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 20),
    text = element_text(size = 20)
    ,
    legend.position = "none"
  ) +
  scale_color_manual(values = pal) +
  labs(list(x = "u", y = expression(bar(chi)[u]))) +
  guides(colour = FALSE)

p4 <- chidf %>%
  filter((alpha == alpha_plt[2]) & (bw == kern_bw_plt)) %>%
  ggplot(aes(
    x = u,
    y = chi,
    group = factor(theta),
    colour = factor(theta)
  )) +
  ylim(0, 1) +
  geom_line() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 20),
    text = element_text(size = 20),
    legend.position = c(0.25, 0.25)
  ) +
  scale_color_manual(values = pal) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  labs(list(x = "u", y = expression(chi[u]))) +
  labs(colour = expression(theta))


# Log-gaussian process basis ---------------------------------------------------
# Load the Monte Carlo samples
load("./data/lnorm_chi.Rdata")
chidfs <- list()
chibardfs <- list()
for (k in 1:length(chilist_all)) {
  chilist <- chilist_all[[k]]
  chid <- chilist[["chid"]]
  chibard <- chilist[["chibard"]]
  gvar <- gppars[k, "gvar"]
  gscl <- gppars[k, "gscl"]
  chidfs[[k]] <- data.frame(chid, gvar = gvar, gscl = gscl)
  chibardfs[[k]] <- data.frame(chibard, gvar = gvar, gscl = gscl)
}
chidf <- do.call(rbind.data.frame, chidfs)
chibardf <- do.call(rbind.data.frame, chibardfs)


# Max-stable (log-GP basis) ----------------------------------------------------
# Make plots
p5 <- chibardf %>%
  filter((alpha == alpha_plt[1]) &
           (gvar == gvar_plt) & gscl == gscl_plt) %>%
  ggplot(aes(
    x = u,
    y = chibar,
    group = factor(theta),
    colour = factor(theta)
  )) +
  ylim(0, 1) +
  geom_line() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 20),
    text = element_text(size = 20)
    ,
    legend.position = "none"
  ) +
  scale_color_manual(values = pal) +
  labs(list(x = "u", y = expression(bar(chi)[u]))) +
  guides(colour = FALSE)

p6 <- chidf %>%
  filter((alpha == alpha_plt[1]) &
           (gvar == gvar_plt) & gscl == gscl_plt) %>%
  ggplot(aes(
    x = u,
    y = chi,
    group = factor(theta),
    colour = factor(theta)
  )) +
  ylim(0, 1) +
  geom_line() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 20),
    text = element_text(size = 20)
    ,
    legend.position = "none"
  ) +
  scale_color_manual(values = pal) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  labs(list(x = "u", y = expression(chi[u]))) +
  labs(colour = expression(theta))


# Max-id (log-GP basis) --------------------------------------------------------
p7 <- chibardf %>%
  filter((alpha == alpha_plt[2]) &
           (gvar == gvar_plt) & gscl == gscl_plt) %>%
  ggplot(aes(
    x = u,
    y = chibar,
    group = factor(theta),
    colour = factor(theta)
  )) +
  ylim(0, 1) +
  geom_line() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 20),
    text = element_text(size = 20)
    ,
    legend.position = "none"
  ) +
  scale_color_manual(values = pal) +
  labs(list(x = "u", y = expression(bar(chi)[u]))) +
  guides(colour = FALSE)

p8 <- chidf %>%
  filter((alpha == alpha_plt[2]) &
           (gvar == gvar_plt) & gscl == gscl_plt) %>%
  ggplot(aes(
    x = u,
    y = chi,
    group = factor(theta),
    colour = factor(theta)
  )) +
  ylim(0, 1) +
  geom_line() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 20),
    text = element_text(size = 20)
  ) +
  scale_color_manual(values = pal) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  labs(list(x = "u", y = expression(chi[u]))) +
  labs(colour = expression(theta)) +
  theme(legend.position = "none")

# combine plots -----------------------------------------------------------
pgrd <- plot_grid(
  p1,  p2,  p3,  p4,
  p5,  p6,  p7,  p8,
  nrow = 2,
  ncol = 4,
  rel_widths = c(1, 1, 1, 1)
)

# Save output
dir.create("./fig/", recursive = T, showWarnings = F)
save_plot(
  pgrd,
  base_aspect_ratio = 2.1,
  base_height = 9,
  filename = "./fig/chi_plot.pdf"
)

