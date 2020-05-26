library(stablemix)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(stablemix)
library(RColorBrewer)

#   -----------------------------------------------------------------------
# Description: Plot the empirical and model estimates of chi_u as a function
#              of distance.
# Output: A figure of chi_u as a function of spatial distance for the log-
#         Gaussian process basis theta >0 model. Figure 5 of the paper.
#         The output figure is called chi_10mi.pdf
#   -----------------------------------------------------------------------


# Load data ---------------------------------------------------------------
load(file.path(data.dir, "mprecip.Rdata"))
load(file.path(samples.dir, "chi_u", "emp_chi.Rdata"))

# Parameters --------------------------------------------------------------
ps <- c(0.25, 0.5, 0.75, 0.9, 0.98)
hmax <- 400 * 1609.34      # convert to miles
hmax_plt <- 300
add_title = FALSE
stub <- "lnorm_maxid"
titles <- "log-GP max-id"
st <- 1

# Plot Spatially ----------------------------------------------------------
emp_chi_plt <-
  data.frame(emp_chi %>% filter(u %in% ps) %>% filter(h < hmax))
pl <- list()
load(file.path(samples.dir, "chi_u", paste0(stub, "_chi_low.Rdata")))
hs <- unique(chimoddf$h)

# Summarize quantiles -----------------------------------------------------
chiqs <- chimoddf %>%
  filter(u %in% ps) %>%
  filter(h < hmax) %>%
  group_by(h, u) %>%
  summarize(ub = quantile(chi, 0.975, na.rm = T),
            lb = quantile(chi, 0.025, na.rm = T))

pal <- brewer.pal(n = 7, name = "Set1")[c(1, 2, 3, 4, 7)]
plt <- chiqs %>%
  ggplot() +
  geom_ribbon(
    aes(
      x = h / 1609.34,
      ymin = lb,
      ymax = ub,
      fill = factor(u)
    ),
    alpha = 0.25,
    show.legend = TRUE,
    colour = NA
  ) +
  geom_line(data = emp_chi_plt,
            aes(
              x = h / 1609.34,
              y = chi,
              colour = factor(u)
            ),
            linetype = 1) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0, size = 1))) +
  ylim(c(0, 1)) +
  xlim(c(0, hmax_plt)) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = expression(h ~ (mi)), y = expression(chi[u])) +
  labs(colour = "u") +
  geom_hline(yintercept = 1 - ps,
             linetype = "twodash",
             colour = "grey") +
  if (add_title) {
    plt <- plt + ggtitle(label = titles[st])
  }

plt$labels$fill <- "u"
legend <- get_legend(plt)
plt <- plt + theme(legend.position = "none")
pl[[st]] <- plt


# Plot for fixed space ----------------------------------------------------
unique(emp_chi$h / 1609.34)
bin_id <- 5
dists = c(20, 100, 180, 260)
bins <- unique(emp_chi$bin[which(abs(emp_chi$h / 1609.34 - dists) < 1)])
emp_chi_plt <- subset(emp_chi, bin %in% bins)
emp_chi_plt$dist <- emp_chi_plt$h / 1609.34
pl_u <- list()
load(file.path(samples.dir, "chi_u", paste0(stub, "_chi_low.Rdata")))

bins <-
  unique(chimoddf$bin[which(abs(chimoddf$h / 1609.34 - dists) < 1)])
chimoddf2 <- subset(chimoddf, bin %in% bins)
chiqs <- chimoddf2 %>%
  group_by(h, u)  %>%
  summarize(ub = quantile(chi, 0.975, na.rm = T),
            lb = quantile(chi, 0.025, na.rm = T))
chiqs$dist <- chiqs$h / 1609.34

pal <- brewer.pal(n = 8, name = "Dark2")[c(1, 2, 4)]
pal <- c(pal, brewer.pal(8, name = "Accent")[5])
plt <- chiqs %>%
  ggplot() +
  geom_ribbon(
    aes(
      x = u,
      ymin = lb,
      ymax = ub,
      fill = factor(dist)
    ),
    alpha = 0.25,
    show.legend = TRUE,
    colour = NA
  )  +
  geom_line(data = emp_chi_plt,
            aes(
              x = u,
              y = chi,
              colour = factor(dist)
            ),
            linetype = 1) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 15,
                                margin = margin(
                                  t = 16,
                                  r = 0,
                                  b = 0,
                                  l = 0
                                )),
    axis.title.y = element_text(size = 15),
    text = element_text(size = 15)
  ) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0, size = 1))) +
  ylim(c(0, 1)) +
  xlim(c(0, 1)) +
  geom_line(
    data = data.frame(p = seq(0, 1, l = 100), y = 1 - seq(0, 1, l = 100)),
    aes(x = p, y = y),
    linetype = "twodash",
    colour = "grey"
  ) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = expression(u), y = expression(chi[u])) +
  labs(colour = "h (mi)")

if (add_title) {
  plt <- plt + ggtitle(label = titles[st])
}

plt$labels$fill <- "h (mi)"
plt <- plt + theme(legend.position = c(0.2, 0.4))
pl_u[[st]] <- plt

# Save plot ---------------------------------------------------------------
pgrd <-
  plot_grid(
    pl_u[[1]],
    pl[[1]],
    legend,
    ncol = 3,
    nrow = 1,
    rel_widths = c(1, 1, 0.12)
  )
save_plot(
  pgrd,
  base_aspect_ratio = 2.12,
  base_height = 6,
  filename = "./fig/chi_10mi.pdf"
)
