setwd("~/Desktop/Publications/alphastableTransient")

# Load necessary library
library(ggplot2)
library(tidyverse)
library(cowplot)


f_prime_flipped <- function(x) {
  x^3  - x + k
}

# Define the flipped equation of motion
find_fixed_points <- function(k) {
  # Define the flipped function U'(-x, k)
  # Find fixed points using the uniroot function for multiple initial guesses
  fixed_points <- numeric()
  for (x0 in seq(-3, 3, by = 0.01)) {
    tryCatch({
      x_fp <- uniroot(f_prime_flipped, c(x0 - 0.5, x0 + 0.5))$root
      if (length(fixed_points) == 0 || !any(abs(fixed_points - x_fp) < 1e-5)) {
        fixed_points <- c(fixed_points, x_fp)
      }
    }, error = function(e) {})
  }
  
  return(sort(fixed_points))
}

# Define range for k values
k_values <- c(seq(-1.1, 1.1, b = 0.01), seq(0.37, 0.40, 0.001), seq(-0.37, -0.40, -0.000001))
fixed_points_data <- data.frame(k = numeric(), x = numeric())

# Compute fixed points for each k
for (k in k_values) {
  fixed_points <- find_fixed_points(k)
  for (fp in fixed_points) {
    fixed_points_data <- rbind(fixed_points_data, data.frame(k = k, x = fp))  # No need to flip y-values here
  }
}

upper_stable = fixed_points_data %>%
  filter(x > 0.6)

under_stable = fixed_points_data %>%
  filter(x < -0.6)

unstable = fixed_points_data %>%
  filter(x > -0.6 & x < 0.6)


fixed_points = data.frame(x = c(-1, 0, 1, 
                                -0.585, 1.16,
                                0.585, -1.16,
                                1.33, -1.33),
                         k = c(0, 0, 0, 
                               -0.385, -0.385, 
                               0.385, 0.385,
                               -1, 1),
                         type = c("stable", "unstable", "stable", "ghost", "stable",  "ghost", "stable", "stable", "stable"))

# Plot the flipped bifurcation diagram
(p1 = ggplot() + theme_classic() +
    geom_line(data = upper_stable, aes(x = k, y = x), size = 1.2) +
    geom_line(data = under_stable, aes(x = k, y = x), size = 1.2) +
    geom_line(data = unstable, aes(x = k, y = x), size = 1.2, linetype = "twodash") +
    geom_vline(xintercept = 0,  color = "black", linewidth = .5) +
    geom_vline(xintercept = -0.385,  color = "#D55E00", linewidth = .5) +
    geom_vline(xintercept = 0.385,  color = "#D55E00", linewidth = .5, linetype = "dashed") +
    geom_vline(xintercept = -1, color = "black", linewidth = .5) +
    geom_vline(xintercept = 1, color = "black", linewidth = .5, linetype = "dashed") +
    geom_point(data = fixed_points, aes(x = k, y = x, shape = type), color = "#0072B2", fill = "#0072B2", size = 3, stroke = 1) +
    scale_shape_manual(values = c("stable" = 19, "unstable" = 1, "ghost" = 10), name = "Fixed points") +
    scale_x_continuous(expand = c(0,0), limits = c(-1.1, 1.1), name = "Bifurcation parameter k") +
    scale_y_continuous(expand = c(0,0), limits = c(-1.5, 1.5), name = "Fixed points X*") +
    theme(axis.text = element_text(size = 15),
          axis.title =  element_text(size = 15),
          legend.background = element_rect(fill='transparent', color = NA),
          legend.box.background = element_rect(fill='transparent', color = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),  
          plot.background = element_rect(fill = "transparent", colour = NA),
          strip.background = element_rect(fill = "transparent", color = NA),
          legend.justification = c(1,0),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15)))



potential = function(x,k) {
  
  U = (0.25*x^4 - 0.5*x^2 + k*x)
  return(U)  
}


df = data.frame(x = seq(-3, 3, length.out = 1000)) %>%
  mutate(`0` = potential(x = x, k = 0),
         `-0.385` = potential(x = x, k = -0.385),
         `0.385` = potential(x = x, k = 0.385),
         `1` = potential(x = x, k = 1),
         `-1` = potential(x = x, k = -1)) %>%
  pivot_longer(cols = c(`0`, `-0.385`, `0.385`, `1`, `-1`), names_to = "k", values_to = "U") 

fixed_points = fixed_points %>%
  mutate(y = potential(x, k))

df$k = factor(df$k, levels = c("-1", "-0.385", "0", "0.385", "1"))

(p2 = ggplot() + theme_classic() +
    geom_line(data = df, aes(x = x, y = U, color = k, linetype = k, linewidth = k), size = .7) +
    geom_point(data = fixed_points, aes(x = x, y = y, shape = type), color = "#0072B2", fill = "#0072B2", size = 3, stroke = 1) +
    scale_shape_manual(values = c("stable" = 19, "unstable" = 1, "ghost" = 10), name = "Fixed points") +
    scale_color_manual(values = c("-1" = "black", "-0.385" =  "#D55E00", "0" = "black", "0.385" =  "#D55E00", "1" = "black"), name = "Bifurcation parameter k") +
    scale_linetype_manual(values = c("-1" = "solid", "-0.385" = "solid", "0" = "solid",  "0.385" = "dashed", "1" = "dashed"), name = "Bifurcation parameter k") +
    scale_linewidth_manual(values = c("-1" = 1.2, "-0.385" = 1, "0" = 1.2, "0.385" = 1,  "1" = 1.2), name = "Bifurcation parameter k") +
    scale_x_continuous(limits = c(-2.2, 2.2), expand = c(0,0), breaks = c(-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5), name = "State X") +
    scale_y_continuous(limits = c(-1.5, .5), expand = c(0,0), name = "Potential U(X)") +
    theme(axis.text = element_text(size = 15),
        axis.title =  element_text(size = 15),
        legend.background = element_rect(fill='transparent', color = NA),
        legend.box.background = element_rect(fill='transparent', color = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA),
        strip.background = element_rect(fill = "transparent", color = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15)) +
    guides(shape = "none"))

##

df_ts = read_csv("data/timeseries_jumping_theory.csv")

fp_ts = data.frame(y = c(-1, 0, 1),
                   type = c("stable", "unstable", "stable"))

(p3 = ggplot() + theme_classic() +
    geom_hline(data = fp_ts, aes(yintercept = y, linetype = type), color = "#0072B2",  size = .75) +
    geom_line(data = df_ts, aes(x = timestep, y = state), linewidth = .01, alpha = .9) +
    geom_point(data = fp_ts, aes(x = -Inf, y = y, shape = type), color = "#0072B2", fill = "#0072B2", size = 3, stroke = 1.2) +
    scale_y_continuous(expand = c(0,0), limits = c(-1.5, 1.5), name = "State X") +
    scale_x_continuous(expand = c(0.005,0), limits = c(1927000, 1985000), name = "Simulation timestep", breaks = c(1920000, 1950000, 1980000)) +
    scale_linetype_manual(values = c("stable" = "solid", "unstable" = "dotdash"), name = "Fixed points") +
    scale_shape_manual(values = c("stable" = 19, "unstable" = 1, "ghost" = 10), name = "Fixed points") +
    theme(axis.text = element_text(size = 15),
          axis.title =  element_text(size = 15),
          legend.background = element_rect(fill='transparent', color = NA),
          legend.box.background = element_rect(fill='transparent', color = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),  
          plot.background = element_rect(fill = "transparent", colour = NA),
          strip.background = element_rect(fill = "transparent", color = NA),
          strip.text = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15)) +
    guides(shape = "none"))

plot_grid(p1, p2, p3, align = "hv", axis = "l", nrow = 1, labels = c("(a)", "(b)", "(c)"))
ggsave("bifurcation_diagram.pdf", width = 14, height = 5, scale = 0.7)

