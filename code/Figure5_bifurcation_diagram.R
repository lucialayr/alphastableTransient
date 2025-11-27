setwd("~/Desktop/Publications/alphastableTransient")

# Load necessary library
library(ggplot2)
library(tidyverse)
library(cowplot)


f_prime_flipped = function(x, k) {
  x^3  - x + k
}


find_fixed_points = function(k) {
  # For cubic equation x^3 + px + q = 0, where p = -1, q = k
  p = -1
  q = k
  
  # Calculate discriminant: Δ = (q/2)^2 + (p/3)^3
  discriminant = (q/2)^2 + (p/3)^3
  
  fixed_points = numeric()
  
  if (discriminant > 1e-10) {
    # One real root - use robust cubic root function
    cubic_root = function(x) {
      if (x >= 0) {
        return(x^(1/3))
      } else {
        return(-(-x)^(1/3))
      }
    }
    u = cubic_root(-q/2 + sqrt(discriminant))
    v = cubic_root(-q/2 - sqrt(discriminant))
    fixed_points = c(fixed_points, u + v)
  } else if (abs(discriminant) <= 1e-10) {
    # Two real roots (one repeated)
    # When discriminant = 0: roots are 3q/p and -3q/(2p)
    root1 = 3*q/p
    root2 = -3*q/(2*p)
    fixed_points = c(fixed_points, root1, root2)
  } else {
    # Three real roots using trigonometric solution
    rho = sqrt(-(p/3)^3)
    theta = acos(-q/(2*rho))
    
    # Three roots: 2*sqrt(-p/3) * cos((theta + 2*pi*n)/3) for n = 0, 1, 2
    sqrt_term = 2*sqrt(-p/3)
    fixed_points = c(fixed_points, sqrt_term * cos(theta/3))
    fixed_points = c(fixed_points, sqrt_term * cos((theta + 2*pi)/3))
    fixed_points = c(fixed_points, sqrt_term * cos((theta + 4*pi)/3))
  }
  
  # Remove duplicates and sort
  fixed_points = unique(round(fixed_points, 10))
  
  return(sort(fixed_points))
}

# Define range for k values
k_values = c(seq(-1.1, 1.1, by = 0.01), seq(0.37, 0.40, 0.001), seq(-0.37, -0.40, -0.000001))

# Pre-allocate lists for better performance
k_list = list()
x_list = list()

# Compute fixed points for each k
for (k in k_values) {
  fixed_points = find_fixed_points(k)
  if (length(fixed_points) > 0) {
    k_list[[length(k_list) + 1]] = rep(k, length(fixed_points))
    x_list[[length(x_list) + 1]] = fixed_points
  }
}

# Efficiently combine all results
fixed_points_data = data.frame(k = unlist(k_list), x = unlist(x_list))

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
    geom_line(data = upper_stable, aes(x = k, y = x), size = 1.2, color = "black") +
    geom_line(data = under_stable, aes(x = k, y = x), size = 1.2, color = "black") +
    geom_line(data = unstable, aes(x = k, y = x), size = 1.2, linetype = "twodash", color = "black") +
    geom_vline(xintercept = 0,  color = "black", linewidth = .5) +
    geom_vline(xintercept = -0.385,  color = "#0072B2", linewidth = .5) +
    geom_vline(xintercept = 0.385,  color = "#0072B2", linewidth = .5, linetype = "dashed") +
    geom_vline(xintercept = -1, color = "black", linewidth = .5) +
    geom_vline(xintercept = 1, color = "black", linewidth = .5, linetype = "dashed") +
    geom_point(data = fixed_points, aes(x = k, y = x, shape = type), color = "#D55E00", fill = "#D55E00", size = 3, stroke = 1) +
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
    geom_point(data = fixed_points, aes(x = x, y = y, shape = type), color = "#D55E00", fill = "#D55E00", size = 3, stroke = 1) +
    scale_shape_manual(values = c("stable" = 19, "unstable" = 1, "ghost" = 10), name = "Fixed points") +
    scale_color_manual(values = c("-1" = "black", "-0.385" =  "#0072B2", "0" = "black", "0.385" =  "#0072B2", "1" = "black"), name = "Bifurcation parameter k") +
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
        legend.position = 'bottom',
        legend.direction = "horizontal",
        legend.justification = "left",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        legend.key.size = unit(1.2, 'cm')) +
    guides(shape = "none"))

##

df_ts = read_csv("data/timeseries_jumping_theory.csv")

fp_ts = data.frame(y = c(-1, 0, 1),
                   type = c("stable", "unstable", "stable"))

(p3 = ggplot() + theme_classic() +
    geom_hline(data = fp_ts, aes(yintercept = y, linetype = type), color = "#D55E00",  size = .75) +
    geom_line(data = df_ts, aes(x = timestep, y = state), linewidth = .01, alpha = .9, color = "black") +
    geom_point(data = fp_ts, aes(x = -Inf, y = y, shape = type), color = "#D55E00", fill = "#D55E00", size = 3, stroke = 1.2) +
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
          legend.position = "none",
          legend.direction = "horizontal",
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15)) +
    guides(shape = "none"))

plot_grid(p1, p2, p3, align = "hv", axis = "l", nrow = 1, labels = c("(a)", "(b)", "(c)"))
ggsave("figures/bifurcation_diagram.pdf", width = 18, height = 7, scale = 0.7)

