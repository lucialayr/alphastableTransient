setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/transient_alphastable")

install.packages("scico")
install.packages("cowplot")

library(tidyverse)
library(scico)
library(cowplot)

theme_set(
  theme_classic() + 
    theme(
      axis.text = element_text(color = "black", size = 15),
      axis.title = element_text(color = "black", size = 15),
      plot.title = element_text(color = "black", size = 15),
      plot.subtitle = element_text(color = "black", size = 15),
      plot.caption = element_text(color = "black", size = 15),
      strip.text = element_text(color = "black", size = 15),
      legend.text = element_text(color = "black", size = 15),
      legend.title = element_text(color = "black", size = 15),
      axis.line = element_line(color = "black"),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
      legend.background = element_rect(fill='transparent', color = NA),
      legend.box.background = element_rect(fill='transparent', color = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),  
      plot.background = element_rect(fill = "transparent", colour = NA),
      strip.background = element_rect(fill = "transparent", color = NA)
    )
)


colors = c("2" = "black", "1.8" =  "#0072B2", "1.5" = "#009E73", "1.3" =  "#56B4E9")

##

estimate_potential = function(x) {
  
  U = -0.5*log(x)
  
  return(U)
}

#### LPJ-GUESS


estimate_density_lpjguess = function(df, t1, t2) {
  df = df %>%
    filter(TG > t1 & TG < t2) %>%
    group_by(Year) %>%
    summarize(
      density_data = list(density(relative, kernel = "gaussian", width = 0.3) %>% 
                            with(data.frame("relative" = x, "density"  = y)))
    ) %>%
    unnest(cols = c(density_data)) %>%
    mutate(potential = estimate_potential(density),
           tas = paste0(t1 - 273, "° C - ", t2 - 273, "° C"))
  
  return(df)
}



potentials_lpj_guess= function() {
  
  df3 = read_csv("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/transient_alphastable/data/processed/lpjguess_2d.csv") %>%
    estimate_density_lpjguess(273, 273 + 2)
  
  df2 = read_csv("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/transient_alphastable/data/processed/lpjguess_2d.csv") %>%
    estimate_density_lpjguess(273 + 3, 273 + 5)
  
  df1 = read_csv("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/transient_alphastable/data/processed/lpjguess_2d.csv") %>%
    estimate_density_lpjguess(273 + 7, 273 + 9)
  
  
  df = bind_rows(df1, df2, df3) %>%
    group_by(relative, tas) %>%
    mutate(mean = mean(potential, na.rm = T),
           model = "LPJ-GUESS")
  
  (p = ggplot() + 
      geom_line(data = df, aes(x = relative, y = potential, color = tas, group = interaction(Year, tas)), linewidth = .2, alpha = .4) +
      geom_line(data = df, aes(x = relative, y = mean, color = tas), linewidth = .8) +
      scale_y_continuous(limits = c(-0.70, 1.7), expand = c(0,0), name = expression("Estimated potential"~hat(U)~"("~chi[BNE]^C~")")) +
      scale_x_continuous(limits = c(-0.25, 1.2), expand = c(0,0), breaks = c(0, 0.5, 1), name = expression(Share~of~needleleaf~trees~chi[BNE]^C)) +
      scale_color_scico_d(palette = "lajolla", begin = .2, end = .8, direction = -1, name = "Temperature range") +
      scale_fill_scico_d(palette = "lajolla", begin = .2, end = .8, direction = -1, name = "Temperature range") +
      facet_wrap(~model) +
      theme(legend.background = element_rect(fill='transparent', color = NA),
            legend.box.background = element_rect(fill='transparent', color = NA),
            legend.position = "bottom",
            legend.title.position = "top",
            legend.direction = "horizontal",
            panel.background = element_rect(fill = "transparent", colour = NA),  
            plot.background = element_rect(fill = "transparent", colour = NA),
            strip.background = element_rect(fill = "transparent", color = NA)) +
      guides(color = guide_legend(nrow = 3)))
  
  
  return(p)
  
}

(p1 = potentials_lpj_guess())

### Potential models



k = -1
a = 2
simulated_densities = function() {
  
  data = list()
  
  for (k in c(-1, -0.39, 0)) {
    
    for (a in c(2, 1.8, 1.5, 1.3)) {
      
      df = read_csv(paste0("data/final_states_a", a, "_k", k, ".csv")) %>%
        filter(abs(final_state) < 3) %>%
        pull(final_state) %>%                
        density(kernel = "gaussian", bw = 0.2) %>%  
        with(data.frame("state" = x, "density"  = y)) %>%
        mutate(potential = estimate_potential(density),
               k_label = paste0("k = ", k),
               alpha = a)
      
      
      data = append(data, list(df))
    }
  }
  
  df = purrr::reduce(data, bind_rows) 
  
  df$alpha = factor(df$alpha, levels = c("2", "1.8",  "1.5", "1.3"))
  df$k_label = factor(df$k_label, levels = c("k = -1", "k = -0.39", "k = 0"))
  
  fixed_points = data.frame(fp = c(-1, 0, 1, 
                                  -0.585, 1.16,
                                  -0.585, 1.33),
                            k_label = c("k = 0", "k = 0", "k = 0",
                                  "k = -0.39", "k = -0.39",
                                  "k = -1", "k = -1"),
                            stability = c("stable", "unstable", "stable", "ghost", "stable",  "ghost", "stable"))
  
  fixed_points$k_label = factor(fixed_points$k_label, levels = c("k = -1", "k = -0.39", "k = 0"))
  
  (p = ggplot() + 
      coord_cartesian(clip = "off") +
      geom_vline(data = fixed_points, aes(xintercept = fp, linetype = stability), color = "#D55E00", linewidth = 0.5) +
      geom_point(data = fixed_points, aes(y = -Inf, x = fp, shape = stability), color = "#D55E00", fill = "#D55E00", size = 3, stroke = 1) +
      geom_line(data = df, aes(x = state, y = potential, color = alpha, linewidth = alpha)) +
      facet_wrap(~k_label, scales = "free_y", nrow = 1) +
      scale_x_continuous(limits = c(-2., 2), expand = c(0,0), name = "State X") +
      scale_y_continuous(expand = c(0,0), limits = c(-.5, 5), name = expression("Estimated potential"~hat(U)~"("~X~")"),) +
      scale_color_manual(values = rev(c("2" = "black", "1.8" =  "#0072B2", "1.5" = "#009E73", "1.3" =  "#56B4E9")),
                         name = expression(alpha~of~noise)) +
      scale_linewidth_manual(values = rev(c("2" = 1, "1.8" =  0.75, "1.5" = 0.75, "1.3" =  0.75)),
                             name = expression(alpha~of~noise)) +
      scale_linetype_manual(values = c("stable" = "solid", "unstable" = "dotted", ghost = "dotdash"), name = "Fixed points") +
      scale_shape_manual(values = c("stable" = 19, "unstable" = 1, "ghost" = 10), name = "Fixed points") +
      theme(legend.position = "bottom",
            legend.title.position = "top",
            legend.direction = "horizontal",
            legend.box = "vertical"))
  
  return(p)
  
}

(p2 = simulated_densities())

plot_grid(p1, p2, axis = "l", align = "hv", rel_widths = c(0.5, 1), labels = c("(a)", "(b)"))

ggsave("figures/estimated_potential.pdf", height = 5.5)

##exploring the denstities
df = read_csv("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/transient_alphastable/data/processed/lpjguess_2d.csv") %>%
  filter(TG > 274.5) %>%
  group_by(Year) %>%
  mutate(bin = cut_width(TG, width = 1, center = 1))

df_potential = df

(p = ggplot() + theme_bw() +
    geom_histogram(data = df, aes(x = relative))) +
  facet_wrap(bin ~ Year, ncol = 11, scales = "free_y")
