setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/patch_analysis_paper")
source("code/utils.R")

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/transient_alphastable")

install.packages("scico")
install.packages("cowplot")
install.packages("ggnewscale")

library(duckdb)
library(purrr)
library(scico)
library(terra)
library(sf)
library(cowplot)
library(tidyverse)
library(MASS)
library(ggnewscale)


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


########### TRANSIENTS SIMPLE

plot_simple_model = function(k, t1 = 176000, t2 = 193000, fp = c(-1, 0, 1), stability = c("stable", "unstable", "stable")) {
  
  data = list()
  
  for (a in c(2, 1.8, 1.5, 1.3)) {
    
    for (k in c(k)) {
      
      df = read_csv(paste0("data/trajectories_a", a, "_k", k, ".csv")) %>%
        mutate(timestep = row_number()) %>%
        pivot_longer(cols = -timestep, names_to = "run", values_to = "state") %>%
        mutate(a = a,
               k = k) %>%
        filter(run == 8)
      
      data = append(data, list(df))
    }
  }
  
  df = purrr::reduce(data, bind_rows) %>%
    filter(timestep > t1 & timestep < t2)
  
  df$a = factor(df$a, levels = rev(c("2", "1.8",  "1.5", "1.3")))
  
 (p1 = ggplot() + 
    coord_cartesian(clip = "off") +
    geom_hline(aes(yintercept = fp, linetype = stability), color = "#D55E00", linewidth = 0.5) +
    geom_line(data = df, aes(x = timestep, y = state, color = a, linewidth = a)) +
    geom_point(aes(x = -Inf, y = fp, shape = stability), color = "#D55E00", fill = "#D55E00", size = 3, stroke = 1) +
    scale_color_manual(values = (c("2" = "black", "1.8" =  "#0072B2", "1.5" = "#009E73", "1.3" =  "#56B4E9")),
                       name = expression(alpha~of~noise)) +
    scale_linetype_manual(values = c("stable" = "solid", "unstable" = "dotted", ghost = "dotdash"), name = "Fixed points") +
    scale_shape_manual(values = c("stable" = 19, "unstable" = 1, "ghost" = 10), name = "Fixed points") +
    scale_linewidth_manual(values = (c("2" = 0.3, "1.8" =  0.3, "1.5" = 0.3, "1.3" =  0.3)),
                           name = expression(alpha~of~noise)) +
    scale_x_continuous(expand = c(0,0), name = "Simulation timestep") +
    scale_y_continuous(limits = c(-2.5, 2), breaks = c(-2, -1, 0, 1, 2), expand = c(0,0), name = "State X") +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          plot.margin = unit(c(0.25, 0, 0, 0), "cm")) +
     guides(color = guide_legend(override.aes = list(linewidth = 2))))

  
  ((p2 = ggplot() + 
      coord_flip(clip = "off") +
      geom_vline(aes(xintercept = fp, linetype = stability), color = "#D55E00", linewidth = 0.5) +
      geom_density(data = df, aes(x = state, color = a, fill = a, linewidth = a), alpha = .05,  bw = .1) +
      geom_point(aes(y = -Inf, x = fp, shape = stability), color = "#D55E00", fill = "#D55E00", size = 3, stroke = 1.2) +
      scale_y_continuous(expand = c(0,0), breaks = c(0), name = expression("Density "~hat(p)~"("~X~")")) +
      scale_x_continuous(limits = c(-2.5, 2), breaks = c(-2, -1, 0, 1, 2),  expand = c(0,0), name =  "") +
      scale_color_manual(values = rev(c("2" = "black", "1.8" =  "#0072B2", "1.5" = "#009E73", "1.3" =  "#56B4E9")),
                         name = expression(alpha~of~noise)) +
      scale_linetype_manual(values = c("stable" = "solid", "unstable" = "dotted", ghost = "dotdash"), name = "Fixed points") +
      scale_shape_manual(values = c("stable" = 19, "unstable" = 1, "ghost" = 10), name = "Fixed points") +
      scale_fill_manual(values = rev(c("2" = "black", "1.8" =  "#0072B2", "1.5" = "#009E73", "1.3" =  "#56B4E9")),
                        name = expression(alpha~of~noise)) +
      scale_linewidth_manual(values = rev(c("2" = 0.75, "1.8" =  0.5, "1.5" = 0.5, "1.3" =  0.5)),
                             name = expression(alpha~of~noise)) +
      theme(legend.position = "none",
            legend.direction = "vertical",
            plot.margin = unit(c(0.25, 0, 0, 0), "cm"))))
  
  p = plot_grid(p1, p2, align = "hv", rel_widths = c(1, 0.3))
  
  return(p)
}

(pA = plot_simple_model(-0.39, fp = c(1.16, -0.585), stability = c("stable", "ghost")))

(pB = plot_simple_model(0, t1 = 150000, t2 = 167000))

########### TRANSIENT VEGGIE

read_vegetation_data = function(variable) {
  scenario = "ssp126"
  
  con = dbConnect(duckdb(), dbdir = "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/patch_analysis_paper/patches2.duckdb", read_only = FALSE) #create the database
  dbListTables(con)
  
  
  #these were created by manually looking at subsets of the data and finding trajectories that were able to fully recover twice in a row (pretty rare in this setting ..)
  #see bottom of the script
  
  long_transients = data.frame(Lon = c(-149.25, -157.25, -156.75, -147.75, -148.75, -149.75, -146.25, -147.75, -148.75, -149.75, -142.25, -141.75),
                               Lat = c(64.75, 66.75, 66.25, 64.75, 64.75, 64.25, 66.75, 66.25, 65.25, 64.25, 67.25, 67.25),
                               PID = c(1, 1, 2, 3, 4, 7, 10, 11, 12, 15, 16, 16)) %>%
    unique()
  
  
  dbWriteTable(con, "long_transient_ssp126_pretty", long_transients, overwrite = T)
  
  df_ts_transient = dbGetQuery(con, paste0("SELECT d.Year, d.Lon, d.Lat, d.PID, d.PFT, d.age, d.cmass, d.ndist FROM '", scenario, "_d150_cmass' 
                                AS d INNER JOIN long_transient_ssp126_pretty AS l ON d.PID = l.PID AND d.Lon = l.Lon AND d.Lat = l.Lat")) %>%
    group_by(Year, Lon, Lat, PID) %>%
    mutate(relative = cmass/sum(cmass))  %>% 
    ungroup() %>%
    filter(PFT %in% c("BNE", "IBS"),
           PID %in% c(1), #1, 4, 7
           Lon < -148 & Lon > -150, 
           Year %in% seq(1850, 2300)) %>%
    mutate(across(everything(), ~ifelse(is.na(.), 0, .))) %>% #if sum(cmass) = 0, this will be NA (can happen in the first years after a disturbance)
    unique() %>%
    mutate(PFT = long_names_pfts(tolower(PFT))) 
  
  df_ts_transient$PFT = factor(df_ts_transient$PFT, levels = c("Pioneering broadleaf", "Needleleaf evergreen"))
  
  return(df_ts_transient)
}

vegetation_ts = function(df) {
  
  (p = ggplot() + 
     geom_line(data = df, aes(x = Year, y = relative, color = PFT, linewidth = PFT), alpha = .9) +
     scale_x_continuous(expand = c(0,0), name = "Simulation year") +
     scale_y_continuous(limits = c(-0.1, 1.1),  breaks = c(0, 0.5, 1), expand = c(0,0), name =  expression(PFT~share~chi[i]^C)) +
     scale_color_manual(name = "Plant functional types (PFTs)", drop = TRUE,
                        values = c("Needleleaf evergreen" = "#0072B2", "Pioneering broadleaf" = "#E69F00")) +
     scale_linewidth_manual(name = "Plant functional types (PFTs)", values = c("Needleleaf evergreen" = 0.75, "Pioneering broadleaf" = 0.5)) +
     theme(legend.position = "bottom",
           legend.direction = "horizontal"))
  
  return(p)
}

vegetation_density = function(df) {
  
  (p = ggplot() +
     coord_flip() +
     geom_density(data = df, aes(x = relative, color = PFT, fill = PFT, linewidth = PFT), alpha = .05,  bw = .1) +
     scale_y_continuous(expand = c(0,0), breaks = c(0), name = expression("Density "~hat(p)~"("~chi[i]^C~")")) +
     scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(-0.1, 1.1), expand = c(0,0), name = "") +
     scale_fill_manual(name = "Plant functional types (PFTs)", drop = TRUE,
                       values = c("Needleleaf evergreen" = "#0072B2", "Pioneering broadleaf" = "#E69F00"),
                       breaks = c( "Needleleaf evergreen", "Pioneering broadleaf")) +
     scale_color_manual(name = "Plant functional types (PFTs)", drop = TRUE,
                        values = c("Needleleaf evergreen" = "#0072B2", "Pioneering broadleaf" = "#E69F00"),
                        breaks = c( "Needleleaf evergreen", "Pioneering broadleaf")) +
     scale_linewidth_manual(name = "Plant functional types (PFTs)", values = c("Needleleaf evergreen" = 0.75, "Pioneering broadleaf" = 0.5)) +
     theme(legend.position = "none",
           legend.direction = "horizontal"))
  
  return(p)
}


df_vegetation = read_vegetation_data("cmass")
(p1 = vegetation_ts(df_vegetation))
(p2 = vegetation_density(df_vegetation))


pC = plot_grid(p1, p2, rel_widths = c(1, 0.3), align = "hv", axis = "l")

########### STICH AND SAVE

plot_grid(pC, pB, pA, ncol = 1, labels = c("(a)", "(b)", "(c)"), vjust = 1)

ggsave("figures/trajectories_transients.pdf", height = 8.5, scale = 1)

##################
#OLD
##################
###looking for long transients
fully_recovered = dbGetQuery(con, paste0("SELECT Lon, Lat, PID FROM '", scenario, "_d150_", variable, "' WHERE age = 300 AND Year > 2100" )) %>%
  unique() 

dbWriteTable(con, "fully_recovered", fully_recovered, overwrite = T)

df_cmass = dbGetQuery(con, paste0("SELECT d.Year, d.PFT, d.PID, d.Lon, d.Lat, d.cmass, d.age FROM '", scenario, "_d150_cmass' 
                                AS d INNER JOIN fully_recovered AS l ON d.PID = l.PID AND d.Lon = l.Lon AND d.Lat = l.Lat"))

df_cmass_test = df_cmass %>%
  filter(PID == 16 & Lon < -140 & Lon > -160) %>%
  group_by(Year, Lon, Lat, PID) %>%
  mutate(relative = cmass/sum(cmass))  %>% 
  ungroup() %>%
  filter(PFT == "BNE") %>%
  mutate(across(everything(), ~ifelse(is.na(.), 0, .))) %>% #if sum(cmass) = 0, this will be NA (can happen in the first years after a disturbance)
  unique()


ggplot() + theme_classic() +
  geom_line(data = df_cmass_test, aes(x = Year, y = relative, group = interaction(Lon, Lat, PID))) +
  geom_point(data = df_cmass_test[df_cmass_test$age > 300,], aes(x = Year, y = relative), color = "red", shape = 4) +
  scale_x_continuous(breaks = c(1900, 2100)) +
  facet_wrap(Lon ~ Lat)


long_transients = data.frame(Lon = c(-149.25, -157.25, -156.75, -147.75, -148.75, -149.75, -146.25, -147.75, -148.75, -149.75, -142.25, -141.75),
                             Lat = c(64.75, 66.75, 66.25, 64.75, 64.75, 64.25, 66.75, 66.25, 65.25, 64.25, 67.25, 67.25),
                             PID = c(1, 1, 2, 3, 4, 7, 10, 11, 12, 15, 16, 16)) %>%
  unique()



hist(df_cmass[df_cmass$Lon < -140 & df_cmass$Lon > -160 & df_cmass$PFT == "BNE" & df_cmass$Year == 2100,]$cmass)


df_cmass_subset = df_cmass %>%
  filter(Lon < -140 & Lon > -160 & Lat > 64 & Lat < 67) %>%
  group_by(Year, Lon, Lat, PID) %>%
  mutate(relative = cmass/sum(cmass))  %>% 
  ungroup() %>%
  filter(PFT == "BNE") %>%
  mutate(across(everything(), ~ifelse(is.na(.), 0, .))) %>% #if sum(cmass) = 0, this will be NA (can happen in the first years after a disturbance)
  unique()

hist(df_cmass_subset$relative, breaks = 30)



