setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/transient_alphastable")
source("code/utils.R")

install.packages("scico")
install.packages("cowplot")

library(duckdb)
library(purrr)
library(scico)
library(terra)
library(sf)
library(cowplot)
library(tidyverse)
library(MASS)

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


#select a subset of the data to plot against temperature
year = 2015
scenario = "ssp126"
variable = "cmass"

con = dbConnect(duckdb(), dbdir = "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/patch_analysis_paper/patches2.duckdb", read_only = FALSE) #create the database
dbListTables(con)

ecoregions = read_csv("data/ext/above_ecoregion_lpjguess_grid.csv") %>%
  filter(ecoregion %in% c("TAIGA PLAIN", "TAIGA CORDILLERA"))

ecoregions_shp = st_read("data/processed/above_ecoregion_lpjguess_grid.shp") %>%
  filter(NA_L2NAME%in% c("TAIGA PLAIN", "TAIGA CORDILLERA")) %>%
  terra::vect()

df_vegetation = dbGetQuery(con, paste0("SELECT Year, Lon, Lat, PID, PFT, age, cmass, ndist FROM '", scenario, "_d300' 
                                 WHERE Year > 2090 AND Year < 2100 AND Lon < -102 AND Lon > -163 AND Lat > 54 AND Lat < 68")) %>%
  right_join(ecoregions) %>%
  group_by(Year, Lon, Lat, PID) %>%
  mutate(relative = cmass/sum(cmass))  %>% 
  ungroup() %>%
  filter(PFT == "BNE") %>%
  mutate(across(everything(), ~ifelse(is.na(.), 0, .))) %>% #if sum(cmass) = 0, this will be NA (can happen in the first years after a disturbance)
  unique() 
    

(p1 = ggplot() + coord_flip() +
  geom_histogram(data = df_vegetation, aes(x = relative, y = after_stat(density)), color = "black", fill = "#0072B2", bins = 40, binwidth = .04) +
  geom_density(data = df_vegetation, aes(x = relative, group = Year), color = "black", alpha = 0, linewidth = 0.25) +
  scale_x_continuous(expand = c(0,0), name = "", position = "top") +
  scale_y_reverse(expand = c(0,0), name = expression("Density "~hat(p)~"("~chi[BNE]^C~")")) +
  theme(legend.position = "None",
        plot.margin = unit(c(0,-0.3,0,0), "cm"),
        legend.direction = "vertical"))




get_2d_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
} #adapted from https://slowkow.com/notes/ggplot2-color-by-density/



rast_climate = terra::rast(paste0("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/patch_analysis_paper/data/covariates/mri-esm2-0_r1i1p1f1_ssp126_tas_daily_inverted_1850_2300_boreal_yearlyavg_growingseason.nc")) %>%
  terra::mask(ecoregions_shp)

rast_climate = rast_climate[[time(rast_climate) > as.Date("2014-01-05") & time(rast_climate) < as.Date("2026-12-31")]]


df_climate = as.data.frame(rast_climate, xy = TRUE)

names(df_climate) = c("Lon", "Lat", seq(2090, 2101))

df_climate = df_climate %>%
  pivot_longer(cols = -c(Lon, Lat), names_to = "Year", values_to = "TG") %>%
  filter(Year != "2101") %>%
  mutate(Year = as.numeric(Year))


df_density = df_climate %>%
  right_join(df_vegetation) %>%
  st_drop_geometry() 


write_csv(df_density, "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/ge96dul2/transient_alphastable/data/processed/lpjguess_2d.csv")


bw_x <- max(sd(df_density$TG) * 0.25, 0.25)  # Set minimum bandwidth
bw_y <- max(sd(df_density$relative) * 0.25, 0.25)  # Set minimum bandwidth

# Kernel density estimation
df_density$density = get_2d_density(df_density$TG, df_density$relative, n = 100, lims = c(range(df_density$TG), range(df_density$relative)), h = c(bw_x, bw_y))

(p2 = ggplot() + 
  geom_point(data = df_density, aes(x = TG, y = relative, color = density), size = .2) +
  scale_color_scico(palette = "lajolla", direction = -1,  name = expression("Density "~hat(p)~"("~chi[BNE]^C~","~k[T[G]]~")"), breaks = c(0, 0.25, 0.5)) +
  scale_x_continuous(expand = c(0,0), breaks = c(272, 274, 276, 278, 280, 282), labels = c(272, 274, 276, 278, 280, 282) - 273,
                      name = expression("Growing season temperature"~k[T[G]]~" in ° C")) +
  scale_y_continuous(expand = c(0.01,0.01), name = expression(atop("Share of ", needleleaf~trees~chi[BNE]^C))) +
  theme(legend.position = "bottom",
        plot.margin = unit(c(.2,0,0,0.2), "cm"),
        axis.text.x = element_text(hjust = 1, size = 15),
        legend.direction = "horizontal"))


plot_grid(p1, p2, ncol = 2, labels = c("(a)", "(b)"), 
          align = "hv", axis = "tblr", rel_widths = c(0.7, 1),
          hjust = c(-3, -1))


ggsave("figures/lpj_bimodal.pdf", height = 5, width = 10)



create_one_scenario_density = function(var_climate, variable,  scenario, year, pft, data){
  
  df_climate = terra::rast(paste0("data/external/mri-esm2-0_r1i1p1f1_", scenario, "_", var_climate, "_", year, "_cropped.nc")) %>%
    terra::as.polygons(aggregate = F) %>%
    st_as_sf() %>%
    st_make_valid() %>%
    mutate(Lon = round(st_coordinates(st_centroid(geometry))[, "X"], 2),
           Lat = round(st_coordinates(st_centroid(geometry))[, "Y"], 2))
  
  if (pft == "otherC") {
    
    df_vegetation = dbGetQuery(con, paste0("SELECT Year, Lon, Lat, PID, PFT, ", variable, " FROM '", scenario, "_d150_", variable, "' WHERE Year = ", year, "")) %>%
      group_by(Lon, Lat, PID) %>%
      mutate(relative = !!rlang::sym(variable)/sum(!!rlang::sym(variable))) %>%
      filter(PFT %in%c("BINE", "TeNE", "BNS")) %>%
      group_by(PFT, Lon, Lat, PID, Year) %>%
      summarize(relative = sum(relative)) %>%
      mutate(relative = replace_na(relative, 0))
    
  } else {
    
    df_vegetation = dbGetQuery(con, paste0("SELECT Year, Lon, Lat, PID, PFT, ", variable, " FROM '", scenario, "_d150_", variable, "' WHERE Year = ", year, 
    " AND Lon < -100 AND Lat > 60 AND Lat < 70")) %>%
      group_by(Lon, Lat, PID) %>%
      mutate(relative = !!rlang::sym(variable)/sum(!!rlang::sym(variable))) %>%
      filter(PFT == pft) %>%
      mutate(relative = replace_na(relative, 0))
    
  }
  
  print(head(df_vegetation))
  
  df = df_climate %>%
    dplyr::select(tas, Lon, Lat) %>%
    right_join(df_vegetation) %>%
    st_drop_geometry()
  
  bw_x <- max(sd(df$tas) * 0.5, 0.5)  # Set minimum bandwidth
  bw_y <- max(sd(df$relative) * 0.5, 0.5)  # Set minimum bandwidth
  
  # Kernel density estimation
  df$density = get_2d_density(df$tas, df$relative, n = 100, lims = c(range(df$tas), range(df$relative)), h = c(bw_x, bw_y))
  
  data = append(data, list(df))
  
  return(data)
}

density_2d_data = function(variable, scenario, pft) {
  data = list()
  
  data = create_one_scenario_density("tas", variable, scenario, 2015, pft, data) #todo look at 1850, climate data needs to be adapted for this
  data = create_one_scenario_density("tas", variable, scenario, 2050, pft, data)
  data = create_one_scenario_density("tas", variable, scenario, 2100, pft, data)
  
  print(data)
  
  df = purrr::reduce(data, bind_rows) %>%
    mutate(tas = tas - 273.15)
  
  write_csv(df, paste0("data/processed/P1.C_", variable, "_", pft, ".csv"))
  
  return(df)
}

density_2d_plot = function(variable, scenario, pft) {
  
  print(pft)
  
  if (file.exists(paste0("data/processed/density_", variable, "_", pft, ".csv"))) {
    print("Loading data")
    df = read_csv(paste0("data/processed/density_", variable, "_", pft, ".csv"), show_col_types = F)
  } else {
    print("Creating data")
    df = density_2d_data(variable, scenario, pft)
  }
  
  print(head(df))
  
  (p = ggplot(df, aes(x=tas, y=relative) ) + 
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = .5, ymax = Inf), fill = "grey95") +
      geom_point(aes(color = density), stroke = 0, shape = 16, size = .5) +
      geom_hline(yintercept = 0.5, color = "grey40", linetype = "dotted") +
      facet_wrap(~Year, ncol = 1) +
      scale_x_continuous(name = "Mean annual temperature in °C", expand = c(0,0)) +
      scale_y_continuous(name = bquote(Share ~ "of" ~ .(pft) ~ "per patch" ~ chi[~.(pft)]), breaks = c(0, 0.5, 1), expand = c(0,0)) +
      scico::scale_color_scico(palette = "lipari", direction = -1, name = "Density", begin = .3) +
      add_common_layout(fontsize = 15) +
      theme(text = element_text(size = 15)) +
      annotate("text", x = 8, y = 0.55, label = "Dominance"))
}

density_2d_plot("cmass", "ssp126", "BNE")


##################
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



