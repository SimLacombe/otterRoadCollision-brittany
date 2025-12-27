library(tidyverse)
library(sf)
library(patchwork)
library(tmaptools)

rm(list = ls())

source("fun/myUtilities.R")
source("fun/plot_functions.R")

dpts <- c(22,29,56)

paths <- list(background_grid = "data/background_grid_22-29-56_res1000m_fit_r2000.geojson",
              spatialData = "data/spatial_data_22-29-56.RData",
              collisionData = "data/CollisionData/synthese_observations_2025-03-24T10_12_28.148Z.geojson")

load(paths$spatialData)

# rm(forests, lakes, landUse, rivers, coastLine)

studyArea <- filter(studyArea, code %in% dpts)
roads <- reduce(roads[as.character(dpts)], rbind)
rivers <- reduce(rivers[as.character(dpts)], rbind)
lakes <- reduce(lakes[as.character(dpts)], rbind)
zhp <- reduce(zhp[as.character(dpts)], rbind)
forests <- reduce(forests[as.character(dpts)], rbind)
landUse$Urban <- reduce(landUse$Urban[as.character(dpts)], rbind)
landUse$Cultures <- reduce(landUse$Cultures[as.character(dpts)], rbind)

background_grid <- read_sf(paths$background_grid) %>%
  rename(otterDensity = otterIntensity,
         altitude = z) %>%
  mutate(`river intersection` = factor(crossing)) %>% 
  rowwise %>%
  mutate(roadImportance = which(c(roadLen.1, roadLen.2, roadLen.3, roadLen.4, roadLen.5) != 0)[1])

otter_coll <- read_sf(paths$collisionData) %>%
  filter(year(date_debut) > 2000) %>%
  st_transform(crs = 2154) %>%
  st_intersection(st_union(background_grid))

## Fig1. Raw data + roads + rivers --------------------------------------------- 

ggplot(background_grid) + 
  geom_sf(aes(fill = otterDensity), col = NA, show.legend = FALSE)+
  geom_sf(data = st_union(filter(rivers, numero_d_ordre > 3)),linewidth = 1, color = "blue") +
  geom_sf(data = st_union(filter(rivers, numero_d_ordre == 3)),linewidth = 0.5, color = "blue") +
  geom_sf(data = st_union(filter(roads, importance %in% 1:2)),linewidth = 1, color = "black") +
  geom_sf(data = st_union(filter(roads, importance == 3)), color ="black") +
  geom_sf(data = otter_coll, col = "yellow") +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange")+
  simsim_theme_map(mode = "dark")

ggsave("out/rawDatDark.png",  width = 29.7, height = 21, unit = "cm")

## Fig2. Collision risk map ---------------------------------------------------- 
                   
ggplot(background_grid) +
  geom_sf(aes(fill = logMu), col = NA) +
  paletteer::scale_fill_paletteer_c("grDevices::Lajolla", name = "log(\u00B5)")+
  simsim_theme_map()
  
ggsave("out/riskMap.png",  width = 29.7, height = 21, unit = "cm")

## Fig3. Effect of covariates --------------------------------------------------

vars_w_cross <- c("riverLen","maxStra")
 
plotlist <- map2(vars_w_cross,  c("gam", "lm"), 
                 ~ simsim_scatter(data = background_grid, xVar = .x, yVar = "logMu",
                                  groupVar = "river intersection", method = .y)) %>%
  map2(vars_w_cross, ~ simsim_title(.x, .y, mycol = "black")) %>%
  map(~.x + ylim(min(background_grid$logMu), max(background_grid$logMu)) + ylab("log(\u00B5)"))

vars_wo_cross <- c("roadLen", "otterDensity", "distUrban", "altitude", "riverSin","waterBodyArea",  
                   "distLake", "roadSideVegetation", "roadImportance","pt.interactions")

names(plotlist) <- vars_w_cross

plotlist_2 <- map2(vars_wo_cross,  c("gam", "gam", "gam", "gam", "gam", "gam","gam", "gam", "lm", "lm"), 
                 ~ simsim_scatter(data = background_grid, xVar = .x, yVar = "logMu",
                                  method = .y, palette_points = "black")) %>%
  map2(vars_wo_cross, ~ simsim_title(.x, .y, mycol = "black")) %>%
  map(~.x + ylim(min(background_grid$logMu), max(background_grid$logMu)) + ylab("log(\u00B5)"))

names(plotlist_2) <- vars_wo_cross

plots <- append(plotlist, plotlist_2)
plots$riverIntersection <- ggplot(background_grid) + 
  geom_boxplot(aes(x = `river intersection`, y = logMu, fill = `river intersection`), show.legend = F)+
  geom_point(aes(x = `river intersection`, y = logMu, col = `river intersection`), fill = NA, pch = 21, show.legend = F)+
  scale_fill_manual(values = simsim_palette_d_1()) + 
  scale_color_manual(values = simsim_palette_d_1(2, +0.25)) + 
  ylab("log(\u00B5)") + xlab("") +
  simsim_theme_scatter()

plots$riverIntersection <- simsim_title(plots$riverIntersection, "road river intersection", mycol = "black")

(plots$roadLen | plots$roadImportance + strip_y | plots$distUrban + strip_y)

ggsave("out/traffic_covs.png",  width = 29.7, height = 7.5, unit = "cm")

(plots$riverIntersection | plots$riverLen + strip_y | plots$maxStra +  strip_y )/
  (plots$waterBodyArea | plots$distLake + strip_y | plots$roadSideVegetation +  strip_y )/
  (plots$altitude | plots$otterDensity + strip_y | plots$pt.interactions +  strip_y )/
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("out/habitat_covs.png",  width = 29.7, height = 21, unit = "cm")

## Fig5. Cumulated collision risk ---------------------------------------------- 

cumulated_risk <- st_intersection(background_grid, otter_coll) %>% 
  st_drop_geometry %>%
  group_by(id) %>%
  summarize(ncoll = n()) %>% 
  right_join(background_grid) %>%
  mutate(ncoll = replace_na(ncoll, 0))%>%
  arrange(desc(mu)) %>% 
  mutate(cum_colls = cumsum(ncoll), n_px = 1:n())

cumulated_risk_bygroup <- cumulated_risk %>%
  group_by(`river intersection`) %>%
  arrange(desc(mu)) %>% 
  mutate(cum_colls = cumsum(ncoll), n_px = 1:n())

breaks <- cumulated_risk %>%
  mutate(bin = cut(cum_colls, breaks = floor(nrow(otter_coll)/4*0:4), right = FALSE)) %>% 
  group_by(bin) %>% 
  filter(row_number() == 1) %>% 
  select(n_px, cum_colls)

ggplot(cumulated_risk) +
  geom_segment(data = breaks, aes(x=n_px/nrow(background_grid) * 100, y = cum_colls, yend = 0),
                                     linewidth = 1, linetype = "dashed")+
  geom_line(data = cumulated_risk_bygroup, 
            aes(x = n_px/nrow(background_grid) * 100, y = cum_colls, color = `river intersection`), linewidth = 1) +
  geom_line(aes(x = n_px/nrow(background_grid) * 100, y = cum_colls), linewidth = 0.75) +
  geom_ribbon(aes(x = n_px/nrow(background_grid) * 100, ymin = 0, ymax = logMu * 71 + 250), alpha = 0.25) +
  scale_color_manual(values = ggsimsim::simsim_palette_d_1()) + 
  scale_fill_manual(values = ggsimsim::simsim_palette_d_1()) + 
  scale_x_continuous(breaks = round(c(0,5,7.5,10, breaks$n_px/nrow(background_grid) * 100),2),
                     limits = c(0,10)) +
  scale_y_continuous(sec.axis = sec_axis( trans= ~(.- 250)/71, name="log(\u00B5)") ,limits = c(0,333))+
  xlab("proportion of the study area (%)") + ylab("cumulated number of collisions") +
  simsim_theme_scatter()

ggsave("out/cummulated_collisions.png",  width = 29.7, height = 21, unit = "cm")

## Fig. 6 Tiny maps of pixels --------------------------------------------------

locs <- read_sf("data/SpatialData/pxls_to_display") %>% 
  st_transform(crs = 2154)

ids <- st_intersection(background_grid, locs)$id

ecdf(background_grid$mu)(background_grid$mu[background_grid$id%in%ids])

masks <- split(background_grid[background_grid$id %in% ids, ],
               background_grid$id[background_grid$id %in% ids])

otter_coll_ <- map(masks, st_intersection, otter_coll) %>%
  map(select) %>% 
  map(st_transform, crs = 3857)

rivers_ <- map(masks, st_intersection, rivers) %>%
  map(select) %>% 
  map(st_transform, crs = 3857)

lakes_ <- map(masks, st_intersection, lakes) %>%
  map(select) %>% 
  map(st_transform, crs = 3857)

bboxs <- map(masks, st_bbox)

satviews <- map(bboxs, read_osm, type = "esri-imagery", zoom = 17) %>%
  map(as.data.frame, xy = T) %>%
  map(pivot_wider, names_from = "band", values_from = "X")

plots <- pmap(list(satviews, otter_coll_, rivers_, lakes_), ~
               ggplot(..1) + 
                geom_raster(aes(x=x,y=y, fill = rgb(red, green, blue, maxColorValue = 255))) + 
                scale_fill_identity()+ 
                geom_sf(data = ..3, color = "blue", linewidth = .5)+
                geom_sf(data = ..4, fill = "blue", color = NA,  alpha = 0.25)+
                geom_sf(data = ..2, color = "red", size = 2)+
                coord_sf() + 
                theme_void())

pxls_locs <- ggplot() +
  geom_sf(data = background_grid, col = NA, aes(fill = logMu), show.legend = FALSE) +
  geom_text(data = st_coordinates(st_centroid(reduce(masks, rbind))),
            aes(x = X, y = Y, label = c("A", "B", "C", "D", "E", "F", "G", "H")),
            col = "black") +
  paletteer::scale_fill_paletteer_c("grDevices::Lajolla", name = "log(\u00B5)")+
  theme_void()

ggpubr::ggarrange(plotlist = c(list(pxls_locs), plots),
                  labels = c("","A", "B", "C", "D", "E", "F", "G", "H"))

ggsave("out/pixelZoom.png",  width = 29.7, height = 29.7, unit = "cm")

# args <- map(masks, ~ map(list(rivers = rivers,
#                               roads = roads,
#                               waterbodies = lakes, vegetation = forests, urban = landUse$Urban),
#                          st_intersection, .x))
# 
# plots <- map2(masks, args, ~ do.call(ggsimsim::simsim_map, c(list(.x), .y))) %>%
#   map(~.x + geom_sf(data = otter_coll,  col = "red", size = 2.5))
# 
# 

## Raw covariates values -------------------------------------------------------

background_grid %>%
  mutate(across(c(roadLen, maxStra, treeCover), ~ c(scale(.x)))) %>%
  pivot_longer(cols = c("roadLen", "maxStra", "treeCover"),
               names_to = "param", values_to = "value") %>% 
  ggplot() + 
  geom_sf(aes(fill = value), col = NA, show.legend = FALSE) + 
  facet_wrap(~param) + 
  hrbrthemes::theme_modern_rc() + 
  theme(strip.text = element_text(size = 10, face = "bold", color = "white"))

ggsave("out/rawCovDark.png",  width = 16, height = 9, unit = "cm")

## Seasonality -----------------------------------------------------------------

ggplot(otter_coll) +
  geom_histogram(aes(x = factor(month(date_debut))), stat = "Count") + 
  xlab("Month") +
  theme_bw()

ggsave("out/collMonth.png",  width = 16, height = 9, unit = "cm")

