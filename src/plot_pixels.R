library(tidyverse)
library(sf)
library(patchwork)
library(tmaptools)

rm(list = ls())

source("fun/myUtilities.R")
source("fun/plot_functions.R")

dpts <- c(22,29,56)

paths <- list(bg = "out/out250m.rds",
              BDTopo = map(dpts, ~ paste0("data/SpatialData/BDT_3-4_GPKG_LAMB93-ED2024-09-15/BDT_3-4_GPKG_LAMB93_D0",
                                          .x, "-ED2024-09-15.gpkg")),
              collisionData = "data/CollisionData/synthese_observations_2025-03-24T10_12_28.148Z.geojson")

rivers <-  map(paths$BDTopo, read_sf, layer = "TRONCON_HYDROGRAPHIQUE") %>%
  reduce(rbind) %>%
  filter(nature %in% c("Ecoulement naturel", "Conduit buse", "Canal", "Marais")) %>%
  mutate(numero_d_ordre = as.numeric(numero_d_ordre))

waterBodies <- map(paths$BDTopo, read_sf, layer = "PLAN_D_EAU") %>%
  reduce(rbind)

lakes <- filter(waterBodies,
                nature %in% c("Lac", "Retenue", "Retenue-barrage", "Réservoir-bassin"),
                as.numeric(st_area(waterBodies)) > 50000)

bg <- readRDS(paths$bg) %>%
  filter(roadImportance > 0) %>%
  mutate(px = 1:n())

otter_coll <- read_sf(paths$collisionData) %>%
  filter(year(date_debut) > 2000) %>%
  st_transform(crs = 2154)

locs <- read_sf("data/SpatialData/pxls_to_display") %>% 
  st_transform(crs = 2154) %>%
  mutate(plotID = 1:n())

bg_sf <- st_as_sf(bg, coords = c("X", "Y"), crs = 2154) 

bg_plot <- st_intersection(bg_sf,
                           st_buffer(locs, 300, endCapStyle = "SQUARE"))

pxls_toplot <- split(st_buffer(bg_plot, 250, endCapStyle = "SQUARE"),
                     bg_plot$plotID) 

masks <- map(pxls_toplot, st_union)

otter_coll_ <- map(masks, st_intersection, otter_coll) %>%
  map(st_transform, crs = 3857)

rivers_ <- map(masks, st_intersection, rivers) %>%
  map(st_transform, crs = 3857)

lakes_ <- map(masks, st_intersection, lakes) %>%
  map(st_transform, crs = 3857)

bboxs <- map(masks, st_bbox)

satviews <- map(bboxs, read_osm, type = "esri-imagery", zoom = 17) %>%
  map(as.data.frame, xy = T) %>%
  map(pivot_wider, names_from = "band", values_from = "X")

plots <- pmap(list(satviews, otter_coll_, rivers_, lakes_, pxls_toplot), ~
                ggplot(..1) + 
                geom_raster(aes(x=x,y=y, fill = rgb(red, green, blue, maxColorValue = 255))) + 
                scale_fill_identity()+ 
                geom_sf(data = ..3, color = "blue", linewidth = .5)+
                geom_sf(data = ..4, fill = "blue", color = NA,  alpha = 0.25)+
                geom_sf(data = ..2, color = "red", size = 2)+
                geom_sf(data = ..5, color = "black", fill = NA)+
                theme_void())

pxls_locs <- ggplot() +
  geom_tile(data = bg, col = NA, aes(x=X, y =Y, fill = log(mu)), show.legend = FALSE) +
  geom_text(data = st_coordinates(locs),
            aes(x = X, y = Y, label = c("A", "B", "C", "D", "E", "F", "G", "H")),
            col = "black") +
  paletteer::scale_fill_paletteer_c("grDevices::Lajolla", name = "log(\u00B5)")+
  theme_void() +
  theme(text = element_text(face = "bold"))+
  coord_fixed()

ggpubr::ggarrange(plotlist = c(list(pxls_locs), plots),
                  labels = c("","A", "B", "C", "D", "E", "F", "G", "H"))

ggsave("out/pixelZoom.png",  width = 29.7, height = 29.7, unit = "cm")
