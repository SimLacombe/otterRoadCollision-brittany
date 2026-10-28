library(tidyverse)
library(sf)
library(terra)

rm(list = ls())

dpts <- c(`22` = 22, `29` =  29, `56` = 56)

paths <- list(
  studyArea = "data/SpatialData/contour-des-departements.geojson",
  BDTopo = map(dpts, ~ paste0("data/SpatialData/BDT_3-4_GPKG_LAMB93-ED2024-09-15/BDT_3-4_GPKG_LAMB93_D0", .x, "-ED2024-09-15.gpkg")),
  CLC = "data/SpatialData/U2012_CLC2006_V2020_20u1.gdb"
)

studyArea <- read_sf(paths$studyArea) %>%
  filter(code %in% dpts) %>%
  select(geometry, code) %>%
  st_transform(2154)

cultures <- read_sf(paths$CLC) %>%
  st_transform(crs = 2154) %>%
  st_intersection(studyArea) %>%
  select(Code_06) %>%
  filter(Code_06 %in% c("211", "212", "213", "221", "222", "223",
                        "231", "241", "242", "243", "244"))

treeCover <-  map(paths$BDTopo[as.character(dpts)], read_sf, layer = "ZONE_DE_VEGETATION") %>%
  reduce(rbind) %>%
  filter(grepl("Bois|Forêt", nature))

ext <- st_bbox(studyArea)

r <- rast(xmin = ext["xmin"], xmax = ext["xmax"],
          ymin = ext["ymin"], ymax = ext["ymax"],
          resolution = 10, crs = "EPSG:2154")

treeCover.r <- rasterize(treeCover, r, field = 1, background = 0)
cultures.r <- rasterize(cultures, r, field = 1, background = 0)
