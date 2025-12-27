library(tidyverse)
library(sf)
library(ppmlasso)

rm(list = ls())

source("fun/fix_ppmlasso.R")

dpts <- c(22,29,56)
res <- 1000

tauLake <- 1500
tauCoast <- 1500

paths <- list(background_grid = "data/background_grid_22-29-56_res1000m_fit_r2000.geojson",
              collisions = "data/CollisionData/synthese_observations_2025-03-24T10_12_28.148Z.geojson")

form <- ~ otterIntensity + roadLen + roadLen.1 + roadLen.2 + roadLen.3 + roadLen.4 + roadLen.5 + 
  + crossing * (z + riverLen + riverSin + maxStra + lakeProx + treeCover + roadSideVegetation + riverSideVegetation + waterBodyArea +
                  Cultures + coastProx + distUrban) 

### 1. LOAD DATA ---------------------------------------------------------------

background_grid <- read_sf(paths$background_grid) %>%
  mutate_all(~replace_na(.x, 0)) %>%
  mutate(lakeProx = exp(-distLake/tauLake),
         coastProx = exp(-distCoast/tauCoast)) 

collisions <- read_sf(paths$collisions) %>%
  st_transform(crs = 2154) %>%
  filter(year(date_debut) > 2000) %>%
  st_intersection(st_union(background_grid))

### 2. PERFORM CROSS VALIDATION

Niters <- 50
nBlocks <- 8
muCV <- matrix(0, nrow(background_grid), Niters)

for(i in 1:Niters){
  ### a. GET VALIDATION BLOCKS ---------------------------------------------------
  
  bg_i <- st_make_grid(st_bbox(background_grid), cellsize = c(5000, 5000)) %>%
    st_as_sf(crs = 2154) %>%
    mutate(CVblock = sample(1:nBlocks, n(), replace = TRUE)) %>%
    st_intersection(st_centroid(background_grid)) %>%
    st_drop_geometry %>%
    select(id, CVblock) %>%
    left_join(background_grid) %>%
    st_as_sf(crs = 2154)
  
  ### b. ARRANGE DATA ------------------------------------------------------------
  
  bgxy <- cbind(bg_i, st_coordinates(st_centroid(bg_i))) %>%
    st_drop_geometry()
  
  spxy <- collisions %>%
    st_intersection(select(bg_i, CVblock), .) %>% 
    cbind(., st_coordinates(.)) %>%
    st_drop_geometry() %>%
    select(CVblock, X, Y)
  
  ### c. FIT AND PREDICT ---------------------------------------------------------
  
  bgxy_allbutCV <- map(unique(bg_i$CVblock),
                       ~ filter(bgxy, CVblock != .x))
  bgxy_onlyCV <- map(unique(bg_i$CVblock),
                     ~ filter(bgxy, CVblock == .x))
  
  spxy_allbutCV <- map(unique(bg_i$CVblock),
                       ~ filter(spxy, CVblock != .x))
  spxy_onlyCV <- map(unique(bg_i$CVblock),
                     ~ filter(spxy, CVblock == .x))
  
  ppmdata <- map2(
    spxy_allbutCV,
    bgxy_allbutCV,
    ppmdat,
    coord = c("X", "Y"),
    sp.scale = 1000,
    writefile = FALSE
  )
  ppmfit <- map(
    ppmdata,
    ~ ppmlasso(
      data = .x,
      formula = form,
      sp.scale = res,
      coord = c("X", "Y"),
      r = 2000,
      family = "area.inter",
      n.fits = 1,
      lamb = 5.77,
      standardise = TRUE,
      writefile = FALSE
    )
  )
  
  ixs <- map(unique(bg_i$CVblock), ~ bgxy$pt.interactions[bgxy$CVblock == .x])
  
  ppmpredict <- pmap(list(
    object = ppmfit,
    newdata = bgxy_onlyCV,
    interactions = ixs
  ),
  predict.ppmlasso)
  
  muCV[, i] <- map2(bgxy_onlyCV, ppmpredict, ~ mutate(.x, muCV = .y[, 1] * 1e+06)) %>%
    reduce(rbind) %>%
    select(id, muCV) %>%
    left_join(bg_i, .) %>%
    .$muCV %>%
    c %>%
    replace_na(0)
}

### 5. PLOT --------------------------------------------------------------------

apply(muCV, 2, cor, background_grid$mu) %>% t.test

# 95 percent confidence interval:
#   0.9330619 0.9483371
# sample estimates:
#   mean of x 
# 0.9406995 

apply(log(muCV), 2, function(x,y){mean(abs((x - y)/mean(y)))}, background_grid$logMu) %>% t.test

# 95 percent confidence interval:
#   0.01285255 0.01334518
# sample estimates:
#   mean of x 
# 0.01309887 

background_grid$avg_muCV <- apply(muCV, 1, mean)

ggplot(background_grid) +
  geom_sf(aes(fill = log(avg_muCV/mu)), col = NA)+
  scale_fill_gradient2(low = "red", high = "green", mid = "white", midpoint = 0,
                       name = "log(\u00B5CV/\u00B5)") +
  theme(panel.grid.major = element_line(color = "lightgrey"),
        axis.ticks = element_blank())

ggsave("out/CVmap.png",  width = 29.7, height = 21, unit = "cm")

### 6. BOYCE INDEX -------------------------------------------------------------

BI <- ecospat::ecospat.boyce(background_grid$mu, st_intersection(collisions, background_grid)$mu)
