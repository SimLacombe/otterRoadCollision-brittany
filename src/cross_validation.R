library(tidyverse)

library(ppmlasso)
library(splines)
library(pROC)

library(sf)
library(terra)
library(spatstat)

rm(list = ls())

source("fun/fix_ppmlasso.R")
source("fun/myUtilities.R")

paths <- list(bg = paste0("data/bg250m.rds"),
              studyArea = "data/SpatialData/contour-des-departements.geojson",
              collisions = "data/CollisionData/synthese_observations_2025-03-24T10_12_28.148Z.geojson")

### 1. LOAD DATA ---------------------------------------------------------------

studyArea <- read_sf(paths$studyArea) %>%
  filter(code %in% c(22,29,56)) %>%
  select(geometry, code) %>%
  st_transform(2154)

bg <- readRDS(paths$bg) %>%
  filter(roadImportance > 0) %>%
  mutate(maxStra = ifelse(crossings3, straCross, maxStra),
         maxStra = ifelse(maxStra == 8, 7, maxStra),
         roadImportance = ifelse(crossings3, roadCross, roadImportance)) %>% 
  mutate(id = 1:n())

collisions <- read_sf(paths$collisions) %>%
  filter(year(date_debut) > 2000) %>%
  st_transform(crs = 2154) %>%
  st_intersection(st_union(studyArea)) %>%
  distinct(geometry)

spxy <- collisions %>%
  st_coordinates %>%
  as.data.frame

myppmdata <- ppmdat(sp.xy = spxy, back.xy = bg, 
                    coord = c("X", "Y"), sp.scale = 500, writefile = FALSE)

form <- ~ (rivers == 0) + roadImportance + wetlands  +  crossings3 +
  roads + otterDensity + altitude + cultures + treeCover + rivers + maxStra + 
  exp(-distLake/1000) + exp(-distCoast/10000) + exp(-distUrban/500)

fit <- ppmlasso(form, data = myppmdata, n.fits = 20, sp.scale = 500, coord=c("X","Y"),
                standardise = TRUE, writefile = FALSE)

mu <- predict.ppmlasso(fit, newdata = bg)* 1e+06

### 2. PERFORM CROSS VALIDATION ------------------------------------------------

blocks <- list(as.numeric(as.factor(cut(bg$X, breaks = 8))),
               as.numeric(as.factor(cut(bg$Y, breaks = 8))),
               as.numeric(as.factor(cut(bg$X, breaks = 5))),
               as.numeric(as.factor(cut(bg$Y, breaks = 5))),
               as.numeric(as.factor(cut(c(bg$Y-min(bg$Y))**2 + c(bg$X-min(bg$X))**2, breaks = 5))),
               as.numeric(as.factor(cut(c(bg$Y-max(bg$Y))**2 + c(bg$X-max(bg$X))**2, breaks = 5))),
               as.numeric(as.factor(paste(cut(bg$X, breaks = 4), cut(bg$Y, breaks = 2)))),
               as.numeric(as.factor(paste(cut(bg$X, breaks = 2), cut(bg$Y, breaks = 4)))),
               as.numeric(as.factor(paste(cut(bg$X, breaks = 3), cut(bg$Y, breaks = 3)))))

muCV <- matrix(0, nrow(bg), length(blocks))

bg_sp <- st_as_sf(bg, coords = c("X", "Y"), crs = 2154) %>%
  st_buffer(125, endCapStyle = "SQUARE")

for(i in 1:length(blocks)){
  cat(i, "\n")
  ### a. GET VALIDATION BLOCKS ---------------------------------------------------
  
  ### b. ARRANGE DATA ------------------------------------------------------------
  bg_i <- mutate(bg, block = blocks[[i]])
  bg_sp_i <- mutate(bg_sp, block = blocks[[i]])
  
  coll_i <- collisions %>%
    st_intersection(select(bg_sp_i, block), .) %>% 
    cbind(., st_coordinates(.)) %>%
    st_drop_geometry() %>%
    select(block, X, Y)
  
  ### c. FIT AND PREDICT ---------------------------------------------------------
  
  bg_fit <- map(unique(bg_i$block),
                  ~ filter(bg_i, block != .x))
  bg_pred <- map(unique(bg_i$block),
                     ~ filter(bg_i, block == .x))
  
  coll_fit <- map(unique(bg_i$block),
                       ~ filter(coll_i, block != .x))
  coll_pred <- map(unique(bg_i$block),
                     ~ filter(coll_i, block == .x))
  
  myppmdata_i <- map2(
    coll_fit,
    bg_fit,
    ppmdat,
    coord = c("X", "Y"),
    sp.scale = 1000,
    writefile = FALSE
  )
  fit_i <- map(
    myppmdata_i,
    ~ ppmlasso(
      data = .x,
      formula = form,
      sp.scale = res,
      coord = c("X", "Y"),
      lamb = fit$lambda,
      standardise = TRUE,
      writefile = FALSE
    )
  )
  
  mus_i <- pmap(list(
    object = fit_i,
    newdata = bg_pred
    ),
    predict.ppmlasso)
  
  muCV[, i] <- map2(bg_pred,
                    mus_i, ~ mutate(.x, muCV = .y[, 1] * 1e+06)) %>%
    reduce(rbind) %>%
    select(id, muCV) %>%
    left_join(bg_i, .) %>%
    .$muCV %>%
    c %>%
    replace_na(0)
}

### 5. SUMMARIZE ---------------------------------------------------------------

jaccard <- function(a, b){
  u <- intersect(a,b)
  length(u) / (length(a) + length(b) - length(u))
}

apply(muCV, 2, function(x,y){
  jaccard(which(x > quantile(x, 0.95)),
          which(y > quantile(y, 0.95)))},
  mu) %>% t.test

# 95 percent confidence interval:
#   0.7820319 0.8866666
# sample estimates:
#   mean of x 
# 0.8343492

apply(muCV, 2, cor, mu, method = "spear") %>% t.test

# 95 percent confidence interval:
#   0.9759670  0.9893847
# sample estimates:
#   mean of x 
# 0.9826758  

apply(log(muCV), 2, function(x,y){mean(abs((x - y)/mean(y)))}, log(mu)) %>% t.test

# 95 percent confidence interval:
#   0.03181962 0.04160178
# sample estimates:
#   mean of x 
# 0.0367107 

coll <- 1 - (st_intersects(bg_sp, collisions) %>% sapply(is_empty))

apply(muCV, 2, function(x, y){
  rocdat <- data.frame(mu = x, coll = y)
  as.numeric(roc_(data = rocdat, "coll", "mu")$auc)
}, coll) %>% t.test

# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   0.9207811 0.9258190
# sample estimates:
#   mean of x 
# 0.9233 
