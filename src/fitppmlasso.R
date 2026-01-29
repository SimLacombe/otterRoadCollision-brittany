library(tidyverse)
library(sf)
library(ppmlasso)

rm(list = ls())

dpts <- c(22,29,56)
res <- 1000
r <- 2000

tauLake <- 2500
tauCoast <- 2500

paths <- list(background_grid = paste0("data/background_grid_", paste(dpts, collapse = "-"), "_res", res, "m.geojson"),
              collisions = "data/CollisionData/synthese_observations_2025-03-24T10_12_28.148Z.geojson",
              out = paste0("data/background_grid_", paste(dpts, collapse = "-"), "_res", res, "m_fit_r", r, ".geojson"))

form <- ~ otterIntensity + roadLen + roadLen.1 + roadLen.2 + roadLen.3 + roadLen.4 + roadLen.5 + 
  +  distUrban + crossing * (z + riverLen + riverSin + maxStra + lakeProx + treeCover + roadSideVegetation + riverSideVegetation + waterBodyArea +
                 Cultures + coastProx) 

### 1. LOAD DATA ---------------------------------------------------------------

background_grid <- read_sf(paths$background_grid) %>%
  mutate(across(c("roadLen.1", "roadLen.2", "roadLen.3", "roadLen.4", "roadLen.5"), ~ .x/roadLen)) %>%
  mutate_all(~replace_na(.x, 0))

collisions <- read_sf(paths$collisions) %>%
  filter(year(date_debut) > 2000) %>%
  st_transform(crs = 2154) %>%
  st_intersection(st_union(background_grid)) 

### 2. Prepare data for point process modelling --------------------------------

bgxy <- cbind(background_grid, st_coordinates(st_centroid(background_grid))) %>%
  st_drop_geometry() %>%
  mutate(lakeProx = exp(-distLake/tauLake),
         coastProx = exp(-distCoast/tauCoast)) 

spxy <- collisions %>%
  st_coordinates %>%
  as.data.frame

myppmdata <- ppmdat(sp.xy = spxy, back.xy = bgxy, coord = c("X", "Y"), sp.scale = res, writefile = FALSE)

### 3. Fit PPP -----------------------------------------------------------------

fit <- ppmlasso(form, data = myppmdata, n.fits = 50, sp.scale = res, coord=c("X","Y"),
                r = r, family = "area.inter", standardise = TRUE, writefile = FALSE)

background_grid$mu <- fit$mu[(nrow(spxy)+1):length(fit$mu)]*1e+06
background_grid$logMu <- log(background_grid$mu)
background_grid$psi <- 1 - exp(-background_grid$mu)

background_grid$pt.interactions <- fit$pt.interactions[(nrow(spxy)+1):length(fit$pt.interactions)]

st_write(background_grid, paths$out, driver = "GeoJSON")
# saveRDS(fit, "out/fit_arinterlasso_22-29-56_1000m_r2000.rds")

### 4. DIAGNOSTICS  ------------------------------------------------------------

## a. Effect of the scale parameter on ll

tauLake_vect <- c(500,1000,1500,2500,5000)
tauLake_diagnose <- map(tauLake_vect,
                         ~ ppmlasso(form, sp.xy = spxy, n.fits=1,env.grid = mutate(bgxy, lakeProx = exp(-distLake/.x)),
                                    sp.scale = res,coord=c("X","Y"), writefile = FALSE, lamb = 5.77)) %>%
  map(~.x$likelihoods) %>%
  reduce(c)

diag_lake <- data.frame(tau_lake = tauLake_vect, lik = na.omit(tauLake_diagnose)) %>%
  ggplot()+
  geom_point(aes(x = tau_lake, y = lik))+
  geom_line(aes(x = tau_lake, y = lik))+
  ggthemes::theme_clean()+
  theme(plot.background = element_blank())

tauCoast_vect <- c(500,1000,1500,2500,5000)
tauCoast_diagnose <- map(tauCoast_vect,
                         ~ ppmlasso(form, sp.xy = spxy, n.fits=1,env.grid = mutate(bgxy, coastProx = exp(-distCoast/.x)),
                                    sp.scale = res,coord=c("X","Y"), writefile = FALSE, lamb = 5.77)) %>%
  map(~.x$likelihoods) %>%
  reduce(c)

diag_coast <- data.frame(tau_coast = tauCoast_vect, lik = na.omit(tauCoast_diagnose)) %>%
  ggplot()+
  geom_point(aes(x = tau_coast, y = lik))+
  geom_line(aes(x = tau_coast, y = lik))+
  ggthemes::theme_clean()+
  theme(plot.background = element_blank())

( diag_lake | diag_coast )
# ggsave("out/diagScaleProx.png",  width = 29.7, height = 7.5, unit = "cm")

## b. Effect of the resolution 

findRes(c(1000,2000, 4000, 6000, 10000),
        formula = form, sp.xy = spxy,env.grid = bgxy, coord = c("X", "Y"))

# ggsave("out/diagRes.png",  width = 16, height = 9, unit = "cm")

## c. Check remaining clustering for various values of r

fit1 <- ppmlasso(form, data = myppmdata, n.fits = 50, sp.scale = res, coord=c("X","Y"),
                standardise = TRUE, writefile = FALSE)
background_grid$mu1 <- fit1$mu[(nrow(spxy)+1):length(fit1$mu)]

fit2 <- ppmlasso(form, data = myppmdata, n.fits = 50, sp.scale = res, coord=c("X","Y"),
                r = 1000, family = "area.inter", standardise = TRUE, writefile = FALSE)
background_grid$mu2 <- fit2$mu[(nrow(spxy)+1):length(fit2$mu)]

fit3 <- ppmlasso(form, data = myppmdata, n.fits = 50, sp.scale = res, coord=c("X","Y"),
                r = 2000, family = "area.inter", standardise = TRUE, writefile = FALSE)
background_grid$mu3 <- fit3$mu[(nrow(spxy)+1):length(fit3$mu)]

fit4 <- ppmlasso(form, data = myppmdata, n.fits = 50, sp.scale = res, coord=c("X","Y"),
                r = 3000, family = "area.inter", standardise = TRUE, writefile = FALSE)
background_grid$mu4 <- fit4$mu[(nrow(spxy)+1):length(fit4$mu)]

fit5 <- ppmlasso(form, data = myppmdata, n.fits = 50, sp.scale = res, coord=c("X","Y"),
                 r = 5000, family = "area.inter", standardise = TRUE, writefile = FALSE)
background_grid$mu5 <- fit5$mu[(nrow(spxy)+1):length(fit5$mu)]


background_grid <- background_grid %>%
  mutate(across(all_of(c("mu1", "mu2", "mu3", "mu4", "mu5")), ~ replace_na(.x, 0)))

r <- terra::rast(xmin = st_bbox(background_grid)[1], xmax = st_bbox(background_grid)[3],
                 ymin = st_bbox(background_grid)[2], ymax = st_bbox(background_grid)[4],
                 resolution = c(1000, 1000))

rs <-  map(c("mu1", "mu2", "mu3", "mu4", "mu5"), ~ terra::rasterize(background_grid, r, field = .x))

extract_mu <- function(i) {
  function(x, y) {
    terra::extract(rs[[i]], cbind(x, y))[,1]
  }
}

true_ppp <- ppp(spxy$X, spxy$Y, window = as.owin(st_union(background_grid)))

env_inhom <- map(1:5, ~ envelope(true_ppp, Kinhom, nsim = 50,
                                 simulate = expression(rpoispp(lambda = extract_mu(.x),
                                                               lmax = max(background_grid$mu1),
                                                               win = as.owin(st_union(background_grid))))))

env_inhom <- map2(env_inhom, c("NA", "1000", "2000", "3000", "5000"),
                  ~ mutate(as.data.frame(.x), interaction.radius = .y)) %>%
  reduce(rbind)

ggplot(env_inhom) +
  geom_line(aes(x = r, y = sqrt(obs) - sqrt(r)), col = "red", linewidth = .75)+
  geom_line(aes(x = r, y = sqrt(mmean) - sqrt(r)), linetype = "dashed", linewidth = .75)+
  geom_ribbon(aes(x = r, ymin = sqrt(lo) - sqrt(r), ymax = sqrt(hi)- sqrt(r)), alpha = .5)+
  xlim(c(0,10000))+ylim(c(0,22000))+
  ylab(expression(sqrt(hat(K)) - sqrt(r)))+
  theme_bw() +
  facet_wrap(~ interaction.radius, nrow = 3, labeller = labeller(interaction.radius = label_both)) +
  theme_bw() +
  theme(plot.background = element_blank())

# ggsave("out/Kenveloppe.png",  width = 29.7, height = 21, unit = "cm")

