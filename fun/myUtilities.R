### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ### 

## General ---------------------------------------------------------------------

my_join <- function(x, y, crs, ...){
  
  left_join(y, x) %>%
    mutate(across(c(...), ~ replace_na(.x, 0))) %>%
    st_as_sf(crs = crs)
}

which_max_by <- function(vals, ids){
  split(vals, ids) %>% 
    map(sum) %>%
    reduce(c)%>%
    which.max %>%
    ids[.]
}

## spatial operations ----------------------------------------------------------

match_to_closest <- function(points, grid, id){
  dmat <- st_distance(points, grid)
  matched <- id[apply(dmat, 1, which.min)]
  
  matched
}

st_line_sinuosity <- function(x){
  dist <- st_length(x)
  straight_dist <- st_line_sample(x, sample = c(0,1)) %>%
    st_cast("LINESTRING") %>%
    st_length
  as.numeric(dist/straight_dist)
}

summarize_shape <- function(grid, shape, myFUN, id = id, name_to, name_sin = NULL,
                            byClass = FALSE, class = NA, sinuosity = FALSE, ...){

  if(sinuosity & identical(myFUN, st_length)) {   
    out <- st_intersection(grid, shape) %>%
      st_cast("MULTILINESTRING") %>%
      st_cast("LINESTRING") %>% 
      mutate(tmp = as.numeric(myFUN(x)),
             sin = st_line_sinuosity(x)) 
  } else {
    out <- st_intersection(grid, shape) %>%
      mutate(tmp = as.numeric(myFUN(x))) 
    }
  
  out <- st_drop_geometry(out)
  
  if(byClass){ out <- group_by(out, {{id}}, {{ class }})}
  else { out <- group_by(out, {{id}})}
  
  if(sinuosity & identical(myFUN, st_length)) {
  out <-  out %>%
    summarize({{ name_to }} := sum(tmp), 
              {{ name_sin }} := weighted.mean(sin, tmp),
              ...) %>%
    my_join(grid, crs = 2154, {{ name_to }}, {{name_sin}})
  } else {
    out <-  out %>%
      summarize({{ name_to }} := sum(tmp),
                ...) %>%
      my_join(grid, crs = 2154, {{ name_to }})
    }
  out 
}

## To subsample otter data -----------------------------------------------------

subsample <- function(geometry, thr){
  
  distances <- st_distance(geometry) %>%
    apply(1, as.numeric)
  
  diag(distances) <- +Inf
  
  ids_toprocess <- seq_along(geometry)
  ids_kept <- integer()
  
  neighboors <- which(distances < thr, arr.ind = TRUE) %>% 
    as.data.frame
  
  while(length(ids_toprocess) > 1){
    id <- sample(ids_toprocess, 1)
    ids_kept <- append(ids_kept, id)
    ids_toprocess <- ids_toprocess[!ids_toprocess ==id]
    
    ids_torm <- neighboors$row[neighboors$col == id]
    ids_toprocess <- ids_toprocess[!ids_toprocess %in% ids_torm]
  }
  append(ids_kept, ids_toprocess)
}
