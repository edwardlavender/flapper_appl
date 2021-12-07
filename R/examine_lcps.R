######################################
######################################
#### examine_lcps.R

#### This code:
# 1) Evaluates the use of Euclidean distances as an approximation
# ... of shortest distances, as movement scales relavent to flapper skate
# ... (given their assumed mobility), within the study area.

#### Steps preceding this code:
# 1) Definition of processed spatial fields and graph for study area
# ... via process_data_raw.R and examine_space_use.R

#### Results
# The simulations show that all Euclidean distances (for paths that do not cross a barrier) (n = 99800)
# ... are below a mobility threshold of 500 m up to a Euclidean distance of 226 m.
# For 95 % of the pairwise connections with a Euclidean distance of less than 481 m
# ... the shortest distances remained below the mobility threshold
# Overall, for 89 % of the pairwise connections less than 500 m apart
# ... the shortest distances remained below the mobility threshold.
# For the segments that crossed a barrier (n = 128), all segments
# ... are below a mobility threshold of 500 m up to a Euclidean distance of 345 m
# For 95 % of the pairwise connections that crossed a barrier
# ... with a Euclidean distance of less than 421 m,
# ... the shortest distances remained below the mobility threshold.
# ... Overall, for 55 % of the pairwise connections less than 500 m apart that crossed a barrier
# ... the shortest distances remained below the mobility threshold.


######################################
######################################
#### Set up

#### Wipe workspace & load flapper & global param
rm(list = ls())
library(flapper)
# source("./R/define_global_param.R")

#### Define data
site_bathy         <- raster::raster("./data/spatial/site_bathy.tif")
site_coast         <- readRDS("./data/spatial/site_coast.rds")
site_coast_barrier <- sf::st_as_sf(site_coast)
site_bathy_graph   <- readRDS("./data/movement/space_use/site_bathy_graph.rds")
mob                <- 500


######################################
######################################
#### Implement simulation

#### Implement flapper::lcp_comp()
# For n_max = 5,000   without parallelisation this takes ~1.12 mins
# For n_max = 10,000  without parallelisation this takes ~2.07 mins
# For n_max = 100,000 without parallelisation this takes ~13.41 minutes
run <- FALSE
if(run){
  set.seed(1)
  out <- lcp_comp(surface = site_bathy,
                  barrier = site_coast_barrier,
                  mobility = mob,
                  n = 10000, n_max = 100000,
                  graph = site_bathy_graph)
  saveRDS(out, "./data/spatial/lcp_comp.rds")
} else out <- readRDS("./data/spatial/lcp_comp.rds")


#### Process results
# Drop pairs for which we could not calculate LCPs
nrow(out)
out <- out %>% dplyr::filter(!is.na(.data$dist_lcp))
nrow(out)

#### Examine results
head(out)
table(out$barrier)

#### Visualise sampled locations on map
raster::plot(site_bathy)
raster::lines(site_coast)
invisible(lapply(split(out, seq_len(nrow(out))), function(d){
  arrows(x0 = d$cell_x0, y0 = d$cell_y0,
         x1 = d$cell_x1, y1 = d$cell_y1, length = 0.01, lwd = 1)
}))

#### Compare Euclidean and shortest distances
# Extract the data for the paths that do not versus do cross a barrier
out_barrier0 <- out[out$barrier == 0, ]
out_barrier1 <- out[out$barrier == 1, ]
# Set up plotting window
pp <- graphics::par(mfrow = c(1, 2), oma = c(3, 3, 3, 3), mar = c(4, 4, 4, 4))
# Results for paths that do not cross a barrier
# ... Visualisation
prettyGraphics::pretty_plot(out_barrier0$dist_euclid, out_barrier0$dist_lcp,
                            xlab = "Distance (Euclidean) [m]",
                            ylab = "Distance (shortest) [m]",
                            pch = ".")
graphics::abline(0, 1, col = "red")
graphics::abline(h = mob, col = "royalblue", lty = 3)
# ... Euclidean distance parameter at which mobility is exceeded
limit0 <- min(out_barrier0$dist_euclid[out_barrier0$dist_lcp > mob]); limit0
graphics::abline(v = limit0, col = "royalblue", lty = 3)
# Results for paths that cross a barrier
# ... Visualisation
prettyGraphics::pretty_plot(out_barrier1$dist_euclid, out_barrier1$dist_lcp,
                            xlab = "Distance (Euclidean) [m]",
                            ylab = "Distance (shortest) [m]",
                            pch = ".")
graphics::abline(0, 1, col = "red")
graphics::abline(h = mob, col = "royalblue", lty = 3)
# ... Euclidean distance parameter at which mobility is exceeded
# ... ... (This is lower than for paths that do not cross a barrier)
limit1 <- min(out_barrier1$dist_euclid[out_barrier1$dist_lcp > mob], na.rm = TRUE); limit1
graphics::abline(v = limit1, col = "royalblue", lty = 3)
graphics::par(pp)

#### Visualise how % connections with shortest distances less than mobility
# ... declines with Euclidean distance

## Derive quantiles
out_barrier_for_quants <- out_barrier1 # out_barrier0
quantiles <-
  data.frame(distance = seq(min(out_barrier_for_quants$dist_euclid), 500, length.out = 10000),
             prop = NA)
for(i in seq_len(nrow(quantiles))){
  # i = 1
  tmp <- out_barrier_for_quants[out_barrier_for_quants$dist_euclid <= quantiles$distance[i], ]
  prop <- (length(which(tmp$dist_lcp <= mob)))/nrow(tmp)
  quantiles$prop[i] <- prop
}

## Examine 95 percentiles
# ... 95 % of connections less than this Euclidean distance
# ... remained below mobility in terms of shortest distance:
q95 <- quantiles$distance[which.min(quantiles$prop - 0.95 > 0)]
q95

## Examine the tail of the distribution
tail(quantiles)
min(quantiles$prop)

## Visualise distribtion
axis_ls <- prettyGraphics::pretty_plot(quantiles$distance, quantiles$prop,
                            pretty_axis_args = list(control_digits = 2),
                            type = "l",
                            xlab = "Distance (Euclidean) [m]",
                            ylab = "Proportion",
                            lwd = 2)
xlim <- axis_ls[[1]]$lim
ylim <- axis_ls[[2]]$lim
lines(xlim, c(0.95, 0.95), col = "royalblue", lty = 3, lwd = 2)
lines(c(q95, q95), ylim, col = "royalblue", lty = 3, lwd = 2)


#### End of code.
######################################
######################################

