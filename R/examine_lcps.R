######################################
######################################
#### examine_lcps.R

#### This code:
# 1) Evaluates the use of Euclidean distances as an approximation
# ... of shortest distances, as movement scales relevant to flapper skate
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
n_trial <- 20
if(run){
  out <-
    pbapply::pblapply(1:n_trial, cl = 5L, function(i){
    set.seed(i)
    out <- lcp_comp(surface = site_bathy,
                    barrier = site_coast_barrier,
                    mobility = mob,
                    n = 10000, n_max = 50000,
                    graph = site_bathy_graph)
    saveRDS(out, paste0("./data/spatial/lcp_comp_", i, ".rds"))
    return(out)
  })
  beepr::beep(10L)
} else {
  out <- lapply(1:n_trial, function(i){
    readRDS(paste0("./data/spatial/lcp_comp_", i, ".rds"))
  })
}

#### Process results
## Implement initial processing
out <-
  out %>%
  dplyr::bind_rows() %>%
  # Drop pairs for which we could not calculate LCPs
  dplyr::filter(!is.na(.data$dist_lcp)) %>%
  # Drop duplicate samples
  dplyr::mutate(key = paste0(cell_id0, cell_id1)) %>%
  dplyr::filter(!duplicated(key))
## Drop samples from within the barrier
# ... (resulting from a slight discrepancy between site_bathy and site_coast at the edges)
# Cell one
xy <- sp::SpatialPoints(out[, c("cell_x0", "cell_y0")],
                          proj4string = raster::crs(site_coast))
xy_over <- sp::over(xy, site_coast)
xy_bool <- !is.na(xy_over)
xy_sums <- rowSums(xy_bool, na.rm = TRUE)
xy_pos    <- which(xy_sums > 1)
xy_over[xy_pos, ]
check <- FALSE
if(check){
  for(i in 1:length(xy_pos)){
    # i <- 1
    xy_i <- xy[xy_pos[i], , drop = FALSE]
    xy_i <- sp::coordinates(xy_i)
    adj <- 100
    raster::plot(site_coast,
                 xlim = c(plyr::round_any(xy_i[1], adj, floor), plyr::round_any(xy_i[1], adj, ceiling)),
                 ylim = c(plyr::round_any(xy_i[2], adj, floor), plyr::round_any(xy_i[2], adj, ceiling)),
                 col = "lightgreen", border = "lightgreen",
                 main = i)
    axis(side = 1); axis(side = 2)
    points(xy[xy_pos[i], ], lwd = 5, col = "red")
    readline(prompt = "Press [Enter] to continue...")
  }
}
out <- out[-xy_pos, ]
# Cell two
xy <- sp::SpatialPoints(out[, c("cell_x1", "cell_y1")],
                        proj4string = raster::crs(site_coast))
xy_over <- sp::over(xy, site_coast)
xy_bool <- !is.na(xy_over)
xy_sums <- rowSums(xy_bool, na.rm = TRUE)
xy_pos  <- which(xy_sums > 1)
xy_over[xy_pos, ]
out <- out[-xy_pos, ]

#### Examine results
head(out)
table(out$barrier)
pos   <- which.min(out$dist_euclid[out$barrier == 1L])
trans <- out[out$barrier == 1L, ][pos, ]
raster::plot(site_coast,
             xlim = range(trans[, c("cell_x0", "cell_x1")]),
             ylim = range(trans[, c("cell_y0", "cell_y1")]))
raster::plot(site_bathy, add = TRUE)
raster::lines(site_coast)
lines(out[out$barrier == 1L, c("cell_x0", "cell_x1")][pos, ],
      out[out$barrier == 1L, c("cell_y0", "cell_y1")][pos, ],
      col = "red", lwd = 10)
head(sort(out$dist_euclid[out$barrier == 1L]))
quantile(out$dist_euclid[out$barrier == 1L], 0.01)
as.character(quantile(out$dist_euclid[out$barrier == 1L], 0.05))

#### Visualise sampled locations on map
map <- FALSE
if(map){
  raster::plot(site_bathy)
  raster::lines(site_coast)
  invisible(lapply(split(out, seq_len(nrow(out))), function(d){
    arrows(x0 = d$cell_x0, y0 = d$cell_y0,
           x1 = d$cell_x1, y1 = d$cell_y1, length = 0.01, lwd = 1)
  }))
}

#### Model
mod <- lm(dist_lcp ~ 0 + dist_euclid:barrier, data = out)
summary(mod)
coef(mod)
summary(mod)[["r.squared"]] * 100


#### Compare Euclidean and shortest distances

# Extract the data for the paths that do not versus do cross a barrier
out_barrier0 <- out[out$barrier == 0, ]
out_barrier1 <- out[out$barrier == 1, ]

# Set up plot
png("./fig/study_site/lcp_comp.png",
    height = 4, width = 8, units = "in", res = 600)
pp <- par(mfrow = c(1, 2), oma = c(2, 2, 1, 1), mar = c(2, 2, 2, 2))

# Results for paths that do not cross a barrier
prettyGraphics::pretty_plot(out_barrier0$dist_euclid, out_barrier0$dist_lcp,
                            xlab = "", ylab = "",
                            pch = ".",
                            col = scales::alpha("grey", 0.75))
x <- seq(0, 500)
nd <- data.frame(barrier = factor(0, levels = c(0, 1)), dist_euclid = x)
ci <- prettyGraphics::list_CIs(predict(mod, nd, se.fit = TRUE))
prettyGraphics::add_error_envelope(x, ci,
                                   add_fit = list(col = "royalblue", lwd = 3))
mtext(side = 3, "A", font = 2, adj = 0, cex = 1.2)

# Results for paths that cross a barrier
prettyGraphics::pretty_plot(out_barrier1$dist_euclid, out_barrier1$dist_lcp,
                            xlab = "", ylab = "",
                            pch = ".",
                            col = scales::alpha("grey", 0.75),
                            add_fit = list(col = "royalblue", lwd = 2))
nd$barrier <- factor(1, levels = c(0, 1))
ci <- prettyGraphics::list_CIs(predict(mod, nd, se.fit = TRUE))
prettyGraphics::add_error_envelope(x, ci,
                                   add_fit = list(col = "royalblue", lwd = 3))
mtext(side = 3, "B", font = 2, adj = 0, cex = 1.2)

mtext(side = 1, "Euclidean distance (m)", outer = TRUE, line = 0)
mtext(side = 2, "Shortest distance (m)", outer = TRUE, line = 0.75)
par(pp)
dev.off()

# Define predictive model based on this analysis
as.character(coef(mod)[1])
as.character(coef(mod)[2])
lcp_predict <- function(distance, barrier){
  1.06450396994437 * distance * (barrier == 0) +
    1.17288343596087 * distance * (barrier == 1)
}


#### End of code.
######################################
######################################
