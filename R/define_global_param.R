######################################
######################################
#### define_global_param.R

#### This code:
# 1) Defines global parameters used for applications of the flapper
# ... family of algorithms within this project.

#### Steps preceding this code:
# 1) NA


######################################
######################################
#### Set up

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(flapper)
library(prettyGraphics)


######################################
######################################
#### Spatial parameters

#### Thesis or manuscript colour scheme
use_man_scheme <- TRUE

#### Define projections
proj_wgs84 <- sp::CRS("+init=epsg:4326")
proj_utm <- sp::CRS("+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs")

#### Plotting param
## Coastline
if(use_man_scheme){
  add_coast <- list(x = readRDS("./data/spatial/site_coast.rds"),
                    col = "dimgrey",
                    border = "dimgrey")
} else {
  add_coast <- list(x = readRDS("./data/spatial/site_coast.rds"),
                    col = "#f1f4c7",
                    border = "#bdbf97")
}

## Bathymetry parameters
bathy_zlim <- c(0, 225)
bathy_col_param <- prettyGraphics::pretty_cols_brewer(bathy_zlim,
                                                      scheme = "Blues",
                                                      n_breaks = max(bathy_zlim))
## Map axes
paa <- list(side = 1:4,
            axis = list(labels = FALSE, tck = -0.01)
            )
## North arrow and scale bar
add_map_elements <- function(){
  arrow_y0 <- 6247000 # 6247700
  arrows(x0 = 714000, y0 = arrow_y0, y1 = arrow_y0 + 2000, length = 0.1, lwd = 2)
  raster::scalebar(d = 2 * 1000,
                   label = "2 km", font = 2,
                   xy = c(712000, 6247000),
                   lonlat = FALSE)
}
## Raster contours
add_contour <-
  function(x, p = 0.5, ext = NULL, lwd = 0.5,...){
    if(!is.null(ext)) x <- raster::crop(x, ext)
    x <- spatialEco::raster.vol(x, p = p, sample = FALSE)
    raster::contour(x, nlevels = 1, drawlabels = FALSE, add = TRUE, lwd = lwd,...)
    return(invisible())
  }
## White out zeros
white_out <- function(x){
  if(is.null(x)) return(NULL)
  x[x == 0] <- NA
  return(x)
}


######################################
######################################
#### Depth time series parameters

#### Depth error model
# (tag error + tidal range + bathy error+ possible pelagic movements)
calc_depth_error <- function(depth) {
  e <- 4.77 + 2.5 + sqrt(0.5 ^2 + (0.013 * depth)^2)
  e <- matrix(c(-(e + 5), e), nrow = 2)
  return(e)
}
calc_depth_error(0)
calc_depth_error(300)
calc_depth_error <- Vectorize(calc_depth_error)
calc_depth_error(c(0, 300))


######################################
######################################
#### Acoustic parameters

#### Receiver detection_range (m)
detection_range <- 750

#### Define detection probability model
# Simulate some data that follow the results of drift testing/Klocker (2019)
drifts <- rbind(data.frame(distance = 0, detection = rbinom(100, 1, 0.97)),
                data.frame(distance = 425, detection = rbinom(100, 1, 0.5)),
                data.frame(distance = detection_range + 1, detection = rbinom(100, 1, 0)))
# Use simulated data to guestimate appropriate coefficients
drifts_mod <- glm(detection ~ distance, family = binomial, data = drifts)
summary(drifts_mod)
# Define detection probability model
calc_dpr <- function(distance) {
  pr <- stats::plogis(4 + distance * -0.01)
  pr[distance > detection_range] <- 0
  return(pr)
}
# Examine detection probability model
calc_dpr(0); calc_dpr(425); calc_dpr(750)
plot(1:detection_range, calc_dpr(1:detection_range),
     ylim = c(0, 1),
     type = "l")

#### Mobility
mobility             <- 500
mobility_from_origin <- mobility * 2

#### Euclidean-distance threshold (from examine_lcps.R)
euclid_distance_barrier_limit <- 265
lcp_predict <- function(distance, barrier){
  1.06450396994437 * distance * (barrier == 0) +
    1.17288343596087 * distance * (barrier == 1)
}

######################################
######################################
#### Particle filtering parameters

#### Movement model
calc_mpr <- function(distance, data){
  if(data$state == 0){
    pr <- stats::plogis(7.5 + distance * -0.5)
  } else if(data$state == 1){
    pr <- stats::plogis(7.5 + distance * -0.05)
  }
  pr[pr < 0.0001] <- 0.0001
  pr[distance > mobility] <- 0
  return(pr)
}
# Plot zoom-in of state 0 model
plot(0:50, calc_mpr(0:50, data.frame(state = 0)), ylim = c(0, 1))
# Plot state 0 and state 1 models
plot(0:mobility, calc_mpr(0:mobility, data.frame(state = 0)),
     xlim = c(0, 500), ylim = c(0, 1), col = "red", type = "l")
lines(0:mobility, calc_mpr(0:mobility, data.frame(state = 1)))
# Check probabilities
which.min(calc_mpr(1:mobility, data.frame(state = 0)))
which.min(calc_mpr(1:mobility, data.frame(state = 1)))

#### Movement model from origin
# This is defined as above, but the maximum distance threshold is relaxed
calc_mpr_from_origin <- function(distance, data){
  if(data$state == 0){
    pr <- stats::plogis(7.5 + distance * -0.5)
  } else if(data$state == 1){
    pr <- stats::plogis(7.5 + distance * -0.05)
  }
  pr[distance > mobility_from_origin] <- 0
  return(pr)
}

#### Number of particles
n_particles <- 1000


######################################
######################################
#### Visualise parameters/models

png("./fig/models.png", height = 3.5, width = 10,
    units = "in", res = 600)
pp <- par(mfrow = c(1, 3), oma = c(1, 1, 1, 4))

#### Detection probability model
x <- seq(0, detection_range + 1e-6, length.out = 1000)
pretty_plot(x, calc_dpr(x),
            xlab = "", ylab = "",
            ylim = c(0, 1), cex.axis = 1.2,
            type = "l")
mtext(side = 1, "Distance from receiver (m)", line = 2.5)
mtext(side = 2, "Detection probability", line = 2.5)
mtext(side = 3, "A", font = 2, adj = -0.15, cex = 2, line = 2)

#### Movement model
# Plot state 0 and state 1 models
pretty_plot(0:mobility, calc_mpr(0:mobility, data.frame(state = 0)),
            xlim = c(0, 500), ylim = c(0, 1),
            xlab = "", ylab = "",
            col = "grey", lwd = 2, type = "l", cex.axis = 1.3)
lines(0:mobility, calc_mpr(0:mobility, data.frame(state = 1)), lwd = 2)
legend("topright", legend = c("Mode 0", "Mode 1"),
       lwd = 2, col = c("grey", "black"), box.lty = 3)
mtext(side = 1, "Distance between locations (m)", line = 2.5)
mtext(side = 2, "Movement probability", line = 2.5)
mtext(side = 3, "B", font = 2, adj = -0.15, cex = 2, line = 2)

#### Depth error model
# Define seabed and observed (individual) depths
depth_seabed  <- 0:max(bathy_zlim)
depth_ind_adj <- calc_depth_error(depth_seabed)
depth_ind     <- rbind(depth_seabed + depth_ind_adj[1, ],
                       depth_seabed + depth_ind_adj[2, ])
depth_ind[depth_ind < min(depth_seabed)] <- min(depth_seabed)
depth_ind[depth_ind > max(depth_seabed)] <- max(depth_seabed)
# Define blank plot
lim <- range(depth_seabed)
plot(depth_seabed, depth_ind[1, ] * -1,
     xlim = lim, ylim = sort(lim *-1),
     axes = FALSE,
     xlab = "", ylab = "",
     type = "n")
# Add background shading for depths
for (i in seq_len((length(bathy_col_param$breaks) - 1))) {
  print(i)
  shallow <- bathy_col_param$breaks[i]
  deeper <- bathy_col_param$breaks[i + 1]
  x <- c(shallow, deeper, deeper, 0, 0, shallow)
  y <- c(0, 0, deeper, deeper, shallow, shallow)
  polygon(x, y*-1, col = bathy_col_param$col[i],
          border = bathy_col_param$col[i])
}
# Add depth error envelope
polygon(c(depth_seabed, rev(depth_seabed)),
        c( depth_ind[1, ] * -1, rev(depth_ind[2, ] * -1)),
        border = FALSE,
        col = scales::alpha("grey", 0.5))
# Add 1:1 line for comparison
lines(c(lim[1], lim[2]), c(lim[1], lim[2] * -1), lty = 3)
# Add axes and titles
bins <- seq(0, lim[2], by = 50)
pretty_axis(side = 3:2,
            lim = list(x = lim, y = sort(lim*-1)),
            axis = list(list(at = bins),
                        list(at = sort(bins * -1),
                             labels = abs(sort(bins * -1)))),
            control_axis = list(cex.axis = 1.3),
            add = TRUE)
mtext(side = 3, "Seabed depth (m)", line = 2.5)
mtext(side = 2, "Observed (individual) depth (m)", line = 2.5)
mtext(side = 3, "C", font = 2, adj = -0.15, cex = 2, line = 2)

par(pp)
dev.off()

#### End of code.
######################################
######################################
