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


######################################
######################################
#### Spatial parameters

#### Define projections
proj_wgs84 <- sp::CRS("+init=epsg:4326")
proj_utm <- sp::CRS("+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs")

#### Plotting param
bathy_zlim <- c(0, 225)
bathy_col_param <- prettyGraphics::pretty_cols_brewer(bathy_zlim,
                                                      scheme = "Blues",
                                                      n_breaks = max(bathy_zlim))
paa <- list(side = 1:4,
            axis = list(labels = FALSE, tck = -0.01)
            )


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

#### Euclidean-distance threshold (from evaluate_lcps.R)
euclid_distance_limit         <- 481
euclid_distance_barrier_limit <- 421
# These parameters support particle processing routines, by
# ... restricting the number of LCP calculations that are required.
# ... For distances 421-481, we will check barrier overlaps
# ... and for the subset of cell pairs that overlap with the coastline
# ... we will calculate LCPs.
# ... We will calculate LCPs for all cell connections exceeding 481 m
# ... in Euclidean distances.
# ... Thus, all distances less than 421 m and distances between 421--481 m
# ... are assumed to be less than mobility. LCPs are only calculated for the
# ... remaining cell pairs.


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



#### End of code.
######################################
######################################
