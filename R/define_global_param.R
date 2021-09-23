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





#### End of code.
######################################
######################################
