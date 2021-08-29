######################################
######################################
#### sim_data.R

#### This code
# Simulates arrays, a movement paths and movement
# ... for evaluating the performance of algorithms
# ... for inferring patterns of space use.

#### Steps preceding this code
# NA.

######################################
######################################
#### Set up

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(flapper)

#### Global parameters
seed <- 20210828



######################################
######################################
#### Simulate arrays

######################################
#### Define area

#### Define general area
coast <- dat_coast
ext   <- update_extent(raster::extent(coast), -1250)
coast <- raster::crop(coast, ext)
sea   <- invert_poly(coast)
raster::plot(coast)
coast_buf <- rgeos::gBuffer(coast, width = 500)
sea_buf <- invert_poly(coast_buf)
raster::lines(sea_buf)

#### Define study area grid
# It is essential that this grid is used for simulation and modelling
blank <- raster::raster(raster::extent(dat_gebco), res = c(75, 75))
grid <- raster::resample(dat_gebco, blank)

#### Define study duration
# We will keep this as small as possible for simulation speed.
study_period <- seq(as.POSIXct("2021-01-01 12:00:00", tz = "UTC"),
                    as.POSIXct("2021-01-05 12:00:00", tz = "UTC"), by = "2 mins")

#### Define n receivers in arrays
n_receivers <- c(100, 10)


######################################
#### Simulate arrays

simulate_array <- FALSE
if(simulate_array){
   #### Array 1
   dat_sim_array_1 <- sim_array(boundaries = ext,
                                coastline = coast_buf,
                                n_receivers = n_receivers[1],
                                arrangement = "regular",
                                add_sea = list(),
                                seed = seed
   )
   dat_sim_array_1$array$xy <- raster::xyFromCell(grid, raster::cellFromXY(sp::coordinates(dat_sim_array_1$array$xy)))
   raster::plot(coast)
   points(dat_sim_array_1$array$xy, col = "red")
   # receivers are ~ 1000 km from each other
   saveRDS(dat_sim_array_1, "./data/algorithms/sim/data/dat_sim_array_1.rds")

   #### Array 2
   dat_sim_array_2 <- sim_array(boundaries = ext,
                                coastline = coast_buf,
                                n_receivers = n_receivers[1],
                                arrangement = "random",
                                seed = seed,
                                add_sea = list()
   )
   dat_sim_array_2$array$xy <- raster::xyFromCell(grid, raster::cellFromXY(sp::coordinates(dat_sim_array_2$array$xy)))
   saveRDS(dat_sim_array_2, "./data/algorithms/sim/data/dat_sim_array_2.rds")

   #### Array 3
   dat_sim_array_3 <- sim_array(boundaries = ext,
                                coastline = coast_buf,
                                n_receivers = n_receivers[2],
                                arrangement = "regular",
                                seed = seed,
                                add_sea = list()
   )
   dat_sim_array_3$array$xy <- raster::xyFromCell(grid, raster::cellFromXY(sp::coordinates(dat_sim_array_3$array$xy)))
   saveRDS(dat_sim_array_3, "./data/algorithms/sim/data/dat_sim_array_3.rds")

   #### Array 4
   dat_sim_array_4 <- sim_array(boundaries = ext,
                                coastline = coast_buf,
                                n_receivers = n_receivers[2],
                                arrangement = "random",
                                seed = seed,
                                add_sea = list()
   )
   dat_sim_array_4$array$xy <- raster::xyFromCell(grid, raster::cellFromXY(sp::coordinates(dat_sim_array_4$array$xy)))
   saveRDS(dat_sim_array_4, "./data/algorithms/sim/data/dat_sim_array_4.rds")

   #### List array designs
   dat_sim_array_list <- list(dat_sim_array_1 = dat_sim_array_1,
                              dat_sim_array_2 = dat_sim_array_2,
                              dat_sim_array_3 = dat_sim_array_3,
                              dat_sim_array_4 = dat_sim_array_4)
   saveRDS(dat_sim_array_list, "./data/algorithms/sim/data/dat_sim_array_list.rds")

} else {
   dat_sim_array_list <- readRDS("./data/algorithms/sim/data/dat_sim_array_list.rds")
}

#### Separate array-specific designs
dat_sim_array_1 <- dat_sim_array_list[[1]]
dat_sim_array_2 <- dat_sim_array_list[[2]]
dat_sim_array_3 <- dat_sim_array_list[[3]]
dat_sim_array_4 <- dat_sim_array_list[[4]]

#### Define moorings dataframes
dat_sim_moorings_list <- lapply(dat_sim_array_list, function(array){
   # array <- dat_sim_array_list[[1]]
   moorings <- data.frame(receiver_id = 1:length(array$array$xy),
                          receiver_start_date = as.Date(min(study_period)),
                          receiver_end_date = as.Date(max(study_period))
   )
   return(moorings)
})
dat_sim_moorings_1 = dat_sim_moorings_list[[1]]
dat_sim_moorings_2 = dat_sim_moorings_list[[2]]
dat_sim_moorings_3 = dat_sim_moorings_list[[3]]
dat_sim_moorings_4 = dat_sim_moorings_list[[4]]


######################################
######################################
#### Simulate movement

#### Define movement probability given distance and use this to simulate step lengths
# We will implement this approach so that the simulated/modelled step lengths are based on the same model
calc_mpr <- function(distance,...) {
   pr <- stats::plogis(10 + distance * -0.06)
   pr[distance > 250] <- 0
   return(pr)
}
plot(1:500, calc_mpr(1:500))
steps <- data.frame(distance = seq(0, 250, length.out = 1e4))
steps$pr <- calc_mpr(steps$distance)
sim_step_every_2_mins <- function(...,data = steps, size = 1) {
   sample(x = data$distance, size = size, prob = data$pr)
}
# prettyGraphics::pretty_hist(replicate(1e4, sim_step_every_2_mins()), xlim = c(0, 250), n = 100)
prettyGraphics::pretty_hist(sim_step_every_2_mins(size = 1e3))

simulate_movement <- FALSE
if(simulate_movement){

   #### Explore possible movement models
   # sim_step_every_2_mins <- function(...) stats::rgamma(1, shape = 8, scale = 8)
   # sim_step_every_2_mins <- function(...) stats::rlnorm(1, meanlog = 1, sdlog = 0.75)
   # sim_step_every_2_mins <- function(...) stats::rexp(1, rate = 0.01)
   # sim_step_every_2_mins <- function(...) stats::rweibull(1, shape = 100, scale = 100)

   #### Simulate starting location
   set.seed(seed)
   p_1 <- sp::coordinates(sp::spsample(sea_buf, n = 1, type = "random"))
   raster::plot(sea_buf)
   points(p_1, col = "red")

   #### Simulate movement in area
   dat_sim_path <- sim_path_sa(n = length(study_period),
                               p_1 = p_1,
                               area = sea_buf,
                               sim_step = sim_step_every_2_mins,
                               seed = seed)
   raster::lines(coast, col = "darkgreen")
   saveRDS(dat_sim_path, "./data/algorithms/sim/data/dat_sim_paths.rds")

} else {
   dat_sim_path <- readRDS("./data/algorithms/sim/data/dat_sim_paths.rds")
}

#### Re-express sampled locations on the grid
# Because of the low resolution of the grid, there is a mismatch.
dat_sim_path$xy_mat_on_grid <- raster::xyFromCell(grid, raster::cellFromXY(grid, dat_sim_path$xy_mat))

#### Check distances along the simulated path
## Original path
utils.add::basic_stats(
   raster::pointDistance(
      dat_sim_path$xy_mat[1:(nrow(dat_sim_path$xy_mat)- 1), ],
      dat_sim_path$xy_mat[2:nrow(dat_sim_path$xy_mat), ],
      lonlat = FALSE)
   )
## Path re-expressed on grid
dat_sim_path$dist_on_grid <-   rep(NA, nrow(dat_sim_path$xy_mat_on_grid))
dat_sim_path$dist_on_grid[2:length(dat_sim_path$dist_on_grid)] <-
   raster::pointDistance(
      dat_sim_path$xy_mat_on_grid[1:(nrow(dat_sim_path$xy_mat_on_grid)- 1), ],
      dat_sim_path$xy_mat_on_grid[2:nrow(dat_sim_path$xy_mat_on_grid), ],
      lonlat = FALSE)
utils.add::basic_stats(dat_sim_path$dist_on_grid, na.rm = TRUE)


######################################
######################################
#### Simulate observed movement time series

#### Define detection probability function based on distance
calc_dpr <-
   function(x){
      ifelse(x <= 750, stats::plogis(2.5 + -0.02 * x), 0)
   }

#### Simulate acoustic timeseries for each array design
simulate_detections <- FALSE
if(simulate_detections){
   dat_sim_detections_by_array <-
      pbapply::pblapply(dat_sim_array_list, function(array){
         ## Simulate detections at receivers (across grid)
         dat_sim_detections <- sim_detections(path = dat_sim_path$xy_mat_on_grid,   # expressed on grid
                                              xy = sp::coordinates(array$array$xy), # expressed on grid
                                              calc_detection_pr = calc_dpr,
                                              seed = seed,
                                              plot = FALSE)
         rownames(dat_sim_detections$det_mat) <- as.character(study_period)
         colnames(dat_sim_detections$det_mat) <- as.character(1:ncol(dat_sim_detections$det_mat))
         # Define 'acoustics' dataframe
         dat_sim_acoustics <- make_df_detections(dat_sim_detections$det_mat,
                                                 only_keep_detections = TRUE,
                                                 set_names = TRUE,
                                                 as_POSIXct = function(x) as.POSIXct(x, tz = "UTC"))
         dat_sim_acoustics$receiver_id <- as.integer(as.character(dat_sim_acoustics$receiver_id))

         return(dat_sim_acoustics)
      })
   saveRDS(dat_sim_detections_by_array, "./data/algorithms/sim/data/dat_sim_detections_by_array.rds")

} else{
   dat_sim_detections_by_array <- readRDS("./data/algorithms/sim/data/dat_sim_detections_by_array.rds")
}

#### Isolate array-specific acoustics dataframes
dat_sim_acoustics_1 <- dat_sim_detections_by_array[[1]]
dat_sim_acoustics_2 <- dat_sim_detections_by_array[[2]]
dat_sim_acoustics_3 <- dat_sim_detections_by_array[[3]]
dat_sim_acoustics_4 <- dat_sim_detections_by_array[[4]]

#### Generate archival time series
simulate_archival <- FALSE
if(simulate_archival){
   # Get depth corresponding to locations
   dat_sim_archival <- data.frame(indivdual_id = 1,
                                  timestamp = study_period,
                                  depth_sim = raster::extract(grid, dat_sim_path$xy_mat_on_grid) # on grid
                                  )
   table(is.na(dat_sim_archival$depth_sim))
   raster::plot(coast)
   points(dat_sim_path$xy_mat[is.na(dat_sim_archival$depth_sim), ], pch=".", col = "red")
   # Simulate depths using a depth error model
   dat_sim_archival$depth <- runif(1:nrow(dat_sim_archival),
                                   dat_sim_archival$depth_sim - 5,
                                   dat_sim_archival$depth_sim + 5)
   saveRDS(dat_sim_archival, "./data/algorithms/sim/data/dat_sim_archival.rds")
} else {
   dat_sim_archival <- readRDS("./data/algorithms/sim/data/dat_sim_archival.rds")
}


######################################
######################################
#### Examine simulated patterns of space use

#### Define habitat/non habitat for UD estimation
habitat <- raster::setValues(grid, 0)
habitat <- raster::crop(habitat, ext)
habitat <- raster::mask(habitat, coast, updatevalue = 1)
habitat <- methods::as(habitat, "SpatialPixelsDataFrame")
# sp::plot(habitat)

#### Estimate UD for simulated data
dat_sim_path_spdf <- sp::SpatialPointsDataFrame(
   dat_sim_path$xy_mat_on_grid,
   data = data.frame(ID = factor(rep(1, nrow(dat_sim_path$xy_mat)))),
   proj4string = raster::crs(coast))
dat_sim_path_ud <- kud_around_coastline(xy = dat_sim_path_spdf, grid = habitat)
dat_sim_path_ud <- raster::raster(dat_sim_path_ud[[1]])
dat_sim_path_ud <- raster::mask(dat_sim_path_ud, invert_poly(coast))

#### Plot UD for simulated data
raster::plot(dat_sim_path_ud)


######################################
######################################
#### Set up algorithms (shared param)

#### Define study area
# Define above.

#### Set numeric algorithm parameters
step             <- 120
det_rng          <- 750      # as defined previously because on the scale of the grid
clock_drift      <- 5
mob_on_grid      <- 250 + 75 # mobility is higher than simulated when expressed at the resolution of the grid

#### Redefine detection probability function on grid
# This is unnecessary because we have re-expressed receiver locations on the grid.

#### Define movement model over grid (slightly relaxed)
calc_mpr_on_grid <- function(distance,...) {
   pr <- stats::plogis(5 + distance * -0.035)
   pr[distance > 325] <- 0
   return(pr)
}
plot(1:500, calc_mpr_on_grid(1:500), type = "l")
lines(1:500, calc_mpr(1:500), col = "red")

#### Prepare movement time series
## Examine frequency of detections
prettyGraphics::pretty_line(dat_sim_acoustics_1$timestamp)
utils.add::basic_stats(Tools4ETS::serial_difference(dat_sim_acoustics_1$timestamp, units = "mins"), na.rm = TRUE)

#### Check mobility estimates for each array
mob_est_by_array <- mapply(dat_sim_detections_by_array, dat_sim_array_list, FUN = function(acoustics, array){
   # array <- dat_sim_array_list[[1]]
   rxy <- array$array$xy
   rxy <- sp::SpatialPointsDataFrame(rxy, data.frame(receiver_id = 1:length(rxy)))
   get_mvt_mobility_from_acoustics(data = acoustics,
                                   moorings = rxy,
                                   detection_range = det_rng,
                                   transmission_interval = 120,
                                   step = 120)
   })
get_mvt_mobility_from_archival(dat_sim_archival)

#### Now proceed to implement 'flapper' algorithms in subsequent scripts e.g., sim_acpf_1_1.R

#### End of code.
######################################
######################################
