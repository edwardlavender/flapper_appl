######################################
######################################
#### workflow_acpf_example_id.R

#### This code:
# 1) Follows an example workflow for the application of the AC/ACPF algorithms
# ... to the acoustic time series for an example individual.

#### Steps preceding this code:
#


######################################
######################################
#### Setup

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(flapper)

#### Define global parameters
# Projections
proj_wgs84 <- sp::CRS("+init=epsg:4326")
proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
                          "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#### Algorithm parameters
# Clock drift parameter
clock_drift <- 5
# Detection range (informed by other data)
det_rng     <- 750


######################################
######################################
#### Set up AC algorithm

#### Step (1) Define study area
# Grid resolution is a key consideration
# ... Resolution needs to balance probability estimates with computation time
# ... For the AC algorithm, resolution can be relatively low.
# ... For other algorithms, see the process_surface() function for further comments.
blank <- raster::raster(raster::extent(dat_gebco), res = c(75, 75))
gebco <- raster::resample(dat_gebco, blank)
raster::crs(gebco) <- proj_utm

#### Step (2) Prepare movement time series for algorithm
## (A) Implement standard processing
# E.g., define unique receiver IDs via process_receiver_id()
# E.g., check for false detections via glatos::false_detections() and flapper::process_false_detections_sf()
## Focus on a specific individual and a short sample of time series
id <- 25
acc <- dat_acoustics[dat_acoustics$individual_id == id, ]
acc <- acc[acc$timestamp >= as.POSIXct("2016-03-17 00:00:00", tz = "UTC") &
             acc$timestamp <= as.POSIXct("2016-03-20 00:00:00", tz = "UTC"), ]
## (B) Check for overlapping detections at non-overlapping receivers
# (i) Get detection centroid overlap(s)
moorings_sp   <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj_wgs84)
moorings_sp   <- sp::spTransform(moorings_sp, proj_utm)
centroids_dets <- get_detection_centroids(xy = moorings_sp,
                                          detection_range = det_rng,
                                          coastline = dat_coast,
                                          byid = TRUE)
centroids_df <- dat_moorings[, c("receiver_id",
                                 "receiver_start_date",
                                 "receiver_end_date")]
row.names(centroids_df) <- names(centroids_dets)
centroids_dets <- sp::SpatialPolygonsDataFrame(centroids_dets, centroids_df)
overlaps <- get_detection_centroids_overlap(centroids =  centroids_dets)
# (ii) Check for overlapping detections at non-overlapping receivers
get_detection_overlaps(acoustics = acc,
                       overlaps = overlaps,
                       clock_drift = clock_drift)
## (C) Examine characteristics of time series
nrow(acc)
difftime(max(acc$timestamp), min(acc$timestamp), units = "s")/120
gaps <- Tools4ETS::serial_difference(acc$timestamp, units = "days")
summary(gaps)

#### Step (3) Define mobility parameter
# (A) Estimates for mobility from acoustic data
moorings_spdf <-
  sp::SpatialPointsDataFrame(moorings_sp,
                             dat_moorings[, c("receiver_id",
                                              "receiver_start_date",
                                              "receiver_end_date")])
get_mvt_mobility_from_acoustics(data = acc,
                                moorings = moorings_spdf,
                                detection_range = det_rng,
                                transmission_interval = 90,
                                step = 120)
# (B) Estimates for mobility from archival data
arc <- dat_archival[dat_archival$individual_id == acc$individual_id[1], ]
get_mvt_mobility_from_archival(arc)
## (C) Examine the assumption of a 'constant' mobility parameter
summary(acs_setup_mobility(arc$depth, mobility = 250))
mob <- 250

#### Step (3) Define the acoustic centroids for the algorithm
## (A) Get suggestions for the number of centroids (time steps)
n_timesteps <-
  acs_setup_n_centroids(detections = acc$timestamp,
                         step = 120,
                         moorings = dat_moorings,
                         mobility = mob,
                         boundaries = gebco)
# We'll focus the algorithms within the area of interest,
# ... so the estimates from Method (3) are most relevant.
# This approach suggests about ~105 centroids, but this does not account for the
# ... slight difference between the detection range and mobility, so we will
# ... use a few additional centroids.
n_timesteps <- 110

## (B) Make the acoustic centroids
centroids_acc <-
  acs_setup_centroids(xy = moorings_spdf,
                       detection_range = det_rng,
                       mobility = mob,
                       n_timesteps = n_timesteps,
                       coastline = dat_coast,
                       cl = parallel::makeCluster(10L)
  )

#### Step (4) Make detection detection probability kernels for algorithm
## (A) Define detection probability function based on distance and detection range
calc_dpr <-
  function(x){
    ifelse(x <= 750, stats::plogis(2.5 + -0.02 * x), 0)
  }
## (B) Get detection kernels
kernels <- acs_setup_detection_kernels(xy = moorings_spdf,
                                        centroids = centroids_acc,
                                        overlaps = overlaps,
                                        calc_detection_pr = calc_dpr,
                                        map = gebco,
                                        coastline = invert_poly(dat_coast))


######################################
######################################
#### Implement AC algorithm

#### Implement AC algorithm
run <- TRUE
if(run){
  # Cluster needs to be specified outside of the function and we need to export the raster package.
  cl <- parallel::makeCluster(2L)
  parallel::clusterEvalQ(cl, library(raster))
  out_ac <- ac(acoustics = acc,
               step = 120,
               bathy = gebco,
               detection_range = det_rng,
               detection_kernels = kernels,
               detection_kernels_overlap = overlaps,
               mobility = mob,
               acc_centroids = centroids_acc,
               save_record_spatial = NULL,
               con = "./data/algorithms/workflows_real/ac/",
               write_record_spatial_for_pf = list(filename = "./data/algorithms/workflows_real/ac/record/", format = "GTiff"),
               cl = cl,
               varlist = "kernels"
               )
}


######################################
######################################
#### Examine AC algorithm outputs

#### Simplify AC algorithm outputs
out_ac_s <- acdc_simplify(out_ac, mask = gebco)

#### Plot record (for an example time step)
acdc_plot_record(out_ac_s,
                 plot = 200,
                 add_coastline = list(x = dat_coast),
                 add_receivers = list(x = moorings_sp, pch = 21, bg = "royalblue", col = "royalblue", cex = 0.5)
                 )

#### Animate record
boundaries <- update_extent(raster::extent(dat_coast), -1000)
param <- list(record = out_ac_s,
              plot = NULL,
              add_coastline = list(x = dat_coast),
              add_receivers = list(x = moorings_sp, pch = 21, bg = "royalblue", col = "royalblue", cex = 2),
              xlim = boundaries[1:2], ylim = boundaries[3:4],
              crop_spatial = TRUE
              )
acdc_animate_record(expr_param = param,
                    html_name = "./fig/algorithms/workflows_real/ac/animations/animation_1.html")

#### Check animations
## Implement-algorithm stepwise
out_ac_stepwise <- ac(acoustics = acc,
                      step = 120,
                      bathy = gebco,
                      detection_range = det_rng,
                      detection_kernels = kernels,
                      detection_kernels_overlap = overlaps,
                      mobility = mob,
                      acc_centroids = centroids_acc,
                      save_record_spatial = NULL
                      )

#### Compare cumulative maps from the two approaches
maps_from_stepwise <- lapply(out_ac_stepwise$archive[[1]]$record, function(record_elm){
  lapply(record_elm$spatial, function(spatial_elm) spatial_elm$map_cumulative)
}) %>% unlist()
maps_from_chunkwise <- lapply(out_ac_s$record, function(record_elm){
  lapply(record_elm$spatial, function(spatial_elm) spatial_elm$map_cumulative)
}) %>% unlist()
length(maps_from_stepwise)
length(maps_from_chunkwise)
table(mapply(all.equal, maps_from_stepwise, maps_from_chunkwise))

#### Examine overall map
scale <- raster::maxValue(out_ac_s$map)
png("./fig/algorithms/workflows_real/ac/overall_map.png",
    height = 8, width = 8, units = "in", res = 600)
prettyGraphics::pretty_map(out_ac_s$map/scale,
                           add_rasters = list(x = out_ac_s$map/scale),
                           add_polys = list(x = dat_coast, lwd = 2),
                           xlim = c(700000,  711500),
                           ylim = c(6250000, 6270000),
                           crop_spatial = TRUE)
dev.off()


######################################
######################################
#### Implement particle filtering

#### List files
files <- pf_setup_record("./data/algorithms/workflows_real/ac/record/")
head(files)
str(files)

#### Implement particle filtering to refine maps
out_pf <- pf(record = files,
             mobility = mob,
             con = "./data/algorithms/workflows_real/acpf/acpf_log.txt",
             n = 100)
# saveRDS(out_pf, "./data/algorithms/workflows_real/acpf/out_pf.rds")



######################################
######################################
#### Examine the results of particle filtering

#### Examine particle histories through time
pf_plot_history(out_pf)
pf_animate_history(
  expr_param = list(archive = out_pf,
                    add_particles = list(cex = 0.75, pch = 21,
                                         col = "black", bg = "black"),
                    add_polys = list(x = dat_coast, lwd = 1),
                    xlim = boundaries[1:2], ylim = boundaries[3:4],
                    crop_spatial = TRUE,
                    prompt = FALSE),
  dir = "./fig/algorithms/workflows_real/acpf/animations/")

#### Generate an overall map of space use based on particle histories
png("./fig/algorithms/workflows_real/acpf/map_overall.png",
    height = 10, width = 10, units = "in", res = 600)
map_from_pf <- pf_plot_map(archive = out_pf,
                           map = gebco,
                           scale = "max",
                           add_points = list(x = moorings_sp, pch = 21, bg = "royalblue", col = "royalblue", cex = 0.2),
                           add_polys = list(x = dat_coast),
                           xlim = boundaries[1:2], ylim = boundaries[3:4],
                           crop_spatial = TRUE)
dev.

#### Compare maps with/without particle filtering
# Without particle filtering, possible locations. This is probably driven by
# ... (a) the lack of a movement limitations in the AC algorithm
# ... (b) the lack of sampling in the AC algorithm
png("./fig/algorithms/workflows_real/acpf/ac_versus_acpf.png",
    height = 8, width = 14, units = "in", res = 600)
pp <- par(mfrow = c(1, 2), oma = c(2, 2, 2, 2), mar = c(3, 3, 3, 3))
paa <- list(side = 1:4,
            pretty = list(n = 3),
            axis = list(list(), list(),
                        list(labels = FALSE),
                        list(labels = FALSE)),
            control_sci_notation = list(magnitude = 16L, digits = 0))
prettyGraphics::pretty_map(add_rasters = list(x = out_ac_s$map/raster::cellStats(out_ac_s$map, "max")),
                           add_points = list(x = moorings_sp, pch = 21, bg = "royalblue", col = "royalblue", cex = 0.2),
                           add_polys = list(x = dat_coast, lwd = 2),
                           xlim = boundaries[1:2], ylim = boundaries[3:4],
                           pretty_axis_args = paa,
                           crop_spatial = TRUE,
                           main = "AC")
prettyGraphics::pretty_map(add_rasters = list(x = map_from_pf),
                           add_points = list(x = moorings_sp, pch = 21, bg = "royalblue", col = "royalblue", cex = 0.2),
                           add_polys = list(x = dat_coast, lwd = 2),
                           xlim = boundaries[1:2], ylim = boundaries[3:4],
                           pretty_axis_args = paa,
                           crop_spatial = TRUE,
                           main = "ACPF")
par(pp)
dev.off()

#### Reconstruct specific movement paths
out_pf_paths <- pf_simplify(archive = out_pf, max_n_copies = 10)
# saveRDS(out_pf_paths, "./data/algorithms/workflows_real/acpf/out_pf_paths.rds")

#### Examine movement path probabilities
# Cell probabilities are a combination of
# ... (a) probabilities from AC (due to detection Pr)
# ... (b) the movement model
## Check probabilities for time 1:
r_1 <- raster::raster(files[1])
raster::extract(r_1, pf_paths$cell_id[1]); pf_paths$cell_pr[1];
## Check probabilities for time 2, which includes movement:
r_2 <- raster::raster(files[2])
raster::extract(r_2, pf_paths$cell_id[2]) * out_pf$args$calc_movement_pr(pf_paths$dist[2])

#### Visualise paths
## Example path
out_pf_paths_eg <- out_pf_paths[out_pf_paths$path_id == 1, ]
n_segments <- nrow(out_pf_paths_eg)
## Limits
path_xlim <- c(min(out_pf_paths$cell_x) - 125, max(out_pf_paths$cell_x) + 125)
path_ylim <- c(min(out_pf_paths$cell_y) - 125, max(out_pf_paths$cell_y) + 125)
## 2d example
png("./fig/algorithms/workflows_real/acpf/path_example_2d.png",
    height = 8, width = 14, units = "in", res = 600)
pf_plot_2d(out_pf_paths[out_pf_paths$path_id == 1, ], gebco,
           add_paths = list(length = 0.05, col = viridis::viridis(n = n_segments), lwd = 2),
           xlim = path_xlim,
           ylim = path_ylim,
           crop_spatial = TRUE,
           prompt = FALSE)
points(out_pf_paths_eg$cell_x[1], out_pf_paths_eg$cell_y[1],
       pch = 4, cex = 2, lwd = 5, col = viridis::viridis(1))
points(out_pf_paths_eg$cell_x[n_segments], out_pf_paths_eg$cell_y[n_segments],
       pch = 1, cex = 2, lwd = 5, col = viridis::viridis(n_segments)[n_segments])
dev.off()
## 3d example
out_pf_paths_eg$cell_z <- raster::extract(gebco, out_pf_paths_eg[, c("cell_x", "cell_y")]) - 10
surface <- raster::crop(gebco, raster::extent(path_xlim, path_ylim))
pf_plot_3d(out_pf_paths_eg, surface,
           add_paths = list(line = list(color = viridis::viridis(n_segments), width = 10)))
