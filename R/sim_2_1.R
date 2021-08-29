######################################
######################################
#### sim_2_1.R

#### This code:
# 1) Implements flapper family algorithms for the second scenario (array 2).
# This code is copied from sim_1_1.R with [[1]] changed to [[2]], 1_1 changed to 2_1 and all other 1_ references checked.

######################################
######################################
#### Implement COA approach

#### Guess a suitable delta_t interval
delta_t <- "12 hours"
pp <- par(mfrow = c(1, 2))
coa_setup_delta_t(dat_sim_acoustics_2, delta_t = delta_t, method = 1:2L)
par(pp)

#### Calculate COAs
dat_sim_acoustics_2$individual_id <- 1L
dat_sim_acoustics_2$receiver_id <- factor(dat_sim_acoustics_2$receiver_id, dat_sim_moorings_2$receiver_id)
dat_sim_acoustics_2_mat <- make_matrix_detections(dat_sim_acoustics_2, delta_t = delta_t)
out_coa_2_1 <- coa(dat_sim_acoustics_2_mat, sp::coordinates(dat_sim_array_2$array$xy))

#### Examine COAs
raster::plot(coast)
points(out_coa_2_1, col ="red")

#### Estimate KUD from COAs
out_coa_2_1_spdf <- sp::SpatialPointsDataFrame(
  out_coa_2_1[, c("x", "y")],
  data = data.frame(ID = factor(rep(1, nrow(out_coa_2_1)))),
  proj4string = raster::crs(coast))
out_coa_2_1_ud <- kud_around_coastline(xy = out_coa_2_1_spdf, grid = habitat)
out_coa_2_1_ud <- raster::raster(out_coa_2_1_ud[[1]])
out_coa_2_1_ud <- raster::mask(out_coa_2_1_ud, invert_poly(coast))
raster::plot(out_coa_2_1_ud)


######################################
######################################
#### Flapper algorithm: array-specific set up

#### Define acoustic centroids (~ 4 minutes)
make_acoustic_centroids <- FALSE
if(make_acoustic_centroids){
  array <- dat_sim_array_list[[2]]
  rxy <- array$array$xy
  rxy <- sp::SpatialPointsDataFrame(rxy, data.frame(receiver_id = 1:length(rxy)))
  acoustic_centroids <- acs_setup_centroids(xy = rxy,
                                            detection_range = det_rng,
                                            mobility = mob_on_grid,
                                            n_timesteps = 250,
                                            coastline = coast,
                                            boundaries = NULL,
                                            plot = FALSE,
                                            cl = parallel::makeCluster(4L),
                                            verbose = TRUE
  )
  saveRDS(acoustic_centroids, "./data/algorithms/sim/array_2/acoustic_centroids.rds")
} else {
  acoustic_centroids <- readRDS("./data/algorithms/sim/array_2/acoustic_centroids.rds")
}

#### Detection centroid overlaps
make_detection_centriods_overlaps <- FALSE
if(make_detection_centriods_overlaps){
  centroids_dets <- get_detection_centroids(xy = dat_sim_array_2$array$xy,
                                            detection_range = det_rng,
                                            coastline = coast,
                                            byid = TRUE)
  centroids_df <- dat_sim_moorings_2
  row.names(centroids_df) <- names(centroids_dets)
  centroids_dets <- sp::SpatialPolygonsDataFrame(centroids_dets, centroids_df)
  overlaps <- get_detection_centroids_overlap(centroids =  centroids_dets)
  saveRDS(overlaps, "./data/algorithms/sim/array_2/detection_centroids_overlaps.rds")
} else {
  overlaps <- readRDS("./data/algorithms/sim/detection_centroids_overlaps.rds")
}

#### Define detection kernels
make_detection_kernels <- FALSE
if(make_detection_kernels){
  array <- dat_sim_array_list[[2]]
  rxy <- array$array$xy
  rxy <- sp::SpatialPointsDataFrame(rxy, dat_sim_moorings_2)
  kernels <- acs_setup_detection_kernels(xy = rxy,
                                         centroids = acoustic_centroids,
                                         overlaps = overlaps,
                                         calc_detection_pr = calc_dpr,
                                         map = grid,
                                         coastline = invert_poly(coast))
  saveRDS(kernels, "./data/algorithms/sim/array_2/detection_kernels.rds")
} else {
  kernels <- readRDS("./data/algorithms/sim/array_2/detection_kernels.rds")
}

#......................# Warning: Detection probability is NA at receiver 6.
#......................# ... Is the 'coastline' mask the right way around (masking the land and not the sea)?

######################################
######################################
#### Implement AC and ACDC algorithms

#### Implement AC algorithm (~10 minutes on one core)
run_ac_2_1 <- FALSE
if(run_ac_2_1){
  out_ac_2_1 <- ac(acoustics = dat_sim_detections_by_array[[2]],
                   step = step,
                   bathy = grid,
                   detection_range = det_rng,
                   detection_kernels = kernels, detection_kernels_overlap = overlaps, detection_time_window = clock_drift,
                   acc_centroids = acoustic_centroids,
                   mobility = mob_on_grid,
                   save_record_spatial = 0,
                   con = "./data/algorithms/sim/array_2/ac_2_1/",
                   write_record_spatial_for_pf = list(filename = "./data/algorithms/sim/array_2/ac_2_1/record/", format = "GTiff")
  )
  saveRDS(out_ac_2_1, "./data/algorithms/sim/array_2/ac_2_1/out_ac_2_1.rds")
} else{
  out_ac_2_1 <- readRDS("./data/algorithms/sim/array_2/ac_2_1/out_ac_2_1.rds")
}

#### Implement ACDC algorithm
run_acdc_2_1 <- FALSE
if(run_acdc_2_1){
  dat_sim_archival_2 <- dat_sim_archival[dat_sim_archival$timestamp >= min(dat_sim_detections_by_array[[2]]$timestamp), ]
  head(dat_sim_archival_2)
  head(dat_sim_detections_by_array[[2]])
  out_acdc_2_1 <- acdc(acoustics = dat_sim_detections_by_array[[2]],
                       archival = dat_sim_archival_2,
                       bathy = grid,
                       detection_range = det_rng,
                       detection_kernels = kernels, detection_kernels_overlap = overlaps, detection_time_window = clock_drift,
                       acc_centroids = acoustic_centroids,
                       mobility = mob_on_grid,
                       calc_depth_error = function(...) matrix(c(-5, 5), nrow = 2),
                       save_record_spatial = 0,
                       con = "./data/algorithms/sim/array_2/acdc_2_1/",
                       write_record_spatial_for_pf = list(filename = "./data/algorithms/sim/array_2/acdc_2_1/record/", format = "GTiff")
  )
  saveRDS(out_acdc_2_1, "./data/algorithms/sim/array_2/acdc_2_1/out_acdc_2_1.rds")
} else{
  out_acdc_2_1 <- readRDS("./data/algorithms/sim/array_2/acdc_2_1/out_acdc_2_1.rds")
}


######################################
######################################
#### Examine AC/ACDC results

#### Plotting window
pp <- par(mfrow = c(1, 2))

#### AC
## Simplify outputs
out_ac_2_1_s <- acdc_simplify(out_ac_2_1)
## Overall map
raster::plot(out_ac_2_1_s$map)
prettyGraphics::add_sp_path(dat_sim_path$xy_mat, length = 0.01)

#### ACDC
## Simplify outputs
out_acdc_2_1_s <- acdc_simplify(out_acdc_2_1)
## Overall map
raster::plot(out_acdc_2_1_s$map)
# prettyGraphics::add_sp_path(dat_sim_path$xy_mat, length = 0.01)
par(pp)


######################################
######################################
#### Implement particle filtering

#### ACPF
## Define record
out_ac_2_1_record <- pf_setup_record("./data/algorithms/sim/array_2/ac_2_1/record/")
## Define movement model
# Defined above.
## Implement algorithm (1)
run_acpf_2_1 <- FALSE
if(run_acpf_2_1){
  out_acpf_2_1 <- pf(record = out_ac_2_1_record,
                     calc_movement_pr = calc_mpr_on_grid,
                     mobility = mob_on_grid,
                     n = 500L,
                     con = "./data/algorithms/sim/array_2/acpf_2_1/acpf_log.txt"
  )
  saveRDS(out_acpf_2_1, "./data/algorithms/sim/array_2/acpf_2_1/out_acpf_2_1.rds")
} else {
  out_acpf_2_1 <- readRDS("./data/algorithms/sim/array_2/acpf_2_1/out_acpf_2_1.rds")
}

#### ACDCPF
## Define record
out_acdc_2_1_record <- pf_setup_record("./data/algorithms/sim/array_2/acdc_2_1/record/")
## Define movement model
# Defined above.
## Implement algorithm (1)
run_acdcpf_2_1 <- FALSE
if(run_acdcpf_2_1){
  out_acdcpf_2_1 <- pf(record = out_acdc_2_1_record,
                       calc_movement_pr = calc_mpr_on_grid, # note relaxed movement model for grid
                       mobility = mob_on_grid,              # note relaxed mobility
                       n = 500L,
                       con = "./data/algorithms/sim/array_2/acdcpf_2_1/acdcpf_log.txt"
  )
  saveRDS(out_acdcpf_2_1, "./data/algorithms/sim/array_2/acdcpf_2_1/out_acdcpf_2_1.rds")
  beepr::beep(10)
} else {
  out_acdcpf_2_1 <- readRDS("./data/algorithms/sim/array_2/acdcpf_2_1/out_acdcpf_2_1.rds")
}

#### Fix paths to record files
out_acpf_2_1$args$record   <- out_ac_2_1_record
out_acdcpf_2_1$args$record <- out_acdc_2_1_record

######################################
######################################
#### Examine particle filtering

#### Examine overall histories for example time steps
pp <- par(mfrow = c(1, 2))
pf_plot_history(out_acpf_2_1, time_steps = 100)
pf_plot_history(out_acdcpf_2_1, time_steps = 100)
par(pp)

#### Assemble particle histories for connected cell pairs
## ACPF
run_pf_simplify <- FALSE
if(run_pf_simplify){
  out_ac_2_1_pairs <- pf_simplify(out_acpf_2_1,
                                  cl = parallel::makeCluster(4L),
                                  return = "archive"
  )
  saveRDS(out_ac_2_1_pairs, "./data/algorithms/sim/array_2/acpf_2_1/out_ac_2_1_pairs.rds")
} else {
  out_ac_2_1_pairs <- readRDS("./data/algorithms/sim/array_2/acpf_2_1/out_ac_2_1_pairs.rds")
}
## ACDCPF
run_pf_simplify <- FALSE
if(run_pf_simplify){
  out_acdc_2_1_pairs <- pf_simplify(out_acdcpf_2_1,
                                    cl = parallel::makeCluster(4L),
                                    return = "archive"
  )
  saveRDS(out_acdc_2_1_pairs, "./data/algorithms/sim/array_2/acdcpf_2_1/out_acdc_2_1_pairs.rds")
} else {
  out_acdc_2_1_pairs <- readRDS("./data/algorithms/sim/array_2/acdcpf_2_1/out_acdc_2_1_pairs.rds")
}

#### Examine overall maps
pf_plot_map(out_ac_2_1_pairs, map = grid, scale = "sum")
pf_plot_map(out_acdc_2_1_pairs, map = grid, scale = "sum")

#### Apply KUD to a sample of sampled particles
out_acpf_2_1_ud   <- pf_kud(out_ac_2_1_pairs,
                            bathy = grid, sample_size = 5000L,
                            estimate_ud = kud_around_coastline, grid = habitat)
out_acdcpf_2_1_ud <- pf_kud(out_acdc_2_1_pairs,
                            bathy = grid, sample_size = 5000L,
                            estimate_ud = kud_around_coastline, grid = habitat)


#### Compare the true path, COA estimates and AC, ACPF and ACDCPF outputs
save_png <- FALSE
if(save_png) png("./fig/algorithms/sim/array_2_eval.png",
                 height = 10, width = 7, units = "in", res = 600)
## Graphical param
pp <- par(mfrow = c(3, 3), oma= c(1, 1, 2, 1), mar = c(1, 0, 1, 0))
add_coastline <- list(x = coast, col = "white")
xlim <- ext[1:2]; ylim <- ext[3:4]
paa <- list(side = 1:4, axis = list(labels = FALSE))
cex_main <- 1.25
adj_1 <- 0.125
adj_2 <- 0.125
spaces <- "        "
## Make plots
# Simulated array
prettyGraphics::pretty_map(add_rasters = list(x = grid, plot_method = raster::plot, legend = FALSE),
                           add_points = list(x = dat_sim_array_2$array$xy, pch = 21, col = "royalblue", bg = "royalblue"),
                           add_polys = add_coastline,
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE
)
mtext(side = 3, "A", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(array)"), adj = adj_2)
# True path and KUD
prettyGraphics::pretty_map(add_rasters = list(x = dat_sim_path_ud, plot_method = raster::plot, legend = FALSE),
                           add_paths = list(x = dat_sim_path$xy_mat, lwd = 0.05, length = 0.01),
                           add_polys = add_coastline,
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "B", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(path & KUD)"), adj = adj_2)
# COAs and COA KUD
prettyGraphics::pretty_map(add_rasters = list(x = out_coa_2_1_ud, plot_method = raster::plot, legend = FALSE),
                           add_points = list(x = out_coa_2_1),
                           add_polys = add_coastline,
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "C", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(COAs & KUD)"), adj = adj_2)
# AC
prettyGraphics::pretty_map(add_rasters = list(x = out_ac_2_1_s$map, plot_method = raster::plot, legend = FALSE),
                           add_polys = add_coastline,
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "D", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(AC)"), adj = adj_2)
# ACPF
pf_plot_map(out_ac_2_1_pairs,
            add_rasters = list(plot_method = raster::plot, legend = FALSE),
            map = grid, scale = "sum",
            add_polys = add_coastline,
            xlim = xlim, ylim = ylim,
            pretty_axis_args = paa,
            crop_spatial = TRUE)
mtext(side = 3, "E", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACPF)"), adj = adj_2)
# ACPF + KUD
prettyGraphics::pretty_map(add_rasters = list(x = out_acpf_2_1_ud, plot_method = raster::plot, legend = FALSE),
                           add_polys = add_coastline,
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "F", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACPF KUD)"), adj = adj_2)
## ACDC
prettyGraphics::pretty_map(add_rasters = list(x = out_acdc_2_1_s$map, plot_method = raster::plot, legend = FALSE),
                           add_polys = add_coastline,
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "G", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDC)"), adj = adj_2)
## ACDCPF
pf_plot_map(out_acdc_2_1_pairs,
            add_rasters = list(plot_method = raster::plot, legend = FALSE),
            map = grid, scale = "sum",
            add_polys = add_coastline,
            xlim = xlim, ylim = ylim,
            pretty_axis_args = paa,
            crop_spatial = TRUE)
mtext(side = 3, "H", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDCPF)"), adj = adj_2)
## ACDCPF + KUD
prettyGraphics::pretty_map(add_rasters = list(x = out_acdcpf_2_1_ud, plot_method = raster::plot, legend = FALSE),
                           add_polys = add_coastline,
                           xlim = xlim, ylicm = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "I", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDCPF KUD)"), adj = adj_2)
par(pp)
if(save_png) dev.off()


#### End of code.
######################################
######################################
