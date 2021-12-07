######################################
######################################
#### examine_space_use.R

#### This code:
# 1) Implements the COA-KUD, ACPF and ACDCPF algorithms to examine
# ... habitat use for an example individual over a restricted time period.
# ... This code is supported by examine_space_use_time_trials.R

#### Steps preceding this code:
# 1) Definition of global parameters via define_global_param.R
# 2) Definition of datasets via process_data_raw.R
# 3) Some of the code follows that laid out in examine_cooccurrences.R


######################################
######################################
#### Set up

#### Wipe workspace & load flapper & global param
rm(list = ls())
source("./R/define_global_param.R")

#### Set raster options
rop <- raster::rasterOptions()
raster::rasterOptions(tmpdir = "./data/tmp/")
raster::tmpDir()

#### Load data
acoustics       <- readRDS("./data/movement/tag/acoustics_eg.rds")
moorings        <- readRDS("./data/movement/generic/moorings.rds")
moorings_xy     <- readRDS("./data/spatial/moorings_xy.rds")
archival        <- readRDS("./data/movement/tag/archival_eg.rds")
site_bathy      <- raster::raster("./data/spatial/site_bathy.tif")
site_coast      <- readRDS("./data/spatial/site_coast.rds")
# Lower resolution 'site habitat' layer for faster KUD estimation with large particle samples
site_habitat_lr <- readRDS("./data/spatial/site_habitat_lr.rds")

#### Process acoustic and archival time stamps
# Processing implemented in process_data_raw.R

#### Define putative 'resting' behaviour (following the threshold in Lavender et al. in review)
# I.e., moments when we suspect horizontal movement to be limited
archival$va      <- Tools4ETS::serial_difference(archival$depth)
archival$va_abs  <- abs(archival$va)
archival$state   <- ifelse(archival$va_abs <= 0.5, 0, 1)
archival$state[nrow(archival)] <- 1
# Define associated colours (for plotting)
cols <- c("gray", "black")
archival$col <- cols[factor(archival$state)]
archival$depth_neg <- archival$depth * - 1
acoustics$col      <- archival$col[match(acoustics$timestamp, archival$timestamp)]
table(is.na(acoustics$col))


######################################
######################################
#### Visualise movement time series

#### Visualise acoustic an archival time series (as in examine_cooccurrences.R)
acc_1 <- acoustics
# acc_1$col <- c("black", "royalblue", "red", "orange")[factor(acoustics$receiver_id)]
arc_1 <- archival
# Set up figure to save
png("./fig/space_use/movement_ts.png",
    height = 5, width = 10, units = "in", res = 600)
# Create blank plot
axis_ls <-
  prettyGraphics::pretty_plot(arc_1$timestamp, arc_1$depth * - 1,
                              pretty_axis_args = list(side = 3:2,
                                                      axis = list(list(format = "%d-%b-%y"), list())),
                              xlim = range(acc_1$timestamp), ylim = c(-225, 0),
                              xlab = "", ylab = "",
                              type = "n")
# Add depth time series
s <- nrow(arc_1)
arrows(x0 = arc_1$timestamp[1:(s-1)],
       x1 = arc_1$timestamp[2:s],
       y0 = arc_1$depth_neg[1:(s-1)],
       y1 = arc_1$depth_neg[2:s],
       col = arc_1$col,
       length = 0, lwd = 1)
# Add acoustic time series
# [plot acc_1 below acc_2 (below) (following pattern of depth time series)]
px <- par(xpd = NA)
prettyGraphics::pretty_line(acc_1$timestamp,
                            pretty_axis_args = list(axis_ls = axis_ls),
                            inherit = 1L, replace = list(pos = -220, labels = FALSE, lwd.ticks = 0),
                            add = TRUE,
                            pch = 21, col = acc_1$col, bg = acc_1$col, lwd = 1)
par(px)
# Add legend
legend(x = acc_1$timestamp[1] + 2.5e4, y = -5,
       lty = c(1, 1), col = c("grey", "black"),
       pch = c(21, 21), pt.bg = c("grey", "black"),
       legend = c("Mode 0", "Mode 1"),
       ncol = 2,
       box.lty = 3)
# Add titles
mtext(side = 3, "Time (dd-mm-yy)", cex = 1, line = 2)
mtext(side = 2, "Depth (m)", cex = 1, line = 2)
# Save
dev.off()

#### Define detection days
acoustics$date <- as.Date(acoustics$timestamp)
detection_days <- get_detection_days(acoustics = acoustics, receiver_id = acoustics$receiver_id)
moorings_xy$detection_days <-
  detection_days$detection_days[match(moorings_xy$receiver_id,
                                      detection_days$receiver_id)]

#### Map receivers with detections, weighted by detection days
prettyGraphics::pretty_map(
  x = site_bathy,
  add_rasters = list(x = site_bathy, zlim = bathy_zlim, col = bathy_col_param$col),
  add_polys = list(x = site_coast, col = "dimgrey"),
  add_points = list(list(x = moorings_xy, pch = 21, col = "red", bg = "red", cex = 0.5),
                    list(x = moorings_xy, cex = moorings_xy$detection_days/4,
                         pch = 21, bg = scales::alpha("brown", 0.5))
  ),
  crop_spatial = TRUE
)


######################################
######################################
#### Implement the COA-KUD algorithm

#### Define delta t
delta_t <-  "48 hours"

#### Check definition of delta_t
pp <- par(mfrow = c(1, 2))
coa_setup_delta_t(acoustics = acoustics,
                  delta_t = delta_t,
                  method = 1:2L,
                  implementation = 1L)
par(pp)

#### Implement COA algorithm
# Define moorings for the acoustic time series
moorings_for_acc    <- moorings[moorings$receiver_id %in% acoustics$receiver_id, ]
moorings_xy_for_acc <- moorings_xy[moorings_xy$receiver_id %in% acoustics$receiver_id, ]
moorings_xy_for_acc <- sp::coordinates(moorings_xy_for_acc)
# Define detection matrix
acc_mat <- make_matrix_detections(acoustics, moorings = moorings_for_acc)
# Check correspondance between detection matrix and mooring locations
colnames(acc_mat)
rownames(moorings_xy_for_acc)
# Implement COA algorithm
out_coa <- coa(mat = acc_mat, xy = moorings_xy_for_acc)
# Plot COAs
raster::plot(site_bathy)
points(out_coa)
points(moorings_xy_for_acc, pch = 4)

#### Implement KUD estimation
out_coa_kud <- NULL
if(nrow(out_coa) >= 5L){
  ## Get COA-KUD surface
  get_coa <- FALSE
  if(get_coa){
    out_coa_spdf <- sp::SpatialPointsDataFrame(
      out_coa[, c("x", "y")],
      data = data.frame(ID = factor(rep(1, nrow(out_coa)))),
      proj4string = raster::crs(site_bathy))
    out_coa_kud <- kud_around_coastline(xy = out_coa_spdf, grid = site_habitat_lr)
    out_coa_kud <- raster::raster(out_coa_kud[[1]])
    out_coa_kud <- raster::resample(out_coa_kud, site_bathy)
    out_coa_kud <- raster::mask(out_coa_kud, site_bathy)
    out_coa_kud <- out_coa_kud/raster::cellStats(out_coa_kud, "sum")
    raster::cellStats(out_coa_kud, "sum")
    # out_coa_kud <- out_coa_kud/raster::cellStats(out_coa_kud, "max")
    raster::writeRaster(out_coa_kud, "./data/movement/space_use/coa/out_coa_kud.tif")
  } else out_coa_kud <- raster::raster("./data/movement/space_use/coa/out_coa_kud.tif")
  ## Visualise surface
  png("./fig/space_use/out_coa_kud.png",
      height = 5, width = 6, res = 600, units = "in")
  out_coa_kud_scaled_to_one <- out_coa_kud/raster::cellStats(out_coa_kud, "max")
  prettyGraphics::pretty_map(
    x = site_bathy,
    add_rasters = list(x = out_coa_kud_scaled_to_one,
                       zlim = c(0, 1),
                       smallplot = c(0.785, 0.825, 0.27, 0.77),
                       axis.args = list(tck = -0.1, mgp = c(2.5, 0.2, 0), cex.axis = 1.75)),
    add_polys = add_coast,
    add_points = list(list(x = moorings_xy, pch = 21, col = "black", bg = "black", cex = 0.5),
                      list(x = moorings_xy, cex = moorings_xy$detection_days/4,
                           pch = 21, col = scales::alpha("brown", 0.5)),
                      list(x = out_coa, pch = 4, col = "red", cex = 1.5, lwd = 0.5)
    ),
    pretty_axis_args = paa
  )
  legend(x = 701700, y = 6248700,
         pch = c(1, 1),
         col = scales::alpha("brown", 0.5),
         pt.cex = c(5, 10)/4,
         legend = c("5", "10"),
         ncol = 2,
         x.intersp = 0.5, y.intersp	= 0.25,
         box.lty = 3, bg = NA,
         cex = 1.5)
  add_map_elements()
  mtext(side = 4, "Intensity", cex = 2, line = 0.5)
  dev.off()
}


######################################
######################################
#### Implement the ACDC algorithms

#### Isolate receivers
## Focus on receivers within the study area
site_poly                 <- as(raster::extent(site_bathy), "SpatialPolygons")
raster::crs(site_poly)    <- proj_utm
moorings_xy$in_study_site <- rgeos::gContains(site_poly, moorings_xy, byid = TRUE)
moorings_xy               <- moorings_xy[which(moorings_xy$in_study_site), ]
moorings                  <- moorings[moorings$receiver_id %in% moorings_xy$receiver_id, ]
## Focus on receivers that were active during the period under consideration
study_interval    <- lubridate::interval(as.Date(min(acoustics$timestamp)), as.Date(max(acoustics$timestamp)))
moorings$interval <- lubridate::interval(moorings$receiver_start_date, moorings$receiver_end_date)
moorings$overlap  <- lubridate::int_overlaps(moorings$interval, study_interval)
moorings          <- moorings[moorings$overlap, ]
moorings_xy       <- moorings_xy[which(moorings_xy$receiver_id %in% moorings$receiver_id), ]
moorings$receiver_id == moorings_xy$receiver_id

#### Define detection centroids
det_centroids <- acs_setup_centroids(moorings_xy,
                                     detection_range = detection_range,
                                     coastline = site_coast,
                                     boundaries = raster::extent(site_bathy),
                                     plot = TRUE)

#### Define detection centroid overlaps
det_centroids_overlaps <-
  get_detection_centroids_overlap(centroids = do.call(raster::bind, plyr::compact(det_centroids)),
                                  services = NULL)

#### Define detection kernels
run <- FALSE
if(run){
  site_sea <- invert_poly(site_coast)
  det_kernels <- acs_setup_detection_kernels(xy = moorings_xy,
                                             services = NULL,
                                             centroids = det_centroids,
                                             overlaps = det_centroids_overlaps,
                                             calc_detection_pr = calc_dpr,
                                             bathy = site_bathy)
  saveRDS(det_kernels, "./data/movement/space_use/det_kernels.rds")
} else det_kernels <- readRDS("./data/movement/space_use/det_kernels.rds")

#### Implement the AC algorithm
run <- FALSE
if(run){
  cl <- parallel::makeCluster(4L)
  parallel::clusterEvalQ(cl,
                         {
                           library(raster)
                           raster::rasterOptions(tmpdir = "./data/tmp/")
                         })
  out_ac <- ac(acoustics = acoustics, # acoustics[1:3, ],
               step = 120,
               bathy = site_bathy,
               detection_centroids = det_centroids,
               detection_kernels = det_kernels,
               detection_kernels_overlap = det_centroids_overlaps,
               mobility = mobility,
               write_record_spatial_for_pf =
                 list(filename = "./data/movement/space_use/ac/record/"),
               con = "./data/movement/space_use/ac/",
               cl = cl, varlist = "det_kernels"
  )
  saveRDS(out_ac, "./data/movement/space_use/ac/out_ac.rds")
  out_ac_s <- acdc_simplify(out_ac, mask = site_bathy)
  saveRDS(out_ac_s, "./data/movement/space_use/ac/out_ac_s.rds")
  raster::writeRaster(out_ac_s$map, "./data/movement/space_use/ac/out_ac_map.tif")
  raster::removeTmpFiles(h = 0)
} else out_ac <- readRDS("./data/movement/space_use/ac/out_ac.rds")


#### Implement the ACDC algorithm
run <- FALSE
if(run){
  cl <- parallel::makeCluster(4L)
  parallel::clusterEvalQ(cl,
                         {
                           library(raster)
                           raster::rasterOptions(tmpdir = "./data/tmp/")
                         })
  out_acdc <- acdc(acoustics = acoustics, # acoustics[1:3, ],
                   archival = archival,
                   bathy = site_bathy,
                   detection_centroids = det_centroids,
                   detection_kernels = det_kernels,
                   detection_kernels_overlap = det_centroids_overlaps,
                   mobility = mobility,
                   calc_depth_error = calc_depth_error,
                   write_record_spatial_for_pf =
                     list(filename = "./data/movement/space_use/acdc/record/"),
                   con = "./data/movement/space_use/acdc/",
                   cl = cl, varlist = "det_kernels")
  saveRDS(out_acdc, "./data/movement/space_use/acdc/out_acdc.rds")
  out_acdc_s <- acdc_simplify(out_acdc, mask = site_bathy)
  saveRDS(out_acdc_s, "./data/movement/space_use/acdc/out_acdc_s.rds")
  raster::writeRaster(out_acdc_s$map, "./data/movement/space_use/acdc/out_acdc_map.tif")
  raster::removeTmpFiles(h = 0)
} else out_acdc <- readRDS("./data/movement/space_use/acdc/out_acdc.rds")


######################################
######################################
#### Implement the ACPF and ACDCPF algorithms

#### Implement ACPF algorithm
run <- FALSE
if(run){

  #### Define data for ACPF
  # Get data
  acpf_data <- acdc_access_dat(acdc_simplify(out_ac))
  acpf_data$depth   <- archival$depth[1:(nrow(archival)-1)]
  acpf_data$va      <- Tools4ETS::serial_difference(acpf_data$depth)
  acpf_data$va_abs  <- abs(acpf_data$va)
  acpf_data$state   <- ifelse(acpf_data$va_abs <= 0.5, 0, 1)
  acpf_data$state[nrow(acpf_data)] <- 1
  ## Get record
  out_ac_record <- pf_setup_record("./data/movement/space_use/ac/record/", pattern = "*.grd")

  #### Time trials for PF
  # See examine_space_use_time_trials.R

  #### Implement algorithm
  sink("./data/movement/space_use/acpf/rgass_log.txt")
  pf_opts <- pf_setup_optimisers(use_calc_distance_euclid_backend_grass = TRUE,
                                 use_grass_dir = "/Applications/GRASS-7.4.4.app/Contents/Resources")
  out_acpf <- pf(record = out_ac_record,
                 data = acpf_data,
                 bathy = site_bathy,
                 calc_movement_pr = calc_mpr,
                 n = n_particles,
                 write_history = list(file = "./data/movement/space_use/acpf/history/"),
                 con = "./data/movement/space_use/acpf/acpf_log.txt",
                 optimisers = pf_opts)
  sink()
  saveRDS(out_acpf, "./data/movement/space_use/acpf/out_acpf.rds")
} else out_acpf <- readRDS("./data/movement/space_use/acpf/out_acpf.rds")

#### Implement ACDCPF algorithm
run <- FALSE
if(run){
  ## Define data for ACDCPF
  acdcpf_data <- acdc_access_dat(acdc_simplify(out_acdc))
  acdcpf_data$depth   <- archival$depth[1:(nrow(archival)-1)]
  acdcpf_data$va      <- Tools4ETS::serial_difference(acdcpf_data$depth)
  acdcpf_data$va_abs  <- abs(acdcpf_data$va)
  acdcpf_data$state   <- ifelse(acdcpf_data$va_abs <= 0.5, 0, 1)
  acdcpf_data$state[nrow(acdcpf_data)] <- 1
  ## Get record
  out_acdc_record <- pf_setup_record("./data/movement/space_use/acdc/record/", pattern = "*.grd")
  ## Implement ACDCPF
  sink("./data/movement/space_use/acdcpf/rgass_log.txt")
  pf_opts <- pf_setup_optimisers(use_calc_distance_euclid_backend_grass = TRUE,
                                 use_grass_dir = "/Applications/GRASS-7.4.4.app/Contents/Resources")
  out_acdcpf <- pf(record = out_acdc_record,
                   data = acdcpf_data,
                   bathy = site_bathy,
                   calc_movement_pr = calc_mpr,
                   n = n_particles, # 10L
                   write_history = list(file = "./data/movement/space_use/acdcpf/history/"),
                   con = "./data/movement/space_use/acdcpf/acdcpf_log.txt",
                   optimisers = pf_opts)
  sink()
  saveRDS(out_acdcpf, "./data/movement/space_use/acdcpf/out_acdcpf.rds")
} else out_acdcpf <- readRDS("./data/movement/space_use/acdcpf/out_acdcpf.rds")


######################################
######################################
#### Particle-based maps

######################################
#### Prepare objects for distance calculations

#### Set mobility_from_origin and mobility param
# (This may improve speed in pf_simplify())
out_acpf$args$mobility_from_origin   <- mobility
out_acpf$args$mobility               <- mobility
out_acdcpf$args$mobility_from_origin <- mobility
out_acdcpf$args$mobility             <- mobility

#### Define movement barrier
site_barrier <- sf::st_as_sf(site_coast)
# plot(site_coast)

#### Define a raster with 'distance from coastline'
# This is used for barrier overlaps in pf_simplify()
# ... to minimise the number of intersections that are required
# ... (see segments_cross_barrier()).
run <- FALSE
if(run){
  ## Step 1: Convert site_coast to a raster
  # This takes ~ 1 minute one minute
  # For other applications, the fasterize package may be faster
  site_coast_grid <- raster::rasterize(site_coast, site_bathy)
  raster::plot(site_coast_grid)
  ## Step 2: Calculate distances from the coast across the grid
  # This takes ~ 5 minutes.
  # For some applications, the fasterRaster::fasterRastDistance() may be faster
  site_coast_distances <- raster::distance(site_coast_grid)
  raster::plot(site_coast_distances)
  ## Step 3: save rasters
  raster::writeRaster(site_coast_grid, "./data/spatial/site_coast_grid.tif")
  raster::writeRaster(site_coast_distances, "./data/spatial/site_coast_distances.tif")
} else {
  site_coast_grid      <- raster::raster("./data/spatial/site_coast_grid.tif")
  site_coast_distances <- raster::raster("./data/spatial/site_coast_distances.tif")
}

#### Define cost surface (~ 2 minutes + time to save)
run <- FALSE
if(run){
  site_bathy_costs <- lcp_costs(site_bathy)
  saveRDS(site_bathy_costs, "./data/movement/space_use/site_bathy_costs.rds")
} else site_bathy_costs <- readRDS("./data/movement/space_use/site_bathy_costs.rds")

#### Define graph (~ 2 minutes + time to save)
run <- FALSE
if(run){
  site_bathy_graph <- lcp_graph_surface(surface = site_bathy,
                                        cost = site_bathy_costs$dist_total
  )
  saveRDS(site_bathy_graph, "./data/movement/space_use/site_bathy_graph.rds")
} else site_bathy_graph <- readRDS("./data/movement/space_use/site_bathy_graph.rds")

#### Graph processing
run <- FALSE
if(run){
  # Get unique cells
  acpf_cells_unq   <- pf_access_particles_unique(out_acpf)
  acdcpf_cells_unq <- pf_access_particles_unique(out_acdcpf)
  cells_unq        <- unique(c(acpf_cells_unq, acdcpf_cells_unq))
  # Simplify graph [duration: a few s]
  t1_a <- Sys.time()
  site_bathy_graph_simp <- cppRouting::cpp_simplify(site_bathy_graph,
                                                    keep = cells_unq,
                                                    rm_loop = TRUE,
                                                    iterate = TRUE,
                                                    silent = FALSE)
  # Confirm that all sampled cells are in the simplified graph
  all(as.character(cells_unq) %in% site_bathy_graph$dict$ref)
  t2_b <- Sys.time()
  difftime(t2_b, t1_a)
  # Contract simplified graph [duration: >3 days --> not implemented]
  # t1_b <- Sys.time()
  # site_bathy_graph_cont <- cppRouting::cpp_contract(site_bathy_graph_simp)
  # t2_b <- Sys.time()
  # difftime(t2_b, t1_b)
  # Save processed graphs
  saveRDS(site_bathy_graph_simp, "./data/movement/space_use/site_bathy_graph_simp.rds")
  # saveRDS(site_bathy_graph_cont, "./data/movement/space_use/site_bathy_graph_cont.rds")
} else {
  site_bathy_graph_simp <- readRDS("./data/movement/space_use/site_bathy_graph_simp.rds")
  # site_bathy_graph_cont <- readRDS("./data/movement/space_use/site_bathy_graph_cont.rds")
}


######################################
#### Define particles for mapping

run <- FALSE
if(run){

  #### ACPF
  sink("./data/movement/space_use/acpf/processing/pf_simplify_1.txt")
  t1 <- Sys.time()
  out_acpf_s <- pf_simplify(out_acpf,
                            calc_distance = "lcp",
                            calc_distance_graph = site_bathy_graph,
                            calc_distance_limit = euclid_distance_limit,
                            calc_distance_barrier = site_barrier,
                            calc_distance_barrier_limit = euclid_distance_barrier_limit,
                            calc_distance_restrict = TRUE,
                            calc_distance_barrier_grid = site_coast_grid,
                            calc_distance_algorithm = "Dijkstra",
                            summarise_pr = max,
                            write_history = list(file = "./data/movement/space_use/acpf/processing/"),
                            cl = 10L,
                            return = "archive")
  t2 <- Sys.time()
  difftime(t2, t1)
  sink()
  saveRDS(out_acpf_s, "./data/movement/space_use/acpf/out_acpf_s.rds")

  #### ACDCPF
  sink("./data/movement/space_use/acdcpf/pf_simplify_1.txt")
  t1 <- Sys.time()
  out_acdcpf_s <- pf_simplify(out_acdcpf,
                              calc_distance = "lcp",
                              calc_distance_graph = site_bathy_graph,
                              calc_distance_limit = euclid_distance_limit,
                              calc_distance_barrier = site_barrier,
                              calc_distance_barrier_limit = euclid_distance_barrier_limit,
                              calc_distance_restrict = TRUE,
                              calc_distance_barrier_grid = site_coast_grid,
                              calc_distance_algorithm = "Dijkstra",
                              summarise_pr = max,
                              write_history = list(file = "./data/movement/space_use/acdcpf/processing/"),
                              cl = 10L,
                              return = "archive")
  t2 <- Sys.time()
  difftime(t2, t1)
  sink()
  saveRDS(out_acdcpf_s, "./data/movement/space_use/acdcpf/out_acdcpf_s.rds")

} else {
  ## ACPF
  out_acpf_s <- readRDS("./data/movement/space_use/acpf/out_acpf_s.rds")
  ## ACDCPF
  out_acdcpf_s <- readRDS("./data/movement/space_use/acdcpf/out_acdcpf_s.rds")
}


######################################
#### POU Maps

run <- FALSE
if(run){

  #### ACPF
  out_acpf_pou <- pf_plot_map(out_acpf_s,
                              map = site_bathy)
  raster::cellStats(out_acpf_pou, "sum")
  raster::writeRaster(out_acpf_pou, "./data/movement/space_use/acpf/out_acpf_pou.tif")

  #### ACDCPF [38 s]
  out_acdcpf_pou <- pf_plot_map(out_acdcpf_s,
                                map = site_bathy)
  raster::cellStats(out_acdcpf_pou, "sum")
  raster::writeRaster(out_acdcpf_pou, "./data/movement/space_use/acdcpf/out_acdcpf_pou.tif")

} else {
  out_acpf_pou   <- raster::raster("./data/movement/space_use/acpf/out_acpf_pou.tif")
  out_acdcpf_pou <- raster::raster("./data/movement/space_use/acdcpf/out_acdcpf_pou.tif")
}

######################################
#### Fit KUDs

run <- FALSE
if(run){

  #### ACPF [6.57 hours]
  t1 <- Sys.time()
  out_acpf_kud <- pf_kud_2(xpf = out_acpf_s,
                           # sample_size = 100L,
                           bathy = site_bathy,
                           estimate_ud = kud_around_coastline,
                           grid = site_habitat_lr,
                           mask = site_bathy)
  t2 <- Sys.time()
  difftime(t2, t1)
  raster::writeRaster(out_acpf_kud, "./data/movement/space_use/acpf/out_acpf_kud.tif")

  ## ACDCPF [4.70 hours]
  t1 <- Sys.time()
  out_acdcpf_kud <- pf_kud_2(xpf = out_acdcpf_s,
                             # sample_size = 100L,
                             bathy = site_bathy,
                             estimate_ud = kud_around_coastline,
                             grid = site_habitat_lr,
                             mask = site_bathy
  )
  t2 <- Sys.time()
  difftime(t2, t1)
  raster::writeRaster(out_acdcpf_kud, "./data/movement/space_use/acdcpf/out_acdcpf_kud.tif")

} else {
  out_acpf_kud   <- raster::raster("./data/movement/space_use/acpf/out_acpf_kud.tif")
  out_acdcpf_kud <- raster::raster("./data/movement/space_use/acdcpf/out_acdcpf_kud.tif")
}


######################################
######################################
#### Visualise maps of space use

#### Process POU layers
# [Not currently implemented for consistency across all plots]
# Set zero areas to NA (white) so that they stand out more clearly
out_acpf_pou_2 <- out_acpf_pou
# out_acpf_pou_2[out_acpf_pou_2 == 0] <- NA
out_acdcpf_pou_2 <- out_acdcpf_pou
# out_acdcpf_pou_2[out_acdcpf_pou_2 == 0] <- NA

#### POU maps
png("./fig/space_use/pou_maps.png",
    height = 4, width = 8, res = 600, units = "in")
pp <- par(mfrow = c(1, 2), oma = c(0, 0, 0, 4), mar = c(1, 2, 1, 2))

legend_param <- list(smallplot = c(0.88, 0.92, 0.27, 0.77),
                     axis.args = list(tck = -0.1, mgp = c(2.5, 0.2, 0), cex.axis = 1.25)
                     )
prettyGraphics::pretty_map(add_rasters = rlist::list.merge(list(x = out_acpf_pou_2), legend_param),
                           add_polys = add_coast,
                           pretty_axis_args = paa)
add_map_elements()
mtext(side = 3, bquote(bold(.("A")) ~ "(ACPF)"), adj = 0.05, line = -2, cex = 1.25)
prettyGraphics::pretty_map(add_rasters = rlist::list.merge(list(x = out_acdcpf_pou_2), legend_param),
                           add_polys = add_coast,
                           pretty_axis_args = paa)
mtext(side = 3, bquote(bold(.("B")) ~ "(ACDCPF)"), adj = 0.05, line = -2, cex = 1.25)
mtext(side = 4, "POU", line = 3.75, cex = 1.25)
add_map_elements()
par(pp)
dev.off()

#### KUD maps [code modified from examine_depth_use.R to match graphical parameters]
## ACPF
png("./fig/space_use/out_acpf_kud.png",
    height = 5, width = 6, res = 600, units = "in")
out_acpf_kud_scaled_to_one <- out_acpf_kud/raster::cellStats(out_acpf_kud, "max")
prettyGraphics::pretty_map(add_rasters = list(x = out_acpf_kud_scaled_to_one,
                                              zlim = c(0, 1),
                                              plot_method = raster::plot,
                                              legend = FALSE),
                           add_polys = add_coast,
                           pretty_axis_args = paa)
add_map_elements()
dev.off()
## ACDCPF
png("./fig/space_use/out_acdcpf_kud.png",
    height = 5, width = 6, res = 600, units = "in")
out_acdcpf_kud_scaled_to_one <- out_acdcpf_kud/raster::cellStats(out_acdcpf_kud, "max")
prettyGraphics::pretty_map(add_rasters = list(x = out_acdcpf_kud_scaled_to_one,
                                              zlim = c(0, 1),
                                              plot_method = raster::plot,
                                              legend = FALSE),
                           add_polys = add_coast,
                           pretty_axis_args = paa)
add_map_elements()
dev.off()


#### End of code.
######################################
######################################
