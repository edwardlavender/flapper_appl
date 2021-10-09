######################################
######################################
#### examine_space_use.R

#### This code:
# 1) Implements the COA-KUD, ACPF and ACDCPF algorithms to examine
# ... habitat use for an example individual over a restricted time period.

#### Steps preceding this code:
# 1)


######################################
######################################
#### Set up

#### Wipe workspace & load flapper & global param
rm(list = ls())
source("./R/define_global_param.R")

#### Set raster options
rop <- raster::rasterOptions()
raster::rasterOptions(tmpdir= "./data/tmp/")

#### Load data
acoustics    <- readRDS("./data/movement/tag/acoustics_eg.rds")
moorings     <- readRDS("./data/movement/generic/moorings.rds")
moorings_xy  <- readRDS("./data/spatial/moorings_xy.rds")
archival     <- readRDS("./data/movement/tag/archival_eg.rds")
site_bathy   <- raster::raster("./data/spatial/site_bathy.tif")
site_coast   <- readRDS("./data/spatial/site_coast.rds")
site_habitat <- readRDS("./data/spatial/site_habitat.rds")

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
    out_coa_ud <- kud_around_coastline(xy = out_coa_spdf, grid = site_habitat)
    out_coa_ud <- raster::raster(out_coa_ud[[1]])
    out_coa_ud <- out_coa_ud/raster::cellStats(out_coa_ud, "max")
    saveRDS(out_coa_ud, "./data/movement/space_use/coa/out_coa_ud.rds")
  } else out_coa_ud <- readRDS("./data/movement/space_use/coa/out_coa_ud.rds")
  ## Visualise surface
  png("./fig/out_coa_kud.png",
      height = 4, width = 5, res = 600, units = "in")
  prettyGraphics::pretty_map(
    x = site_bathy,
    add_rasters = list(x = out_coa_ud,
                       smallplot = c(0.75, 0.78, 0.3, 0.75),
                       axis.args = list(tck = -0.1, mgp = c(2.5, 0.2, 0))),
    add_polys = list(x = site_coast, col = "dimgrey"),
    add_points = list(list(x = moorings_xy, pch = 21, col = "black", bg = "black", cex = 0.5),
                      list(x = moorings_xy, cex = moorings_xy$detection_days/4,
                           pch = 21, col = scales::alpha("brown", 0.5)),
                      list(x = out_coa, pch = 4, col = "red", cex = 1.5, lwd = 0.5)
                      ),
    pretty_axis_args = paa
  )
  legend(x = 701700, y = 6248600,
         pch = c(1, 1),
         col = scales::alpha("brown", 0.5),
         pt.cex = c(5, 10)/4,
         legend = c("5", "10"),
         ncol = 2,
         x.intersp = 0.75, y.intersp	= 0.25,
         box.lty = 3, bg = NA)
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
                                     plot = TRUE
                                     )

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
                                             map = site_bathy,
                                             coastline = site_sea,
                                             boundaries = raster::extent(site_bathy))
  saveRDS(det_kernels, "./data/movement/space_use/det_kernels.rds")
} else det_kernels <- readRDS("./data/movement/space_use/det_kernels.rds")

#### Implement the AC algorithm
run <- FALSE
if(run){
  cl <- parallel::makeCluster(4L)
  parallel::clusterEvalQ(cl, library(raster))
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
  raster::removeTmpFiles(h = 0)
} else out_ac <- readRDS("./data/movement/space_use/ac/out_ac.rds")


#### Implement the ACDC algorithm
run <- FALSE
if(run){
  cl <- parallel::makeCluster(4L)
  parallel::clusterEvalQ(cl, library(raster))
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
  raster::removeTmpFiles(h = 0)
} else out_acdc <- readRDS("./data/movement/space_use/acdc/out_acdc.rds")


######################################
######################################
#### Implement the ACPF and ACDCPF algorithms

#### Implement ACPF algorithm
run <- FALSE
if(run){
  ## Define data for ACPF
  acpf_data <- acdc_access_dat(acdc_simplify(out_ac))
  acpf_data$depth   <- acpf_data$archival_depth
  acpf_data$va      <- Tools4ETS::serial_difference(acpf_data$depth)
  acpf_data$va_abs  <- abs(acpf_data$va)
  acpf_data$state   <- ifelse(acpf_data$va_abs <= 0.5, 0, 1)
  acpf_data$state[nrow(acpf_data)] <- 1
  ## Get record
  out_ac_record <- pf_setup_record("./data/movement/space_use/ac/record/")
  ## Implement ACPF
  out_acpf <- pf(record = out_ac_record,
                 data = acpf_data,
                 bathy = site_bathy,
                 calc_movement_pr = calc_mpr,
                 mobility = mobility,
                 n = n_particles, # 10L
                 con = "./data/movement/space_use/acpf/acpf_log.txt")
  saveRDS(out_acpf, "./data/movement/space_use/acpf/out_acpf.rds")
} else out_acpf <- readRDS("./data/movement/space_use/acpf/out_acpf.rds")

#### Implement ACDCPF algorithm
run <- FALSE
if(run){
  ## Define data for ACDCPF
  acdcpf_data <- acdc_access_dat(acdc_simplify(out_acdc))
  acdcpf_data$depth   <- acdcpf_data$archival_depth
  acdcpf_data$va      <- Tools4ETS::serial_difference(acdcpf_data$depth)
  acdcpf_data$va_abs  <- abs(acdcpf_data$va)
  acdcpf_data$state   <- ifelse(acdcpf_data$va_abs <= 0.5, 0, 1)
  acdcpf_data$state[nrow(acdcpf_data)] <- 1
  ## Get record
  out_acdc_record <- pf_setup_record("./data/movement/space_use/acdc/record/")
  ## Implement ACDCPF
  out_acdcpf <- pf(record = out_acdc_record,
                   data = acdcpf_data,
                   bathy = site_bathy,
                   calc_movement_pr = calc_mpr,
                   mobility = mobility,
                   n = n_particles, # 10L
                   con = "./data/movement/space_use/acdcpf/acdcpf_log.txt")
  saveRDS(out_acdcpf, "./data/movement/space_use/acdcpf/out_acdcpf.rds")
} else out_acdcpf <- readRDS("./data/movement/space_use/acdcpf/out_acdcpf.rds")



######################################
######################################
#### Particle-based maps

#### Define particles for mapping
run <- FALSE
if(run){
  ## ACPF
  out_acpf_s_1 <- pf_simplify(out_acpf, return = "archive")
  out_acpf_s_2 <- pf_simplify(out_acpf_s_1, return = "archive", summarise_pr = max)
  saveRDS(out_acpf_s_1, "./data/movement/space_use/acpf/out_acpf_s_1.rds")
  saveRDS(out_acpf_s_2, "./data/movement/space_use/acpf/out_acpf_s_2.rds")
  ## acdcpf
  out_acdcpf_s_1 <- pf_simplify(out_acdcpf, return = "archive")
  out_acdcpf_s_2 <- pf_simplify(out_acdcpf_s_1, return = "archive", summarise_pr = max)
  saveRDS(out_acdcpf_s_1, "./data/movement/space_use/acdcpf/out_acdcpf_s_1.rds")
  saveRDS(out_acdcpf_s_2, "./data/movement/space_use/acdcpf/out_acdcpf_s_2.rds")
} else {
  ## ACPF
  out_acpf_s_1 <- readRDS("./data/movement/space_use/acpf/out_acpf_s_1.rds")
  out_acpf_s_2 <- readRDS("./data/movement/space_use/acpf/out_acpf_s_2.rds")
  ## ACDCPF
  out_acdcpf_s_1 <- readRDS("./data/movement/space_use/acdcpf/out_acdcpf_s_1.rds")
  out_acdcpf_s_2 <- readRDS("./data/movement/space_use/acdcpf/out_acdcpf_s_2.rds")
}

#### Fit KUDs
run <- FALSE
if(run){
  ## ACPF
  out_acpf_kud <- pf_kud_2(xpf = out_acpf_s_2,
                           bathy = site_bathy,
                           estimate_ud = kud_around_coastline,
                           grid = site_habitat,
                           mask = site_bathy)
  saveRDS(out_acpf_kud, "./data/movement/space_use/acpf/out_acpf_kud.rds")
  ## ACDCPF
  out_acdcpf_kud <- pf_kud_2(xpf = out_acdcpf_s_2,
                             bathy = site_bathy,
                             estimate_ud = kud_around_coastline,
                             grid = site_habitat,
                             mask = site_bathy)
  saveRDS(out_acdcpf_kud, "./data/movement/space_use/acdcpf/out_acdcpf_kud.rds")
} else {
  out_acpf_kud <- readRDS("./data/movement/space_use/acpf/out_acpf_kud.rds")
  out_acdcpf_kud <- readRDS("./data/movement/space_use/acdcpf/out_acdcpf_kud.rds")
}


######################################
######################################
#### Path-based maps (accounting for LCPs)

#### Build paths
max_n_paths <- 1000L
out_acpf_paths <- pf_simplify(out_acpf_1, return = "path", max_n_paths = max_n_paths)

#### Interpolate LCPs

#### Only consider paths for which we have been able to construct LCPS

#### Process paths in line with LCPs

#### Apply KUD estimation to retained paths


######################################
######################################
#### Visually compare maps of space use

#### Define plotting param

#### Map detection days

#### Map COA-KUD approach

#### Map selected ACPF-KUD (from particles or paths)

#### Map selected ACDCPF-KUD (from particles or paths)

#### Save plot


#### End of code.
######################################
######################################
