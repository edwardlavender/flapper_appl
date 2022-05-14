######################################
######################################
#### examine_cooccurrences.R

#### This code:
# 1) Implements the ACDCPF algorithm to reconstruct fine-scale movement
# ... paths of a pair of individuals regularly detected at the same time
# ... to examine possible interactions and fine-scale spatial partitioning.

#### Steps preceding this code:
# 1) Define global parameters via define_global_param.R
# 2) Process raw data via process_data_raw.R


######################################
######################################
#### Set up

#### Wipe workspace & load flapper & global param
rm(list = ls())
source("./R/define_global_param.R")

#### Load data
# Acoustic and archival time series
acc_1    <- readRDS("./data/movement/cooccurrences/acc_1.rds")
acc_2    <- readRDS("./data/movement/cooccurrences/acc_2.rds")
arc_1    <- readRDS("./data/movement/cooccurrences/arc_1.rds")
arc_2    <- readRDS("./data/movement/cooccurrences/arc_2.rds")
# Moorings data
moorings     <- readRDS("./data/movement/generic/moorings.rds")
moorings_xy  <- readRDS("./data/spatial/moorings_xy.rds")
# Study site fields
site_bathy   <- raster::raster("./data/spatial/site_bathy.tif")
site_coast   <- readRDS("./data/spatial/site_coast.rds")


######################################
######################################
#### Visualise time series

#### Define 'resting' behaviour as in examine_post_release_paths.R
# Define colors
cols <- c("gray", "black")
# arc_1
arc_1$va      <- Tools4ETS::serial_difference(arc_1$depth)
arc_1$va_abs  <- abs(arc_1$va)
arc_1$state <- ifelse(arc_1$va_abs <= 0.5, 0, 1)
arc_1$state[nrow(arc_1)] <- 1
arc_1$col <- cols[factor(arc_1$state)]
arc_1$depth_neg <- arc_1$depth * - 1
# arc_2
arc_2$va      <- Tools4ETS::serial_difference(arc_2$depth)
arc_2$va_abs  <- abs(arc_2$va)
arc_2$state <- ifelse(arc_2$va_abs <= 0.5, 0, 1)
arc_2$state[nrow(arc_2)] <- 1
arc_2$depth_neg <- arc_2$depth * - 1
arc_2$col <- cols[factor(arc_2$state)]
# Copy colours to acoustics dataframe
acc_1$col <- arc_1$col[match(acc_1$timestamp, arc_1$timestamp)]
acc_2$col <- arc_2$col[match(acc_2$timestamp, arc_2$timestamp)]

#### Plot acoustic and archival time series
# Set up plot to save
png("./fig/cooccurrences/movement_ts.png",
    height = 5, width = 5, units = "in", res = 600)
# Create blank plot
axis_ls <-
  prettyGraphics::pretty_plot(arc_1$timestamp, arc_1$depth * - 1,
                              pretty_axis_args = list(side = 3:2,
                                                      axis = list(list(format = "%H:%M"), list())),
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
       length = 0, lwd = 2, lty = 3)
arrows(x0 = arc_2$timestamp[1:(s-1)],
       x1 = arc_2$timestamp[2:s],
       y0 = arc_2$depth_neg[1:(s-1)],
       y1 = arc_2$depth_neg[2:s],
       col = arc_2$col,
       length = 0, lwd = 2)
# Add acoustic time series
# [plot acc_1 below acc_2 (below) (following pattern of depth time series)]
px <- par(xpd = NA)
prettyGraphics::pretty_line(acc_1$timestamp,
                            pretty_axis_args = list(axis_ls = axis_ls),
                            inherit = 1L, replace = list(pos = -220, labels = FALSE, lwd.ticks = 0),
                            add = TRUE,
                            pch = 21, col = acc_1$col, lwd = 2)
prettyGraphics::pretty_line(acc_2$timestamp,
                            pretty_axis_args = list(axis_ls = axis_ls),
                            inherit = 1L, replace = list(pos = -200, labels = FALSE, lwd.ticks = 0),
                            add = TRUE,
                            pch = 21, col = acc_2$col, bg = acc_2$col)
par(px)

# Add legend
legend(x = acc_1$timestamp[1] + 50, y = -5,
       lty = c(1, 3), lwd = c(2, 2),
       pch = c(21, 21), pt.bg = c("black", NA),
       legend = c(as.character(acc_2$acc_id[1]), as.character(acc_1$acc_id[1])),
       ncol = 1,
       box.lty = 3)
# Add titles
mtext(side = 3, "Time (hh:mm)", cex = 1, line = 1.75)
mtext(side = 2, "Depth (m)", cex = 1, line = 2.5)
# Save
dev.off()


######################################
######################################
#### Define study area

#### Define a buffer around the receivers at which individuals were detected
# ... The selected individuals were detected at a single receiver almost continuously
# ... So we can focus on this region (plus a small buffer) to speed up the algorithms
rxy <- moorings_xy[moorings_xy$receiver_id %in%
                     unique(c(acc_1$receiver_id, acc_2$receiver_id)), ]
max_gap <- max(Tools4ETS::serial_difference(acc_1$timestamp, units = "mins", na.rm = TRUE),
               Tools4ETS::serial_difference(acc_2$timestamp, units = "mins", na.rm = TRUE))
max_gap <- plyr::round_any(as.numeric(max_gap), 2, f = ceiling)
rxy_buf <- rgeos::gBuffer(rxy,
                          width = detection_range + mobility * max_gap,
                          quadsegs = 1e4)

#### Define the study site
site_bathy <- raster::crop(site_bathy, rxy_buf)
site_coast <- raster::crop(site_coast, rxy_buf)
site_habitat <- kud_habitat(site_bathy)
raster::writeRaster(site_bathy, "./data/movement/cooccurrences/site_bathy.tif")
# Visualise the study site
raster::plot(site_bathy)
raster::lines(site_coast, lwd = 2)


######################################
######################################
#### Set up the ACDC algorithm

#### Isolate receivers (as in examine_space_use.R)
## Focus on receivers within the study area
site_poly                 <- as(raster::extent(site_bathy), "SpatialPolygons")
raster::crs(site_poly)    <- proj_utm
moorings_xy$in_study_site <- rgeos::gContains(site_poly, moorings_xy, byid = TRUE)
moorings_xy               <- moorings_xy[which(moorings_xy$in_study_site), ]
moorings                  <- moorings[moorings$receiver_id %in% moorings_xy$receiver_id, ]
## Focus on receivers that were active during the period under consideration
study_interval    <- lubridate::interval(as.Date(min(c(acc_1$timestamp, acc_2$timestamp))),
                                         as.Date(max(c(acc_1$timestamp, acc_2$timestamp))))
moorings$interval <- lubridate::interval(moorings$receiver_start_date, moorings$receiver_end_date)
moorings$overlap  <- lubridate::int_overlaps(moorings$interval, study_interval)
moorings          <- moorings[moorings$overlap, ]
moorings_xy       <- moorings_xy[which(moorings_xy$receiver_id %in% moorings$receiver_id), ]
moorings$receiver_id == moorings_xy$receiver_id

#### Define detection centroids (as in examine_space_use.R)
det_centroids <- acs_setup_centroids(moorings_xy,
                                     detection_range = detection_range,
                                     coastline = site_coast,
                                     boundaries = raster::extent(site_bathy),
                                     plot = TRUE
                                     )

#### Define detection centroid overlaps (as in examine_space_use.R)
det_centroids_overlaps <-
  get_detection_centroids_overlap(centroids = do.call(raster::bind, plyr::compact(det_centroids)),
                                  services = NULL)

#### Define detection kernels (as in examine_space_use.R)
run <- FALSE
if(run){
  det_kernels <- acs_setup_detection_kernels(xy = moorings_xy,
                                             services = NULL,
                                             centroids = det_centroids,
                                             overlaps = det_centroids_overlaps,
                                             calc_detection_pr = calc_dpr,
                                             bathy = site_bathy)
  saveRDS(det_kernels, "./data/movement/cooccurrences/det_kernels.rds")
} else det_kernels <- readRDS("./data/movement/cooccurrences/det_kernels.rds")


######################################
######################################
#### Implement the ACDC algorithm

#### Select individual and movement time series
id <- 1:2
id <- id[2]
if(id == 1){
  acoustics <- acc_1; acc_1$acc_id[1]
  archival  <- arc_1
} else if(id == 2){
  acoustics <- acc_2
  archival  <- arc_2
} else stop("ID unrecognised.")
acoustics$index <- 1:nrow(acoustics)

#### Process acoustic and archival time stamps (as in examine_space_use.R)
# Implemented in process_data_raw.R

#### Implement the ACDC algorithm (for the selected individual)
# This takes ~20 s on one core for the selected time series.
run <- TRUE
if(run){
  out_acdc <- acdc(acoustics = acoustics,
                   archival = archival,
                   bathy = site_bathy,
                   detection_centroids = det_centroids,
                   detection_kernels = det_kernels,
                   detection_kernels_overlap = det_centroids_overlaps,
                   mobility = mobility,
                   calc_depth_error = calc_depth_error,
                   write_record_spatial_for_pf =
                     list(filename = paste0("./data/movement/cooccurrences/", id, "/acdc/record/"),
                          format = "GTiff"),
                   con = paste0("./data/movement/cooccurrences/", id, "/acdc/"))
  saveRDS(out_acdc, paste0("./data/movement/cooccurrences/", id, "/acdc/out_acdc.rds"))
} else out_acdc <- readRDS(paste0("./data/movement/cooccurrences/", id, "/acdc/out_acdc.rds"))


######################################
######################################
#### Implement PF

#### Define data for PF
pf_data <- acdc_access_dat(acdc_simplify(out_acdc))
pf_data$depth   <- pf_data$archival_depth
pf_data$va      <- Tools4ETS::serial_difference(pf_data$depth)
pf_data$va_abs  <- abs(pf_data$va)
pf_data$state   <- ifelse(pf_data$va_abs <= 0.5, 0, 1)
pf_data$state[nrow(pf_data)] <- 1

#### Implement ACDCPF algorithm
# This takes ~ 25 minutes per individual.
run <- TRUE
if(run){
  sink(paste0("./data/movement/cooccurrences/", id, "/acdcpf/rgass_log.txt"))
  pf_opts <- pf_setup_optimisers(use_calc_distance_euclid_backend_grass = TRUE,
                                 use_grass_dir = "/Applications/GRASS-7.4.4.app/Contents/Resources")
  out_acdc_record <- pf_setup_record(paste0("./data/movement/cooccurrences/", id, "/acdc/record/"))
  out_acdcpf <- pf(record = out_acdc_record,
                   data = pf_data,
                   bathy = site_bathy,
                   calc_movement_pr = calc_mpr,
                   mobility = mobility,
                   n = n_particles, # 10L
                   con = paste0("./data/movement/cooccurrences/", id, "/acdcpf/acdcpf_log.txt"),
                   optimisers = pf_opts)
  sink()
  saveRDS(out_acdcpf, paste0("./data/movement/cooccurrences/", id, "/acdcpf/out_adcpf.rds"))
} else out_acdcpf <- readRDS(paste0("./data/movement/cooccurrences/", id, "/acdcpf/out_adcpf.rds"))


######################################
######################################
#### Build paths

#### Build paths (as in examine_post_release_paths.R)
# This takes 0.83 minutes for ID 1.
# This takes 0.42 minutes for ID 2.
run_pf_simplify <- TRUE
if(run_pf_simplify){
  out_pf_paths <- pf_simplify(archive = out_acdcpf,
                              return = "path",
                              max_n_paths = 1000L)
  saveRDS(out_pf_paths, paste0("./data/movement/cooccurrences/", id, "/acdcpf/out_pf_paths.rds"))
} else out_pf_paths <- readRDS(paste0("./data/movement/cooccurrences/", id, "/acdcpf/out_pf_paths.rds"))


######################################
######################################
#### Examine paths

#### The number of reconstructed paths (as in examine_post_release_paths.R)
max(out_pf_paths$path_id)

#### Check distances using LCPs (as in examine_post_release_paths.R)
# This takes 1.63 minutes for ID 1 with 1,000 paths
# This takes 1.88 minutes for ID 2 with 1,000 paths
run_lcp_interp <- TRUE
if(run_lcp_interp){
  out_pf_lcps <- lcp_interp(paths = out_pf_paths, surface = site_bathy)
  saveRDS(out_pf_lcps,  paste0("./data/movement/cooccurrences/", id, "/acdcpf/out_pf_lcps.rds"))
} else out_pf_lcps <- readRDS(paste0("./data/movement/cooccurrences/", id, "/acdcpf/out_pf_lcps.rds"))

#### Check LCPs
head(out_pf_lcps$path_lcp)
head(out_pf_lcps$dist_lcp)

#### Only consider paths for which we have been able to construct LCPS (as in examine_post_release_paths.R)
out_pf_paths$dist_lcp <- out_pf_lcps$dist_lcp$dist

#### Process paths in line with LCPs (similar to in examine_post_release_paths.R)
run_process_paths <- TRUE
if(run_process_paths){
  ## Examine distances
  range(out_pf_paths$dist, na.rm = TRUE)
  range(out_pf_paths$dist_lcp, na.rm = TRUE)
  ## Drop any paths with overly large steps
  # Drop any paths with distances beyond mobility
  path_id_impossible <-
    unique(
      out_pf_paths$path_id[which(out_pf_paths$dist_lcp[out_pf_paths$timestep != 1] > mobility)]
    )
  length(path_id_impossible)
  out_pf_paths <- out_pf_paths[!(out_pf_paths$path_id %in% path_id_impossible), ]
  ## Recalculate probabilities based on LCPs
  out_pf_paths$cell_pr_from_euclid <- out_pf_paths$cell_pr
  for(i in 1:nrow(out_pf_paths)){
    if(out_pf_paths$timestep[i] >= 1){
      out_pf_paths$cell_pr[i] <-
        calc_mpr(out_pf_paths$dist_lcp[i],
                 pf_data[out_pf_paths$timestep[i], ])
    }
  }
  range(out_pf_paths$cell_pr_from_euclid - out_pf_paths$cell_pr, na.rm = TRUE)
  saveRDS(out_pf_paths, paste0("./data/movement/cooccurrences/", id, "/acdcpf/out_pf_paths_pro.rds"))

} else {

  #### Load processed paths for both individuals
  ## Processed paths
  out_pf_paths_1 <- readRDS(paste0("./data/movement/cooccurrences/", 1, "/acdcpf/out_pf_paths_pro.rds"))
  out_pf_paths_2 <- readRDS(paste0("./data/movement/cooccurrences/", 2, "/acdcpf/out_pf_paths_pro.rds"))
  ## LCPs
  out_pf_lcps_1 <- readRDS(paste0("./data/movement/cooccurrences/", 1, "/acdcpf/out_pf_lcps.rds"))
  out_pf_lcps_2 <- readRDS(paste0("./data/movement/cooccurrences/", 2, "/acdcpf/out_pf_lcps.rds"))

}


######################################
######################################
#### Visualise paths for both individuals

#### Define the likelihood of paths
out_pf_paths_ll_1 <- pf_loglik(out_pf_paths_1)
out_pf_paths_ll_2 <- pf_loglik(out_pf_paths_2)

#### Select a sample of paths for plotting
## ID 1
path_id_egs_1      <- out_pf_paths_ll_1$path_id[1]
out_pf_paths_1_sbt <- out_pf_paths_1[out_pf_paths_1$path_id %in% path_id_egs_1, ]
head(out_pf_paths_1_sbt)
out_pf_lcps_1_sbt  <- out_pf_lcps_1$path_lcp[out_pf_lcps_1$path_lcp$path_id %in% path_id_egs_1, ]
out_pf_lcps_1_sbt  <- out_pf_lcps_1_sbt %>% dplyr::arrange(timestep)
head(out_pf_lcps_1_sbt)
## ID 2
path_id_egs_2      <- out_pf_paths_ll_2$path_id[1]
out_pf_paths_2_sbt <- out_pf_paths_2[out_pf_paths_2$path_id %in% path_id_egs_2, ]
head(out_pf_paths_2_sbt)
out_pf_lcps_2_sbt  <- out_pf_lcps_2$path_lcp[out_pf_lcps_2$path_lcp$path_id %in% path_id_egs_2, ]
out_pf_lcps_2_sbt  <- out_pf_lcps_2_sbt %>% dplyr::arrange(timestep)
head(out_pf_lcps_2_sbt)
## Compare time series
# The time series are slightly different lengths
nrow(out_pf_paths_1_sbt); nrow(out_pf_paths_2_sbt)

#### Visualise paths (2d) [use Euclidean paths for clarity] (modified from examine_post_release_paths.R)
## Set up plot to save
png("./fig/cooccurrences/out_pf_2d.png",
    height = 10, width = 12, units = "in", res = 800)
## Define graphical parameters
pp <- par(mfrow = c(1, 2))
round_fct <- 1000
xlim <- range(c(out_pf_lcps_1_sbt$cell_x, out_pf_lcps_2_sbt$cell_x)); xlim
xlim[1] <- plyr::round_any(xlim[1], accuracy = round_fct, f = floor)
xlim[2] <- plyr::round_any(xlim[2], accuracy = round_fct, f = ceiling)
ylim <- range(c(out_pf_lcps_1_sbt$cell_y, out_pf_lcps_2_sbt$cell_y)); ylim
ylim[1] <- plyr::round_any(ylim[1], accuracy = round_fct, f = floor)
ylim[2] <- plyr::round_any(ylim[2], accuracy = round_fct, f = ceiling)
## Make plots
lapply(list(out_pf_paths_1_sbt, out_pf_paths_2_sbt), function(out_pf_paths_sbt){
  # out_pf_paths_sbt <- out_pf_paths_1_sbt
  path_cols <- data.frame(timestep = 1:max(out_pf_paths_sbt$timestep),
                          col = viridis::viridis(max(out_pf_paths_sbt$timestep)))
  out_pf_paths_sbt$col <- path_cols$col[match(out_pf_paths_sbt$timestep, path_cols$timestep)]
  pf_plot_2d(paths = out_pf_paths_sbt,
             bathy = site_bathy,
             xlim = xlim, ylim = ylim,
             add_bathy = list(col = bathy_col_param$col, zlim = bathy_zlim),
             add_paths = list(length = 0.05, col = out_pf_paths_sbt$col, lwd = 2),
             add_polys = list(x = site_coast, col = "dimgrey"),
             crop_spatial = TRUE)
}) %>% invisible()
dev.off()

#### Visualise paths (3d)
# Define graphical parameters
shift     <- 25
stretch   <- -7.5
line_lwd  <- 15
site_bathy_sbt <- raster::crop(site_bathy, raster::extent(xlim, ylim))

# Define receiver coordinates
rxyz <- data.frame(x = sp::coordinates(rxy)[, 1],
                   y = sp::coordinates(rxy)[, 2],
                   z = rxy$receiver_depth)
# Plot the path for ID 1
out_pf_3d <- pf_plot_3d(paths = out_pf_lcps_1_sbt,
                        add_paths = list(line = list(color = viridis::viridis(nrow(out_pf_lcps_1_sbt)),
                                                     width = line_lwd)),
                        bathy = site_bathy_sbt,
                        # aggregate = list(fact = 2),
                        xlim = xlim, ylim = ylim, zlim = bathy_zlim,
                        shift = shift, stretch = stretch,
                        add_surface = list(colorscale = bathy_col_param$col),
                        add_markers = list(x = rxyz$x, y = rxyz$y, z = rxyz$z,
                                           marker = list(color = "black", size = 20))
                        )
# Add the path for ID 2
out_pf_3d <-
  out_pf_3d %>%
  plotly::add_paths(x = out_pf_lcps_2_sbt$cell_x,
                    y = out_pf_lcps_2_sbt$cell_y,
                    z = (out_pf_lcps_2_sbt$cell_z * stretch + shift),
                    line = list(color = viridis::viridis(nrow(out_pf_lcps_2_sbt)),
                                width = line_lwd)
                    ) %>%
  plotly::layout(showlegend = FALSE)
# Visualise plot
out_pf_3d

#### Plot colour bars
png(paste0("./fig/cooccurrences/colour_bar_time.png"),
    height = 5, width = 3, units = "in", res = 600)
pp <- par(oma = c(1, 1, 1, 4))
fields::image.plot(site_bathy,
                   zlim = c(0, nrow(out_pf_lcps_1_sbt)),
                   legend.only = TRUE,
                   col = viridis::viridis(nrow(out_pf_lcps_1_sbt) + 1))
mtext(side = 4, "Time (steps)", cex = 1.25, line = 2.5)
dev.off()

#### Results
# 1) One individual probably resting in the deep channel
# 2) One individual was moving along depth contours at the top of the channel and then moved out of the channel
# ... For this pair of individuals, for this time window, there is no evidence of interaction
# ... ... or even behaving in the same way

#### End of code.
######################################
######################################

