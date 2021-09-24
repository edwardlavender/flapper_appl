######################################
######################################
#### examine_post_release_paths.R

#### This code:
# 1) Implements the DCPF algorithm for example individuals that
# ... appear to exhibit irregular post-release behaviour by
# ... Lavender et al (in review) to reconstruct possible
# ... post-release movement paths.

#### Steps preceding this code:
# 1) Processing of raw data   ... via process_data_raw.R
# 2) Define global parameters ... via define_global_param.R


######################################
######################################
#### Set up

#### Wipe workspace & load flapper & global param
rm(list = ls())
source("./R/define_global_param.R")

#### Define individual
id <- c(1507, 1558)
id <- id[1]

#### Load data
root <- "./data/movement/post_release_paths/"
archival_pr <- readRDS(paste0(root, id, "/archival_pr.rds"))
xy_release  <- readRDS(paste0(root, id, "/xy_release.rds"))
xy_release  <- sp::coordinates(xy_release)
site_coast  <- readRDS(paste0(root, id, "/site_coast.rds"))
site_bathy  <- raster::raster(paste0(root, id, "/site_bathy.tif"))


######################################
######################################
#### Implement algorithm(s)

#### Define putative 'resting' behaviour (following the threshold in Lavender et al. in review)
# I.e., moments when we suspect horizontal movement to be limited
archival_pr$va      <- Tools4ETS::serial_difference(archival_pr$depth)
# archival_pr$va[nrow(archival_pr)] <- 0
archival_pr$va_abs  <- abs(archival_pr$va)
archival_pr$state <- ifelse(archival_pr$va_abs <= 0.5, 0, 1)
archival_pr$state[nrow(archival_pr)] <- 1

#### Filter archival_pr to focus on movements following return to the seabed
archival_pr_from_seabed <- archival_pr[archival_pr$include, ]

#### Examine depth time series
prettyGraphics::pretty_plot(archival_pr_from_seabed$timestamp, archival_pr_from_seabed$depth * - 1,
                            col = factor(archival_pr_from_seabed$state))

#### Implement DC algorithm (for selected individual)
# This takes ~ 14 s on three cores or ~ 1 minute on one core
run_dc <- FALSE
if(run_dc){
  n_cores <- 3L
  split <- floor(nrow(archival_pr_from_seabed)/n_cores)
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterEvalQ(cl, library(raster))
  out_dc <- dc(archival = archival_pr_from_seabed,
               bathy = site_bathy,
               calc_depth_error = calc_depth_error,
               write_record_spatial_for_pf = list(filename = paste0(root, id, "/dc/record/"), format = "GTiff", overwrite = TRUE),
               con = paste0(root, id, "/dc/dc_log.txt"),
               split = split,
               cl = cl)
  saveRDS(out_dc, paste0(root, id, "/dc/out_dc.rds"))
} else out_dc <- readRDS(paste0(root, id, "/dc/out_dc.rds"))


#### Process out_dc and check availability of depth contours at each time step
out_dc_s <- acdc_simplify(out_dc, type = "dc", mask = site_bathy)
out_dc_dat <- acdc_access_dat(out_dc_s)
table(out_dc_dat$availability)

#### Setup particle filtering
## Define record for particle filtering
out_dc_record <- pf_setup_record(paste0(root, id, "/dc/record/"), type = "dc")

#### Implement particle filtering
# For 1507, PF will run using 10 particles in 1 minute
# For 1558, larger numbers of particles are required (1000 particles --> 16 minutes)
# For consistency, for both individuals, n = 1000 particles are used.
run_pf <- FALSE
if(run_pf){
  out_pf <- pf(record = out_dc_record,
               data = archival_pr_from_seabed,
               origin = xy_release,
               bathy = site_bathy,
               calc_movement_pr = calc_mpr,
               calc_movement_pr_from_origin = calc_mpr_from_origin,
               mobility = mobility,
               mobility_from_origin = mobility_from_origin,
               n = n_particles,
               con = paste0(root, id, "/pf/pf_log.txt")
  )
  saveRDS(out_pf, paste0(root, id, "/pf/out_pf.rds"))
} else out_pf <- readRDS(paste0(root, id, "/pf/out_pf.rds"))

#### Process paths
# Step 1 (Stepping through time steps to join coordinate pairs) takes 18 s on one core or ~12 s on four cores
run_pf_simplify <- FALSE
if(run_pf_simplify){
  out_pf_paths <- pf_simplify(archive = out_pf,
                              cl = parallel::makeCluster(4L),
                              varlist = c("mobility", "mobility_from_origin"),
                              return = "path",
                              max_n_paths = 1000L)
  saveRDS(out_pf_paths, paste0(root, id, "/pf/out_pf_paths.rds"))
} else out_pf_paths <- readRDS(paste0(root, id, "/pf/out_pf_paths.rds"))


######################################
######################################
#### Examine paths

#### The number of reconstructed paths
max(out_pf_paths$path_id)

#### Check distances using LCPs
# This takes ~15 minute with 1,000 particles
run_lcp_interp <- FALSE
if(run_lcp_interp){
  out_pf_lcps <- lcp_interp(paths = out_pf_paths, surface = site_bathy)
  saveRDS(out_pf_lcps,  paste0(root, id, "/pf/out_pf_lcps.rds"))
} else out_pf_lcps <- readRDS(paste0(root, id, "/pf/out_pf_lcps.rds"))

#### Only consider paths for which we have been able to construct LCPS
out_pf_paths$dist_lcp <- out_pf_lcps$dist_lcp$dist

#### Process paths in line with LCPs
## Examine distances
range(out_pf_paths$dist, na.rm = TRUE)
range(out_pf_paths$dist_lcp, na.rm = TRUE)
## Drop any paths with overly large steps
# Drop any paths with distances beyond mobility_from_origin/mobility
path_id_impossible <-
  unique(
    c(out_pf_paths$path_id[which(out_pf_paths$dist_lcp[out_pf_paths$timestep == 1] > mobility_from_origin)],
      out_pf_paths$path_id[which(out_pf_paths$dist_lcp[out_pf_paths$timestep != 1] > mobility)])
  )
length(path_id_impossible)
out_pf_paths <- out_pf_paths[!(out_pf_paths$path_id %in% path_id_impossible), ]
## Recalculate probabilities based on LCPs
out_pf_paths$cell_pr_from_euclid <- out_pf_paths$cell_pr
for(i in 1:nrow(out_pf_paths)){
  if(out_pf_paths$timestep[i] == 1) {
    out_pf_paths$cell_pr[i] <-
      calc_mpr_from_origin(out_pf_paths$dist_lcp[i],
                           archival_pr_from_seabed[out_pf_paths$timestep[i],])
  } else if(out_pf_paths$timestep[i] > 1){
    out_pf_paths$cell_pr[i] <-
      calc_mpr(out_pf_paths$dist_lcp[i],
               archival_pr_from_seabed[out_pf_paths$timestep[i], ])
  }
}
range(out_pf_paths$cell_pr_from_euclid - out_pf_paths$cell_pr, na.rm = TRUE)

#### Define the likelihood of paths
head(out_pf_paths)
out_pf_paths_ll <- pf_loglik(out_pf_paths)

#### Select a sample of paths for plotting
path_id_egs <- out_pf_paths_ll$path_id[1:10]
out_pf_paths_sbt <- out_pf_paths[out_pf_paths$path_id %in% path_id_egs, ]
out_pf_lcps_sbt  <- out_pf_lcps$path_lcp[out_pf_lcps$path_lcp$path_id %in% path_id_egs, ]

#### Visualise paths (1d)
pf_plot_1d(paths = out_pf_paths_sbt, archival = archival_pr_from_seabed)

#### Visualise paths (2d) [use Euclidean paths for clarity]
## Define graphical param
path_cols <- data.frame(timestep = 1:max(out_pf_paths_sbt$timestep),
                        col = viridis::viridis(max(out_pf_paths_sbt$timestep)))
out_pf_paths_sbt$col <- path_cols$col[match(out_pf_paths_sbt$timestep, path_cols$timestep)]
## Wide angle
png(paste0("./fig/post_release_paths/", id, "/out_pf_2d.png"),
    height = 8, width = 10, units = "in", res = 800)
pf_plot_2d(paths = out_pf_paths_sbt,
           bathy = site_bathy,
           add_bathy = list(col = bathy_col_param$col, zlim = bathy_zlim),
           add_paths = list(length = 0.05, col = out_pf_paths_sbt$col),
           add_polys = list(x = site_coast, col = "dimgrey"))
points(xy_release, pch = 21, col = "black", bg = "black", cex = 2)
dev.off()
## Zoom
xlim <- range(pretty(range(out_pf_paths_sbt$cell_x)))
ylim <- range(pretty(range(out_pf_paths_sbt$cell_y)))
site_bathy_sbt <- raster::crop(site_bathy, raster::extent(xlim, ylim))
site_coast_sbt <- raster::crop(site_coast, raster::extent(xlim, ylim))
png(paste0("./fig/post_release_paths/", id, "/out_pf_2d_zoom.png"),
    height = 8, width = 10, units = "in", res = 800)
pf_plot_2d(paths = out_pf_paths_sbt,
           bathy = site_bathy_sbt,
           add_bathy = list(col = bathy_col_param$col, zlim = bathy_zlim),
           add_paths = list(length = 0.05, col = out_pf_paths_sbt$col, lwd = 2),
           add_polys = list(x = site_coast_sbt, col = "dimgrey"))
points(xy_release, pch = 21, col = "black", bg = "black", cex = 2)
dev.off()

#### Visualise paths (3d) [use interpolated LCPs]
xyz_release <- cbind(xy_release, 0)
colnames(xyz_release) <- c("x", "y", "z")
out_pf_3d <- pf_plot_3d(paths = out_pf_lcps_sbt,
                        bathy = site_bathy,
                        aggregate = list(fact = 2),
                        zlim = bathy_zlim, shift = 50,
                        add_surface = list(colorscale = bathy_col_param$col),
                        coastline = site_coast,
                        add_markers = list(x = xyz_release[1], y = xyz_release[2], z = xyz_release[3],
                                           marker = list(color = "black", width = 10))
                        )
out_pf_3d <- out_pf_3d %>% plotly::layout(showlegend = FALSE)

#### Visualise paths [use interpolated LCPs] [zoom in]
out_pf_3d_zoom <- pf_plot_3d(paths = out_pf_lcps_sbt,
                             bathy = site_bathy_sbt,
                             # aggregate = list(fact = 2),
                             xlim = xlim, ylim = ylim,
                             zlim = bathy_zlim, shift = 50,
                             add_surface = list(colorscale = bathy_col_param$col),
                             coastline = site_coast_sbt,
                             add_markers = list(x = xyz_release[1], y = xyz_release[2], z = xyz_release[3],
                                                marker = list(color = "black", width = 10))
                             )
out_pf_3d_zoom <- out_pf_3d_zoom %>% plotly::layout(showlegend = FALSE)
# Save plotly plots manually.


#### End of code.
######################################
######################################
