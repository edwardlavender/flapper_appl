######################################
######################################
#### examine_post_release_paths.R

#### This code:
# 1) Implements the DCPF algorithm for example individuals that
# ... appear to exhibit irregular post-release behaviour by
# ... Lavender et al (2022) to reconstruct possible
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
id <- id[2]

#### Load data
root <- "./data/movement/post_release_paths/"
archival_pr <- readRDS(paste0(root, id, "/archival_pr.rds"))
xy_release  <- readRDS(paste0(root, id, "/xy_release.rds"))
xy_release  <- sp::coordinates(xy_release)
site_coast  <- readRDS(paste0(root, id, "/site_coast.rds"))
site_bathy  <- raster::raster(paste0(root, id, "/site_bathy.tif"))


######################################
######################################
#### Visualise post-release depth time series

#### Load individual time series
arc_1 <- readRDS(paste0(root, "1507/archival_pr.rds"))
arc_2 <- readRDS(paste0(root, "1558/archival_pr.rds"))

#### Assign resting behaviour (as implemented below)
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

#### Plot depth/behaviour time series
# Set up plot to save
png("./fig/post_release_paths/depth_ts.png",
    height = 5, width = 5, units = "in", res = 600)
# Create blank plot
ylim <- c(-200, 0)
prettyGraphics::pretty_plot(arc_1$time_index, arc_1$depth * - 1,
                            pretty_axis_args = list(side = 3:2),
                            xlab = "", ylab = "",
                            type = "n")
# Add depth time series (as points)
points(arc_1$time_index, arc_1$depth_neg, cex = 0.75, pch = 21, col = arc_1$col)
points(arc_2$time_index, arc_2$depth_neg, cex = 0.75, pch = 21, col = arc_2$col, bg = arc_2$col)
# Mark the start of the 'post-release' period
spos <- which.max(arc_1$include)
arrows(x0 = spos, x1 = spos,
       y0 = arc_1$depth_neg[spos] + 40,
       y1 = arc_1$depth_neg[spos] + 10,
       length = 0.04)
spos <- which.max(arc_2$include)
arrows(x0 = spos, x1 = spos,
       y0 = arc_2$depth_neg[spos] + 40,
       y1 = arc_2$depth_neg[spos] + 10,
       length = 0.04)
# lines(rep(which.max(arc_1$include), 2), ylim, lty = 3, lwd = 0.5)
# lines(rep(which.max(arc_2$include), 2), ylim, lwd = 0.5)
# Add depth time series (as lines)
s <- nrow(arc_1)
arrows(x0 = arc_1$time_index[1:(s-1)],
       x1 = arc_1$time_index[2:s],
       y0 = arc_1$depth_neg[1:(s-1)],
       y1 = arc_1$depth_neg[2:s],
       col = arc_1$col,
       length = 0, lwd = 2, lty = 3)
arrows(x0 = arc_2$time_index[1:(s-1)],
       x1 = arc_2$time_index[2:s],
       y0 = arc_2$depth_neg[1:(s-1)],
       y1 = arc_2$depth_neg[2:s],
       col = arc_2$col,
       length = 0, lwd = 2)
# Add legend
legend(x = 0.5, y = -175,
       lty = c(3, 1), lwd = c(2, 2),
       pch = c(21, 21), pt.bg = c(NA, "black"),
       legend = c(as.character(arc_1$dst_id[1]), as.character(arc_2$dst_id[1])),
       ncol = 2, x.intersp = 0.75,
       box.lty = 3)
# Add titles
mtext(side = 3, "Time (index)", cex = 1, line = 1.75)
mtext(side = 2, "Depth (m)", cex = 1, line = 2.5)
# Save
dev.off()


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
nrow(archival_pr_from_seabed)

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
               normalise = TRUE,
               write_record_spatial_for_pf = list(filename = paste0(root, id, "/dc/record/"), format = "GTiff", overwrite = TRUE),
               con = paste0(root, id, "/dc/dc_log.txt"),
               split = split,
               cl = cl)
  saveRDS(out_dc, paste0(root, id, "/dc/out_dc.rds"))
} else out_dc <- readRDS(paste0(root, id, "/dc/out_dc.rds"))


#### Process out_dc and check availability of depth contours at each time step
# (normalisation = FALSE for legacy reasons)
out_dc_s <- acdc_simplify(out_dc, type = "dc", mask = site_bathy, normalisation = FALSE)
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
  raster::cellStats(raster::raster(out_dc_record[1]), "sum")
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
# For 1507, this takes 0.36 minutes.
# For 1558, this takes 0.39 minutes.
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
# For 1507, this takes 2.82 minutes with 1,000 paths.
# For 1558, this takes 3.05 minutes with 1,000 paths.
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
pt_size <- 15
xyz_release <- cbind(xy_release, 0)
colnames(xyz_release) <- c("x", "y", "z")
out_pf_3d <- pf_plot_3d(paths = out_pf_lcps_sbt,
                        bathy = site_bathy,
                        aggregate = list(fact = 2),
                        zlim = bathy_zlim, shift = 50,
                        add_surface = list(colorscale = bathy_col_param$col),
                        coastline = site_coast,
                        add_markers = list(x = xyz_release[1], y = xyz_release[2], z = xyz_release[3],
                                           marker = list(color = "black", size = pt_size)),

                        )
out_pf_3d <- out_pf_3d %>% plotly::layout(showlegend = FALSE)
out_pf_3d

#### Visualise paths [use interpolated LCPs] [zoom in]
# [out_pf_lcps_sbt$path_id == unique(out_pf_lcps_sbt$path_id)[5], ]
out_pf_3d_zoom <- pf_plot_3d(paths = out_pf_lcps_sbt ,
                             bathy = site_bathy_sbt,
                             # aggregate = list(fact = 2),
                             xlim = xlim, ylim = ylim,
                             zlim = bathy_zlim, shift = 50,
                             add_surface = list(colorscale = bathy_col_param$col),
                             # coastline = site_coast_sbt,
                             # coastline_paths = list(line = list(width = 20, color = "black")),
                             add_paths = list(line = list(width = 15)),
                             add_markers = list(x = xyz_release[1], y = xyz_release[2], z = xyz_release[3],
                                                marker = list(color = "black", size = 30))
                             )
out_pf_3d_zoom <- out_pf_3d_zoom %>% plotly::layout(showlegend = FALSE)
out_pf_3d_zoom
# Save plotly plots manually.


#### End of code.
######################################
######################################
