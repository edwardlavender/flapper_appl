######################################
######################################
#### examine_depth_use.R

#### This code:
# 1) Implements the DC algorithm for an example individual to examine habitat
# ... representation in terms of depth use.

#### Steps preceding this code:
# 1) Processing of raw data   ... via process_data_raw.R
# 2) Define global parameters ... via define_global_param.R


######################################
######################################
#### Set up

#### Wipe workspace & load flapper & global param
rm(list = ls())
source("./R/define_global_param.R")

#### Load data
archival       <- readRDS("./data/movement/tag/archival_eg.rds")
site_bathy     <- raster::raster("./data/spatial/site_bathy.tif")
site_coast     <- readRDS("./data/spatial/site_coast.rds")
site_sediments <- readRDS("./data/spatial/site_sediments.rds")


######################################
######################################
#### Implement DC algorithm

#### Implement algorithm
# This takes 1 hr 40 mins with 4 cores (3 hr 38 m with write_record_spatial_for_pf)
run_dc <- FALSE
if(run_dc){

  ## Define algorithm param
  n_cores <- 4L
  split <- floor(nrow(archival)/n_cores)
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterEvalQ(cl, library(raster))

  ## Implement algorithm
  out_dc <- dc(archival = archival,
               bathy = site_bathy,
               calc_depth_error = calc_depth_error,
               check_availability = TRUE,
               write_record_spatial_for_pf = list(filename = "./data/movement/depth_use/dc/record/"),
               con = "./data/movement/depth_use/dc/dc_log.txt",
               split = split,
               cl = cl)
  saveRDS(out_dc, "./data/movement/depth_use/dc/out_dc.rds")

} else out_dc <- readRDS("./data/movement/depth_use/dc/out_dc.rds")


######################################
######################################
#### Implement DC processing

run_dc_processing <- FALSE
if(run_dc_processing){
  #### Process out_dc
  # (normalisation implemented below & set to false here for legacy reasons)
  out_dc_s <- acdc_simplify(out_dc, type = "dc", mask = site_bathy, normalise = FALSE)

  #### Check availability of depth contours at each time step [all TRUE]
  out_dc_dat <- acdc_access_dat(out_dc_s)
  table(out_dc_dat$availability)

  #### Process map (for plotting)
  out_dc_map <- out_dc_s$map
  raster::writeRaster(out_dc_map, "./data/movement/depth_use/dc/out_dc_map.tif")

} else out_dc_map <- raster::raster("./data/movement/depth_use/dc/out_dc_map.tif")

#### Map processing
out_dc_map[is.na(out_dc_map)] <- 0
out_dc_map_pc <- (out_dc_map/nrow(archival))
raster::cellStats(out_dc_map_pc, "sum")
out_dc_map_pc_scaled <- out_dc_map_pc/raster::cellStats(out_dc_map_pc, "max")
raster::plot(out_dc_map)
raster::plot(out_dc_map > 0)
raster::lines(site_coast)


######################################
######################################
#### Map results

#### Plot map
png("./fig/depth_use/out_dc_map.png",
    height = 5, width = 6, res = 600, units = "in")
prettyGraphics::pretty_map(add_rasters = list(x = white_out(out_dc_map_pc_scaled),
                                              smallplot = c(0.785, 0.825, 0.27, 0.77),
                                              axis.args = list(tck = -0.1, mgp = c(2.5, 0.2, 0), cex.axis = 1.75),
                                              zlim = c(0, 1)),
                           add_polys = add_coast,
                           pretty_axis_args = paa)
mtext(side = 4, "POU", cex = 2, line = 1.25)
add_contour(out_dc_map_pc)
add_map_elements()
dev.off()


######################################
######################################
#### Explore home/core ranges

#### Map home and core ranges
## full range
out_dc_map_full <- out_dc_map
out_dc_map_full[out_dc_map_full > 0] <- 1
raster::plot(out_dc_map_full); raster::lines(site_coast, lwd = 2)
## 'home' range
out_dc_map_home <- out_dc_map
home <- raster::quantile(out_dc_map_home, 0.5)
out_dc_map_home[out_dc_map_home < home]  <- 0
out_dc_map_home[out_dc_map_home >= home] <- 1
prettyGraphics::pretty_map(add_rasters = list(x = out_dc_map_home),
                           add_polys = list(x = site_coast, col = "dimgrey", border = "dimgrey"))
## 'core' range
out_dc_map_core <- out_dc_map
core <- raster::quantile(out_dc_map_core, 0.95)
out_dc_map_core[out_dc_map_core < core]  <- 0
out_dc_map_core[out_dc_map_core >= core] <- 1
prettyGraphics::pretty_map(add_rasters = list(x = out_dc_map_core),
                           add_polys = list(x = site_coast, col = "dimgrey", border = "dimgrey"))


#### Relate putative patterns of depth use to sediment availability
## 'background' versus 'used' sediment availability in region examined
run_sediments <- FALSE
if(run_sediments){
  core_xy <- raster::rasterToPoints(out_dc_map_core, fun = function(x) x==1)
  core_xy <- sp::SpatialPoints(core_xy[, 1:2], proj4string = proj_utm)
  core_sediments <- sp::over(core_xy, site_sediments)
  saveRDS(core_sediments, "./data/movement/depth_use/dc/core_sediments.rds")
} else {
  core_sediments <- readRDS("./data/movement/depth_use/dc/core_sediments.rds")
  background_sediments <- readRDS("./data/spatial/site_sediments_background.rds")
}
## Examine putative sediment use versus availability
# Get sediment availability versus use proportions
dat_sediments <- table(background_sediments$seabed_sub)/sum(table(background_sediments$seabed_sub))
dat_sediments <- data.frame(sediment = factor(names(dat_sediments)),
                            availability = as.numeric(dat_sediments))

dat_sediments_core <- table(core_sediments$seabed_sub)/sum(table(core_sediments$seabed_sub))
dat_sediments$core <- as.numeric(dat_sediments_core)[match(as.character(dat_sediments$sediment),
                                                           names(dat_sediments_core))]
dat_sediments$core[is.na(dat_sediments$core)] <- 0
dat_sediments <- dat_sediments %>% dplyr::arrange(availability)
# Plot differences in sediment availability versus use
pp <- par(mfrow = c(1, 2))
barplot(dat_sediments$availability)
barplot(dat_sediments$core)
par(pp)


#### End of code.
######################################
######################################
