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
archival   <- readRDS("./data/movement/tag/archival_eg.rds")
site_bathy <- raster::raster("./data/spatial/site_bathy.tif")
site_coast <- readRDS("./data/spatial/site_coast.rds")


######################################
######################################
#### Implement DC algorithm

#### Define algorithm param
n_cores <- 10L
split <- floor(nrow(archival)/n_cores)
cl <- parallel::makeCluster(n_cores)
parallel::clusterEvalQ(cl, library(raster))

#### Implement algorithm
out_dc <- dc(archival = archival[1:50, ],
             bathy = site_bathy,
             calc_depth_error = calc_depth_error,
             check_availability = TRUE,
             # write_record_spatial_for_pf = list(filename = "./data/movement/depth_use/dc/record/", format = "GTiff"),
             con = "./data/movement/depth_use/dc/dc_log.txt",
             split = split,
             cl = cl)

#### Process out_dc
out_dc_s <- acdc_simplify(out_dc, type = "dc", mask = site_bathy)

#### Check availability of depth contours at each time step
out_dc_dat <- acdc_access_dat(out_dc_s)
table(out_dc_dat$availability)

#### Process map (for plotting)
out_dc_map <- out_dc_s$out_dc_map
out_dc_map[is.na(out_dc_map)] <- 0

#### Plot map
prettyGraphics::pretty_map(add_rasters = list(x = out_dc_map),
                           add_polys = list(x = site_coast, col = "dimgrey", border = "dimgrey"))

#### End of code.
######################################
######################################
