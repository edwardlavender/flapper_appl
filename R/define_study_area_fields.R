######################################
######################################
#### define_study_area_fields.R

#### This code:
# 1) Defines study area spatial fields (e.g., study site boundaries)
# ... for mapping the study area via the QGIS projecct in ./fig/study_site/

#### Steps preceding this code:
# 1) Definition of spatial layers in preparatory and analytical scripts.


######################################
######################################
#### Set up

#### Wipe workspace & load flapper & global param
rm(list = ls())
source("./R/define_global_param.R")


######################################
######################################
#### Make fields

######################################
#### Bathymetry

#### Copy bathymetry file to QGIS project
copy_bathy <- FALSE
if(copy_bathy) file.copy("./data-raw/spatial/bathy/bathy_res_full_ext_full_abs.tif",
                         './fig/study_site/bathy.tif')


######################################
#### Study area boundaries

#### Get study area rasters
r_1 <- raster::raster("./data/spatial/site_bathy.tif")
r_2 <- raster::raster("./data/movement/post_release_paths/1507/site_bathy.tif")
r_3 <- raster::raster("./data/movement/post_release_paths/1558/site_bathy.tif")
r_4 <- raster::raster("./data/movement/cooccurrences/site_bathy.tif")

#### Examine areas
pp <- par(mfrow = c(2, 2))
raster::plot(r_1, main = "r_1")
raster::plot(r_2, main = "r_2")
raster::plot(r_3, main = "r_3")
raster::plot(r_4, main = "r_4")
par(pp)

#### Get boundaries of study area
ext_by_site <- lapply(list(r_1, r_2, r_3, r_4), function(r) raster::extent(r)[1:4])
ext <- do.call(rbind, ext_by_site)
range(ext[, 1:2])
# 698954.7 714429.6
range(ext[, 3:4])
# 6241369 6260064

#### Get study area polygons
poly_by_site <-
  lapply(1:length(ext_by_site), function(i){
    ext <- raster::extent(ext_by_site[[i]])
    poly <- as(ext, "SpatialPolygons")
    poly <- sp::SpatialPolygonsDataFrame(poly, data.frame(ID = 1L))
    raster::crs(poly) <- proj_utm
    rgdal::writeOGR(poly,
                    dsn = "./fig/study_site/",
                    layer = paste0("study_site_boundaries_", i),
                    driver ="ESRI Shapefile",
                    overwrite_layer = TRUE)
})


######################################
#### Tag deployment locations

#### Define selected individuals
# 540 (tag deployment)  [35]
# 1507 (tag deployment) [38]
# 1558 (tag deployment) [13]
# 542 and 560 (tag deployment locations) [29 and 28]

#### Get tag deployment locations
skateids <- readRDS("./data/movement/generic/skateids.rds")
tag_deployment_xy <-
  skateids %>%
  dplyr::filter(individual_id %in% c(35, 38, 13, 29, 28)) %>%
  dplyr::select(x = tag_x, y = tag_y)

#### Save locations
write.csv(tag_deployment_xy,
          "./fig/study_site/tag_deployment_xy.csv",
          quote = FALSE,
          row.names = FALSE)


#### End of code
######################################
######################################
