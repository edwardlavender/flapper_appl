######################################
######################################
#### examine_habitat_preferences.R

#### This code:
# 1) Examines the implications of different maps of space use
# ... for an example flapper skate, as derived from difference algorithms,
# ... for 'habitat preferences', focusing on sediment use.
# ... A simple, visual, method is used to compare the distribution of available
# ... sediment types (within the study area) to the distribution of 'used'
# ... sediment types found within the core range (50 % contour) of the UD
# ... derived from each algorithm.

#### Steps preceding this code:
# 1) Definition of maps of depth/space use from:
# ... examine_depth_use.R
# ... examine_space_use.R


######################################
######################################
#### Set up

#### Wipe workspace & load flapper & global param
rm(list = ls())
source("./R/define_global_param.R")

#### Load data
out_dc_map         <- raster::raster("./data/movement/depth_use/dc/out_dc_map.tif")
out_ac_map         <- raster::raster("./data/movement/space_use/ac/out_ac_map.tif")
out_acdc_map       <- raster::raster("./data/movement/space_use/acdc/out_acdc_map.tif")
out_acpf_kud       <- raster::raster("./data/movement/space_use/acpf/out_acpf_kud.tif")
out_acdcpf_kud     <- raster::raster("./data/movement/space_use/acdcpf/out_acdcpf_kud.tif")
out_coa_kud        <- raster::raster("./data/movement/space_use/coa/out_coa_kud.tif")
site_bathy         <- raster::raster("./data/spatial/site_bathy.tif")
site_coast         <- readRDS("./data/spatial/site_coast.rds")
site_sediments     <- readRDS("./data/spatial/site_sediments.rds")
site_sediments_bkg <- readRDS("./data/spatial/site_sediments_background.rds")


######################################
######################################
#### Map processing

#### Scale AC/DC/ACDC maps for comparability with UDs
## DC
out_dc_map <- out_dc_map/raster::cellStats(out_dc_map, "sum")
raster::cellStats(out_dc_map, "sum")
## AC
out_ac_map <- out_ac_map/raster::cellStats(out_ac_map, "sum")
raster::cellStats(out_ac_map, "sum")
## ACDC
out_acdc_map <- out_acdc_map/raster::cellStats(out_acdc_map, "sum")
raster::cellStats(out_acdc_map, "sum")

#### Scale KUD maps
raster::cellStats(out_coa_kud, "sum")
raster::cellStats(out_acpf_kud, "sum")
raster::cellStats(out_acdcpf_kud, "sum")

#### Confirm raster properties are identical
raster::compareRaster(out_dc_map, out_ac_map, out_acdc_map,
                      out_coa_kud, out_acpf_kud, out_acdcpf_kud)

#### Define maps list (select all maps or KUD maps)
## All maps
maps <- list(out_dc_map = out_dc_map,
             out_ac_map = out_ac_map,
             out_acdc_map = out_acdc_map,
             out_coa_kud = out_coa_kud,
             out_acpf_kud = out_acpf_kud,
             out_acdcpf_kud = out_acdcpf_kud)
## KUD maps
maps <- list(out_coa_kud = out_coa_kud,
             out_acpf_kud = out_acpf_kud,
             out_acdcpf_kud = out_acdcpf_kud)
maps_names <- names(maps)

#### Force zero/NA consistency (if necessary)
# NA: 'impossible' areas (land)
# 0 : 'unused' areas
maps <- cl_lapply(maps, function(map){
  map[is.na(map)] <- 0
  map <- raster::mask(map, site_bathy)
  return(map)
})
names(maps) <- maps_names

#### Get core ranges from each algorithm
## Notes
# a) The KUD maps are scaled below to have a maximum value of one.
# ... This means we can use the same scale bar for all plots. These needs to be
# ... added manually afterwards.
# b) The core-range maps are derived from get_hr_core().
# ... Within this function, zero cells (i.e., unused areas) are
# ... ignored by default, so the home ranges are derived based on quantiles
# ... in 'used' areas.
## Set up figure to save
png("./fig/habitat_preferences/core_ranges.png",
    height = 8, width = 6, units = "in", res = 600)
pp <- par(mfcol = c(length(maps), 2),
          oma = c(1, 1, 1, 6), mar = rep(0.5, 4),
          xaxs = "i", yaxs = "i")
## Define raster plotting options
add_map <- list(plot_method = raster::plot, legend = FALSE)
## Loop over each map and plot the KUD and the core range
core_hr_by_map <- cl_lapply(1:length(maps), fun = function(i){
  # Define tidy map titles
  map           <- maps[[i]]
  alg           <- names(maps)[i]
  alg_names     <- data.frame(name = c("(COA)", "(ACPF)", "(ACDCPF)"),
                              code = c("out_coa_kud", "out_acpf_kud", "out_acdcpf_kud"))
  alg           <- alg_names$name[match(alg, alg_names$code)]
  paa_for_map_1 <- paa_for_map_2 <- paa
  title_1       <- LETTERS[i*2-1]
  title_2       <- LETTERS[i*2]
  title_1       <- bquote(bold(.(title_1)) ~ .(alg))
  title_2       <- bquote(bold(.(title_2)) ~ .(alg))
  # Option to add x and y axis labels
  use_detailed_axes <- FALSE
  if(use_detailed_axes){
    if(title_1 %in% c("A", "C", "E")){
      paa_for_map_1 <- list(side = 1:4,
                            axis = list(list(labels = FALSE),
                                        list(),
                                        list(labels = FALSE),
                                        list(labels = FALSE)),
                            control_axis = list(tck = -0.01, las = TRUE),
                            control_sci_notation = list(magnitude = 16L, digits = 0))
      if(title_1 == "E") paa_for_map_1$axis[[1]]$labels <- NULL
    }
    if(title_2 == "F") {
      paa_for_map_2 <- list(side = 1:4,
                            axis = list(list(),
                                        list(labels = FALSE),
                                        list(labels = FALSE),
                                        list(labels = FALSE)),
                            control_axis = list(tck = -0.01, las = TRUE),
                            control_sci_notation = list(magnitude = 16L, digits = 0))
    }
  }
  # Process KUD map onto common scale, with a maximum value of one
  map_scaled_to_one <- map/raster::cellStats(map, "max")
  print(raster::cellStats(map_scaled_to_one, "max"))
  # Plot KUD map
  prettyGraphics::pretty_map(add_rasters =
                               rlist::list.merge(list(x = map_scaled_to_one),
                                                 add_map,
                                                 list(zlim = c(0, 1))),
                             add_polys = add_coast,
                             pretty_axis_args = paa_for_map_1)
  add_map_elements()
  mtext(side = 3, title_1, cex = 1.25, font = 2, adj = 0.025, line = -2.25)
  # Plot core range
  get_hr_core(map,
              add_raster = add_map,
              add_polys = add_coast,
              pretty_axis_args = paa_for_map_2)
  add_map_elements()
  mtext(side = 3, title_2, cex = 1.25, font = 2, adj = 0.025, line = -2.25)
})
add_map_elements()
names(core_hr_by_map) <- names(maps)
par(pp)
dev.off()

## Save legend
png("./fig/habitat_preferences/core_ranges_legend.png",
    height = 8, width = 3, units = "in", res = 600)
# pp <- par(oma = c(2, 2, 2, 2))
fields::image.plot(zlim = c(0, 1),
                   col = rev(terrain.colors(255)),
                   axis.args = list(cex.axis = 1.25),
                   legend.only = TRUE)
mtext(side = 4, "Intensity", line = 1, cex = 1.25)
# par(pp)
dev.off()


######################################
######################################
#### Sediments analysis

#### Summarise sediment 'availability'
dat_sediments_bkg <- table(site_sediments_bkg$seabed_sub)/sum(table(site_sediments_bkg$seabed_sub))
dat_sediments_bkg <- data.frame(sediment = factor(names(dat_sediments_bkg)),
                                availability = as.numeric(dat_sediments_bkg))

#### Get sediment 'use', according to each algorithm [~ 15 minutes]
run_sediments <- FALSE
if(run_sediments){
  core_sediments_by_map <-
    cl_lapply(1:length(core_hr_by_map), fun = function(i){
      # i = 1
      name <- names(core_hr_by_map)[i]
      map  <- core_hr_by_map[[i]]
      con_1 <- paste0("./data/movement/habitat_preferences/sediments/", name, ".rds")
      if(!file.exists(con_1)){
        core_xy <- raster::rasterToPoints(map, fun = function(x) x==1)
        core_xy <- sp::SpatialPoints(core_xy[, 1:2], proj4string = proj_utm)
        core_sediments <- sp::over(core_xy, site_sediments)
        saveRDS(core_sediments, con_1)
      }
      con_2 <- paste0("./data/movement/habitat_preferences/sediments/", name, "_summary.rds")
      if(!file.exists(con_2)){
        dat_sediments <- dat_sediments_bkg
        dat_sediments_core <- table(core_sediments$seabed_sub)/sum(table(core_sediments$seabed_sub))
        dat_sediments$core <- as.numeric(dat_sediments_core)[match(as.character(dat_sediments$sediment),
                                                                   names(dat_sediments_core))]
        dat_sediments$core[is.na(dat_sediments$core)] <- 0
        saveRDS(dat_sediments, con_2)
      } else {
        dat_sediments <- readRDS(con_2)
      }
      return(dat_sediments[, "core", drop = FALSE])
    })
} else {
  core_sediments_by_map <-
    cl_lapply(names(core_hr_by_map), fun = function(file){
      readRDS(paste0("./data/movement/habitat_preferences/sediments/", file, "_summary.rds"))
    })
}
core_sediments_by_map <- lapply(core_sediments_by_map, function(d) d[, "core", drop = FALSE])
names(core_sediments_by_map) <- names(maps)

#### Process 'dat_sediments' for visualisation
dat_sediments <-
  list(dat_sediments_bkg, core_sediments_by_map) %>%
  dplyr::bind_cols() %>%
  dplyr::arrange(availability)
colnames(dat_sediments)[3:ncol(dat_sediments)] <- names(maps)


######################################
######################################
#### Visualisation

#### Define matrix for barplot
# Define matrix
rownames(dat_sediments) <- dat_sediments$sediment
dat_sediments_mat <- as.matrix(dat_sediments[2:ncol(dat_sediments)])
dat_sediments_mat <- t(dat_sediments_mat)

#### Tidy row and column names for axis/legend labels
## Algorithm names
# Names for all or selected COA maps
if(nrow(dat_sediments_mat) == 7L){
  rownames(dat_sediments_mat) <-
    c("AV", "DC", "AC", "ACDC",
      "COA", "ACPF", "ACDCPF")
} else if(nrow(dat_sediments_mat) == 4L){
  rownames(dat_sediments_mat) <-
    c("AV", "COA", "ACPF", "ACDCPF")
} else {
  message("Unable to assign rownames.")
}
## Sediment categories
# Define basic labels
xlab <- colnames(dat_sediments_mat)
# Use sentence case
xlab <-
  stringr::str_to_sentence(xlab)
# Wrap the longer labels
xlab <-
  stringr::str_replace(xlab, " and ", "\n & ")
# Hack to ensure even spacing from the y axis for all labels
xlab[!stringr::str_detect(xlab, "&")] <-
  paste0(xlab[!stringr::str_detect(xlab, "&")], "\n", " ")
xlab

#### Define graphical param,
# Colour for each algorithm
gp <- data.frame(legend = rownames(dat_sediments_mat))
                  # col = c("black", "dimgrey", "grey", "lightgrey"))
# Bar width and spacing
width <- 1
space <- c(0, 1)

#### Define barplot
png("./fig/habitat_preferences/habitat_preferences_sediments.png",
    height = 5, width = 12, units = "in", res = 600)
pp <- par(oma = c(2, 2, 2, 2))
b <- barplot(dat_sediments_mat,
             # Define bar placement
             beside = TRUE, width = width, space = c(0, 1),
             # Define bar colours
             # col = gp$col, border = NA,
             # Suppress axes and labelling
             axes = FALSE, names.arg = rep("", nrow(dat_sediments)),
             # Add legend
             legend.text = rownames(dat_sediments_mat),
             # args.legend = list(x = 8, y = 0.4, bty = "n", border = NA)) # legend placement for all maps
             args.legend = list(x = 5, y = 0.4, bty = "n", border = NA))

#### Add axes
cex.axis <- 0.8
xat <- apply(b, 2, mean)
pm <- par(mgp = c(3, 1.5, 0))
axis(side = 1, at = c(0, max(xat) + width), labels = FALSE, lwd.ticks = 0, pos = 0)
axis(side = 1, at = xat, labels = xlab, pos = 0, cex.axis = cex.axis)
par(pm)
axis(side = 2, las = TRUE, pos = 0, cex.axis = cex.axis)

#### Add titles
mtext(side = 1, "Sediment", line = 3)
mtext(side = 2, "Proportion", line = 2.5)
dev.off()


#### End of code.
######################################
######################################
