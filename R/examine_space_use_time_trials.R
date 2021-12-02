######################################
######################################
#### examine_space_use_time_trials.R

#### This code:
# 1) Uses a trial-and-error approach to examine ways to decrease the wall time
# ... for pf(), pf_simplify() and pf_kud_2() in examine_space use.R

#### Steps preceding this code:
# 1) For the 'pf() exploration' part of the code
# ... run initial set-up code for pf() in examine_space_use.R
# ... to just before the implementation of pf().
# 2) For the 'pf_simplify()' part of the code
# ... continue running the remaining code in examine_space_use.R
# ... up to just before the section where pf_simplify() is implemented
# ... ('Define particles for mapping')
# 3) For the 'kud_around_coastline() part of the code,
# ... continue running the remaining code in examine_space_use.R
# ... up to just before the section where pf_kud_2() is implemented
# ... ('Fit KUDs')


######################################
######################################
#### pf() exploration

run_pf_exploration <- FALSE
if(run_pf_exploration){

   #### Profile code
   pf(record = out_ac_record[1:5],
      data = acpf_data[1:5, ],
      bathy = site_bathy,
      calc_movement_pr = calc_mpr,
      # mobility = mobility,
      n = n_particles, # 10L
      write_history = list(file = "./data/movement/space_use/acpf/history/"),
      con = "./data/movement/space_use/acpf/acpf_log.txt")

   #### Speed trial 1 with mobility
   system.time(pf(record = out_ac_record[1:5],
                  data = acpf_data[1:5, ],
                  bathy = site_bathy,
                  calc_movement_pr = calc_mpr,
                  mobility = mobility,
                  n = n_particles, # 10L
                  write_history = list(file = "./data/movement/space_use/acpf/history/"),
                  con = "./data/movement/space_use/acpf/acpf_log.txt")
   )
   # user  system elapsed
   # 161.053   7.663 169.632
   ((169.632/5) * 21158)/(60*60*24) # Expected code duration: 8 days

   #### Speed trial 2 without mobility: MUCH FASTER
   system.time(pf(record = out_ac_record[1:5],
                  data = acpf_data[1:5, ],
                  bathy = site_bathy,
                  calc_movement_pr = calc_mpr,
                  # mobility = mobility,
                  n = n_particles, # 10L
                  write_history = list(file = "./data/movement/space_use/acpf/history/"),
                  con = "./data/movement/space_use/acpf/acpf_log.txt")
   )
   # user  system elapsed
   # 53.324   1.674  55.149
   ((50.297/5) * 21158)/(60*60*24) # Expected code duration: 2.5 days

   #### Time trial 3 without verbose: NO DIFFERENCE
   system.time(pf(record = out_ac_record[1:5],
                  data = acpf_data[1:5, ],
                  bathy = site_bathy,
                  calc_movement_pr = calc_mpr,
                  # mobility = mobility,
                  n = n_particles, # 10L
                  write_history = list(file = "./data/movement/space_use/acpf/history/"),
                  verbose = FALSE)
   )
   # user  system elapsed
   # 54.078   1.595  55.915

   #### Time trial 4 adjusting raster memory

   ## Doubling memory is slower:
   raster::rasterOptions(maxmemory = 10e9)
   system.time(pf(record = out_ac_record[1:5],
                  data = acpf_data[1:5, ],
                  bathy = site_bathy,
                  calc_movement_pr = calc_mpr,
                  # mobility = mobility,
                  n = n_particles, # 10L
                  write_history = list(file = "./data/movement/space_use/acpf/history/"),
                  verbose = FALSE)
   )
   raster::rasterOptions(rop)
   # user  system elapsed
   # 67.000   4.612  70.938

   ## Doubling the memory and reducing chunk size is similar
   raster::rasterOptions(chunksize = 1e5, maxmemory = 10e9)
   system.time(pf(record = out_ac_record[1:5],
                  data = acpf_data[1:5, ],
                  bathy = site_bathy,
                  calc_movement_pr = calc_mpr,
                  # mobility = mobility,
                  n = n_particles, # 10L
                  write_history = list(file = "./data/movement/space_use/acpf/history/"),
                  verbose = FALSE)
   )
   raster::rasterOptions(rop)
   # user  system elapsed
   # 55.189   2.569  58.640

   ## Increasing chunk size alone is similar
   raster::rasterOptions(chunksize = 1e+09)
   system.time(pf(record = out_ac_record[1:5],
                  data = acpf_data[1:5, ],
                  bathy = site_bathy,
                  calc_movement_pr = calc_mpr,
                  # mobility = mobility,
                  n = n_particles, # 10L
                  write_history = list(file = "./data/movement/space_use/acpf/history/"),
                  verbose = FALSE)
   )
   raster::rasterOptions(rop)

   #### Time trial 5 with GRASS: MARGINALLY FASTER
   sink("./data/movement/space_use/acpf/rgass_log.txt")
   pf_opts <- pf_setup_optimisers(use_calc_distance_euclid_backend_grass = TRUE,
                                  use_grass_dir = "/Applications/GRASS-7.4.4.app/Contents/Resources")
   system.time(pf(record = out_ac_record[1:5],
                  data = acpf_data[1:5, ],
                  bathy = site_bathy,
                  calc_movement_pr = calc_mpr,
                  # mobility = mobility,
                  n = n_particles, # 10L
                  write_history = list(file = "./data/movement/space_use/acpf/history/"),
                  con = "./data/movement/space_use/acpf/acpf_log.txt",
                  optimisers = pf_opts)
   )
   sink()
   # Time record without sink:
   # user  system elapsed
   # 19.863   4.552  25.562
   # Time record with sink()
   # user  system elapsed
   # 19.890   4.557  25.565
   ((39.437/5) * 21158)/(60*60*24) # Expected code duration: 2.7 days

   #### Time trail 6 with use_raster_options = FALSE: NO DIFFERENCE
   sink("./data/movement/space_use/acpf/rgass_log.txt")
   pf_opts <- pf_setup_optimisers(use_calc_distance_euclid_backend_grass = TRUE,
                                  use_grass_dir = "/Applications/GRASS-7.4.4.app/Contents/Resources",
                                  use_raster_operations = FALSE)
   system.time(pf(record = out_ac_record[1:5],
                  data = acpf_data[1:5, ],
                  bathy = site_bathy,
                  calc_movement_pr = calc_mpr,
                  # mobility = mobility,
                  n = n_particles, # 10L
                  write_history = list(file = "./data/movement/space_use/acpf/history/"),
                  con = "./data/movement/space_use/acpf/acpf_log.txt",
                  optimisers = pf_opts)
   )
   sink()
   # user  system elapsed
   # 19.418   4.620  25.173

   #### Time trials with fewer particles
   ## Time trial with default options
   system.time(pf(record = out_ac_record[1:5],
                  data = acpf_data[1:5, ],
                  bathy = site_bathy,
                  calc_movement_pr = calc_mpr,
                  # mobility = mobility,
                  n = 100L,
                  write_history = list(file = "./data/movement/space_use/acpf/history/"),
                  con = "./data/movement/space_use/acpf/acpf_log.txt")
   )
   # user  system elapsed
   # 10.169   1.344  14.520
   ((14.520/5) * 21158)/(60*60*24) # Expected code duration: 0.7 days
   ## Time trial with GRASS: SLOWER with fewer particles
   # ... So this approach is slower at small numbers of particles but seems to scale
   # ... more effectively with large numbers of particles
   sink("./data/movement/space_use/acpf/rgass_log.txt")
   pf_opts <- pf_setup_optimisers(use_calc_distance_euclid_backend_grass = TRUE,
                                  use_grass_dir = "/Applications/GRASS-7.4.4.app/Contents/Resources")
   system.time(pf(record = out_ac_record[1:5],
                  data = acpf_data[1:5, ],
                  bathy = site_bathy,
                  calc_movement_pr = calc_mpr,
                  # mobility = mobility,
                  n = 100L,
                  write_history = list(file = "./data/movement/space_use/acpf/history/"),
                  con = "./data/movement/space_use/acpf/acpf_log.txt",
                  optimisers = pf_opts)
   )
   sink()
   # user  system elapsed
   # 20.170   4.438  27.852

}


######################################
######################################
#### pf_simplify() exploration

run_pf_simplify_exploration <- FALSE
if(run_pf_simplify_exploration){

   #### Define testing data
   out_acdcpf_testr         <- out_acdcpf
   out_acdcpf_testr$history <- out_acdcpf_testr$history[1:100]
   message(length(out_acdcpf_testr$history))

   #### Test the duration of Euclidean distance calculations on the first and second time steps [v. fast]
   tmp <- flapper:::.planedist2(raster::xyFromCell(site_bathy, out_acdcpf_testr$history[[1]]$id_current),
                                raster::xyFromCell(site_bathy, out_acdcpf_testr$history[[2]]$id_current))
   tmp <- data.frame(id_previous = out_acdcpf_testr$history[[1]]$id_current[as.vector(row(tmp))],
                     id_current = out_acdcpf_testr$history[[2]]$id_current[as.vector(col(tmp))],
                     dist_current = as.vector(tmp))
   tmp <- tmp %>% dplyr::filter(.data$dist_current <= mobility)
   nrow(tmp)

   #### Test the duration of segment-cross-barrier calculations
   ## Define index for barrier calculations
   bathy_xy <- raster::coordinates(site_bathy)
   index_for_barrier <- which(tmp$dist_current > calc_distance_barrier_limit & tmp$dist_current <= calc_distance_limit)
   length(index_for_barrier)
   ## Scenario (A): all points are away from the barrier [v. fast]
   system.time(
      {
         index_in_barrier <-
            which(segments_cross_barrier(bathy_xy[tmp$id_previous[index_for_barrier], , drop = FALSE],
                                         bathy_xy[tmp$id_current[index_for_barrier], , drop = FALSE],
                                         site_coast))
         index_in_barrier <- index_for_barrier[index_in_barrier]
      }
   )
   ## Scenario (B): points within barrier [~one minute]
   # We can check the speed of internal routines when cell pairs need to be checked for intersections
   # ... by extracting the key bits of the internal code from segments_cross_barrier()).
   # Depending on (i) how many LCP calculations this prevents and (ii) the duration of LCP algorithms,
   # ... it may be better to ignore this option.
   # In this example, the implementation of barrier-overlap routines would prevent
   # ... 75,518 LCP calculations, but at a cost of 49 s. But this will often be faster
   # ... due to the incorporation of the 'distance' argument in segments_cross_barrier().
   system.time(
      {
         # Define variables
         start   <- bathy_xy[tmp$id_previous[index_for_barrier], , drop = FALSE]
         end     <- bathy_xy[tmp$id_current[index_for_barrier], , drop = FALSE]
         barrier <- site_coast
         # Define point matrices as dataframes
         start <- data.frame(start[, 1:2])
         end   <- data.frame(end[, 1:2])
         colnames(start) <- colnames(end) <- c("x", "y")
         # Assign line IDs
         start$linestring_id <- end$linestring_id <- seq_len(nrow(start))
         # Define lines
         lines <-
            dplyr::bind_rows(start, end) %>%
            dplyr::arrange(.data$linestring_id)
         lines <- sfheaders::sf_linestring(lines,
                                           x = "x", y = "y",
                                           linestring_id = "linestring_id") # fast step
         # Get intersection
         int <- sf::st_intersects(lines, barrier, sparse = FALSE) # slow step
         print(table(int))
      }
   )
   # For 430,958 rows [nrow(tmp)], necessitating 75,518 intersections [length(index_for_barrier)]:
   # user  system elapsed
   # 45.706   0.342  46.208

   #### Examine the number of time steps for which point pairs COULD cross a barrier for ACDCPF
   # There are 3252 time steps for which at least some intersections are required in
   # ... flapper::segments_cross_barrier(), which can be slow (see above). But the number
   # ... of intersections at each time step may be relatively low.
   run <- FALSE
   if(run){
      # [Duration: ~ 5 minutes]
      # [This code was adapted from flapper::segments_cross_barrier()]
      ## Define barrier
      barrier <- sf::st_as_sf(dat_coast)
      sf::st_crs(barrier) <- NA
      # Set up cluster
      cl <- parallel::makeCluster(10L)
      parallel::clusterExport(cl, c("site_bathy", "out_acdcpf", "barrier", "site_coast_distances", "mobility"))
      # Loop over pairs of time steps and work out whether or not the bounding polygon
      # ... around the points intersects with the coastline. If yes, then
      # ... we may need to examine potential intersections more closely. This provides
      # ... an indication of the extra time potentially required during particle processing.
      out_acdcpf_barrier_chk <-
         pbapply::pblapply(1:(length(out_acdcpf$history)-1), cl = cl, function(t){
            ## Define coordinates for sampled cells at the current and next time step
            # t = 1
            ht1 <- out_acdcpf$history[[t]]
            ht2 <- out_acdcpf$history[[t+1]]
            # pairs <- expand.grid(id_previous = ht1$id_current, id_current = ht2$id_current)
            # start <- raster::xyFromCell(site_bathy, pairs$id_previous)
            # end   <- raster::xyFromCell(site_bathy, pairs$id_current)
            start <- raster::xyFromCell(site_bathy, ht1$id_current)
            end   <- raster::xyFromCell(site_bathy, ht2$id_current)
            ## Define the bounding box around those coordinates
            xlim <- range(c(start[, 1], end[, 1]))
            ylim <- range(c(start[, 2], end[, 2]))
            boundaries <- matrix(c(xlim[1], ylim[1],
                                   xlim[2], ylim[1],
                                   xlim[2], ylim[2],
                                   xlim[1], ylim[2]), ncol = 2, byrow = TRUE)
            boundaries <- rbind(boundaries, boundaries[1, , drop = FALSE])
            boundaries <- sf::st_polygon(list(boundaries))
            ## Identify whether or not the boundary box intersects with the barrier
            int <- sf::st_intersects(boundaries, barrier, sparse = FALSE)
            if(!int){
               return(int)
            } else {
               pairs <- expand.grid(id_previous = ht1$id_current, id_current = ht2$id_current)
               dat <- data.frame(dist_1 = raster::extract(site_coast_distances, pairs$id_previous),
                                 dist_2 = raster::extract(site_coast_distances, pairs$id_current))
               return(matrix(any((dat$dist_1 < mobility) | (dat$dist_2 < mobility)), ncol = 1))
            }
         })
      out_acdcpf_barrier_chk <- unlist(out_acdcpf_barrier_chk)
      saveRDS(out_acdcpf_barrier_chk, "./data/movement/space_use/acdcpf/out_acdcpf_barrier_chk.rds")
   } else {
      out_acdcpf_barrier_chk <-
         readRDS("./data/movement/space_use/acdcpf/out_acdcpf_barrier_chk.rds")
   }
   # The number of occasions when we would need to examine barrier overlaps for cell pairs:
   table(out_acdcpf_barrier_chk)
   # FALSE  TRUE
   # 17905  3252

   #### Test the LCP algorithm speeds, including
   # ... alternative algorithms: "Dijkstra", "bi", "A*", "NBA"
   # ... alternative graphs (default, simplified, simplified + contracted)
   # ... parallelisation versus non parallelisation
   # Define dataframe with alternative options
   # ... This includes 'graph', 'algorithm' and 'allcores', which can be changed.
   # ... 'n' is the number of location pairs for which these test LCP calculations will be made
   # ... 'time' will record the time for each combination.
   dat_lcp_op <-
      expand.grid(
         graph = c("default", "simplified"), # c("default", "simplified", "contracted"),
         algorithm = c("Dijkstra", "bi", "A*", "NBA"),
         allcores = c(FALSE, TRUE),
         n = 1000L,
         time = NA
      )
   # Determine the duration of each option
   dat_lcp_by_op <- lapply(split(dat_lcp_op, 1:nrow(dat_lcp_op)), function(op){
      if(op$graph == "default"){
         gph <- site_bathy_graph
      } else if(op$graph == "simplified"){
         gph <- site_bathy_graph_simp
      } else if(op$graph == "contracted"){
         gph <- site_bathy_graph_cont
      }
      t1 <- Sys.time()
      lcp_for_op <- cppRouting::get_distance_pair(Graph = gph,
                                                  from = tmp$id_previous[seq_len(op$n)],
                                                  to = tmp$id_current[seq_len(op$n)],
                                                  algorithm = op$algorithm,
                                                  allcores = op$allcores)
      t2 <- Sys.time()
      t3 <- as.numeric(difftime(t2, t1, units = "mins"))
      op$time <- t3
      return(op)
   })
   # Examine the results
   dat_lcp_by_op <-
      dat_lcp_by_op %>%
      dplyr::bind_rows() %>%
      dplyr::tibble() %>%
      dplyr::arrange(.data$time)

   #### Estimate pf_simplify() wall times for different implementations
   # --> specifically, test the benefits of parallelising time steps relative to LCP calculations (above)

   # A) Baseline (Euclidean distances)
   # ~ 3.2 s for 3 time steps
   # ~ 15.1 s for 10 time steps
   # ~ 120 s for 100 time steps
   t1_a <- Sys.time()
   out_acdcpf_testr_s_b <- pf_simplify(out_acdcpf_testr,
                                       write_history = list(file = "./data/movement/space_use/acdcpf/processing/testr/"),
                                       return = "archive")

   t2_a <- Sys.time()
   difftime(t2_a, t1_a)

   # B) Baseline (shortest distances)
   # ~ 6..4 for 3 time steps
   # ~30 s for 10 time steps
   # ~327 s for 100 time steps
   t1_b <- Sys.time()
   out_acdcpf_testr_s_b <- pf_simplify(out_acdcpf_testr,
                                       calc_distance = "lcp",
                                       calc_distance_graph = site_bathy_graph,
                                       calc_distance_limit = euclid_distance_limit,
                                       # mobility = 481,
                                       calc_distance_barrier = site_barrier,
                                       calc_distance_barrier_limit = euclid_distance_barrier_limit,
                                       calc_distance_barrier_grid = site_coast_grid,
                                       calc_distance_restrict = TRUE,
                                       calc_distance_algorithm = "Dijkstra",
                                       write_history = list(file = "./data/movement/space_use/acdcpf/processing/testr/"),
                                       return = "archive")
   t2_b <- Sys.time()
   difftime(t2_b, t1_b)

   # C) Baseline shortest distances with reduced mobility
   # ~6.4 s for 3 time steps
   # ~ 37.7 s for 10 time steps
   t1_c <- Sys.time()
   out_acdcpf_testr_s_c <- pf_simplify(out_acdcpf_testr,
                                       calc_distance = "lcp",
                                       calc_distance_graph = site_bathy_graph,
                                       # calc_distance_limit = euclid_distance_limit,
                                       mobility = euclid_distance_limit,
                                       calc_distance_barrier = site_barrier,
                                       # calc_distance_barrier_limit = euclid_distance_barrier_limit,
                                       calc_distance_barrier_grid = site_coast_grid,
                                       calc_distance_restrict = TRUE,
                                       calc_distance_algorithm = "Dijkstra",
                                       write_history = list(file = "./data/movement/space_use/acdcpf/processing/testr/"),
                                       return = "archive")
   t2_c <- Sys.time()
   difftime(t2_c, t1_c)

   # D) Parallelise time steps (using the fastest one-core LCP calculation method)
   # With pre-chunking method (parallelising over time steps) & 10 cores:
   # ... [10 time steps]  45.7 s (FORK) versus 94.8 s (SOCKET)
   # ... [100 time steps] 145 s  (FORK) versus 840 s  (SOCKET)
   # With chunking method (parallelising over chunks, looping over time steps in serial within each chunk):
   # ... [10 time steps]  31 s versus 93 s (SOCKET)
   # ... [100 time steps] 84.5 s (FORK) (or 88.7 s if euclid_distance_limit supplied)
   # With chunking method and fewer cores
   # ... [100 time steps] 94 s (FORK) (cl = 5L)
   # ... [100 time steps] 151 s (FORK) (cl = 2L)
   # With chunking method and fewer (5L) cores and no progress bar
   # ... [100 time steps] 88 s (FORK)
   pb_op <- pbapply::pboptions(type = "none")
   t1_d <- Sys.time()
   out_acdcpf_testr_s_d <- pf_simplify(out_acdcpf_testr,
                                       calc_distance = "lcp",
                                       calc_distance_graph = site_bathy_graph,
                                       calc_distance_limit = euclid_distance_limit, # slightly slower but required
                                       calc_distance_barrier = site_barrier,
                                       calc_distance_barrier_limit = euclid_distance_barrier_limit,
                                       calc_distance_restrict = TRUE,
                                       calc_distance_barrier_grid = site_coast_grid,
                                       calc_distance_algorithm = "Dijkstra",
                                       write_history = list(file = "./data/movement/space_use/acdcpf/processing/testr/"),
                                       cl = 10L,
                                       # cl = parallel::makeCluster(10L), varlist = "mobility",
                                       return = "archive")
   t2_d <- Sys.time()
   difftime(t2_d, t1_d)
   pbapply::pboptions(pb_op)

   # E) Parallelise LCP calculations (using the fastest multi-core method)
   t1_e <- Sys.time()
   out_acdcpf_testr_s_e <- pf_simplify(out_acdcpf_testr,
                                     calc_distance = "lcp",
                                     calc_distance_graph = site_bathy_graph,
                                     calc_distance_barrier = site_barrier,
                                     calc_distance_barrier_limit = euclid_distance_barrier_limit,
                                     calc_distance_barrier_grid = site_coast_grid,
                                     calc_distance_restrict = TRUE,
                                     calc_distance_algorithm = "Dijkstra",
                                     write_history = list(file = "./data/movement/space_use/acdcpf/processing/testr/"),
                                     use_all_cores = TRUE,
                                     return = "archive")
   t2_e <- Sys.time()
   difftime(t2_e, t1_e)

}


######################################
######################################
#### kud_around_coastline() exploration

run_pf_simplify_exploration <- FALSE
if(run_pf_simplify_exploration){

   #### Load particle samples
   out_acdcpf_s <- readRDS("./data/movement/space_use/acdcpf/out_acdcpf_s.rds")
   site_bathy   <- raster::raster("./data/spatial/site_bathy.tif")

   #### Examine the duration of POU maps [38 s]
   t1 <- Sys.time()
   out_acdcpf_pou <- pf_plot_map(out_acdcpf_s, bathy)
   t2 <- Sys.time()
   difftime(t2, t1)

   #### Define site_habitat for KUD estimation
   # The resolution of site_habitat significantly impacts the speed
   # ... of KUD estimation, as shown below. Here, we reduce the resolution of
   # ... bathy to define a 20 x 20 m grid.
   site_grid <- raster::aggregate(site_bathy, fact = 4, fun = mean, na.rm = TRUE)
   raster::res(site_grid)
   raster::plot(site_grid)
   site_habitat <- kud_habitat(site_grid)

   #### Define objects to test the duration kud_around_coastline()
   # We will use internal code in pf_kud_2() to estimate the duration of
   # ... kud_around_coastline() given the resolution of site_habitat and
   # ... the number of particles sampled, so we will define objects
   # ... required by that code here.
   xpf         <- out_acdcpf_s
   bathy       <- site_bathy
   crs         <- raster::crs(bathy)
   estimate_ud <-  kud_around_coastline
   grid        <- site_habitat

   #### Get cell probabilities and coordinates
   # Processed particles for each time step
   particles_by_t <-
      lapply(1:length(xpf$history), function(t) {
         elm <- xpf$history[[t]]
         history_for_t <- elm[, c("id_current", "pr_current"), drop = FALSE]
         colnames(history_for_t) <- c("cell_id", "cell_pr")
         history_for_t$timestep <- t
         return(history_for_t)
      })
   # A dataframe with all particles
   p <- do.call(rbind, particles_by_t)
   p[, c("cell_x", "cell_y")] <- raster::xyFromCell(bathy, p$cell_id)
   nrow(p) # 19328309

   #### Estimate the duration of KUD estimation for different (small) particle samples
   n_seq <- c(10, 100, 1000, 10000, 20000, 50000, 100000)
   kud_duration_by_n <-
      pbapply::pblapply(n_seq, function(sample_size){
         ## Get the full set of particles
         particles <- p
         ## Sample locations according to their probability
         if(!is.null(sample_size)){
            particles <-
               particles %>%
               dplyr::slice_sample(n = sample_size, weight_by = .data$cell_pr, replace = TRUE)
         }
         ## Make SPDF
         # cat_to_console("... Implementing KUD estimation...")
         particles_spdf <- sp::SpatialPointsDataFrame(
            particles[, c("cell_x", "cell_y")],
            data = data.frame(ID = factor(rep(1, nrow(particles)))),
            proj4string = crs)
         ## Estimate the time of KUD estimation
         duration <- system.time(ud <- estimate_ud(xy = particles_spdf, grid = grid))
         raster::plot(raster::raster(ud[[1]]))
         return(duration)
      })

   #### Define a dataframe of estimation durations
   kud_duration <- data.frame(n = n_seq,
                              time = sapply(kud_duration_by_n, function(x) x[3]))

   #### Visualise estimation duration with increasing numbers of particles
   plot(kud_duration$n, kud_duration$time)
   mod <- lm(time ~ n, data = kud_duration)
   abline(mod)

   #### Estimate the duration (hours) of KUD estimation
   # ... for all particles at each trialled grid resolution
   predict(mod, data.frame(n = nrow(p)))/60/60
   # ~1.2 hours with 50 x 50 m grid
   # ~4.3 hours with 25 x 25 m grid
   # ~6.8 hours with 20 x 20 m grid --> selected.
   # ~27 hours with 10 x 10 m grid
}


#### End of code.
######################################
######################################
