######################################
######################################
#### examine_space_use_time_trials.R

#### This code:
# 1) Uses a trial-and-error approach to examine ways to decrease the wall time
# ... for pf() in examine_space use.

#### Steps preceding this code:
# 1) Run initial set-up code for pf() in examine_space_use, to just before the implementation of pf().


######################################
######################################
#### Exploration

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

#### End of code.
######################################
######################################
