# script 1/3
# Dean Koch, Feb 2023
# Executes a NetLogo simulation of the territory model


# This is slow, takes about 2 hours on a single core (I haven't implemented
# parallel processing). At the end a uniquely named geotiff file is written
# to `proj_dir` containing the estimated territory sizes by pixel.

# The projection info is baked into to the geotiff but I think this info
# is wrong. See the process_netlogo_results.R script for details

# run this script multiple times to get multiple repetitions. There should
# be no issues with older files getting overwritten

library(RNetLogo)
library(terra)

# paths to: Sells et al supplement, directory for results, helper functions
supplement_zip = 'D:/wolves/refs/data_supplement/Model.zip'
proj_dir = 'D:/wolves/bias_simulations'
helper_path = 'D:/wolves/R/0_netlogo_helpers.R'

# this defines `run_netlogo()`
source(helper_path)

# unzip supplement to get project folder and a place to store output
utils::unzip(supplement_zip, exdir=proj_dir)

# run the simulations in a loop
n_rep = 3L
for(i in seq(n_rep)) {
  
  run_netlogo(proj_dir) 
  gc()
}
