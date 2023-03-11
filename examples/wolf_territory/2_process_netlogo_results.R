# script 2/3
# Dean Koch, Feb 2023
# Opens the simulated territory data and runs analysis

# This writes its analysis output - a dataframe and several lists of rasters -
# to a binary (.rds) file in the project directory `proj_dir` (defined below).

# You can run this script at any point after the first simulation result has
# been written to `proj_dir`. All raster (.tif) files in `proj_dir` are loaded,
# so don't copy any other raster files into this directory.


# rasters
library(terra)
library(snapKrig)

# data frames
library(dplyr)

# helper functions
source('D:/wolves/R/0_netlogo_helpers.R')

# directory containing the output files (and no other .tif files!) 
proj_dir = 'D:/wolves/bias_simulations'

# the analysis results file name
analysis_name = 'territory_analysis.rds'


## Load the data

# load the outputs as multiband SpatRaster
output_files = list.files(proj_dir, '.tif$')
output_rast = terra::rast(file.path(proj_dir, output_files))


## Run the bias experiment

# define true occupancy grid and number of test grids
g_source = sk(output_rast) > 0
n_test = 30L
n_rep = ncol(g_source[])
cat('\nprocessing', n_rep, 'simulation output files...\n')

# compute bias at all test resolutions (takes about 2 minutes)
bias_results = run_analysis(g_source, n_test)

# write the results to disk
saveRDS(bias_results, file.path(proj_dir, analysis_name))

