# script 3/3
# Dean Koch, Feb 2023
# make some plots summarizing the NetLogo analysis

# This makes three plots illustrating results from the process_netlogo_results.R
# workflow. Run that script first so that a file named `analysis_name` can be
# found in the project directory `proj_dir`. The script creates two graphs from
# the results data and writes them as png files to `proj_dir`

# rasters
library(terra)
library(snapKrig)

# data frames
library(dplyr)

# helper functions
source('D:/wolves/R/0_netlogo_helpers.R')

# directory containing the output files (and no other .tif files!) 
proj_dir = 'D:/wolves/bias_simulations'

# the analysis results file and output image files
analysis_name = 'territory_analysis.rds'
png_chart = 'bias_chart.png'
png_chart_v2 = 'bias_chart_v2.png'
png_raster = 'wolf_territory.png'

## Load the data

# a list containing results in a data frame and nested lists of sk objects
bias_results = file.path(proj_dir,analysis_name) |> readRDS()

# this data frame summarizes all results
bias_df = bias_results[['df']]

# specific resolutions highlighted
res_plot = c(16L, 100L, 576L)

# make a line chart
png(file.path(proj_dir, png_chart), width=800, height=400, unit='px', pointsize=20)
plot_bias(bias_df,
          res_plot,
          ylab='overestimate',
          xlab='cell area',
          main='bias in territory size estimates')
dev.off()

# make another line chart with adjusted axes showing abundance estimate
png(file.path(proj_dir, png_chart_v2), width=800, height=400, unit='px', pointsize=20)
true_area = sum(bias_results[['g_fine']][[1]][[1]] == 2L)
bias_df |>
  dplyr::mutate(percent_over = 1.15 * (1 + units::drop_units(percent_over)/100) * true_area * 5 / 447) |>
  plot_bias(res_plot,
            ylab='abundance estimate',
            xlab='cell area',
            main='abundance estimates with pack size 5 (truth is 730)')
dev.off()

# make a raster plot showing the example resolutions
png(file.path(proj_dir, png_raster), width=1200, height=400, unit='px', pointsize=18)
plot_territory(bias_results,
               res_plot,
               cex.main=1.3,
               cex.z=1.2,
               leg_just=0.3)
dev.off()
