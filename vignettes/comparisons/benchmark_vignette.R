#' ---
#' title: "benchmark_vignette"
#' author: "Dean Koch"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'
#' This code generated the results plotted in Figure 8 and 9 of the paper
#'
#' It creates a set of test data-sets, writes them to disk, then runs an ordinary
#' kriging workflow in a loop over test data-sets, using 5 different packages.
#'
#'
#'
#'
#' ## Important parameters

# project directory
dir_project = 'D:/snapKrig/vignettes/comparisons'

# directory to write input and output files (about 2-4GB total)
dir_storage = file.path(dir_project, 'data')

# number of repetitions for measuring evaluation time
n_rep = 5

# maximum allowed time for kriging prediction or variance (over all repetitions, in seconds)
timeout = n_rep * 10 * 60

# proportion of the observed data to use for cross validation based MSPE (0 to skip)
p_cv = 0

# number of points to sample from treed for irregular points example
n_pts_treed = 1e3

# factors by which to up-scale the treed raster (see Input Files section)
up_test = 2^(7:0)

# factors by which to down-scale the treed raster
down_test = 2^(1:2)

# packages tested (pick any subset of the 'snapKrig', 'fields', 'RandomFields', 'geoR', 'gstat')
pkg_test = c('snapKrig', 'fields', 'RandomFields', 'geoR', 'gstat')


#'
#' This document walks us through the code for making the benchmarking results plotted in the
#' paper. It records evaluation time for ordinary kriging on several datasets, and compares five
#' different R packages. If you're going to run the code yourself, make sure to first change
#' `dir_storage` to something that works on your system.
#'
#' ## Dependencies
#'

# essential GIS packages
library(terra)
library(sf)

# data sources
library(sp)
library(rasterbc)

# snapKrig dev version
library(devtools)
load_all()

# other GP packages to compare
library(fields)
library(RandomFields)
library(geoR)
library(gstat)

# timing and timeout
library(microbenchmark)
library(R.utils)

# plotting results as we go
library(ggplot2)
library(gridExtra)

#'
#' ## Helper functions
#'
#' The code is broken up into several files containing helper function definitions
#'
#' * `get_data.R`: loads source data (`get_ozone`, `get_meuse`, `get_mpb`)
#' * `make_inputs.R`: prepares input data, writing input files (`make_inputs`)
#' * `run_benchmark_methods.R`: defines individual kriging methods for different packages ( `{pkg}_bench_fit`, and `{pkg}_bench_OK`)
#' * `run_benchmark.R`: executes the long benchmark loop, writing output files (`run_benchmark`)
#'
#' These files must be found in your `dir_project` directory

source(file.path(dir_project, 'get_data.R'))
source(file.path(dir_project, 'make_inputs.R'))
source(file.path(dir_project, 'run_benchmark_methods.R'))
source(file.path(dir_project, 'run_benchmark.R'))


#' ## Example Data
#'
#' The example data comes from the `fields`, `sp`, and `rasterbc` packages. The first two are included
#' with the package installation, and third requires downloading two files (13MB total).
#'
#' These rasterbc downloads fo in the sub-folders "pine" and "dem" of `dir_storage` (keep them there
#' to skip the download in future). A cropped version of each one (with no NAs) is also saved to
#' `dir_storage`  in the files "treed.tif" and "dem.tif".

# create the storage directory as needed and define output paths
if(!dir.exists(dir_storage)) dir.create(dir_storage)
path_treed = file.path(dir_storage, 'treed.tif')
path_treed_dem = file.path(dir_storage, 'treed_dem.tif')

# Ozone concentration points (fields::ChicagoO3)
ozone = get_ozone(chi=TRUE)

# Meuse river points (from `sp`)
meuse = get_meuse(dfMaxLength=NA)

# download (or re-load) forest density raster and matching DEM (uses `rasterbc`)
treed = get_mpb(out='treed', rasterbc_storage_path=dir_storage)
treed_dem = get_mpb(out='dem', rasterbc_storage_path=dir_storage)

# save the treed and dem rasters as geotiff
terra::writeRaster(treed, path_treed, overwrite=TRUE)
terra::writeRaster(treed_dem, path_treed_dem, overwrite=TRUE)

# select a subset of the points in treed to serve as irregular example
idx_sub = sort(sample.int(prod(dim(treed)), n_pts_treed))
treed_pts = sk_coords(treed, out='sf')[idx_sub,]


#' ## Input Files
#'
#' Next we define the input and output point locations for each example. The output is always a grid,
#' and we test a range of resolutions. To make comparisons easier, we created a set of test resolutions
#' by multiplying or dividing the dimensions of the treed source data by powers of 2. These powers are
#' defined the `down_test` and `up_test` vectors, respectively.
#'
#' For the treed example itself, we use `sk_rescale` to generate up-scaled versions where points are
#' sampled at regular intervals, forming a sub-grid that we will treat as observed. We also down-scale
#' (but don't impute missing values) to make at least one grid the has finer resolution than treed.

# use snapKrig to up-scale and down-scale
treed_up = lapply(up_test, function(up) sk_rescale(treed, up=up))
treed_down = lapply(down_test, function(down) sk_rescale(treed, down=down))

# export to a big list of terra SpatRasters
eg_rast = lapply(c(treed_up, treed_down), sk_export)

# gather all non-gridded example datasets in a list for storage
eg_pts = list(ozone=ozone, meuse=meuse[['soils']]['log_zinc'], treed=treed_pts)

#' Pass the results to `make_inputs` to write input files corresponding to each example dataset, and
#' each output resolution - a total of 48 test cases. This writes the 48 input files to the "inputs"
#' sub-folder of `dir_storage`, with file-names of the form "{example name}_{n_in}_{n_out}.rds"
#' (where the "_{n_in}_" dropped for the irregular point datasets). It also writes a CSV with details
#' about each of these files to the file "inputs.csv" in `dir_storage`.
#'
#' This should take about a minute to complete and needs about 51 MB of disk space.

# make the input files (writes files and returns path to a csv with more info)
#path_inputs_df = make_inputs(eg_pts, eg_rast, dir_storage, p_cv, n_rep)
path_inputs_df = file.path(dir_storage, 'inputs.csv')

# load the resulting file info data-frame
inputs_df = read.csv(path_inputs_df)

# remove the big objects from workspace and run garbage collection
rm(ozone, meuse, treed_pts, eg_pts, treed, treed_dem, treed_up, treed_down, eg_rast)
gc(verbose=FALSE)

#'
#' ## Testing loop
#'
#' Now we can run the benchmarking loop. This opens each input file and runs the OK workflow
#' from each package, timing results and writing results to a new rds file in the "outputs"
#' sub-directory of `dir_storage`. For each input file, `length(pkg_test)` such output files
#' are created (a total of 240). When the function call is finishing, it writes a data-frame
#' with details and results from each test case to the file "outputs.csv" in `dir_storage`.
#'
#' This takes a long time (12-24 hours) and requires about 2.5GB of disk space.
#'

# inputs_df = inputs_df[44:48,]
# pkg_test = 'fields'
# path_outputs_df = run_benchmark(inputs_df, dir_storage, pkg_test, n_rep, timeout, make_plots=TRUE)

# maximum number of prediction points allowed for each package in the various test datasets
# this helps to limit total compute time if we have a set of bounds in mind for the plot already
n_max = list(

  # irregular point examples
  ozone = list(fit=c(snapKrig=Inf, fields=Inf, RandomFields=Inf, geoR=Inf, gstat=Inf),
               pred=c(snapKrig=Inf, fields=Inf, RandomFields=Inf, geoR=Inf, gstat=10^7),
               var=c(snapKrig=Inf, fields=Inf, RandomFields=Inf, geoR=Inf, gstat=10^7)),

  meuse = list(fit=c(snapKrig=Inf, fields=Inf, RandomFields=Inf, geoR=Inf, gstat=Inf),
               pred=c(snapKrig=Inf, fields=Inf, RandomFields=10^7, geoR=10^7, gstat=10^7),
               var=c(snapKrig=10^7, fields=10^7, RandomFields=10^7, geoR=10^7, gstat=10^7)),

  treed = list(fit=c(snapKrig=Inf, fields=Inf, RandomFields=Inf, geoR=Inf, gstat=Inf),
               pred=c(snapKrig=Inf, fields=10^7, RandomFields=10^6, geoR=10^6, gstat=10^6),
               var=c(snapKrig=10^6.5, fields=10^6, RandomFields=10^6, geoR=10^6, gstat=10^6)),

  # raster examples

  treed_88 = list(fit=c(snapKrig=Inf, fields=Inf, RandomFields=Inf, geoR=Inf, gstat=Inf),
                  pred=c(snapKrig=Inf, fields=Inf, RandomFields=10^7, geoR=10^7, gstat=10^7),
                  var=c(snapKrig=10^7, fields=10^7, RandomFields=10^7, geoR=10^7, gstat=10^7)),

  treed_352 = list(fit=c(snapKrig=Inf, fields=Inf, RandomFields=Inf, geoR=Inf, gstat=Inf),
                   pred=c(snapKrig=Inf, fields=10^7, RandomFields=10^6.5, geoR=10^6.5, gstat=10^6.5),
                   var=c(snapKrig=10^7, fields=10^6.5, RandomFields=10^6.5, geoR=10^6.5, gstat=10^6.5)),

  treed_1376 = list(fit=c(snapKrig=Inf, fields=Inf, RandomFields=Inf, geoR=Inf, gstat=Inf),
                    pred=c(snapKrig=Inf, fields=10^6.5, RandomFields=10^5, geoR=10^5, gstat=10^5),
                    var=c(snapKrig=10^6.5, fields=10^5, RandomFields=10^5, geoR=10^5, gstat=10^5)),

  # these are extras where we only want the fitting times

  treed_5440 = list(fit=c(snapKrig=Inf, fields=0, RandomFields=0, geoR=0, gstat=0),
                    pred=c(snapKrig=0, fields=0, RandomFields=0, geoR=0, gstat=0),
                    var=c(snapKrig=0, fields=0, RandomFields=0, geoR=0, gstat=0)),

  treed_21632 = list(fit=c(snapKrig=Inf, fields=0, RandomFields=0, geoR=0, gstat=0),
                     pred=c(snapKrig=0, fields=0, RandomFields=0, geoR=0, gstat=0),
                     var=c(snapKrig=0, fields=0, RandomFields=0, geoR=0, gstat=0))
)

# run the big loop
path_outputs_df = run_benchmark(inputs_df, dir_storage, pkg_test, n_rep, timeout, n_max, make_plots=TRUE)
#outputs_df = read.csv(path_outputs_df)

# open all the output files
outputs_fname = list.files(file.path(dir_storage, 'outputs'))
outputs_df_list = lapply(outputs_fname, function(f) readRDS(file.path(dir_storage, 'outputs', f))[['results_ij']])
outputs_df = do.call(rbind, outputs_df_list)
results_df = outputs_df
results_df = results_df[order(results_df$name_in, results_df$n_out),]
outputs_csv_path = file.path(dir_storage, 'outputs.csv')
write.csv(results_df, file=outputs_csv_path, quote=FALSE)



# extract all data for this example so far
name_in = 'treed_1376'
{
results_df_batch = outputs_df[outputs_df[['name_in']] == name_in,]
results_df_batch[['log_10_n_out']] = log(results_df_batch[['n_out']], base=10)
col_names = c('black', 'blue', 'violet', 'red')
plot_colors = col_names[ match(results_df_batch[['pkg']], pkg_test) ]

# set zero eval time for situations where variance was skipped or doesn't apply
results_df_batch[['var_s']][ results_df_batch[['pkg']] %in% c('gstat', 'geoR', 'RandomFields') ] = 0

# create variance + prediction eval time column
results_df_batch[['all_s']] = results_df_batch[['pred_s']] + results_df_batch[['var_s']]
results_df_batch[['with_var']] = 'yes'
results_df_batch[['with_var']][results_df_batch[['pkg']] == 'RandomFields'] = 'no'

# add prediction time rows for two packages that allowed separate variance calculation
idx_pred_only = results_df_batch[['pkg']] %in% c('fields', 'snapKrig')
results_df_batch_new = results_df_batch[idx_pred_only,]
results_df_batch_new[['all_s']] = results_df_batch_new[['pred_s']]
results_df_batch_new[['with_var']] = 'no'
results_df_batch_combo = rbind(results_df_batch, results_df_batch_new)

# create prediction + variance time plot
gg1 = ggplot(results_df_batch_combo) +
  aes(x=n_out, y=all_s, color=pkg, lty=with_var) +
  geom_point() +
  geom_line() +
  geom_line() +
  xlab('prediction points') +
  ylab('time (seconds)') +
  labs(color='R package',
       lty='with variance') +
  scale_x_log10(
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x))
  ) +
  theme_bw() +
  theme(text=element_text(size=8),
        strip.text.x=element_text(face='bold'),
        strip.text.y=element_text(face='bold'))

print(gg1)
}
#inputs_df


#'
#' ## Files written
#'
#' Here is a log of all files written by the script

print(path_inputs_df)
print(path_treed)
print(path_treed_dem)
print(inputs_df[['path_in']])

print(path_outputs_df)
print(outputs_df[['path_out']])

