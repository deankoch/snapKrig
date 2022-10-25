#' ---
#' title: "Benchmarking with snapKrig"
#' author: "Dean Koch"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'
#'
#' This document walks us through the code for making the benchmarking results plotted in the
#' paper. It records evaluation time for ordinary kriging on several datasets, and compares four
#' different R packages. If you're going to run the code yourself, make sure `dir_storage` (below)
#' is set to something sensible for your PC.
#'
#' The entire script takes about 6 hours to finish with default parameters and it writes 3.8GB.
#' The main output file is summary of results in a CSV file (`path_main_result`, see below)
#'
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
library(geoR)
library(gstat)

# timing and timeout
library(microbenchmark)
library(R.utils)

# for plotting results
library(ggplot2)
library(gridExtra)

# for managing working directories
library(here)


#'
#' The goal is to generate the results plotted in Figure 8 and 9 of the paper.
#' The script will create a set of test data-sets, then (in a loop) run an ordinary
#' kriging workflow using four different packages.
#'
#' We use helper functions to load the data and to time operations in different
#' packages. To keep this document tidy, I defined these helper functions in separate
#' .R files, which are sourced below. Make sure you have these files in your working
#' directory, so that `here` points to the right place:
#'

# set project directory and find two helper functions
dir_project = here('vignettes/comparisons')
path_deps = file.path(dir_project, c('helper_get_data.R', 'helper_benchmark_methods.R'))
invisible(sapply(path_deps, source))

#'
#' The first group of helper functions (`get_ozone`, `get_meuse`, `get_mpb`) is for
#' loading example data. The second group (`bench_fit`, `bench_likelihood`, `bench_kriging`)
#' is for executing various steps in the modeling workflow, and timing computations.
#' `make_all_results` and `make_results_plot` are called during the big testing loop below
#' to assemble results and display them as a line chart.
#'
#' Output paths are defined next. Note that `dir_bench_results` should point to a
#' location with at least 3.8GB of free storage space. This does not need to be
#' a sub-directory of `dir_project`.
#'

# directory for file output (CAUTION: THE SCRIPT MAY OVERWRITE EXISTING FILES IN THIS DIRECTORY!)
dir_storage = file.path(dir_project, 'data')

# define output paths: treed raster data, and matching DEM, benchmarking outputs folder
path_treed = file.path(dir_storage, 'treed.tif')
path_treed_dem = file.path(dir_storage, 'treed_dem.tif')
dir_bench_results = file.path(dir_storage, 'benchmarks')

# define the main output file
path_main_result = file.path(dir_storage, 'bench_results.csv')

#'
#' Results will be slightly different each time the script is run because of the
#' random point selection for the treed dataset (below), and because other
#' running processes might slow the computer in unpredictable ways during the
#' testing loop.
#'
#' However the big picture should always look the same - on large
#' sample examples we expect 1-2 orders of magnitude speedup in prediction with
#' snapKrig and even better speedups in likelihood in the complete grid case.
#'
#' ## Important parameters
#'

# number of repetitions for measuring evaluation time
n_rep = 5

# maximum allowed time for each operation (over all repetitions, in seconds)
timeout = n_rep * 180

# number of points to sample for irregular points example "treed" (maximum 1377329)
n_pts_treed = 1e3

# factors by which to up-scale the treed raster (see Input Files section)
up_test = 2^(7:0)

# factors by which to down-scale the treed raster
down_test = 2^(1:2)

# packages
pkg_test = c('snapKrig', 'fields', 'geoR', 'gstat')

# flag indicating to display existing results
overwrite = FALSE


#'
#' ## Example Data
#'
#' The example data comes from the `fields`, `sp`, and `rasterbc` packages. The first two are included
#' with the package installation, and third requires downloading two files (13MB total).
#'
#' These rasterbc downloads fo in the sub-folders "pine" and "dem" of `dir_storage` (keep them there
#' to skip the download in future). A cropped version of each one (with no NAs) is also saved to
#' `dir_storage`  in the files "treed.tif" and "dem.tif".

# create the storage directory as needed
if(!dir.exists(dir_storage)) dir.create(dir_storage)

# Ozone concentration points (fields::ChicagoO3)
ozone = get_ozone(chi=TRUE)

# Meuse river points (from `sp`)
meuse = get_meuse(dfMaxLength=NA)

# download (or re-load) forest density raster and matching DEM (uses `rasterbc`)
treed = get_mpb(out='treed', rasterbc_storage_path=dir_storage)
treed_dem = get_mpb(out='dem', rasterbc_storage_path=dir_storage)

# save the treed and dem rasters as geotiff
if(overwrite)
{
  terra::writeRaster(treed, path_treed, overwrite=TRUE)
  terra::writeRaster(treed_dem, path_treed_dem, overwrite=TRUE)
}

# select a subset of the points in treed to serve as irregular example
idx_sub = sort(sample.int(prod(dim(treed)), n_pts_treed))
treed_pts = sk_coords(treed, out='sf')[idx_sub,]


#' ## Inputs and Outputs
#'
#' Next we define the input and output point locations for each example. The output is always a grid,
#' and we test a range of resolutions. To make comparisons easier, we establish test resolutions
#' by multiplying or dividing the dimensions of the treed source data by powers of 2. These powers are
#' defined the `down_test` and `up_test` vectors.
#'
#' For the treed example itself, we use `sk_rescale` to generate up-scaled versions where points are
#' sampled at regular intervals, forming a sub-grid of observed data. We also down-scale
#' (but don't impute missing values) so that at least one grid has finer resolution than treed.

# use snapKrig to up-scale and down-scale
treed_up = lapply(up_test, function(up) sk_rescale(treed, up=up))
treed_down = lapply(down_test, function(down) sk_rescale(treed, down=down))
treed_grids = c(treed_up, treed_down)
names(treed_grids) = paste0('treed_', sapply(treed_grids, length))

# collect output grid resolutions to test
gdim_out = lapply(treed_grids, dim)
ny_out = sapply(gdim_out, function(x) x[1])
nx_out = sapply(gdim_out, function(x) x[2])
print(gdim_out)

# omit the down-scaled (larger) grids which are replicates of treed_1377329
treed_grids = treed_grids[seq_along(treed_up)]

# print the results
invisible(sapply(treed_grids, print))

#' Next we bundle all of the example datasets-  ozone, meuse, treed, and the 11
#' additional up-scaled versions of the complete treed raster.
#'

# gather all datasets in a list for storage
eg_pts = list(ozone=ozone, meuse=meuse[['soils']]['log_zinc'], treed=treed_pts)
eg_all = c(eg_pts, treed_grids)

# set up paths to fitting results files
path_fit_results = file.path(dir_storage, paste0(names(eg_pts), '_fit_results.rds'))
names(path_fit_results) = names(eg_pts)

# set up file name structure and location for benchmark results files
if( !dir.exists(dir_bench_results) ) dir.create(dir_bench_results)
suffix_basename = paste0(names(eg_all), '_bench_results')
full_basename = paste0(pkg_test, '_', rep(suffix_basename, each=length(pkg_test)))

#' Results are bundled by package and example. For each of of the examples in
#' `eg_all` the loop writes its results to files (one per package) in `dir_bench_results` with
#' names beginning with the package name and ending with `suffix_basename`

# paths to outputs from kriging benchmarks, as table and as detailed list of objects
path_bench_results_csv = file.path(dir_bench_results, paste0(full_basename, '.csv'))
path_bench_results_rds = file.path(dir_bench_results, paste0(full_basename, '.rds'))

# catch unintended execution
if(!overwrite) stop('end of script (run again with overwrite=TRUE to continue past this line)')


#'
#' ## Fit models
#'
#' We start by fitting models to the three irregular points datasets. snapKrig is fitted
#' first, and its fitted parameters are used as starting values in the other packages
#' (where required). The result are bundled by example dataset, and saved to a file.
#'
#' This loops over the datasets, with an inner loop fitting with each package. It
#' uses the `bench_fit` helper function defined in run_benchmark_methods.R. Each
#' iteration of the outer loop writes one of the files in `path_fit_results`.
#'
#' This exercise is about execution times, so we are not particularly in the quality of
#' the model fit as long as the fitted parameters are plausible, and they don't lead to
#' errors later on when we call kriging methods.
#'

# look over all irregular point examples
for(eg_i in seq_along(path_fit_results))
{
  # storage of outputs
  results_list = vector(mode='list', length=length(pkg_test))
  names(results_list) = pkg_test

  # snap to coarse resolution grid and fit with snapKrig
  pts = eg_pts[[eg_i]]
  g_train = sk_snap(pts)
  cat(paste('fitting with snapKrig', '...\n'))
  results_list[['snapKrig']] = bench_fit(g_train, pkg='snapKrig', n_rep)
  pars_ini = results_list[['snapKrig']][['result']]

  # fit with the other packages (without snapping)
  for(pkg in pkg_test[-1])
  {
    cat(paste('fitting with', pkg, '...\n'))

    # suppress console output with sink
    sink(nullfile())
      # passes snapKrig fitted values as initial values for parameters
      results_list[[pkg]] = bench_fit(pts, pkg, pars_ini, n_rep)
    sink()
  }

  # save all relevant objects
  saveRDS(list(pts=pts, g_train=g_train, results=results_list), file=path_fit_results[eg_i])
}



#'
#' ## Time likelihood and kriging
#'
#' Up next is a much slower loop over example datasets, where we include the up-scaled
#' raster data examples in addition to the irregular ones fitted above. We re-use the
#' same fitted values from above in all the examples below. This helps speed things up.
#'
#' The next code chunk defines data frames to hold the results generated in the big
#' testing loop that comes next.
#'

# specify which fitted values go with which example dataset
nm_fitted = c(names(eg_pts), rep('treed', length(eg_all) - length(eg_pts)))

# initialize data frame for likelihood outputs from each example dataset
results_lik_df = data.frame(name=names(eg_all),
                            n_in=c(sapply(eg_pts, nrow), sapply(treed_grids, length)),
                            complete = c(rep(FALSE, length(eg_pts)), rep(TRUE, length(treed_grids))),
                            teval_lik = NA)

# replicate rows for each package
results_lik_df = results_lik_df[rep(seq(length(eg_all)), each=length(pkg_test)),]
results_lik_df[['pkg']] = rep(pkg_test, length(eg_all))


# define a template data frame for output from each package
results_df_template = data.frame(n_out = stats::setNames(nm=ny_out*nx_out),
                                 ny_out = ny_out,
                                 nx_out = nx_out,
                                 teval_pred = NA,
                                 teval_var = NA,
                                 teval_both = NA)



#' ## The Big Testing Loop
#'
#' An outer loop over examples, with two nested inner loops over packages and
#' and output size. In each iteration, evaluation times are recorded for likelihood,
#' prediction, and variance.
#'
#' A linear regression on existing results is used to estimate computation times
#' for the current operation. If the estimated time (of all `n_rep` repetitions) exceeds
#' `timeout`, the operation is skipped.
#'
#' WARNING: THIS TAKES ABOUT 6 HOURS TO COMPLETE

# look over examples
for(eg_i in seq_along(eg_all))
{
  # print console output
  eg_nm = names(eg_all)[eg_i]
  cat(paste('\n--------------\nprocessing', eg_nm, 'example', '...\n'))

  # load the model fitting results file
  path_prev = path_fit_results[nm_fitted[eg_i]]
  cat(paste('loading previously fitted model from', path_prev, '\n'))
  fit_prev = readRDS(path_prev)

  # loop over packages
  for(pkg_j in seq_along(pkg_test))
  {
    # copy info about this package
    pkg_nm = pkg_test[pkg_j]
    cat(paste0('\n__', pkg_nm, '__\n'))
    pars = fit_prev[['results']][[pkg_nm]][['result']]
    is_combo = pkg_nm %in% c('gstat', 'geoR')

    # find current row number in likelihood results data frame and copy some attributes
    idx_lik = pkg_j + length(pkg_test) * ( eg_i - 1 )
    is_complete = results_lik_df[['complete']][idx_lik]
    n_in = results_lik_df[['n_in']][idx_lik]

    # set path and file name to write
    path_csv_j = path_bench_results_csv[idx_lik]
    path_rds_j = path_bench_results_rds[idx_lik]

    ## estimate computation time ahead based on previous results

    # find subset of comparable likelihood results from earlier
    sub_lik_df = subset(results_lik_df, ( pkg == pkg_nm ) & ( complete == is_complete ))

    # linear model for computation times against n_in (proceed if we have at least one time recorded)
    if( all(is.na(sub_lik_df[['teval_lik']])) ) { teval_lik_est = -1 } else {

      # fit model on log-log scale and suppress warnings in prediction (about single obs case)
      lm_lik = lm(log(teval_lik)~log(n_in), sub_lik_df, weights=log(sub_lik_df[['n_in']]))
      teval_lik_est = n_rep * suppressWarnings( exp(predict(lm_lik, newdata=data.frame(n_in))) )
    }

    ## Likelihood

    # proceed with likelihood only if estimated time is reasonable
    if( teval_lik_est > timeout )
    {
      result_lik = list(result=NA, teval=NA, error='estimated time exceeded timeout')

    } else {

      # time likelihood, using sink to suppress console printouts
      cat(paste0('timing likelihood (', round(teval_lik_est, 3), 's estimated)...\n'))
      sink(nullfile())
        result_lik = bench_likelihood(pars, eg_all[[eg_i]], pkg=pkg_nm, n_rep, 2*timeout)
      sink()
    }

    # copy result to likelihood timing data frame
    results_lik_df[['teval_lik']][idx_lik] = ifelse(pkg_nm=='gstat', NA, result_lik[['teval']])

    # create storage lists for outputs created below
    results_df = cbind(results_lik_df[idx_lik,], results_df_template, row.names=NULL)
    results_pred_list = vector(mode='list', length=nrow(results_df))
    results_var_list = vector(mode='list', length=nrow(results_df))
    names(results_pred_list) = results_df[['n_out']]
    names(results_var_list) = results_df[['n_out']]

    ## Test a range of output sizes

    # loop over output dimensions
    cat('timing kriging prediction and variance...\n')
    #for(i in 1:4)
    for(i in seq(nrow(results_df)))
    {
      # requested output dimensions
      gdim = unlist(results_df[i, c('ny_out', 'nx_out')])
      names(gdim) = c('y', 'x')
      n_out = prod(gdim)
      cat(paste('', n_out,  'output points\n'))

      ## estimate computation time ahead based on previous results

      # linear model for computation times against n_out (proceed if we have at least one time recorded)
      if( all(is.na(results_df[['teval_both']])) ) { teval_both_est = -1 } else {

        # fit model on log-log scale and suppress warnings in prediction (about single obs case)
        lm_both = lm(log(teval_both)~log(n_out), results_df, weights=log(results_df[['n_out']]))
        teval_both_est = n_rep * suppressWarnings( exp(predict(lm_both, newdata=data.frame(n_out))) )
      }

      # linear model for prediction only (proceed if we have at least one time recorded)
      if( all(is.na(results_df[['teval_pred']])) ) { teval_pred_est = teval_both_est } else {

        lm_pred = lm(log(teval_pred)~log(n_out), results_df, weights=log(results_df[['n_out']]))
        teval_pred_est = n_rep * suppressWarnings( exp(predict(lm_pred, newdata=data.frame(n_out))) )
      }

      # linear model for variance only (proceed if we have at least one time recorded)
      if( all(is.na(results_df[['teval_var']])) ) { teval_var_est = teval_both_est } else {

        lm_var = lm(log(teval_var)~log(n_out), results_df, weights=log(results_df[['n_out']]))
        teval_var_est = n_rep * suppressWarnings( exp(predict(lm_var, newdata=data.frame(n_out))) )
      }

      # exception to skip large observed n examples with gstat, geoR, fields
      if( ( (pkg_nm %in% c('gstat', 'geoR')) & (n_in > 2e3) ) | ( (pkg_nm == 'fields') & (n_in > 6e3) ) )
      {
        teval_pred_est = timeout + 1
        teval_var_est = timeout + 1
      }

      ## Kriging Prediction and Variance

      # proceed with prediction only if estimated time is reasonable
      if( teval_pred_est > timeout )
      {
        results_pred_list[[i]] = list(result=NA, teval=NA, error='estimated time exceeded timeout')

      } else {

        # time kriging prediction
        cat(paste0('  prediction (', round(teval_pred_est, 3), 's estimated)...'))
        sink(nullfile())
          results_pred_list[[i]] = bench_kriging(eg_all[[eg_i]], gdim, pars, pkg_nm, out='p', n_rep, 2*timeout)
        sink()
      }

      # proceed with variance only if estimated time is reasonable (and skip when unsupported)
      if( (teval_var_est > timeout) | (teval_pred_est > timeout) | (pkg_nm %in% c('gstat', 'geoR')) )
      {
        results_var_list[[i]] = list(result=NA, teval=NA, error='estimated time exceeded timeout')

      } else {

        # time kriging variance
        cat(paste0('\n  variance (', round(teval_var_est, 3), 's estimated)...'))
        sink(nullfile())
        results_var_list[[i]] = bench_kriging(eg_all[[eg_i]], gdim, pars, pkg_nm, out='v', n_rep, 2*timeout)
        sink()
      }

      # copy time results to data frame storage
      teval_add = c(results_pred_list[[i]][['teval']], results_var_list[[i]][['teval']])
      all_skipped = all(is.na(teval_add))
      results_df[i, c('teval_pred', 'teval_var')] = teval_add
      results_df[i, c('teval_both')] = ifelse(all_skipped, NA, sum(teval_add, na.rm=TRUE))

      # these packages don't support separate prediction, variance (both returned in one call)
      if(is_combo) results_df[i, c('teval_pred', 'teval_var')] = NA

      # if variance was skipped for these packages, set combo time to NA
      if(is.na(teval_add[2]) & !is_combo) results_df[i, c('teval_both')] = NA
      cat(ifelse(all_skipped, 'skipped\n', '\n'))
    }

    # bundle detailed results into a big list
    write_list = list(df=results_df,  list=list(lik = result_lik,
                                                pred = results_pred_list,
                                                var = results_var_list))

    # write to disk
    cat(paste('writing results to', dir_bench_results, '\n'))
    cat(paste(' ', basename(path_csv_j), '...\n'))
    write.csv(results_df, file=path_csv_j, quote=FALSE)
    cat(paste(' ', basename(path_rds_j), '...\n'))
    saveRDS(write_list, file=path_rds_j)
    cat('done\n\n')

    print(results_df)
    cat('\n\n')

    # print a plot of results so far
    make_results_plot(make_all_results(path_bench_results_csv), eg_nm)
  }
}


#'
#' ## Compile results
#'
#' If the above loop was successful, then information we want to analyze is spread out
#' over 44 csv small files. This last chunk opens all the files and merges them into a
#' single table, which it saves to a csv file

all_results_df = make_all_results(path_bench_results_csv)
write.csv(all_results_df, path_main_result)
