#' Helper function for creating benchmarking input files
#'
#' This takes the example data from the main vignette, (optionally) splits it into
#' test and training sets, then snaps the data to an output grid. The results are then
#' written to a file.
#'
#' This is done for a variety of output grid sizes, which are defined through `eg_rast`.
#' A sequence of output grid resolutions is selected to match the grid dimensions of the
#' rasters in `eg_rast`. So for each example dataset in `eg_pts`, the function makes a
#' sequence of `length(eg_rast)` test data-sets, each one having a different output grid
#' size specified.
#'
#' `eg_rast` should be a list of SpatRasters with the property that entry in the list
#' is a sub-grid of all list entries that follow. `eg_rast` is also used as example data
#' for testing the complete data case. For these examples, instead of withholding a test
#' set, we just use the super-grids (ie larger grids) as test data (exclude the sub-grid
#' of observed points).
#'
#' For testing, a random sample of the observed data (proportion `p_cv`) can be selected
#' to be the training points, and the remainder are test points used for computing MSPE.
#' The training data (and not the testing data) are snapped to the output grid.
#'
#' A training grid is also defined for each case. This is the same as the output grid
#' (with snapped training data) but it usually has smaller dimensions. The size of this
#' grid is determined as the smallest available size such that neither dimension has fewer
#' grid-lines than there are training points.
#
make_inputs = function(eg_pts, eg_rast, dir_storage, p_cv, n_rep)
{
  # eg_pts: list of simple features collections of points
  # eg_rast: list containing sequence of increasingly large SpatRasters
  # dir_storage: directory to write output files
  # p_cv: proportion of the observed data to use for testing (should be within 0-0.9)
  # number of repetitions for measuring evaluation time

  # define directory to write to and path to CSV with info on all files written
  dir_inputs = file.path(dir_storage, 'inputs')
  if(!dir.exists(dir_inputs)) dir.create(dir_inputs)
  inputs_csv_path = file.path(dir_storage, 'inputs.csv')

  # extract grid dimensions from rasters
  output_gdim_test = lapply(eg_rast, function(r) dim(r)[1:2])
  output_n_test = sapply(output_gdim_test, prod)

  # initialize a data-frame to store file info
  inputs_df = data.frame()

  # write the input files in a loop
  cat('\nwriting point datasets...')
  for(i in seq(eg_pts)) {

    # copy points dataset and select a test set
    pts = eg_pts[[i]]
    name_in = names(eg_pts)[i]
    cat(paste0(name_in, '...'))
    n_in = nrow(pts)
    n_test = round(p_cv*n_in)
    idx_test = sort(sample.int(n_in, n_test))

    # copy test and training sets, set a common variable name
    if( n_test == 0 )
    {
      training_pts = pts
      testing_pts = NULL

    } else {

      training_pts = pts[-idx_test,]
      testing_pts = pts[idx_test,]
      names(testing_pts)[1] = 'gval'
    }

    # name the files to write and prepare metadata
    path_in = file.path(dir_inputs, paste0(name_in, '_', output_n_test, '.rds'))
    i_df = data.frame(name_in, n_in, n_train=n_in-n_test, gdim_y=NA, gdim_x=NA, n_out=output_n_test, path_in)
    inputs_df = rbind(inputs_df, i_df)

    # snap training data to high resolution grid (for use by snapKrig)
    idx_res_fit = which( sapply(output_gdim_test, function(d) all(d > n_in)) )[1]
    training_grid = bk_snap(training_pts, output_gdim_test[[idx_res_fit]], quiet=TRUE)

    # write the input files for each of the output grid sizes
    for(j in seq(output_n_test))
    {
      # define output grid (with test data included) to get big enough extent
      output_grid = bk_snap(pts, output_gdim_test[[j]], quiet=TRUE)

      # overwrite with only the training points
      output_grid = bk_snap(training_pts, output_grid, quiet=TRUE)

      # save results to disk
      saveRDS(list(idx_test = idx_test,
                   testing_pts = testing_pts,
                   training_pts = training_pts,
                   output_grid = output_grid,
                   training_grid = training_grid), file=path_in[j])

    }
  }

  cat('done.\nwriting raster datasets...')
  # use a more limited selection of inputs and outputs for rasters
  idx_in = seq(sum(output_n_test < 3e4))
  #is_incomplete = sapply(eg_rast, function(x) global(x, fun='isNA') > 0)
  #idx_out = seq_along(sum(!is_incomplete))
  idx_out = seq_along(output_n_test)

  # write the input files for treed examples in a loop
  for(i in idx_in) {

    # copy training sub-grid
    training_grid = bk_grid(eg_rast[[i]])
    training_pts = bk_coords(training_grid, out='sf', quiet=TRUE)
    n_in = prod(training_grid[['gdim']])
    name_in = paste0('treed_', n_in)
    cat(paste0(name_in, '...'))

    # select testing grids (sub-grids larger than the training grid)
    #idx_out_i = idx_out[idx_out > i]
    idx_out_i = idx_out
    n_out = output_n_test[idx_out_i]

    # name the files to write and prepare metadata
    path_in = file.path(dir_inputs, paste0(name_in, '_', n_out, '.rds'))
    gy = as.numeric(training_grid[['gdim']][1])
    gx = as.numeric(training_grid[['gdim']][2])
    i_df = data.frame(name_in, n_in, n_train=n_in, n_out, gdim_y=gy, gdim_x=gx, path_in)
    inputs_df = rbind(inputs_df, i_df)

    # write the input files for each of the output grid sizes
    for(j in seq_along(idx_out_i))
    {
      # define output grid (with test data included) to get big enough extent
      output_grid_complete = bk_grid( eg_rast[[ idx_out_i[j] ]] )

      # copy with only the training points snapped
      output_grid = bk_snap(training_pts, output_grid_complete, quiet=TRUE)

      # skip test points in this case
      if( p_cv == 0 )
      {
        testing_pts = NULL
        idx_test = NULL

      } else {

        # (excludes points shared with training set)
        idx_test = is.na(output_grid[['gval']])
        testing_pts = bk_coords(output_grid_complete, out='sf', quiet=TRUE)[idx_test,]
      }

      # save results to disk
      saveRDS(list(idx_test = idx_test,
                   testing_pts = testing_pts,
                   training_pts = training_pts,
                   output_grid = output_grid,
                   training_grid = training_grid), file=path_in[j])
    }
  }

  cat('done.')

  #' The resulting files, along with the corresponding problem size in terms of inputs and outputs,
  #' are listed in the dataframe `inputs_df`, which is also saved to disk below:

  # save inputs_df to disk
  write.csv(inputs_df, file=inputs_csv_path, quote=FALSE, row.names=FALSE)
  return(inputs_csv_path)
}



