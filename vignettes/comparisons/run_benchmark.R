# Helper function for running long benchmarking test loop
#
# This takes the input data files created in the main vignette and opens each one
# in a loop, running the ordinary kriging workflow (creating predictions, and if possible,
# prediction variance) for the specified output grid size.
#
# An inner loop repeats this for all the packages listed in `pkg_test`, which must start with 'snapKrig',
# and include any subset of 'fields', 'RandomFields', 'geoR', 'gstat'. For each of these packages there
# should be a corresponding pair of helper functions defined (in 'run_methods.R') with the names
# `{pkg}_bench_fit` and `{pkg}_bench_OK` (where '{pkg}' is the package name). These fit the model
# and run the kriging workflow, timing results.
#
run_benchmark = function(inputs_df, dir_storage, pkg_test, n_rep, timeout, n_max, make_plots=FALSE)
{
  # inputs_df: data-frame containing info on each input test file (and its path)
  # pkg_test: names of packages to test
  # make_plots: logical, indicates to plot the results as the loop progresses

  # define directory to write to and path to CSV with info on all files written
  dir_outputs = file.path(dir_storage, 'outputs')
  if(!dir.exists(dir_outputs)) dir.create(dir_outputs)
  outputs_csv_path = file.path(dir_storage, 'outputs.csv')

  # begin long loop over input files to fill results_df and write outputs
  results_df = data.frame()

  #for(i in 20:80)
  for(i in seq(nrow(inputs_df)))
  {
    # index of all output grids for this example
    name_in = inputs_df[i, 'name_in']
    i_all = which(inputs_df[['name_in']] == name_in)
    i_eg = paste(', output', which(i_all==i), 'of', length(i_all))

    # number of output points
    n_out =  inputs_df[i, 'n_out']

    # print info about the example
    cat(paste0('\n\nEXAMPLE ', i, ' OF ', nrow(inputs_df), ' (', name_in, i_eg, ')'))
    cat(paste('\nnumber of training points:', inputs_df[i, 'n_train']))
    cat(paste('\nnumber of output points:', n_out))
    cat(paste('\nloading input file:', inputs_df[i, 'path_in']))

    # load the input data file and unpack contents
    path_in = inputs_df[i, 'path_in']
    data_in = readRDS(path_in)
    testing_pts = data_in[['testing_pts']]
    training_pts = data_in[['training_pts']]
    output_grid = data_in[['output_grid']]
    training_grid = data_in[['training_grid']]

    # initial values for the problem (assigned later)
    initial = NULL

    # run workflow loop over packages
    for(j in seq_along(pkg_test))
    {
      # copy functions to use in kriging
      package_name = pkg_test[j]
      cat(paste('\n\nrunning workflow with', package_name, 'package...'))
      fit_fun = get(paste0(package_name, '_bench_fit'))
      pred_fun = get(paste0(package_name, '_bench_OK'))

      # define path to write output
      file_prefix = substr(basename(path_in), 1, nchar(basename(path_in))-4)
      file_out = paste0(file_prefix, paste0('_', package_name, '.Rdata'))
      path_out = file.path(dir_outputs, file_out)

      # only fit once to each dataset
      is_first = (i == i_all[1])
      if(is_first)
      {
        # copy snapKrig's fitted parameters to use as initial values in geoR and fields
        if(package_name != 'snapKrig')
        {
          # load fitted parameters from snapKrig
          is_snapKrig = ( results_df[['pkg']] == 'snapKrig' ) & ( results_df[['name_in']] == name_in )
          if(!any(is_snapKrig)) stop('snapKrig model fitting results expected but not found')
          prev_result = readRDS(results_df[ which(is_snapKrig)[1], 'path_out'])
          sk_pars = prev_result[['fit_result']][['result']][['pars']]
          rm(prev_result)

          # copy to snapKrig fitted values to initial values vector
          range_initial = sk_pars[['y']][['kp']]
          psill_initial = sk_pars[['psill']]
          initial = c(psill=psill_initial, range_initial)
        }

        cat('fitting...')

        # geoR and fields receive a copy of snapKrig's fitted values to use as initials
        initial_pass = initial
        if(package_name == 'fields')
        {
          # fields produced singular matrix errors for the supplied initial values in
          # some cases, so we substitute values that will work manually here.
          # NULL means to use fields' defaults
          if(name_in == 'treed') initial_pass = NULL
          if(name_in == 'treed_88') initial_pass = list(rho=1e4)
          if(name_in == 'treed_352') initial_pass = NULL
          if(name_in == 'treed_1376') initial_pass = NULL
        }

        # skip model fitting if above maximum set in n_max
        if( n_out < n_max[[name_in]][['fit']][package_name] )
        {
          # sink suppresses messy console printouts
          sink(nullfile())

          # no timeout on fit calls - user interrupt fails to halt them in some cases anyway
          fit_result = fit_fun(training_grid, training_pts, n_rep, initial_pass, timeout=Inf)
          sink()

        } else {

          fit_result = list(error='skipped')
        }

      } else {

        # load previously fitted parameters
        is_prev = ( results_df[['pkg']] == package_name ) & ( results_df[['name_in']] == name_in )
        if(!any(is_prev)) stop('previous model fitting results expected but not found')
        prev_result = readRDS(results_df[ which(is_prev)[1], 'path_out'])
        fit_result = prev_result[['fit_result']]
        rm(prev_result)
      }

      # initialize ordinary kriging output
      OK_pred_result = OK_var_result = list()
      OK_err = NA
      error_string_pred = NA
      error_string_var = NA
      if( !is.null(fit_result[['error']]) )
      {
        cat('failed. Skipping prediction, MSPE, and variance.')

      } else {

        # run kriging prediction with timeout
        cat('prediction...')
        if( n_out < n_max[[name_in]][['pred']][package_name] )
        {
          # sink suppresses messy console printouts
          sink(nullfile())
          OK_pred_result = pred_fun(output_grid, training_pts, fit_result[['result']],
                                    n_rep, 'p', timeout)

          sink()

        } else {

          OK_pred_result = list(error='skipped')
        }

        # convert output to terra object for MSPE computation
        if( !is.null(OK_pred_result[['error']]) )
        {
          cat('failed. Skipping MSPE and variance.')

          # copy class of error
          error_string_pred = class(OK_pred_result[['error']])[1]

        } else {

          # compute MSPE when there is a test set
          if(!is.null(testing_pts))
          {
            # extract predictions and calculate MSPE
            OK_rast = sk_export(OK_pred_result[['result']])
            OK_diff = testing_pts[['gval']] - terra::extract(OK_rast, testing_pts)[,2]
            OK_err = mean(OK_diff^2, na.rm=TRUE)
          }

          # run kriging variance with timeout
          cat('variance...')
          if( n_out < n_max[[name_in]][['var']][package_name] )
          {
            sink(nullfile())
            OK_var_result = pred_fun(output_grid, training_pts, fit_result[['result']],
                                     n_rep, 'v', timeout)
            sink()

          } else {

            OK_var_result = list(error='skipped')
          }


          # copy class of error (if any)
          if( is.null(OK_var_result[['error']]) ) { cat('done.') } else {

            error_string_var = class(OK_var_result[['error']])[1]
          }
        }
      }

      # assign NAs when execution failed
      if( is.null(fit_result[['teval']]) ) fit_result[['teval']] = NA
      if( is.null(OK_pred_result[['teval']]) ) OK_pred_result[['teval']] = NA
      if( is.null(OK_var_result[['teval']]) ) OK_var_result[['teval']] = NA
      if( is.null(OK_err) ) OK_err = NA

      # gather into data-frame
      results_ij = cbind(inputs_df[i,], data.frame(pkg = package_name,
                                                   fit_s = fit_result[['teval']],
                                                   pred_s = OK_pred_result[['teval']],
                                                   var_s = OK_var_result[['teval']],
                                                   mspe = OK_err,
                                                   path_out = path_out,
                                                   pred_failure = error_string_pred,
                                                   var_failure = error_string_var))

      # copy results to storage data-frame
      results_df = rbind(results_df, results_ij)

      # save important objects to disk
      cat(paste('\nwriting results to', path_out))
      saveRDS(list(fit_fun = fit_fun,
                   pred_fun = pred_fun,
                   path_in = path_in,
                   fit_result = fit_result,
                   OK_pred_result = OK_pred_result,
                   OK_var_result = OK_var_result,
                   results_ij = results_ij), file=path_out)
    }

    # optionally make a plot to show results so far for this example
    if(make_plots)
    {
      # extract all data for this example so far
      results_df_batch = results_df[results_df[['name_in']] == name_in,]
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

      # # create MSPE plot
      # gg2 = ggplot(results_df_batch) +
      #   aes(x=n_out, y=mspe, color=pkg) +
      #   geom_point() +
      #   geom_line() +
      #   geom_line() +
      #   xlab('prediction points') +
      #   ylab('MSPE') +
      #   labs(color='R package',
      #        lty='with variance') +
      #   scale_x_log10(
      #     breaks = scales::trans_breaks('log10', function(x) 10^x),
      #     labels = scales::trans_format('log10', scales::math_format(10^.x))
      #   ) +
      #   scale_y_log10(
      #     breaks = scales::trans_breaks('log10', function(x) 10^x),
      #     labels = scales::trans_format('log10', scales::math_format(10^.x))
      #   ) +
      #   theme_bw() +
      #   theme(text=element_text(size=8),
      #         strip.text.x=element_text(face='bold'),
      #         strip.text.y=element_text(face='bold'))

      # # draw the plots
      # suppressMessages(grid.arrange(gg1, gg2, nrow=2))
    }
  }

  # save results_df to disk as a csv
  write.csv(results_df, file=outputs_csv_path, quote=FALSE)
  return(outputs_csv_path)

}
