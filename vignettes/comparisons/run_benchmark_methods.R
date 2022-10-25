# this is a wrapper for `fields::Exp.cov` to fix a bug in `fields::spatialProcess`
my_Expcov = function(x1, x2 = NULL, aRange = 1, p = 1, distMat = NA, C = NA,
                     marginal = FALSE, onlyUpper = FALSE, theta = NULL,
                     Distance=NULL, Dist.args=NULL)
{
  # this function simply accepts and ignores the problematic arguments Distance, Dist.args
  Exp.cov(x1, x2=x2, aRange = aRange, p=p, distMat = distMat,
          C = C, marginal = marginal, onlyUpper=onlyUpper, theta=theta)


}


# time model-fitting with snapKrig, geoR, fields, gstat
bench_fit = function(g, pkg, pars_initial, n_rep=1, timeout=5*60)
{
  # g: sk grid or (if not using snapKrig) sf POINTS data frame containing the observed data
  # pkg: one of 'snapKrig', 'geoR', 'fields', 'gstat'
  # pars_initial: fitted parameters list from snapKrig (ignored if pkg=='snapKrig')
  # n_rep: number of repetitions
  # timeout: maximum time (seconds) allowed for all repetitions

  # catch and return any errors
  error_result = tryCatch({

    # interrupt evaluation when timeout is exceeded
    mb_fit = withTimeout({

      # extract non-NA points to data-frame
      if(pkg != 'snapKrig')
      {
        if(inherits(g, 'sk')) g = sk_coords(g, out='sf', na_omit=TRUE, quiet=TRUE)
        g_coords = sf::st_coordinates(g)
        g_vals = sf::st_drop_geometry(g)

        # extract parameters from sk fitted model to use as starting values for the other packages
        p_ini = c(unlist(pars_initial[c('psill', 'eps')]), pars_initial[['y']][['kp']]['rho'])
      }

      # prepare arguments for gstat
      if(pkg == 'gstat')
      {
        pts_gstat = data.frame(g_vals, g_coords)
        names(pts_gstat) = c('z', 'x', 'y')
      }

      # prepare arguments for geoR
      if(pkg=='geoR')
      {
        # export observed data to geodata list and copy starting parameters
        pts_geoR = geoR::as.geodata(cbind(g_coords, g_vals))
        pars_geoR = p_ini[c('psill', 'rho')]
      }

      # prepare arguments for fields
      if(pkg=='fields')
      {
        # fields specifies scaled partial sill directly
        p_ini['lambda'] = p_ini['eps']/p_ini['psill']
        pars_fields = list(aRange=p_ini['rho'], lambda=p_ini['lambda'])
      }

      # repeat evaluations several times to get a median time
      microbenchmark({

        # different likelihood calls for different packages
        fit_result = switch(pkg,

                            # fast and direct
                            'snapKrig' = sk_fit(g, iso=TRUE, quiet=TRUE),

                            # fits a variogram by least squares (much faster, but not likelihood based)
                            'gstat' =  fit.variogram(variogram(z~1, locations=~x+y, data=pts_gstat),
                                                     model=vgm('Gau')),

                            # this does profile likelihood on 1 test value (to get one LL evaluation)
                            'geoR' = geoR::likfit(pts_geoR,
                                                  ini.cov.pars = pars_geoR,
                                                  cov.model = 'gaussian',
                                                  fix.nugget = TRUE,
                                                  nugget = p_ini['eps']),

                            # more direct, but syntax is a bit complicated
                            'fields' = spatialProcess(g_coords, g_vals,
                                                      cov.function = 'my_Expcov',
                                                      cov.args = list(p=2),
                                                      mKrig.args = list(m=1),
                                                      cov.params.start = pars_fields,
                                                      Dist.args = list(compact=FALSE))
        ) # end switch

      }, times=n_rep) # end microbenchmark repetitions

    }, timeout=timeout, onTimeout='error') # end timeout

    NULL # return value from trycatch if the execution succeeds

  }, error = identity) # end tryCatch

  # return from failed execution
  if( !is.null(error_result) ) return( list(result=NA, teval=NA, error=error_result) )

  # return all results in list (convert teval to units of seconds)
  return( list(result = fit_result,
               teval = median(mb_fit[['time']]) / 1e9,
               error = error_result) )
}

# time a likelihood function evaluation with snapKrig, geoR, fields
bench_likelihood = function(pars, g, pkg, n_rep=10, timeout=5*60)
{
  # pars: the (package-specific) results object returned from bench_fit
  # g: sk grid or (if not using skapKrig) sf POINTS data frame containing the observed data
  # pkg: one of 'snapKrig', 'geoR', 'fields'
  # n_rep: number of repetitions
  # timeout: maximum time (seconds) allowed for all repetitions

  # catch and return any errors
  error_result = tryCatch({

    # interrupt evaluation when timeout is exceeded
    mb_lik = withTimeout({

      # reshape data as needed
      if(pkg == 'snapKrig')
      {
        # if a points dataframe is passed with snapKrig, snap it to a grid
        if(!inherits(g, 'sk')) g = sk_snap(g)

      } else {

        # extract non-NA points to data-frame
        if(inherits(g, 'sk')) g = sk_coords(g, out='sf', na_omit=TRUE, quiet=TRUE)
        g_coords = sf::st_coordinates(g)
        g_vals = sf::st_drop_geometry(g)

      }

      # prepare arguments for geoR
      if(pkg=='geoR')
      {
        # make a geodata list
        g_geoR = geoR::as.geodata(cbind(g_coords, g_vals))

        # pars must be a likfit class object matching the input data so we fit again (1 iteration)
        pars_geoR = geoR::likfit(g_geoR, cov.model='gaussian',
                                 ini.cov.pars = c(pars[['sigmasq']], sqrt(pars[['phi']])),
                                 control=list(maxit=1))
      }

      # prepare arguments for fields
      if(pkg=='fields')
      {
        # fields accepts covariance parameters in arguments - copy fitted values
        pars_fields = list(lambda = pars[['MLEInfo']][['pars.MLE']]['lambda'],
                           aRange = pars[['MLEInfo']][['pars.MLE']]['aRange'],
                           sigma2 = pars[['rhohat']])
      }

      # repeat evaluations several times to get a median time
      microbenchmark({

        # different likelihood calls for different packages
        lik_result = switch(pkg,

                            # fast and direct
                            'snapKrig' = sk_LL(pars, g, quiet=TRUE),

                            # this does profile likelihood on 1 test value (to get one LL evaluation)
                            'geoR' = geoR::proflik(pars_geoR, g_geoR, nugget.values=pars[['nugget']]),

                            # more direct, but syntax is a bit complicated
                            'fields' = fields::mKrigMLEJoint(g_coords, g_vals,
                                                             cov.function = 'my_Expcov',
                                                             mKrig.args = list(m=1),
                                                             cov.args = pars_fields)
        ) # end switch

      }, times=n_rep) # end microbenchmark repetitions

    }, timeout=timeout, onTimeout='error') # end timeout

    NULL # return value from trycatch if the execution succeeds

  }, error = identity) # end tryCatch

  # return from failed execution
  if( !is.null(error_result) ) return( list(result=NA, teval=NA, error=error_result) )

  # return all results in list (convert teval to units of seconds)
  return( list(result = lik_result,
               teval = median(mb_lik[['time']]) / 1e9,
               error = error_result) )
}

# time kriging prediction and variance with snapKrig, geoR, fields, gstat
bench_kriging = function(g, gdim, pars, pkg, out='p', n_rep=10, timeout=5*60)
{
  # g: sk grid or sf POINTS data frame containing the observed data
  # gdim: the desired output resolution (points in g are snapped to a grid of this size)
  # pars: the (package-specific) results object returned from bench_fit
  # pkg: one of 'snapKrig', 'geoR', 'fields', 'gstat'
  # out: either 'p' (prediction) or 'v' (variance)
  # n_rep: number of repetitions
  # timeout: maximum time (seconds) allowed for all repetitions

  # Note: only snapKrig and fields allow variance to be computed separately from the
  # predictions. For the others, we compute both in prediction call (out='p'), and in
  # variance calls (out='v') we compute nothing and return an error.

  # catch and return any errors
  error_result = tryCatch({

    # handle gridded versus irregular point inputs
    if(inherits(g, 'sk'))
    {
      # find the required up/down-scaling factor
      scale_factor = gdim / dim(g)

      # if scaling factor is 1, 1 output is same grid as input
      if( all(scale_factor == 1) ) g_out = g
      if( any(scale_factor < 1) ) g_out = sk_rescale(g, up=round(1/scale_factor))
      if( any(scale_factor > 1) ) g_out = sk_rescale(g, down=round(scale_factor))

    } else {

      # irregular point case: snap to a grid of the desired resolution
      g_out = sk_snap(g, gdim)
    }

    # output grid lines
    xy = g_out[['gyx']][c('x', 'y')]

    # prepare arguments for other packages
    if(pkg != 'snapKrig')
    {
      if(inherits(g, 'sk')) g = sk_coords(g, out='sf', na_omit=TRUE, quiet=TRUE)
      g_coords = sf::st_coordinates(g)
      g_vals = sf::st_drop_geometry(g)
    }

    # prepare arguments for gstat
    if(pkg == 'gstat')
    {
      # make a data frame out of points, then combine with pars to export to gstat object
      df_gstat = data.frame(g_vals, g_coords)
      names(df_gstat) = c('z', 'x', 'y')
      pts_gstat = gstat::gstat(formula=z~1, locations=~x+y, data=df_gstat, model=pars)
    }

    # prepare arguments for geoR
    if(pkg=='geoR')
    {
      # export observed data to geodata list and copy starting parameters
      pts_geoR = geoR::as.geodata(cbind(g_coords, g_vals))
    }

    # interrupt evaluation when timeout is exceeded
    mb_krig = withTimeout({

      # repeat evaluations several times to get a median time
      microbenchmark({

        # different likelihood calls for different packages
        krig_result = switch(pkg,

                             # fast and direct
                             'snapKrig' = sk_cmean(g_out, pars, what=out, quiet=TRUE),

                             # prediction and variance computed in the same call
                             'gstat' = terra::interpolate(sk_export(g_out), model=pts_gstat),

                             # prediction and variance computed in the same call
                             'geoR' = geoR::krige.conv(pts_geoR,
                                                       loc = expand.grid(xy),
                                                       krige = krige.control(obj.model=pars)),

                             # fields allows separate calls but they are handled by different functions
                             'fields' = if(out=='p') {

                               fields::predictSurface(pars, grid.list=xy, extrap=TRUE)

                             } else {

                               fields::predictSurfaceSE(pars, grid.list=xy, extrap=TRUE)
                             }
        ) # end switch

      }, times=n_rep) # end microbenchmark repetitions

    }, timeout=timeout, onTimeout='error') # end timeout

    NULL # return value from trycatch if the execution succeeds

  }, error = identity) # end tryCatch

  # return from failed execution
  if( !is.null(error_result) ) return( list(result=NA, teval=NA, error=error_result) )
  if( is.null(krig_result) ) return( list(result=NA, teval=NA, error='NULL result') )

  # export non-snapKrig results to SpatRaster (gstat produces a SpatRaster by default)
  if(pkg %in% c('fields', 'geoR'))
  {
    # extract the data vector (which is named differently in the two packages)
    result_nm = c(fields='z', geoR='predict')
    z = krig_result[[result_nm[pkg]]]

    # flip to correct for different vectorization ordering
    krig_result = terra::flip(terra::t(rast(matrix(c(z), sapply(xy, length)), crs=g_out[['crs']])))
    ext(krig_result) = c(range(xy[[1]]), range(xy[[2]]))
  }

  # return all results in list (convert teval to units of seconds)
  return( list(result = sk(krig_result),
               teval = median(mb_krig[['time']]) / 1e9,
               error = error_result) )
}

# define a function to open all results files and combine them
make_all_results = function(csv_paths)
{
  csv_paths_ready = csv_paths[file.exists(csv_paths)]
  csv_list = lapply(csv_paths_ready, read.csv)
  do.call(rbind, csv_list)[,-1]
}

# define a function to create a plot of results from a given example
make_results_plot = function(all_results, eg_nm)
{
  # log-log plot of all likelihood timing results
  gg0 = ggplot(all_results) +
    aes(x=n_in, y=teval_lik, color=pkg, lty=complete) +
    geom_point() +
    geom_line() +
    geom_line() +
    xlab('input points') +
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

  # log-log plot of prediction and variance timing results for example eg_nm
  results_plot_df = subset(all_results, name==eg_nm)

  # make a plotting data frame with single column for both times
  n_plot = nrow(results_plot_df)
  results_plot_df = results_plot_df[rep(seq(n_plot), 2),]
  results_plot_df[['with_var']] = rep(c(TRUE, FALSE), each=n_plot)
  results_plot_df[['teval']] = results_plot_df[['teval_pred']]
  results_plot_df[['teval']][seq(n_plot)] = results_plot_df[['teval_both']][seq(n_plot)]

  # create prediction + variance time plot
  gg1 = ggplot(results_plot_df) +
    aes(x=n_out, y=teval, color=pkg, lty=with_var) +
    geom_point() +
    geom_line() +
    geom_line() +
    xlab('output points') +
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

  gridExtra::grid.arrange(gg0, gg1, nrow=1)
}
