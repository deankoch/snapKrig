
# fits a model to grid g; returns parameters and execution time, or any errors
# (pts and initial argument ignored)
sk_bench_fit = function(g, pts, n_rep=10, initial=NULL, timeout=60)
{
  # pts argument ignored

  # catch and return any errors (otherwise return NA)
  error_result = tryCatch({

    # microbenchmark for median execution time
    mb_fit = withTimeout({

      microbenchmark({

        # model fitting
        sk_OK_result = sk_fit(g, quiet=TRUE, iso=TRUE)

      }, times=n_rep)

    }, timeout=timeout, onTimeout='error')
    NULL

  }, error = identity)

  # return from failed executions
  if( !is.null(error_result) ) return( list(result=NULL, teval=NULL, error=error_result) )

  # return in list
  return( list(result = sk_OK_result,
               teval = median(mb_fit[['time']]) / 1e9,
               error = error_result) )
}

# runs ordinary kriging on g with parameters pars; returns execution time, or any errors
# (pts argument ignored)
sk_bench_OK = function(g, pts, pars, n_rep=10, out='p', timeout=60)
{
  # out = 'p' for prediction or 'v' for variance

  # catch and return any errors (otherwise return NA)
  error_result = tryCatch({

    # model fitting
    mb_pred = withTimeout({

      microbenchmark({

        g_predicted = sk_cmean(g, pars[['pars']], out=out, quiet=TRUE)

      }, times=n_rep)

    }, timeout=timeout, onTimeout='error')
    NULL

  }, error = identity)

  # return from failed executions
  if( !is.null(error_result) ) return( list(result=NULL, teval=NULL, error=error_result) )

  # output as terra SpatRaster object in list
  return( list(result = g_predicted,
               teval = median(mb_pred[['time']]) / 1e9,
               error = error_result) )
}


#' ## geoR implementation
#'

# fits a model to pts ; returns parameters and execution time, or any errors
geoR_bench_fit = function(g, pts, n_rep=10, initial=NULL, timeout=60)
{
  # the dimensions and resolution of g are used to estimate initial parameter values.
  # Argument initial can be supplied instead to override these defaults

  # geoR has likelihood-based fitting but requires initial values from the user
  if( is.null(initial) ) initial = sk_bds(sk_pars(1), g)[c('psill', 'y.rho'), 'initial']

  # catch and return any errors (otherwise return NA)
  error_result = tryCatch({

    # unpack observed points data and coordinates and create a geodata list
    coords = st_coordinates(pts)
    vals = st_drop_geometry(pts)
    pts_geoR = as.geodata(cbind(coords, vals))

      # microbenchmark for median execution time
      mb_fit = withTimeout({

        microbenchmark({

          # model fitting
          geoR_OK_result = geoR::likfit(pts_geoR, ini.cov.pars = initial, cov.model = 'gaussian')

        }, times=n_rep)

    }, timeout=timeout, onTimeout='error')
    NULL

  }, error = identity)

  if( !is.null(error_result) ) return( list(result=NULL, teval=NULL, error=error_result) )

  # return in list
  return( list(result = geoR_OK_result,
               teval = median(mb_fit[['time']]) / 1e9,
               error = error_result) )
}

# runs ordinary kriging from pts onto g with model pars; returns execution time, or any errors
# (out ignored)
geoR_bench_OK = function(g, pts, pars, n_rep=10, out='p', timeout=60)
{
  # catch and return any errors (otherwise return NA)
  error_result = tryCatch({

    # supply prediction surface grid lines in expected order for fields
    xy = g[['gyx']][c('x', 'y')]
    geoR_loc = expand.grid(xy)

    # unpack observed points data and coordinates and create a geodata list
    coords = st_coordinates(pts)
    vals = st_drop_geometry(pts)
    pts_geoR = as.geodata(cbind(coords, vals))

    # model fitting
    mb_pred = withTimeout({

      microbenchmark({

        # predictions AND VARIANCE happen in the one call - we only want to do this once
        if(out == 'v') stop('see output for prediction (includes compute time for variance)')
        geoR_z = krige.conv(pts_geoR, loc = geoR_loc, krige = krige.control(obj.model=pars))

      }, times=n_rep)

    }, timeout=timeout, onTimeout='error')
    NULL

  }, error = identity)
  if( !is.null(error_result) ) return( list(result=NULL, teval=NULL, error=error_result) )

  # convert to terra object
  z = rast(matrix(c(geoR_z[['predict']]), rev(g[['gdim']])), crs=g[['crs']]) |> t() |> flip()
  ext(z) = c(range(xy[[1]]), range(xy[[2]]))

  # output as snapKrig object in list
  return( list(result = sk(z),
               teval = median(mb_pred[['time']]) / 1e9,
               error = error_result) )
}

#' ## gstat implementation
#'

# fits a model to pts; returns parameters and execution time, or any errors
# (g and initial arguments ignored)
gstat_bench_fit = function(g, pts, n_rep=10, initial=NULL, timeout=60)
{
  # catch and return any errors (otherwise return NA)
  error_result = tryCatch({

    # simplify sf object to get same data in data.frame (use lowercase x, y)
    pts_df = data.frame(st_drop_geometry(pts), st_coordinates(pts))
    names(pts_df) = c('z', 'x', 'y')

    # microbenchmark for median execution time
    mb_fit = withTimeout({

        microbenchmark({

        # fit a variogram by least squares with the default initial values
        pts_vgm = variogram(z~1, locations=~x+y, data=pts_df)
        gstat_pars_df = fit.variogram(pts_vgm, model=vgm('Gau'))

      }, times=n_rep)

    }, timeout=timeout, onTimeout='error')
    NULL

  }, error = identity)
  if( !is.null(error_result) ) return( list(result=NULL, teval=NULL, error=error_result) )

  # return in list
  return( list(result = gstat_pars_df,
               teval = median(mb_fit[['time']]) / 1e9,
               error = error_result) )
}

# runs ordinary kriging from pts onto g with model pars; returns execution time, or any errors
gstat_bench_OK = function(g, pts, pars, n_rep=10, out='p', timeout=60)
{
  # catch and return any errors (otherwise return NA)
  error_result = tryCatch({

    # simplify sf object to get same data in data.frame (use lowercase x, y)
    pts_df = data.frame(st_drop_geometry(pts), st_coordinates(pts))
    names(pts_df) = c('z', 'x', 'y')
    pts_gstat = gstat(formula=z~1, locations=~x+y, data=pts_df, model=pars)

    # benchmarking
    mb_pred = withTimeout({

      microbenchmark({

        # predictions AND VARIANCE happen in the one call - we only want to do this once
        if(out == 'v') stop('see output for prediction (includes compute time for variance)')
        gstat_pred_rast = interpolate(sk_export(g), model=pts_gstat)

      }, times=n_rep)

    }, timeout=timeout, onTimeout='error')
    NULL

  }, error = identity)
  if( !is.null(error_result) ) return( list(result=NULL, teval=NULL, error=error_result) )

  # output as terra SpatRaster object in list - NOTE we omit variance from return
  return( list(result = sk(gstat_pred_rast[[1]]),
               teval = median(mb_pred[['time']]) / 1e9,
               error = error_result) )
}


#' ## fields implementation
#'

# this is a wrapper for `fields::Exp.cov` to fix a bug in `fields::spatialProcess`
my_Expcov = function(x1, x2 = NULL, aRange = 1, p = 1, distMat = NA, C = NA,
                     marginal = FALSE, onlyUpper = FALSE, theta = NULL,
                     Distance=NULL, Dist.args=NULL)
{

  # this function simply accepts and ignores the problematic arguments Distance, Dist.args
  Exp.cov(x1, x2=x2, aRange = aRange, p=p, distMat = distMat,
          C = C, marginal = marginal, onlyUpper=onlyUpper, theta=theta)


}

# fits a model to pts; returns parameters and execution time, or any errors
fields_bench_fit = function(g, pts, n_rep=10, initial=NULL, timeout=60)
{
  # argument g is ignored

  # catch and return any errors (otherwise return NA)
  error_result = tryCatch({

    # simplify sf object to get same data in data.frame (use lowercase x, y)
    coords = st_coordinates(pts)
    vals = st_drop_geometry(pts)

    # microbenchmark for median execution time
    mb_fit = withTimeout({

      microbenchmark({

        # no initial values supplied
        if( is.null(initial) )
        {
          # fit the model by REML
          fields_result = spatialProcess(x=coords, y=vals,
                                         cov.function='my_Expcov',
                                         cov.args=list(p=2),
                                         Dist.args=list(compact=FALSE))
        } else {

          # pass range from initial values
          fields_result = spatialProcess(x=coords, y=vals,
                                         cov.function='my_Expcov',
                                         cov.args=list(p=2),
                                         cov.params.start=list(aRange=initial['rho']),
                                         Dist.args=list(compact=FALSE))
        }

      }, times=n_rep)

    }, timeout=timeout, onTimeout='error')
    NULL

  }, error = identity)
  if( !is.null(error_result) ) return( list(result=NULL, teval=NULL, error=error_result) )

  # return in list
  return( list(result = fields_result,
               teval = median(mb_fit[['time']]) / 1e9,
               error = error_result) )
}

# runs ordinary kriging from pts onto g with model pars; returns execution time, or any errors
fields_bench_OK = function(g, pts, pars, n_rep=10, out='p', timeout=60)
{
  # argument pts is ignored (training data info is in pars)

  # catch and return any errors (otherwise return NA)
  error_result = tryCatch({

    # supply prediction surface grid lines in expected order for fields
    xy = g[['gyx']][c('x', 'y')]

    # benchmarking
    mb_pred = withTimeout({

      microbenchmark({

        # separate prediction and variance calls
        if(out=='p') fields_gridlist = fields::predictSurface(pars, grid.list=xy, extrap=TRUE)
        if(out=='v') fields_gridlist = fields::predictSurfaceSE(pars, grid.list=xy, extrap=TRUE)

      }, times=n_rep)

    }, timeout=timeout, onTimeout='error')
    NULL

  }, error = identity)
  if( !is.null(error_result) ) return( list(result=NULL, teval=NULL, error=error_result) )

  # convert to terra object
  z = rast(matrix(c(fields_gridlist[['z']]), sapply(xy, length)), crs=g[['crs']]) |> t() |> flip()
  ext(z) = c(range(xy[[1]]), range(xy[[2]]))

  # output as terra SpatRaster object in list
  return( list(result = sk(z),
               teval = median(mb_pred[['time']]) / 1e9,
               error = error_result))
}


# fits a model to pts; returns parameters and execution time, or any errors
RandomFields_bench_fit = function(g, pts, n_rep=10, initial=NULL, timeout=60)
{
  # argument g is ignored

  # catch and return any errors (otherwise return NA)
  error_result = tryCatch({

    # simplify sf object to get a copy of data in dataframe
    rf_data = cbind(st_coordinates(pts), st_drop_geometry(pts))

    # microbenchmark for median execution time
    mb_fit = withTimeout({

      microbenchmark({

        # define the model then fit it by MLE
        rf_model = RMgauss(var=NA, scale=NA) + RMnugget(var=NA) + RMtrend(mean=NA)
        rf_result = suppressMessages(RFfit(rf_model, data=rf_data))

      }, times=n_rep)

    }, timeout=timeout, onTimeout='error')
    NULL

  }, error = identity)
  if( !is.null(error_result) ) return( list(result=NULL, teval=NULL, error=error_result) )

  # return in list
  return( list(result = rf_result,
               teval = median(mb_fit[['time']]) / 1e9,
               error = error_result) )
}

# runs ordinary kriging from pts onto g with model pars; returns execution time, or any errors
RandomFields_bench_OK = function(g, pts, pars, n_rep=10, out='p', timeout=60)
{
  # argument pts is ignored (training data info is in pars)

  # catch and return any errors (otherwise return NA)
  error_result = tryCatch({

    # get coordinates of points on output grid
    xy = g[['gyx']][c('x', 'y')]
    xy_all = expand.grid(xy)

    # simplify sf object to get a copy of data in dataframe
    rf_data = cbind(st_coordinates(pts), st_drop_geometry(pts))

    # benchmarking
    mb_pred = withTimeout({

      microbenchmark({

        # option return_variance=TRUE currently causes 'Error in predictGauss...' so we skip it here
        if(out=='v') stop('variance not supported at this time')
        rf_output = suppressMessages(RFinterpolate(pars, x=xy_all[['x']], y=xy_all[['y']], data=rf_data))

      }, times=n_rep)

    }, timeout=timeout, onTimeout='error')
    NULL

  }, error = identity)
  if( !is.null(error_result) ) return( list(result=NULL, teval=NULL, error=error_result) )

  # convert to terra object
  z = rast(matrix(c(rf_output[[1]]), sapply(xy, length)), crs=g[['crs']]) |> t() |> flip()
  ext(z) = c(range(xy[[1]]), range(xy[[2]]))

  # output as terra SpatRaster object in list
  return( list(result = sk(z),
               teval = median(mb_pred[['time']]) / 1e9,
               error = error_result))
}
