

# copy of the fitted parameters for Meuse, to display them at the beginning
meuse_OK_fitted_pars_iso = list(y=list(k='gau', kp=c(rho=590.4317)),
                                x=list(k='gau', kp=c(rho=590.4317)),
                                eps=0.1192025,
                                psill=0.9024344)


# load the Meuse data into a convenient format
get_meuse = function(dfMaxLength = units::set_units(50, m))
{
  # Note: dfMaxLength sets the interval used to sample line geometries of the river
  # using Voronoi tiles. This is a fussy and not well-tested algorithm for finding the
  # centre line of a river polygon, but it seems to work well enough for the example here

  # EPSG code for the coordinate system
  epsg_meuse = 28992

  # open river location data
  utils::data(meuse.riv)
  crs_meuse = sf::st_crs(epsg_meuse)[['wkt']]

  # reshape the river (edge) point data as a more densely segmented polygon
  colnames(meuse.riv) = c('x', 'y')
  meuse_river_points = sf::st_as_sf(as.data.frame(meuse.riv), coords=c('x', 'y'), crs=crs_meuse)
  meuse_river_seg = sf::st_cast(sf::st_combine(meuse_river_points), 'LINESTRING')
  meuse_river_poly = sf::st_cast(st_segmentize(meuse_river_seg, dfMaxLength), 'POLYGON')

  # skeletonization trick to get a single linestring at center of the river
  meuse_river_voronoi = sf::st_cast(sf::st_voronoi(meuse_river_poly, bOnlyEdges=TRUE), 'POINT')
  meuse_river_skele = sf::st_intersection(meuse_river_voronoi, meuse_river_poly)
  n_skele = length(meuse_river_skele)

  # compute distance matrix
  dmat_skele = units::drop_units(sf::st_distance(meuse_river_skele))

  # re-order to start from northernmost point
  idx_first = which.max(st_coordinates(meuse_river_skele)[,2])
  idx_reorder = c(idx_first, integer(n_skele-1L))
  for(idx_skele in seq(n_skele-1L))
  {
    # find least distance match
    idx_tocheck = seq(n_skele) != idx_first
    idx_first = which(idx_tocheck)[ which.min(dmat_skele[idx_tocheck, idx_first]) ]
    idx_reorder[1L+idx_skele] = idx_first

    # modify distance matrix so the matching point is not selected again
    dmat_skele[idx_first, ] = Inf
  }

  # connect the points to get the spine
  meuse_river = sf::st_cast(sf::st_combine(meuse_river_skele[idx_reorder]), 'LINESTRING')

  # load soil points data
  utils::data(meuse)
  meuse_soils = sf::st_as_sf(meuse, coords=c('x', 'y'), crs=epsg_meuse)

  # add 'distance' (to river) and 'logzinc' columns
  meuse_soils[['distance']] = units::drop_units( sf::st_distance(meuse_soils, meuse_river))
  meuse_soils[['log_zinc']] = log(meuse_soils[['zinc']])

  # crop the river objects to buffered bounding box of soils data
  bbox_padded = st_buffer(sf::st_as_sfc(sf::st_bbox(meuse_soils)), units::set_units(500, m))
  meuse_river_poly = sf::st_crop(meuse_river_poly, bbox_padded)
  meuse_river = sf::st_crop(meuse_river, bbox_padded)

  # return three geometry objects in a list
  return( list(soils=meuse_soils, river_poly=meuse_river_poly, river_line=meuse_river) )
}

# run a short CV analysis on Meuse: this randomly selects folds and runs model-fitting and prediction
run_cv = function(g, g_X, n_fold=5, n_rep=25)
{
  # g and X should be sk grid list objects, containing the response data and covariates matrix.
  # workflows tested: OK with defaults, UK with defaults, UK with anisotropy, UK with mat x mat
  # The function splits the data into n_fold folds (n_rep times) and runs kriging on each to
  # compute MSPE error on the test set

  # count covariates and describe models to fit
  X = g_X[]
  n_lm = ncol(X)
  out_df = data.frame(name = paste0('fit_result_', c('ok', 'uk', 'uk_gau', 'uk_mat')),
                      covariance = c(rep('gau', 3), 'mat'),
                      covariates = c(0, rep(n_lm, 3)),
                      isotropic = c(TRUE, TRUE, FALSE, FALSE),
                      parameters = 1 + 2 + c(1, 1, 2, 4),
                      MSDR = 0,
                      MSPE = 0,
                      MSPEb = 0,
                      AIC = 0,
                      BIC = 0)

  # count observations and number per fold
  is_obs = !is.na(g)
  n_obs = sum(is_obs)
  n_per = floor(n_obs/n_fold)

  # build a list of test sets
  idx_test_list = do.call(c, lapply(seq(n_rep), function(ii) {

    # randomly shuffle observed points
    idx_obs = which(is_obs)[sample(n_obs)]
    lapply(seq(n_fold), function(jj) idx_obs[n_per*(jj-1) + seq(n_per)])

  }))

  # run folds loop to fill the list
  results_list = rep(list(out_df), n_fold * n_rep)
  pb = utils::txtProgressBar(0, n_fold * n_rep, style=3)
  for(i in seq(n_fold * n_rep))
  {
    setTxtProgressBar(pb, i)

    # copy training objects
    is_obs_train = is_obs
    g_train = g

    # wipe test data
    idx_test = idx_test_list[[i]]
    is_obs_train[idx_test] = FALSE
    g_train[idx_test] = NA

    # loop over models
    for(j in seq(nrow(out_df)))
    {
      # model setup
      pars = out_df[j, 'covariance']
      iso = out_df[j, 'isotropic']
      X_train = X_test = NA
      if(out_df[j, 'covariates'] > 0)
      {
        X_train = X[is_obs_train,]
        X_test = X
      }

      # fit and compute information scores on training set
      fit_result = sk_fit(g_train, pars=pars, X=X_train, iso=iso, quiet=TRUE)
      aic_train = sk_LL(fit_result[['pars']], g_train, X=X_train, quiet=TRUE, out='a')
      bic_train = sk_LL(fit_result[['pars']], g_train, X=X_train, quiet=TRUE, out='b')

      # predict and compute variance, transform to original scale
      g_pred = sk_cmean(g, pars=fit_result[['pars']], X=X_test)
      g_var = sk_cmean(g, pars=fit_result[['pars']], X=X_test, what='v', quiet=TRUE)
      g_pred_b = exp(g_pred + g_var/2)

      # standardized residuals on test set
      g_res = g_pred - g
      z_res = g_res[idx_test]
      z_res_std = g_res[idx_test] / sqrt(g_var[idx_test])
      z_res_b = g_pred_b[idx_test] - exp(g[idx_test])

      # compute mean squares store results
      results_list[[i]][j, 'MSPE'] = mean(z_res^2)
      results_list[[i]][j, 'MSDR'] = mean(z_res_std^2)
      results_list[[i]][j, 'MSPEb'] = mean(z_res_b^2)
      results_list[[i]][j, 'AIC'] = aic_train
      results_list[[i]][j, 'BIC'] = bic_train
    }
  }

  # take means of the means to get the mean
  out_df['MSDR'] = rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSDR'])))
  out_df['MSPE'] = rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSPE'])))
  out_df['MSPEb'] = rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSPEb'])))
  out_df['AIC'] = rowMeans(do.call(cbind, lapply(results_list, function(x) x['AIC'])))
  out_df['BIC'] = rowMeans(do.call(cbind, lapply(results_list, function(x) x['BIC'])))

  # take square roots
  out_df['rMSDR'] = sqrt(rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSDR']))))
  out_df['rMSPE'] = sqrt(rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSPE']))))
  out_df['rMSPEb'] = sqrt(rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSPEb']))))

  # format model names and add number of parameters
  out_df['covariance'] = sapply(out_df['covariance'], function(k) paste(k, 'x', k))
  out_df['parameters'] = out_df['covariates'] + out_df['parameters']

  close(pb)
  return(out_df)
}

# load the results from the benchmarking vignette


