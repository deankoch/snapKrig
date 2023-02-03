#' ---
#' title: "Cross validation with snapKrig"
#' author: "Dean Koch"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#'
#' This performs a simple cross-validation experiment on the Meuse dataset and
#' saves the results to a CSV table, "cv_results.csv".
#'
#' On my PC this takes about 10-20 minutes to complete.
#'
#' ## Dependencies
#'
# essential GIS packages
library(terra)
library(sf)

# data sources
library(sp)

# snapKrig
library(snapKrig)

# for managing working directories
library(here)

# define the main output file
dir_project = here('rjarticle')
dir_storage = file.path(dir_project, 'data')
path_main_result = file.path(dir_storage, 'cv_results.csv')

# load helpers
source(file.path(dir_storage, 'data_helpers.R'))


#'
#' ## Repeat the data prep workflow in the rjarticle markdown script
#'
#' This opens the Meuse log zinc data, snaps it to a grid, then computes
#' distance to river for each grid point.

# get Meuse data in a list
meuse_list = get_meuse(dfMaxLength=NA)

# copy the log zinc points data
pts = meuse_list[['soils']]['log_zinc']

# snap log zinc data to grid of specified resolution
g = sk_snap(pts, g=list(gres=c(y=50, x=50)))

# measure distances for every point in the grid
river_dist = sf::st_distance(sk_coords(g, out='sf'), meuse_list[['river_line']])
river_dist = units::drop_units(river_dist)

# include both distance and its square root
g_X = g
g_X[] = scale(river_dist)


# (optionally) check likelihood function
#
# fit_result_ok = sk_fit(g, X=NA, iso=TRUE, quiet=TRUE)
# fit_result_uk = sk_fit(g, X=g_X, iso=TRUE, quiet=TRUE)
#
# sk_LL(fit_result_ok, g, X=NA, quiet=TRUE, out='m')
# sk_LL(fit_result_uk, g, X=g_X, quiet=TRUE, out='m')
# sk_LL(fit_result_uk, g, X=g_X, quiet=TRUE, out='b')

#'
#' ## Define the cross-validation script
#'
#' This randomly selects folds and runs model-fitting and prediction. It tests the following
#' covariance models: OK with defaults (Gaussian isotropic), UK with defaults (Gaussian isotropic),
#' UK with anisotropic Gaussian, UK with mat x mat covariance (anisotropic).
#'
#' The function splits the data into n_fold folds (n_rep times) and runs each of the kriging
#' workflows, computing MSPE and MSDR (residuals scaled by variance) on the test set, as well as
#' the AIC and BIC values on the training set. It returns the averages over folds in a data frame,
#' as well as square roots (prefix 'r') and values after transformation to the original scale
#' (suffix 'b'), where applicable.
#'

# define a function to execute cross validation
run_cv = function(g, g_X, n_fold=5, n_rep=25)
{
  # g and g_X should be sk grid list objects, containing the response data and covariates matrix.

  # correlograms
  cg_nm = c('gau', 'mat', 'gxp', 'sph', 'exp')

  # count covariates and describe models to fit
  n_lm = ncol(g_X[])
  out_df = data.frame(name = paste0('fit_result_', c('ok_iso', 'uk_iso', paste0('uk_', cg_nm))),
                      covariance = c(rep('gau', 2), cg_nm),
                      covariates = c(0, rep(n_lm, 1+length(cg_nm))),
                      isotropic = c(TRUE, TRUE, rep(FALSE, length(cg_nm))),
                      parameters = 3L + c(1, 1, 2, 4, 4, 2, 2),
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
      X = NA
      if(out_df[j, 'covariates'] > 0) X = g_X

      # fit and compute information scores on training set
      fit_result = sk_fit(g_train, pars=pars, X=X, iso=iso, quiet=TRUE)
      aic_train = sk_LL(fit_result, g_train, X=X, quiet=TRUE, out='a')
      bic_train = sk_LL(fit_result, g_train, X=X, quiet=TRUE, out='b')

      # predict and compute variance, transform to original scale
      g_pred = sk_cmean(g, pars=fit_result, X=X)
      g_var = sk_cmean(g, pars=fit_result, X=X, what='v', quiet=TRUE)
      g_pred_b = exp(g_pred + g_var/2)

      # standardized residuals on test set
      g_res = g_pred - g
      z_res = g_res[idx_test]
      z_res_std = g_res[idx_test] / sqrt(g_var[idx_test])
      z_res_b = g_pred_b[idx_test] - exp(g[idx_test])

      # compute mean squares, store results
      results_list[[i]][j, 'MSPE'] = mean(z_res^2)
      results_list[[i]][j, 'MSDR'] = mean(z_res_std^2)
      results_list[[i]][j, 'MSPEb'] = mean(z_res_b^2)
      results_list[[i]][j, 'AIC'] = aic_train
      results_list[[i]][j, 'BIC'] = bic_train
    }
  }

  # take means of the means
  out_df['MSDR'] = rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSDR'])))
  out_df['MSPE'] = rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSPE'])))
  out_df['MSPEb'] = rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSPEb'])))
  out_df['AIC'] = rowMeans(do.call(cbind, lapply(results_list, function(x) x['AIC'])))
  out_df['BIC'] = rowMeans(do.call(cbind, lapply(results_list, function(x) x['BIC'])))

  # take square roots to get rMS... statistics
  out_df['rMSDR'] = sqrt(rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSDR']))))
  out_df['rMSPE'] = sqrt(rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSPE']))))
  out_df['rMSPEb'] = sqrt(rowMeans(do.call(cbind, lapply(results_list, function(x) x['MSPEb']))))

  # format model names and add number of parameters
  out_df['covariance'] = sapply(out_df['covariance'], function(k) paste(k, 'x', k))
  out_df['parameters'] = out_df['covariates'] + out_df['parameters']

  close(pb)
  return(out_df)
}


#'
#' ## Run the cross-validation script and write results to file
#'
#' This part is slow. Expect about 20 minutes wait before the script completes.
#' Note that the "cv_results.csv" is overwritten by write.csv

# set random seed for reproducibility
seed_int = sum(as.integer(charToRaw('snapKrig_CV_experiment')))
set.seed(seed=seed_int)

# call the function
cv_results = run_cv(g, g_X, n_fold=5, n_rep=25)
print(cv_results)

# write results to file (this overwrites old results!)
write.csv(cv_results, path_main_result)





