# bk_model.R
# Dean Koch, 2022
# Functions for modeling with blitzKrig

#' Likelihood of covariance model `pars` given data `g_obs`
#'
#' Computes the log-likelihood for the Gaussian process covariance model `pars`,
#' given 2-dimensional grid data `g_obs`, and, optionally, linear trend data in `X`.
#'
#' The function evaluates:
#'
#' `-log( 2 * pi ) - ( 1/2 ) * ( log_det + quad_form )`,
#'
#' where `log_det` is the logarithm of the determinant of covariance matrix V, and
#' `quad_form` is z^T V^{-1} z, for the observed response vector z. This z is constructed
#' by subtracting the trend specified in `X` (if any) from the non-NA values in `g_obs$gval`.
#'
#' If the trend is known, it can be supplied in argument `X` as a numeric scalar or vector of
#' length equal to the number of non-NA values in `g_obs$gval`, in matching order. Equivalently,
#' users can simply subtract the trend from `g_obs` beforehand and set `X=0` (the default).
#' If the trend is unknown, the function optionally estimates it by GLS using the model
#' `pars`. To estimate a spatially constant mean, set `X=NA`. To estimate a spatially variable
#' mean, supply linear predictors as columns of a matrix argument to `X` (see `bk_GLS`).
#'
#' `fac_method` specifies how to factorize V, either using the Cholesky factor ('chol')
#' or eigen-decomposition ('eigen'). A pre-computed factorization `fac` can be supplied by
#' first calling `bk_var(..., scaled=TRUE)` (in which case `fac_method` is ignored).
#'
#' When `more=TRUE`, the function returns a list of diagnostics: a count of the
#' number of observations, the likelihood function value, and its two major components; the
#' log-determinant `log_det`, and the quadratic form `quad_form`.
#'
#'
#' @param pars list of form returned by `bk_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param g_obs list of form returned by `bk_grid` (with entries 'gdim', 'gres', 'gval')
#' @param X numeric, vector, matrix, or NA, a fixed mean value, or matrix of linear predictors
#' @param fac_method character, the factorization to use: 'chol' (default) or 'eigen'
#' @param fac matrix or list, (optional) pre-computed covariance factorization
#' @param more logical, indicates to return list with likelihood components
#' @param quiet logical indicating to suppress console output
#'
#' @return numeric, the likelihood of `pars` given `g_obs` and `X`, or list (if `more=TRUE`)
#' @export
#'
#' @examples
#' # set up example grid, covariance parameters
#' gdim = c(25, 12)
#' n = prod(gdim)
#' g_all = bk_grid(gdim)
#' pars = modifyList(bk_pars(g_all, 'gau'), list(psill=0.7, eps=5e-2))
#'
#' # generate some covariates and complete data
#' n_betas = 3
#' betas = rnorm(n_betas)
#' X_all = cbind(1, matrix(rnorm(n*(n_betas-1)), n))
#' z = as.vector( bk_sim(g_all) + (X_all %*% betas) )
#' g_all[['gval']] = z
#'
#' # two methods for likelihood
#' LL_chol = bk_LL(pars, g_all, fac_method='chol')
#' LL_eigen = bk_LL(pars, g_all, fac_method='eigen')
#'
#' # compare to working directly with matrix inverse
#' V = bk_var(g_all, pars, fac_method='none', sep=FALSE)
#' V_inv = chol2inv(chol(V))
#' quad_form = as.numeric( t(z) %*% crossprod(V_inv, z) )
#' log_det = as.numeric( determinant(V, logarithm=TRUE) )[1]
#' LL_direct = (-1/2) * ( n * log( 2 * pi ) + log_det + quad_form )
#'
#' # relative errors
#' abs( LL_direct - LL_chol ) / max(LL_direct, LL_chol)
#' abs( LL_direct - LL_eigen ) / max(LL_direct, LL_eigen)
#'
#' # repeat with pre-computed variance factorization
#' fac_eigen = bk_var(g_all, pars, fac_method='eigen', sep=TRUE)
#' bk_LL(pars, g_all, fac=fac_eigen) - LL_eigen
#'
#' # repeat with most data missing
#' n_obs = 50
#' idx_obs = sort(sample.int(n, n_obs))
#' z_obs = g_all$gval[idx_obs]
#' g_obs = modifyList(g_all, list(gval=rep(NA, n)))
#' g_obs[['gval']][idx_obs] = z_obs
#' LL_chol_obs = bk_LL(pars, g_obs, fac_method='chol')
#' LL_eigen_obs = bk_LL(pars, g_obs, fac_method='eigen')
#'
#' # working directly with matrix inverse
#' V_obs = bk_var(g_obs, pars, fac_method='none')
#' V_obs_inv = chol2inv(chol(V_obs))
#' quad_form_obs = as.numeric( t(z_obs) %*% crossprod(V_obs_inv, z_obs) )
#' log_det_obs = as.numeric( determinant(V_obs, logarithm=TRUE) )[1]
#' LL_direct_obs = (-1/2) * ( n_obs * log( 2 * pi ) + log_det_obs + quad_form_obs )
#' abs( LL_direct_obs - LL_chol_obs ) / max(LL_direct_obs, LL_chol_obs)
#' abs( LL_direct_obs - LL_eigen_obs ) / max(LL_direct_obs, LL_eigen_obs)
#'
#' # again using a pre-computed variance factorization
#' fac_chol_obs = bk_var(g_obs, pars, fac_method='chol', scaled=TRUE)
#' fac_eigen_obs = bk_var(g_obs, pars, fac_method='eigen', scaled=TRUE)
#' bk_LL(pars, g_obs, fac=fac_chol_obs) - LL_chol_obs
#' bk_LL(pars, g_obs, fac=fac_eigen_obs) - LL_eigen_obs
#'
#' # copy covariates (don't pass the intercept column in X)
#' X = X_all[idx_obs, -1]
#'
#' # use GLS to de-trend, with and without covariatea
#' g_detrend_obs = g_detrend_obs_X = g_obs
#' g_detrend_obs[['gval']][idx_obs] = z_obs - bk_GLS(g_obs, pars)
#' g_detrend_obs_X[['gval']][idx_obs] = z_obs - bk_GLS(g_obs, pars, X, out='z')
#'
#' # pass X (or NA) to bk_LL to do this automatically
#' LL_detrend_obs = bk_LL(pars, g_detrend_obs)
#' LL_detrend_obs_X = bk_LL(pars, g_detrend_obs_X)
#' LL_detrend_obs - bk_LL(pars, g_obs, X=NA)
#' LL_detrend_obs_X - bk_LL(pars, g_obs, X=X)
#'
#' # equivalent sparse input specification
#' idx_grid = match(seq(n), idx_obs)
#' g_sparse = modifyList(g_all, list(gval=matrix(z_obs, ncol=1), idx_grid=idx_grid))
#' LL_chol_obs - bk_LL(pars, g_sparse)
#' LL_eigen_obs - bk_LL(pars, g_sparse)
#' LL_detrend_obs - bk_LL(pars, g_sparse, X=NA)
#' LL_detrend_obs_X - bk_LL(pars, g_sparse, X=X)
#'
#' # repeat with complete data
#'
#' # (don't pass the intercept column in X)
#' X = X_all[,-1]
#' LL_X_chol = bk_LL(pars, g_all, X=X)
#' LL_X_eigen = bk_LL(pars, g_all, fac_method='eigen', X=X)
#' z_obs = g_all$gval[!is.na(g_all$gval)]
#' #z_mat = matrix(z_obs, ncol=1)
#' V = bk_var(g_all, pars, sep=FALSE)
#' V_inv = chol2inv(chol(V))
#' X_tilde_inv = chol2inv(chol( crossprod(crossprod(V_inv, X_all), X_all) ))
#' betas_gls = X_tilde_inv %*% crossprod(X_all, (V_inv %*% z_obs))
#' z_gls = z_obs - (X_all %*% betas_gls)
#' z_gls_trans = crossprod(V_inv, z_gls)
#' quad_form = as.numeric( t(z_gls) %*% z_gls_trans )
#' log_det = as.numeric( determinant(V, logarithm=TRUE) )[1]
#' LL_direct = (-1/2) * ( n * log( 2 * pi ) + log_det + quad_form )
#' abs( LL_direct - LL_X_chol ) / max(LL_direct, LL_X_chol)
#' abs( LL_direct - LL_X_eigen ) / max(LL_direct, LL_X_eigen)
#'
#' # return components of likelihood with more=TRUE
#' LL_result = bk_LL(pars, g_all, X=X, more=TRUE)
#' LL_result$LL - LL_X_chol
#' LL_result$q - quad_form
#' LL_result$d - log_det
#' LL_result$n - n
#'
bk_LL = function(pars, g_obs, X=0, fac_method='chol', fac=NULL, quiet=TRUE, more=FALSE)
{
  # set flag for GLS and default one layer count
  if(is.data.frame(X)) X = as.matrix(X)
  use_GLS = is.matrix(X) | anyNA(X)

  # multi-layer support
  n_layer = 1
  is_multi = !is.null(g_obs[['idx_grid']])
  if(is_multi)
  {
    # coerce vector to 1-column matrix and identify non-NA points
    if( !is.matrix(g_obs[['gval']]) ) g_obs[['gval']] = matrix(g_obs[['gval']], ncol=1L)
    is_obs = !is.na(g_obs[['idx_grid']])

    # reorder the non-NA data matrix to grid-vectorized order
    reorder_z = g_obs[['idx_grid']][is_obs]
    n_layer = ncol(g_obs[['gval']])

    # matrix(.., ncol) avoids R simplifying to vector in 1 column case
    z = matrix(g_obs[['gval']][reorder_z,], ncol=n_layer)
    n_obs = nrow(z)

  } else {

    # single layer mode - copy non-NA data
    is_obs = as.vector(!is.na(g_obs[['gval']]))
    z = matrix(g_obs[['gval']][is_obs], ncol=1L)
    if( is.null(z) ) stop('No non-NA values found in g_obs')
    n_obs = length(z)
  }

  # complete data case triggers separability methods, which require eigen
  is_sep = all(is_obs)
  if(is_sep) fac_method = 'eigen'

  # compute factorization (scaled=TRUE means partial sill is factored out)
  if( is.null(fac) ) fac = bk_var(g_obs, pars, scaled=TRUE, fac_method=fac_method, sep=is_sep)

  # detect supplied factorization type
  fac_method = ifelse(is.matrix(fac), 'chol', 'eigen')

  # GLS estimate of mean based on predictors in X
  if( use_GLS ) X = bk_GLS(g_obs, pars, X=X, fac=fac, out='z')

  # matricize scalar and vector input to X
  if( !is.matrix(X) ) X = matrix(X, ncol=1L)
  if( nrow(X) == 1 ) X = matrix(rep(X, n_obs), ncol=ncol(X))

  # reorder X to match z
  if(is_multi) X = X[reorder_z,]

  # matrix of de-trended Gaussian random vectors to evaluate
  z_centered = matrix(z-X, ncol=n_layer)

  # Cholesky factor method is fastest
  if( fac_method == 'chol' )
  {
    # check for bad fac input (it should be matrix t(C), where C is the output of chol)
    if( !is.matrix(fac) ) stop('Cholesky factor (matrix) not found in fac')

    # get quadratic form of de-trended variable(s)
    quad_form = apply(z_centered, 2, function(z_i) {
      as.numeric(bk_var_mult(z_i, pars, fac=fac, quad=TRUE))
    })

    # determinant is the squared product of the diagonal elements in Cholesky factor t(C)
    log_det_corr = sum(log(diag(fac)))

    # partial sill was factored out with scaled=TRUE, so it is multiplied back in here
    log_det = n_obs * log(pars[['psill']]) + 2 * log_det_corr
  }

  # eigen-decomposition method
  if( fac_method == 'eigen' )
  {
    # check for bad fac input (it should be the list output of eigen)
    if( !is.list(fac) ) stop('eigen-decomposition (list) not found in fac')

    # get quadratic form of de-trended variable(s)
    quad_form = apply(z_centered, 2, function(z_i) {
      as.numeric(bk_var_mult(z_i, pars, fac=fac, quad=TRUE))
    })

    # determinant is product of the eigenvalues
    if( !is_sep )
    {
      # partial sill was factored out earlier so we multiply it back in (on log scale)
      if( !('values' %in% names(fac)) ) stop('non-separable eigen-decomposition not found in fac')
      log_det = n_obs * log(pars[['psill']]) + sum(log(fac[['values']]))

    } else {

      # separable case: eigenvalues are a kronecker product plus diagonal nugget effect
      if( !all(c('y', 'x') %in% names(fac)) ) stop('separable eigen-decomposition not found in fac')
      ev_scaled = kronecker(fac[['y']][['values']], fac[['x']][['values']])
      log_det = sum(log(pars[['eps']] + pars[['psill']] * ev_scaled))
    }
  }

  # compute log likelihood, print to console then return
  log_likelihood = (-1/2) * ( n_layer*( n_obs * log( 2 * pi ) + log_det ) + sum(quad_form) )
  if( !quiet ) cat( paste(round(log_likelihood, 5), '\n') )
  if(more) return(list(LL=log_likelihood, q=quad_form, d=log_det, n_obs=n_obs))
  return(log_likelihood)
}


#' Negative log-likelihood for parameter vector `p`
#'
#' Returns the negative log-likelihood of parameter vector `p` for the covariance
#' model `pars_fix`, given data grid `g_obs`.
#'
#' This is a wrapper for `-bk_LL()` allowing parameters to be passed as a numeric
#' vector instead of a list (for use in optimization etc). Parameters in `p` are copied
#' to `pars_fix` and passed to the likelihood computer.
#'
#' `p` is the vector of covariance parameters to test. Names in `p` are ignored; Its length
#' and order should correspond with the pattern of NAs in `pars_fix`. Users should check that
#' the desired parameter list is being constructed correctly by testing with:
#' `bk_pars_update(pars_fix, p, iso=iso, na_omit=TRUE)`.
#'
#' @param p numeric vector of covariance parameters accepted by `bk_pars_update`
#' @param g_obs list of form returned by `bk_grid` (with entries 'gdim', 'gres', 'gval')
#' @param pars_fix list of form returned by `bk_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param X numeric, vector, matrix, or NA, the mean or its linear predictors, passed to `bk_LL`
#' @param iso logical, indicates to use identical kernels for x and y (`pars$x` is ignored)
#' @param quiet logical indicating to suppress console output
#'
#' @return numeric, the negative log-likelihood of `p` given `g_obs`
#' @export
#'
#' @examples
#' # set up example grid and data
#' g_obs = bk_grid(10)
#' g_obs$gval = rnorm(10^2)
#'
#' # get some default parameters and vectorize them
#' pars = bk_pars(g_obs, 'gau')
#' p = bk_pars_update(pars)
#' bk_nLL(p, g_obs, pars)
#'
#' # change a parameter and re-evaluate
#' p_compare = p
#' p_compare[1] = 2*p_compare[1]
#' bk_nLL(p_compare, g_obs, pars)
#'
#' # repeat by calling bk_LL directly
#' pars_compare = pars
#' pars_compare$eps = 2*pars_compare$eps
#' -bk_LL(pars_compare, g_obs)
#'
#' # set up a subset of parameters for fitting
#' pars_fix = pars
#' pars_fix$eps = NA
#' pars_fix$y$kp = NA
#'
#' # names in p_fit are for illustration only (only the order matters)
#' p_fit = c(eps=1, y.rho=1)
#' bk_nLL(p_fit, g_obs, pars_fix)
#'
#' # equivalently:
#' pars_fit = pars
#' pars_fit$eps = p_fit[1]
#' pars_fit$y$kp = p_fit[2]
#' -bk_LL(pars_fit, g_obs)
#'
#' # check an input specification
#' bk_pars_update(pars_fix, p_fit, na_omit=TRUE)
#' pars_fit
#'
#'
bk_nLL = function(p, g_obs, pars_fix, X=0, iso=FALSE, quiet=TRUE, log_scale=FALSE)
{
  # transform back from log scale
  if(log_scale)
  {
    # parameter vector and parameter list
    p = exp(p)
    pars_fix = bk_pars_update(pars_fix, exp(bk_pars_update(pars_fix)) )
  }

  # update covariance parameter list with new values then vectorize it
  pars_complete = bk_pars_update(pars_fix, p, iso=iso, na_omit=TRUE)
  p_complete = bk_pars_update(pars_complete)

  # print complete parameters vector then compute likelihood
  if(!quiet) cat( paste(paste(format(p_complete), collapse=', '), ' :: LL = ') )
  return(-bk_LL(pars=pars_complete, g_obs=g_obs, X=X, quiet=quiet))
}


#' Random draw from multivariate normal distribution for grids
#'
#' Generates a random draw from the multivariate Gaussian distribution for the
#' covariance model `pars` on grid `g`, with mean zero.
#'
#' `pars` and `g` define the model's covariance matrix V. This function uses `base::rnorm`
#' to get a vector of independent standard normal variates, which it multiplies by the square
#' root of the covariance matrix, V, for the desired model (as defined by `pars` and `g`). The
#' result has a multivariate normal distribution with mean zero and covariance V.
#'
#' Multiple independent draws can be computed more efficiently by reusing the factorization
#' of V. This can be pre-computed with `bk_var` and supplied in `fac`, or a multi-layer
#' `g` can be supplied (see examples).
#'
#' @param g any object accepted or returned by `bk_grid`
#' @param pars list, covariance parameters in form returned by `bk_pars`
#' @param fac list, optional pre-computed factorization of component correlation matrices
#'
#' @return numeric vector, the vectorized grid data
#' @export
#'
#' @examples
#'
#' # example grid and covariance parameters
#' gdim = c(100, 200)
#' g = bk_grid(gdim)
#' pars_gau = bk_pars(g)
#'
#' # this example has a large nugget effect
#' gval = bk_sim(g, pars=pars_gau)
#' bk_plot(matrix(gval, gdim))
#'
#' # plot with yx coordinates
#' g_sim = modifyList(g, list(gval=gval))
#' bk_plot(g_sim)
#'
#' # repeat with smaller nugget effect for less noisy data
#' pars_smooth = modifyList(pars_gau, list(eps=1e-2))
#' gval_smooth = bk_sim(g, pars_smooth)
#' g_sim_smooth = modifyList(g, list(gval=gval_smooth))
#' bk_plot(g_sim_smooth)
#'
#' # the nugget effect can be very small, but users should avoid eps=0
#' pars_smoother = modifyList(pars_gau, list(eps=1e-12))
#' gval_smoother = bk_sim(g, pars_smoother)
#' g_sim_smoother = modifyList(g, list(gval=gval_smoother))
#' bk_plot(g_sim_smoother)
#'
#' # multi-layer example
#' n_pt = prod(gdim)
#' n_layer = 3
#' g_multi = bk_grid(list(gdim=gdim, gval=matrix(NA, n_pt, n_layer)))
#' gval_multi = bk_sim(g_multi, pars_smoother)
#' g_sim_multi = modifyList(g, list(gval=gval_multi))
#' bk_plot(g_sim_multi, layer=1)
#' bk_plot(g_sim_multi, layer=2)
#' bk_plot(g_sim_multi, layer=3)
#'
#'
bk_sim = function(g, pars=bk_pars(g), fac=NULL)
{
  # extract grid dimensions
  gdim = g[['gdim']]
  n = prod(gdim)
  n_layer = ifelse(is.matrix(g[['gval']]), ncol(g[['gval']]), 1L)

  # eigen-decompositions of separable components of full grid correlation matrix
  g_empty = modifyList(g, list(gval=NULL))
  if(is.null(fac)) fac = bk_var(g_empty, pars, fac_method='eigen', sep=TRUE)

  # report eigenvalue problems
  is_ev_negative = lapply(fac, function(eig) !( (pars$eps + eig[['values']]) > 0 ) )
  if( any(unlist(is_ev_negative)) ) stop('component correlation matrix has negative eigenvalue(s)')

  # multiply random iid normal vector(s) by covariance matrix square root
  seed_noise = matrix(rnorm(n*n_layer), ncol=n_layer)
  sim_gval = bk_var_mult(seed_noise, pars, fac=fac, fac_method='eigen', p=1/2)
  if(n_layer==1) { return(as.vector(sim_gval)) } else { return(sim_gval) }
}

