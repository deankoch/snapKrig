# sk_model.R
# Dean Koch, 2022
# Functions for modeling with snapKrig

#' Likelihood of covariance model `pars` given the data in sk grid `g`
#'
#' This computes the log-likelihood for the Gaussian process covariance model `pars`,
#' given 2-dimensional grid data `g`, and, optionally, linear trend data in `X`.
#'
#' The function evaluates:
#'
#' `-log( 2 * pi ) - ( 1/2 ) * ( log_det + quad_form )`,
#'
#' where `log_det` is the logarithm of the determinant of covariance matrix V, and
#' `quad_form` is z^T V^{-1} z, for the observed response vector z, which is constructed
#' by subtracting the trend specified in `X` (if any) from the non-NA values in `g`.
#'
#' If the trend spatially uniform and known, it can be supplied in argument `X` as a
#' numeric scalar. The default is zero-mean model `X=0`, which assumes users has
#' subtracted the trend from `g` beforehand.
#'
#' If the trend is unknown, the function will automatically use GLS to estimate it.
#' This is profile likelihood on the covariance function parameters (not REML). To
#' estimate a spatially constant mean, set `X=NA`. To estimate a spatially variable mean,
#' supply linear predictors as columns of a matrix argument to `X` (see `sk_GLS`). Users
#' can also pass a multi-layer bk grid `X` with covariates in layers.
#'
#' `fac_method` specifies how to factorize V; either by using the Cholesky factor ('chol')
#' or eigen-decomposition ('eigen'). A pre-computed factorization `fac` can be supplied by
#' first calling `sk_var(..., scaled=TRUE)` (in which case `fac_method` is ignored).
#'
#' When `out='a'`, the function instead returns the AIC value and when `out='b'` it returns
#' the BIC (see `stats::AIC`). This adjusts the likelihood for the number of covariance
#' and trend parameters (and in the case of BIC, the sample size), producing an index that
#' can be used for model comparison (lower is better).
#'
#' When `out='more'`, the function returns a list containing the log-likelihood and both
#' information criteria, along with several diagnostics: the number of observations, the
#' number of parameters, the log-determinant `log_det`, and the quadratic form `quad_form`.
#'
#' @param pars list of form returned by `sk_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param g an sk grid (or list with entries 'gdim', 'gres', 'gval' and/or 'idx_grid')
#' @param X numeric, vector, matrix, or NA, a fixed mean value, or matrix of linear predictors
#' @param fac_method character, the factorization to use: 'chol' (default) or 'eigen'
#' @param fac matrix or list, (optional) pre-computed covariance factorization
#' @param out character, either 'l' (likelihood), 'a' (AIC), 'b' (BIC), or 'more' (see details)
#' @param quiet logical indicating to suppress console output
#'
#' @return numeric, the likelihood of `pars` given `g_obs` and `X`, or list (if `more=TRUE`)
#' @export
#'
#' @family likelihood functions
#' @family variance-related functions
#' @seealso sk sk_GLS sk_var stats::AIC
#'
#' @examples
#' # set up example grid, covariance parameters
#' gdim = c(25, 12)
#' n = prod(gdim)
#' g_empty = g_lm = sk(gdim)
#' pars = modifyList(sk_pars(g_empty, 'gau'), list(psill=0.7, eps=5e-2))
#'
#' # generate some coefficients
#' n_betas = 3
#' betas = rnorm(n_betas)
#'
#' # generate covariates and complete data in grid and vector form
#' g_X = sk_sim(g_empty, pars, n_layer=n_betas-1L)
#' X = g_X[]
#' X_all = cbind(1, X)
#' g_lm = g_empty
#' g_lm[] = c(X_all %*% betas)
#'
#' # add some noise
#' g_all = sk_sim(g_empty, pars) + g_lm
#' z = g_all[]
#'
#' # two methods for likelihood
#' LL_chol = sk_LL(pars, g_all, fac_method='chol')
#' LL_eigen = sk_LL(pars, g_all, fac_method='eigen')
#'
#' # compare to working directly with matrix inverse
#' V = sk_var(g_all, pars, fac_method='none', sep=FALSE)
#' V_inv = chol2inv(chol(V))
#' quad_form = as.numeric( t(z) %*% crossprod(V_inv, z) )
#' log_det = as.numeric( determinant(V, logarithm=TRUE) )[1]
#' LL_direct = (-1/2) * ( n * log( 2 * pi ) + log_det + quad_form )
#'
#' # relative errors
#' abs( LL_direct - LL_chol ) / max(LL_direct, LL_chol)
#' abs( LL_direct - LL_eigen ) / max(LL_direct, LL_eigen)
#'
#' # get AIC or BIC directly
#' sk_LL(pars, g_all, out='a')
#' sk_LL(pars, g_all, out='b')
#'
#' # repeat with pre-computed variance factorization
#' fac_eigen = sk_var(g_all, pars, fac_method='eigen', sep=TRUE)
#' sk_LL(pars, g_all, fac=fac_eigen) - LL_eigen
#'
#' # repeat with multi-layer example
#' n_layer = 10
#' g_noise_multi = sk_sim(g_empty, pars, n_layer)
#' g_multi = g_lm + g_noise_multi
#' LL_chol = sk_LL(pars, g_multi, fac_method='chol')
#' LL_eigen = sk_LL(pars, g_multi, fac_method='eigen')
#' LL_direct = sum(sapply(seq(n_layer), function(j) {
#'  quad_form = as.numeric( t(g_multi[,j]) %*% crossprod(V_inv, g_multi[,j]) )
#'  (-1/2) * ( n * log( 2 * pi ) + log_det + quad_form )
#' }))
#'
#' # relative errors
#' abs( LL_direct - LL_chol ) / max(LL_direct, LL_chol)
#' abs( LL_direct - LL_eigen ) / max(LL_direct, LL_eigen)
#'
#' # repeat with most data missing
#' is_obs = seq(n) %in% sort(sample.int(n, 50))
#' n_obs = sum(is_obs)
#' g_obs = g_empty
#' z_obs = g_all[is_obs]
#' g_obs[is_obs] = z_obs
#'
#' # take subsets of covariates
#' g_X_obs = g_X
#' g_X_obs[!is_obs,] = NA
#' X_obs = X[is_obs,]
#'
#' LL_chol_obs = sk_LL(pars, g_obs, fac_method='chol')
#' LL_eigen_obs = sk_LL(pars, g_obs, fac_method='eigen')
#'
#' # working directly with matrix inverse
#' V_obs = sk_var(g_obs, pars, fac_method='none')
#' V_obs_inv = chol2inv(chol(V_obs))
#' quad_form_obs = as.numeric( t(z_obs) %*% crossprod(V_obs_inv, z_obs) )
#' log_det_obs = as.numeric( determinant(V_obs, logarithm=TRUE) )[1]
#' LL_direct_obs = (-1/2) * ( n_obs * log( 2 * pi ) + log_det_obs + quad_form_obs )
#' abs( LL_direct_obs - LL_chol_obs ) / max(LL_direct_obs, LL_chol_obs)
#' abs( LL_direct_obs - LL_eigen_obs ) / max(LL_direct_obs, LL_eigen_obs)
#'
#' # again using a pre-computed variance factorization
#' fac_chol_obs = sk_var(g_obs, pars, fac_method='chol', scaled=TRUE)
#' fac_eigen_obs = sk_var(g_obs, pars, fac_method='eigen', scaled=TRUE)
#' sk_LL(pars, g_obs, fac=fac_chol_obs) - LL_chol_obs
#' sk_LL(pars, g_obs, fac=fac_eigen_obs) - LL_eigen_obs
#'
#' # detrend the data by hand, with and without covariates then compute likelihood
#' g_obs_dtr = g_obs - sk_GLS(g_obs, pars)
#' g_obs_X_dtr = g_obs - sk_GLS(g_obs, pars, g_X)
#' LL_dtr = sk_LL(pars, g_obs_dtr, X=0)
#' LL_X_dtr = sk_LL(pars, g_obs_X_dtr, X=0)
#'
#' # or pass a covariates grid (or matrix) to de-trend automatically
#' LL_dtr - sk_LL(pars, g_obs, X=NA)
#' LL_X_dtr - sk_LL(pars, g_obs, X=g_X)
#'
#' # note that this introduce new unknown parameter(s), so AIC and BIC increase (worsen)
#' sk_LL(pars, g_obs, X=NA, out='a') > sk_LL(pars, g_obs_dtr, X=0, out='a')
#' sk_LL(pars, g_obs, X=g_X, out='a') > sk_LL(pars, g_obs_X_dtr, X=0, out='a')
#'
#' # X can be the observed subset, or the full grid (as sk grid or as matrix)
#' sk_LL(pars, g_obs, X=X)
#' sk_LL(pars, g_obs, X=X_obs)
#' sk_LL(pars, g_obs, X=g_X)
#' sk_LL(pars, g_obs, X=g_X_obs)
#'
#' # equivalent sparse input specification
#' g_sparse = g_all
#' g_sparse[] = matrix(g_obs[], ncol=1)
#' g_sparse = sk(gval=matrix(g_obs[], ncol=1), gdim=gdim)
#' LL_chol_obs - sk_LL(pars, g_sparse)
#' LL_eigen_obs - sk_LL(pars, g_sparse)
#' LL_dtr - sk_LL(pars, g_sparse, X=NA)
#' LL_X_dtr - sk_LL(pars, g_sparse, X=g_X)
#'
#' ## repeat with complete data
#'
#' # the easy way to get likelihood
#' LL_X_chol = sk_LL(pars, g_all, g_X)
#' LL_X_eigen = sk_LL(pars, g_all, g_X, fac_method='eigen')
#'
#' # the hard way
#' V = sk_var(g_all, pars, sep=FALSE)
#' V_inv = chol2inv(chol(V))
#' X_tilde_inv = chol2inv(chol( crossprod(crossprod(V_inv, X_all), X_all) ))
#' betas_gls = X_tilde_inv %*% crossprod(X_all, (V_inv %*% z))
#' z_gls = z - (X_all %*% betas_gls)
#' z_gls_trans = crossprod(V_inv, z_gls)
#' quad_form = as.numeric( t(z_gls) %*% z_gls_trans )
#' log_det = as.numeric( determinant(V, logarithm=TRUE) )[1]
#' LL_direct = (-1/2) * ( n * log( 2 * pi ) + log_det + quad_form )
#' abs( LL_direct - LL_X_chol ) / max(LL_direct, LL_X_chol)
#' abs( LL_direct - LL_X_eigen ) / max(LL_direct, LL_X_eigen)
#'
#' # return detailed list of components with out='more'
#' LL_result = sk_LL(pars, g_all, X=X, out='more')
#' LL_result[['LL']] - LL_X_chol
#' LL_result[['quad_form']] - quad_form
#' LL_result[['log_det']] - log_det
#' LL_result[['n_obs']] - n
#'
sk_LL = function(pars, g, X=0, fac_method='chol', fac=NULL, quiet=TRUE, out='l')
{
  # count layers
  n_layer = ifelse(is.matrix(g[['gval']]), ncol(g[['gval']]), 1L)

  # extract data from input list object g
  if(inherits(g, 'sk'))
  {
    # generics for sk objects to extract non-NA data
    is_obs = !is.na(g)
    z = matrix(g[is_obs], ncol=n_layer)

  } else {

    # handle non-sparse input
    if( is.null(g[['idx_grid']]) )
    {
      if( is.matrix(g[['gval']]) ) stop('gval was a matrix (expected a vector)')

      # identify observed data points and copy their index
      is_obs = if(is.null(g[['is_na']])) !is.na(g[['gval']]) else !g[['is_na']]

      # extract from vector to represent as 1-column matrix
      z = matrix(g[['gval']][is_obs], ncol=1L)

    } else {

      # copy matrix
      if( !is.matrix(g[['gval']]) ) stop('gval was a vector (expected a matrix)')
      is_obs = if(is.null(g[['is_na']])) !is.na(g[['idx_grid']]) else !g[['is_na']]
      z = g[['gval']]
    }
  }

  # complete data case triggers separability methods, which require eigen
  n_obs = sum(is_obs)
  is_sep = all(is_obs)
  if(is_sep) fac_method = 'eigen'

  # compute factorization (scaled=TRUE means partial sill is factored out)
  if( is.null(fac) ) fac = sk_var(g, pars, scaled=TRUE, fac_method=fac_method, sep=is_sep)

  # detect factorization method based on class of fac (in case user supplied it)
  if( is_sep )
  {
    # in separable case there are factorizations in a list(check only the first)
    fac_method = ifelse(is.matrix(fac[[1]]), 'chol', 'eigen')

  } else { fac_method = ifelse(is.matrix(fac), 'chol', 'eigen') }

  # # unpack X as needed
  # if( inherits(X, 'sk') )
  # {
  #   # extract only the required points
  #   if( is.matrix(X[['gval']]) ) { X = X[is_obs,] } else {X = matrix(X[is_obs], ncol=1L)}
  # }
  # TESTING

  # unpack X as matrix
  if(inherits(X, 'sk')) X = X[is_obs]

  # GLS estimate of mean based on predictors in X
  n_X = 0L
  if(is.data.frame(X)) X = as.matrix(X)
  use_GLS = is.matrix(X) | anyNA(X)
  if( use_GLS )
  {
    # count trend parameters (one for mean + one per covariate)
    n_X = 1L
    if( is.matrix(X) ) n_X = n_X + ncol(X)
    X = sk_GLS(g, pars, X=X, fac=fac, out='z')
  }

  # turn scalar and vector input into matrix X
  if( !is.matrix(X) ) X = matrix(X, ncol=1L)
  if( nrow(X) == 1 ) X = matrix(rep(X, n_obs), ncol=1L)

  # take subset of observed data as needed
  if( nrow(X) > n_obs ) X = X[is_obs]

  # matrix of de-trended Gaussian random vectors to evaluate
  z_centered = sweep(z, 1, X)

  # Cholesky factor method is fastest
  if( fac_method == 'chol' )
  {
    # check for bad fac input (it should be matrix t(C), where C is the output of chol)
    if( !is.matrix(fac) ) stop('Cholesky factor (matrix) not found in fac')

    # get quadratic form of de-trended variable(s)
    quad_form = apply(z_centered, 2, function(z_i) {
      as.numeric(sk_var_mult(z_i, pars, fac=fac, quad=TRUE))
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
      as.numeric(sk_var_mult(z_i, pars, fac=fac, quad=TRUE))
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

  # compute log likelihood and optionally print to console
  log_likelihood = (-1/2) * ( n_layer * ( n_obs * log( 2 * pi ) + log_det ) + sum(quad_form) )
  if( !quiet ) cat( paste(round(log_likelihood, 5), '\n') )

  # count number of parameters and compute AIC, BIC
  n_pars  = ( length(unlist(pars)) - 2 ) + n_X
  aic_out = ( -2 * log_likelihood ) + ( 2 * n_pars )
  bic_out = ( -2 * log_likelihood ) + ( log(n_obs * n_layer) * n_pars  )

  # return requested info
  out = tolower(out)
  if(startsWith(out, 'l')) return(log_likelihood)
  if(startsWith(out, 'a')) return(aic_out)
  if(startsWith(out, 'b')) return(bic_out)
  if(startsWith(out, 'm')) return(list(LL = log_likelihood,
                                       AIC = aic_out,
                                       BIC = bic_out,
                                       quad_form = quad_form,
                                       log_det = n_layer * log_det,
                                       n_obs = n_layer * n_obs,
                                       n_pars = n_pars,
                                       n_layer = n_layer))

  stop('argument out must start with one of the letters: l, a, b, or m')
}


#' Negative log-likelihood for parameter vector `p`
#'
#' Returns the negative log-likelihood of covariance model `pars_fix`, given the observations
#' in data grid `g_obs`. Parameter values are copied from the first argument, vector `p`, so
#' that the function can be passed to numerical optimizers (etc).
#'
#' This is a wrapper for `sk_LL` (times -1) that allows parameters to be passed as a numeric
#' vector instead of a list. Parameters in `p` are copied to `pars_fix` and passed to the
#' likelihood computer.
#'
#' `p` is the vector of covariance parameters to test. Names in `p` are ignored; Its length
#' and order should correspond with the pattern of NAs in `pars_fix`. Users should check that
#' the desired parameter list is being constructed correctly by testing with:
#' `sk_pars_update(pars_fix, p, iso=iso, na_omit=TRUE)`.
#'
#' @param p numeric vector of covariance parameters accepted by `sk_pars_update`
#' @param g_obs sk object or list with entries 'gdim', 'gres', 'gval'
#' @param pars_fix list of form returned by `sk_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param X numeric, vector, matrix, or NA, the mean or its linear predictors, passed to `sk_LL`
#' @param iso logical, indicates to use identical kernels for x and y (`pars$x` is ignored)
#' @param quiet logical indicating to suppress console output
#'
#' @return numeric, the negative log-likelihood of `p` given data `g_obs`
#' @export
#'
#' @family likelihood functions
#' @family variance-related functions
#' @seealso sk sk_GLS sk_var sk_pars_update
#'
#' @examples
#' # set up example grid and data
#' g = sk(gdim=10, gval=rnorm(10^2))
#'
#' # get some default parameters and vectorize them
#' pars = sk_pars(g, 'gau')
#' p = sk_pars_update(pars)
#' sk_nLL(p, g, pars)
#'
#' # change a parameter in the numeric vector and re-evaluate
#' p_compare = p
#' p_compare[1] = 2 * p_compare[1]
#' sk_nLL(p_compare, g, pars)
#'
#' # repeat by calling sk_LL directly with modified parameters list
#' pars_compare = pars
#' pars_compare[['eps']] = 2 * pars_compare[['eps']]
#' -sk_LL(pars_compare, g)
#'
#' # set up a subset of parameters to replace - eg when fitting those parameters
#' pars_fix = pars
#' pars_fix[['eps']] = NA
#' pars_fix[['y']][['kp']] = NA
#'
#' # names in p_fit are for illustration only (only the order matters)
#' p_fit = c(eps=1, y.rho=1)
#'
#' # replace NA parameter values in pars_fix to get completed parameters list
#' sk_pars_update(pars_fix, p_fit, na_omit=TRUE)
#'
#' # make the replacement and evaluate likelihood in one call
#' sk_nLL(p_fit, g, pars_fix)
#'
#' # equivalently:
#' pars_fit = pars
#' pars_fit[['eps']] = p_fit[1]
#' pars_fit[['y']][['kp']] = p_fit[2]
#' -sk_LL(pars_fit, g)
#'
sk_nLL = function(p, g_obs, pars_fix, X=0, iso=FALSE, quiet=TRUE, log_scale=FALSE)
{
  # note the grid object argument name is "g_obs" and not "g" (like in other functions)
  # because there seems to be a name conflict somewhere in the use of ... with
  # `stats::optim` that causes it to fail to copy arguments named "g"

  # transform back from log scale
  if(log_scale)
  {
    # parameter vector and parameter list
    p = exp(p)
    pars_fix = sk_pars_update(pars_fix, exp(sk_pars_update(pars_fix)) )
  }

  # update covariance parameter list with new values then vectorize it
  pars_complete = sk_pars_update(pars_fix, p, iso=iso, na_omit=TRUE)
  p_complete = sk_pars_update(pars_complete)

  # print complete parameters vector then compute likelihood
  if(!quiet) cat( paste(paste(format(p_complete), collapse=', '), ' :: LL = ') )
  return(-sk_LL(pars=pars_complete, g=g_obs, X=X, quiet=quiet))
}


#' Random draw from multivariate normal distribution for sk grids
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
#' of V. This can be pre-computed with `sk_var` and supplied in `fac`, or users can set
#' `n_layer` and the function will do this automatically.
#'
#' @param g an sk object or any grid object accepted by `sk`
#' @param pars list, covariance parameters in form returned by `sk_pars`
#' @param fac list, optional pre-computed factorization of component correlation matrices
#' @param n_layer positive integer, the number of draws to return
#'
#' @return sk grid or its vectorized form (vector for single-layer case, matrix for multi-layer case)
#' @export
#'
#' @family variance-related functions
#' @seealso sk sk_pars base::rnorm
#'
#' @examples
#'
#' # example grid and covariance parameters
#' gdim = c(100, 200)
#' g = sk(gdim)
#' pars_gau = sk_pars(g)
#'
#' # this example has a large nugget effect
#' g_sim = sk_sim(g, pars_gau)
#' plot(g_sim)
#'
#' # repeat with smaller nugget effect for less noisy data
#' pars_smooth = modifyList(pars_gau, list(eps=1e-2))
#' g_sim = sk_sim(g, pars_smooth)
#' plot(g_sim)
#'
#' # the nugget effect can be very small, but users should avoid eps=0
#' pars_smoother = modifyList(pars_gau, list(eps=1e-12))
#' g_sim = sk_sim(g, pars_smoother)
#' plot(g_sim)
#'
#' # multi-layer example
#' g_sim_multi = sk_sim(g, pars_smoother, n_layer=3)
#' plot(g_sim_multi, layer=1)
#' plot(g_sim_multi, layer=2)
#' plot(g_sim_multi, layer=3)
#'
sk_sim = function(g, pars=sk_pars(g), n_layer=1, fac=NULL, out='sk')
{
  # extract as sk object if not already a list
  if( !is.list(g) ) g = sk(g)

  # extract grid dimensions
  gdim = g[['gdim']]
  n = prod(gdim)

  # get eigen-decompositions of separable components of full grid correlation matrix
  if(is.null(fac)) fac = sk_var(g[c('gdim', 'gres')], pars, fac_method='eigen', sep=TRUE)

  # report any eigenvalue problems
  is_ev_negative = lapply(fac, function(eig) !( (pars[['eps']] + eig[['values']]) > 0 ) )
  if( any(unlist(is_ev_negative)) ) stop('component correlation matrix has negative eigenvalue(s)')

  # multiply random iid normal vector(s) by covariance matrix square root
  seed_noise = matrix(rnorm(n*n_layer), ncol=n_layer)
  sim_gval = sk_var_mult(seed_noise, pars, fac=fac, fac_method='eigen', p=1/2)

  # build the sk object to return
  if(n_layer==1) sim_gval = as.vector(sim_gval)
  if(out != 'sk') return(sim_gval)
  return( sk(modifyList(g, list(idx_grid=NULL, gval=sim_gval))) )
}

