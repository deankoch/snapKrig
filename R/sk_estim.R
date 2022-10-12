# sk_estim.R
# Dean Koch, 2022
# Functions for parameter inference and spatial prediction

#' Generalized least squares (GLS) with Kronecker covariances
#'
#' Computes coefficients b of the linear model E(Z) = Xb using the GLS equation
#' for sk grid `g` and covariance model `pars`. By default the function returns the
#' linear predictor as an sk object
#'
#' This is the maximum likelihood estimator for the linear trend Xb if we
#' assume the covariance parameters (in `pars`) are specified correctly.
#'
#' The GLS solution is: b = ( X^T V^{-1} X )^{-1} X^T V^{-1} z,
#'
#' where V is the covariance matrix for data vector z (which is `g[!is.na(g)]`), and X
#' is a matrix of covariates. V is generated from the covariance model `pars` with grid
#' layout `g`.
#'
#' Operations with V^{-1} are computed using the factorization `fac` (see `sk_var`), or
#' else as specified in `fac_method`.
#'
#' Argument `X` can be an sk grid (matching `g`) with covariates in layers; or it can be
#' a matrix of covariates. DO NOT include an intercept layer (all 1's) in argument `X` or
#' you will get collinearity errors. Matrix `X` should have independent columns, and its
#' rows should match the order of `g[]` or `g[!is.na(g)]`.
#'
#' Use `X=NA` to specify an intercept-only model; ie to fit a spatially constant mean. This
#' replaces X in the GLS equation by a vector of 1's.
#'
#' By default `out='s'` returns the linear predictor in an sk grid object. Change this to
#' `'z'` to return it as a vector, or `'b'` to get the GLS coefficients only. Set it to `'a'`
#' to get the second two return types (in a list) along with matrix `X` and its factorization.
#'
#' The length of the vector output for `out='z'` will match the number of rows in `X`.
#' This means that if `NA` grid points are excluded from `X`, they will not appear in
#' the output (and vice versa). In the `X=NA` case, the length is equal to the number of
#' non-`NA` points in `g`. Note that if a point is observed in `g`, the function will expect
#' its covariates to be included `X` (ie `X` should have no `NA`s corresponding to non-`NA`
#' points in `g`).
#'
#' If `g[]` is a matrix (a multi-layer grid), the covariates in `X` are recycled
#' for each layer. Layers are assumed mutually independent and the GLS equation is evaluated
#' using the corresponding block-diagonal V. This is equivalent to (but faster than) calling
#' `sk_GLS` separately on each layer with the same `X` and averaging the resulting b estimates.
#'
#'
#' @param g a sk grid object (or list with entries 'gdim', 'gres', 'gval')
#' @param pars list of form returned by `sk_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param X sk grid, matrix or NA, the linear predictors (in columns) excluding intercept
#' @param out character, either 'b' (coefficients), 'z' or 's' (Xb), or 'a' (all)
#' @param fac_method character, factorization method: 'eigen' (default) or 'chol' (see `sk_var`)
#' @param fac matrix or list, (optional) pre-computed covariance matrix factorization
#'
#' @return linear predictor Xb as an sk grid, or numeric vector, or coefficients (see details)
#'
#' @export
#' @seealso sk
#' @family estimators
#' @family variance-related functions
#'
#' @examples
#' # set up example grid and covariance parameters
#' gdim = c(45, 31)
#' g_empty = sk(gdim)
#' n = length(g_empty)
#' pars = modifyList(sk_pars(g_empty, 'gau'), list(psill=2))
#'
#' # generate spatial noise
#' g_noise = sk_sim(g_empty, pars)
#' plot(g_noise)
#'
#' # generate more spatial noise to use as covariates
#' n_betas = 3
#' betas = rnorm(n_betas, 0, 10)
#' g_X = sk_sim(g_empty, pars, n_layer=n_betas-1L)
#' X = g_X[]
#' X_all = cbind(1, X)
#' g_lm = g_empty
#' g_lm[] = as.vector(X_all %*% betas)
#' plot(g_lm)
#'
#' # combine with noise to make "observed" data
#' g_obs = g_lm + g_noise
#' plot(g_obs)
#'
#' # By default (out='s') the function returns the linear predictor
#' g_lm_est = sk_GLS(g_obs, pars, g_X, out='s')
#' g_lm_est
#' plot(g_lm_est)
#'
#' # equivalent, but slightly faster to get vector output
#' max(abs( sk_GLS(g_obs, pars, g_X, out='z') - g_lm_est[] ))
#'
#' # repeat with matrix X
#' max(abs( sk_GLS(g_obs, pars, g_X[], out='z') - g_lm_est[] ))
#'
#' # return the GLS coefficients
#' betas_est = sk_GLS(g_obs, pars, g_X, out='b')
#' print(betas_est)
#' print(betas)
#'
#' # compute trend manually as product of betas with X and intercept
#' lm_est = X_all %*% betas_est
#' max( abs(lm_est - g_lm_est[] ) )
#'
#' # de-trend observations by subtracting linear predictor
#' plot(g_obs - g_lm_est)
#'
#' # repeat with pre-computed eigen factorization (same result but faster)
#' fac_eigen = sk_var(g_obs, pars, fac_method='eigen', sep=TRUE)
#' betas_est_compare = sk_GLS(g_obs, pars, g_X, fac=fac_eigen, out='b')
#' max( abs( betas_est_compare - betas_est ) )
#'
#' # missing data example
#' n_obs = 10
#' g_miss = g_obs
#' idx_miss = sort(sample.int(n, n-n_obs))
#' g_miss[idx_miss] = NA
#' is_obs = !is.na(g_miss)
#' plot(g_miss)
#'
#' # coefficient estimates are still unbiased but less precise
#' betas_est = sk_GLS(g_miss, pars, g_X, out='b')
#' print(betas_est)
#' print(betas)
#'
#' # set X to NA to estimate the spatially constant trend
#' b0 = sk_GLS(g_miss, pars, X=NA, out='b')
#'
#' # matrix X does not need to include unobserved points, but output is filled to match X
#' X_obs = X[is_obs,]
#' sk_GLS(g_miss, pars, X=X_obs)
#' sk_GLS(g_miss, pars, X=X)
#'
#' # generate some extra noise for 10-layer example
#' g_noise_multi = sk_sim(g_empty, pars, n_layer=10)
#' g_multi = g_lm + g_noise_multi
#' betas_complete = sk_GLS(g_multi, pars, g_X, out='b')
#' print(betas_complete)
#' print(betas)
#'
#' # multi-layer input shares covariates matrix X, and output is to a single layer
#' summary(sk_GLS(g_multi, pars, g_X))
#' summary(sk_GLS(g_multi, pars, X))
#'
#' # note that X cannot be missing data where `g` is observed
#' \dontrun{
#' summary(sk_GLS(g_multi, pars, X_obs))
#' }
#'
#' # repeat with missing data
#' g_multi[!is_obs,] = NA
#' g_X_obs = g_X
#' g_X_obs[!is_obs,] = NA
#' betas_sparse = sk_GLS(g_multi, pars, X, out='b')
#' print(betas_sparse)
#' print(betas)
#' summary(sk_GLS(g_multi, pars, g_X))
#' summary(sk_GLS(g_multi, pars, X))
#' summary(sk_GLS(g_multi, pars, g_X_obs))
#' summary(sk_GLS(g_multi, pars, X_obs))
#'
sk_GLS = function(g, pars, X=NA, out='s', fac_method='eigen', fac=NULL)
{
  # multi-layer support
  is_multi = !is.null(g[['idx_grid']])
  if(is_multi)
  {
    # identify non-NA points and extract them as matrix
    is_obs = !is.na(g[['idx_grid']])
    n_layer = ncol(g[['gval']])
    z = g[['gval']]

  } else {

    # copy non-NA data as a 1-column matrix
    is_obs = as.vector(!is.na(g[['gval']]))
    z = matrix(g[['gval']][is_obs], ncol=1L)
  }

  # set default factorization method
  if( sum(is_obs) < 2 ) stop('Not enough non-NA values in g')
  is_sep = all(is_obs)
  if( is_sep & (fac_method=='chol') ) stop('eigen method is required for complete grids')

  # compute variance factorization (scaled=TRUE -> partial sill is factored out)
  if( is.null(fac) ) fac = sk_var(g, pars, scaled=TRUE, fac_method=fac_method, sep=is_sep)

  # build matrix of covariate values from intercept column and (optionally) X
  n = length(is_obs)
  n_obs = sum(is_obs)
  if( anyNA(X) & !inherits(X, 'sk') )
  {
    # if X is a matrix with NAs, discard it with a warning
    if( is.matrix(X) ) warning('X contained NA value(s). Setting X=NA')
    X = X_obs = matrix(1L, nrow=n_obs)
    is_X_full = FALSE

  } else {

    # unpack sk grid X
    if(inherits(X, 'sk')) X = X[!is.na(X)]

    # check for invalid input to X
    if( !is.matrix(X) ) stop('X must be a matrix of covariate values')

    # check for incorrect length in X
    msg_expected = ifelse(n==n_obs, n, paste(n, 'or', n_obs))
    msg_got = paste('expected', msg_expected, 'but got', nrow(X))
    if( !( nrow(X) %in% c(n, n_obs) ) ) stop(paste('incorrect number of rows in X:', msg_got))

    # append intercept column and take observed subset of X if the full matrix was supplied
    is_X_full = nrow(X) > n_obs
    X = X_obs = cbind(1L, X)
    if(is_X_full) X_obs = matrix(X[is_obs,], ncol=ncol(X))
  }

  # find the factorization of quadratic form with X (scaling by V inverse)
  fac_X = sk_var(g, pars, X=X_obs, scaled=TRUE, fac=fac, fac_method='eigen')

  # compute GLS coefficients using whitened observation data
  z_trans = sk_var_mult(z, pars, fac=fac)
  betas_gls = sk_var_mult(t(crossprod(z_trans, X_obs)), pars, fac=fac_X)
  if(is_multi) betas_gls = rowMeans(betas_gls)

  # for partial matching out argument
  out = tolower(out)

  # return betas by default
  if(startsWith(out, 'b')) return(as.numeric(betas_gls))

  # or return the linear predictor as vector
  z_gls = as.numeric( tcrossprod(X, t(betas_gls)) )
  if(startsWith(out, 'z')) return(z_gls)

  # copy linear predictor to grid list
  if(is_X_full) { g[['gval']] = z_gls } else {

    # for multi-layer input, the output is a single layer (shared)
    if(is_multi) g[['gval']] = g[['idx_grid']]

    # only the observed points get new data assigned
    g[['gval']][is_obs] = z_gls
  }

  # return as single-layer sk grid object
  g[['idx_grid']] = NULL
  if(startsWith(out, 's')) return(sk(g))

  # or return a list with everything
  if(startsWith(out, 'a')) return(list(s = sk(g),
                                       z = z_gls,
                                       b = as.numeric(betas_gls),
                                       x = X_obs,
                                       fac_X = fac_X))
}


#' Compute ordinary kriging predictor (or variance) for data on a grid
#'
#' Evaluates the ordinary kriging equations in section 3 of Cressie (1993) over the
#' grid defined in `g_obs`. These are the predicted values minimizing mean squared
#' prediction error under the covariance model specified by `pars`.
#'
#' Set `makev=TRUE` to return the pointwise kriging variance. This takes approximately
#' n_obs times longer to evaluate than `makev=FALSE`. A progress bar will be printed to
#' console unless `quiet=TRUE`.
#'
#' The covariance factorization `fac` can be pre-computed using `sk_var(..., scaled=TRUE)`
#' to speed up repeated calls where only the observed data values change (ie same covariance
#' structure `pars`, and same NA structure in the data). Note that the kriging variance does
#' not change in this case and only needs to be computed once.
#'
#' @param g_obs list of form returned by `sk` (with entries 'gdim', 'gres', 'gval')
#' @param pars list of form returned by `sk_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param X numeric, vector, matrix, or NA: the mean, or its linear predictors
#' @param out character, the return value, one of 'predictor', 'variance', or 'm'
#' @param fac (optional) pre-computed factorization of covariance matrix scaled by partial sill
#' @param quiet logical indicating to suppress console output
#'
#' @return numeric matrix, the predicted values (or their variance)
#' @export
#'
#' @examples
#' # make example grid and covariance parameters
#' g = sk_sim(100)
#' pars = sk_pars(g)
#' g_pred = sk_cmean(g, pars)
#'
#'
#' g_var = sk_cmean(g_obs, pars, makev=TRUE, quiet=TRUE)
#' #g_obs |> sk_plot()
#' #g_obs |> modifyList(list(gval=g_pred)) |> sk_plot()
#' #g_obs |> modifyList(list(gval=g_var)) |> sk_plot()
#'
sk_cmean = function(g_obs, pars, X=NA, fac=NULL, out='p', fac_method='chol', quiet=FALSE)
{
  # check for expected objects in list in g_obs
  nm_expect = c('gdim', 'gres', 'gval')
  msg_expect = paste('g_obs must be a list with entries named', paste(nm_expect, collapse=', '))
  if( !is.list(g_obs) | !all( nm_expect %in% names(g_obs) ) ) stop(msg_expect)
  gdim = g_obs[['gdim']]
  gres = g_obs[['gres']]

  # identify observed data points and copy their index
  is_obs = is_obs_src = !is.na(g_obs[['gval']])
  idx_obs = which( is_obs )

  # copy non-NA data
  z = g_obs[['gval']][is_obs]
  if( is.null(z) ) stop('No non-NA values found in g_obs')
  n_obs = length(z)
  n = prod(gdim)

  # initialize parameters and check if we are estimating linear predictor for the model
  use_GLS = anyNA(X) | is.matrix(X)
  if( !is.matrix(X) ) mu_GLS = X
  m = rep(0, n)

  # construct first rows of the symmetric Toeplitz correlation matrices for y and x
  cy = cy_obs = sqrt(pars[['psill']]) * sk_corr_mat(pars[['y']], gdim[['y']], gres[['y']], i=1)
  cx = cx_obs = sqrt(pars[['psill']]) * sk_corr_mat(pars[['x']], gdim[['x']], gres[['x']], i=1)

  # set up factorization when it's not supplied. Variance mode forces eigen-decomposition
  if( startsWith(out, 'v') ) fac_method = 'eigen'

  # check for completeness and separability
  is_sep = n == n_obs

  # check for sub-grid structure and substitute simpler equivalent problem if possible
  sg = sk_sub_find(is_obs, g_obs[['gdim']])
  is_sg = !is.null(sg)
  if( is_sg )
  {
    # extract sub-grid layout and find separable covariance eigen-decomposition
    fac_method = 'eigen'
    g_obs = list(gval=z, gdim=sg[['gdim']], gres=g_obs[['gres']] * sg[['res_scale']])
    if( is.null(fac) ) fac = sk_var(g_obs, pars, scaled=TRUE, fac_method=fac_method, sep=TRUE)
    cy_obs = cy[ sg[['ij']][['y']] ]
    cx_obs = cx[ sg[['ij']][['x']] ]

    # g_obs should have no missing (NA) data points now
    is_obs = rep(TRUE, n_obs)
    # (is_obs_src still contains original observed index)

  } else {

    # compute factorization (scaled=TRUE means partial sill is factored out)
    if( is.null(fac) ) fac = sk_var(g_obs, pars, scaled=TRUE, fac_method=fac_method, sep=is_sep)
  }

  # transform the observed data by left-multiplying with inverse covariance
  z_tilde = sk_var_mult(g_obs[['gval']][is_obs], pars, fac=fac)

  # left-multiply by cross-covariance to get simple kriging predictor
  z_p_tilde = sk_toep_mult(cy, z_tilde, cx, idx_obs, gdim) |> as.numeric()
  if( !use_GLS & startsWith(out, 'p') ) return(as.numeric(X) + z_p_tilde)

  # compute GLS coefficients and resulting adjustments to predictor as needed
  if(use_GLS)
  {
    # make a copy of the observed locations in X
    if( !anyNA(X) ) { X_obs = matrix(X[is_obs_src,], n_obs) } else {
      X = NULL
      X_obs = NA
    }

    # find betas, predictor, and the data matrix with an intercept column
    gls = sk_GLS(g_obs, pars, X=X_obs, fac=fac, fac_method=fac_method, out='a')
    fac_X = gls[['fac_X']]

    # compute bias adjustment due to estimation of linear predictor
    X_p_tilde = gls[['x']] |>
      sk_var_mult(pars, fac=fac) |>
      apply(2, \(x) sk_toep_mult(cy, x, cx, idx_obs, gdim))

    # uncomment to get exact interpolator (and discontinuities at observations)
    #X_p_tilde[idx_obs,] = gls[['x']]

    # compute trend and 'm'
    X_adj = cbind(rep(1, n), X) - X_p_tilde
    mu_GLS = tcrossprod(X_adj, t(gls[['b']])) |> as.numeric()
    if( startsWith(out, 'm') ) return( as.numeric(sk_var_mult( t(X_adj), pars, fac=fac_X)) )
  }

  # universal kriging predictor
  if( startsWith(out, 'p') ) return( mu_GLS + z_p_tilde )

  # compute variance contribution from GLS (0 when not using GLS)
  v_gls = numeric(n)
  if(use_GLS)
  {
    # small loop over eigen-values in fac_X, adding each contribution to running total in v_gls
    for(idx in seq_along(fac_X[['values']]))
    {
      # eigen-values of inverse scaled by psill
      ev = 1 / ( pars[['psill']] * fac_X[['values']][idx] )
      v_gls[] = v_gls[] + ev * tcrossprod(fac_X[['vectors']][,idx], X_adj)^2
    }
  }

  # use a more efficient method when observed points form a sub-grid
  idx_ev = seq(n_obs)
  if(is_sg)
  {
    # check that the correct factorization (componentwise, in a list) was supplied
    if( !all( c('y', 'x') %in% names(fac) ) ) stop('supplied factorization was not separable')

    # compute eigenvalues for observed covariance matrix inverse
    ev_corr = kronecker(fac[['x']][['values']], fac[['y']][['values']])
    ev = 1 / ( (pars[['psill']] * ev_corr) + pars[['eps']] )

    # directly build and multiply the relatively small component covariance matrices
    c_cross_y = sk_corr_mat(pars[['y']], gdim[['y']], gres[['y']], j=sg[['ij']][['y']])
    c_cross_x = sk_corr_mat(pars[['x']], gdim[['x']], gres[['x']], j=sg[['ij']][['x']])
    add_y2 = ( c_cross_y %*% fac[['y']][['vectors']] )
    add_x2 = ( c_cross_x %*% fac[['x']][['vectors']] )

    # a different ordering for the loop below (largest eigenvalues first)
    idx_ev = order(ev)
  }

  # large loop over eigen-values of covariance matrix, iteratively adding to v_rem
  v_rem = numeric(n)
  if(!quiet) pb = utils::txtProgressBar(max=n_obs, style=3)

  for(idx in seq(n_obs))
  {
    # update progress bar then change to reordered index
    if(!quiet) utils::setTxtProgressBar(pb, idx)
    idx = idx_ev[idx]

    # non-separable case first
    if( !is_sg )
    {
      # eigen-values of inverse are scaled by psill, then slow multiplication with cross covariance
      ev = 1 / ( pars[['psill']] * fac[['values']][idx] )
      v_add = ev * sk_toep_mult(cy, fac[['vectors']][, idx], cx, idx_obs, gdim)^2

    } else {

      # column indices in component correlation matrices corresponding to eigen-value idx
      idx_yx = sk_vec2mat(idx, sg[['gdim']], out='list')

      # kronecker product of component columns to get large column vector
      add_yx = ( pars[['psill']] * kronecker(add_x2[, idx_yx[['j']]], add_y2[, idx_yx[['i']]]) )^2
      v_add = ev[idx] * add_yx
    }

    # add to total
    v_rem = v_rem + as.vector(v_add)
  }
  if(!quiet) close(pb)
  return(pars[['psill']] + pars[['eps']] + as.numeric(v_gls) - v_rem)
}




#' Fit a covariance model to the data by maximum likelihood
#'
#' An automated model fitting procedure
#'
#' documentation unfinished
#'
#' @param g_obs todo
#' @param pars todo
#' @param X todo
#'
#' @return sdfsdfdsfs
#' @export
#'
#' @examples
#'
#' # define a grid
#' gdim = c(50, 53)
#' g_empty = sk(gdim)
#' pars_src = sk_pars(g_empty)
#' pars_src = modifyList(pars_src, list(eps=runif(1, 0, 1e1), psill=runif(1, 0, 1e2)))
#' pars_src[['y']][['kp']] = pars_src[['x']][['kp']] = runif(1, 1, 50)
#'
#' # generate example data and fit to it
#' g_obs = sk_sim(g_empty, pars_src)
#' sk_plot(g_obs)
#' fit_result = sk_fit(g_obs, pars='gau')
#'
#' fit_result$pars |> sk_pars_update()
#' pars_src |> sk_pars_update()
#'
#' # check sequence of other psill values
#' pars_out = fit_result$pars
#' psill_test = ( 2^(seq(5) - 3) ) * pars_out[['psill']]
#' LL_test = sapply(psill_test, function(s) sk_LL(modifyList(pars_out, list(psill=s)), g_obs) )
#' plot(psill_test, LL_test)
#' lines(psill_test, LL_test)
#' print(data.frame(psill=psill_test, likelihood=LL_test))
#'
#' # repeat with most data missing
#' n = prod(gdim)
#' n_obs = 200
#' gval = sk_sim(g_empty, pars_src, quiet=TRUE)
#' g_obs = modifyList(g_empty, list(gval=gval))
#' idx_obs = sample.int(prod(gdim), n_obs)
#' g_miss = modifyList(g_obs, list(gval=rep(NA, n)))
#' g_miss$gval[idx_obs] = g_obs$gval[idx_obs]
#' sk_plot(g_miss)
#'
#'
sk_fit = function(g_obs, pars=NULL, X=NA, iso=TRUE, initial=NULL, quiet=FALSE,
                     lower=NULL, upper=NULL, n_max=1e3)
{
  # unpack vectorized grid as list
  g_obs = sk(g_obs)
  gdim = g_obs[['gdim']]

  # check for missingness and count observations
  is_obs = !is.na(g_obs[['gval']])
  if( !any(is_obs) ) stop('no data found in g_obs')
  n_all = prod(gdim)
  n_obs = sum(is_obs)

  # problem size sanity check
  if( (n_obs < n_all) & (n_obs > n_max) ) stop('number of observed points exceeded n_max')

  # # handle sparse indexing
  # if( is.null(g[['idx_grid']]) )
  # {
  #   g[['idx_grid']]
  #
  #
  # }

  # coerce to matrix of appropriate dimensions
  # n_obs = sum(is_obs)
  # if( is.matrix(X) )
  # {
  #   # pad with NAs as needed
  #   #idx_pad = match(seq(n_all), which(is_obs))
  #   # if( is.vector(X) ) if( length(X) == n_obs ) X = X[idx_pad]
  #   # if( is.matrix(X) ) if( nrow(X) == n_obs ) X = apply(X, 2, function(x) x[idx_pad])
  #   # X = matrix(X, nrow=n_all)
  #
  #   if( is.vector(X) ) if( length(X) == n_obs ) X = X[idx_pad]
  #   X = matrix(X, nrow=n_obs)
  #
  #   if( all(is.na(X)) ) X = NA
  # }

  # # set flags for multi-layer support
  # is_indexed = !is.null(g_obs[['idx_grid']])
  # is_multi = is.matrix(g_obs[['gval']]) | is_indexed
  # if(is_multi) { g_obs[['gval']] = as.matrix(g_obs[['gval']]) } else {
  #
  #   # substitute equivalent sub-grid problem if possible
  #   sub_result = sk_sub_find(g_obs=g)
  #   if( !is.null(sub_result) )
  #   {
  #     # skip when full grid is non-NA
  #     if( any(sub_result[['gdim']] < gdim) )
  #     {
  #       # compute new grid configuration and copy to g
  #       gdim = sub_result[['gdim']]
  #       gres = g[['gres']] * sub_result[['res_scale']]
  #       gyx = Map(function(yx, idx) yx[idx], yx=g[['gyx']], idx=sub_result[['ij']])
  #
  #
  #       # TODO: test and fix this for sparse and matrix-valued cases
  #       g = modifyList(g, list(gdim=gdim, gres=gres, gyx=gyx, gval=g[['gval']][is_obs]))
  #       is_obs = rep(TRUE, prod(gdim))
  #     }
  #   }
  # }

  # set covariance parameter defaults
  if( !is.list(pars) ) pars = sk_pars(g_obs, ifelse(is.null(pars), 'gau', pars))
  p_fixed = sk_pars_update(pars, iso=iso)
  is_fitted = is.na(p_fixed)
  if( !any(is_fitted) ) is_fitted[] = TRUE

  # set initial value defaults
  nm_fitted = names(is_fitted)[is_fitted]
  nm_fixed = names(is_fitted)[!is_fitted]
  if( is.null(initial) ) initial = sk_bds(pars, g_obs)[nm_fitted, 'initial']
  pars = sk_pars_update(pars, p_fixed, iso=iso)

  # fit the model
  #v = var(g[['gval']], na.rm=TRUE)
  result_optim = sk_optim(g_obs, pars, X=X, iso=iso, quiet=quiet,
                             lower=lower, initial=initial, upper=upper)
  pars_out = result_optim[['pars']]

  # # check sequence of likely psill substitutions
  # psill_test = ( 2^(seq(5) - 3) ) * pars_out[['psill']]
  # pars_test = lapply(psill_test, function(s) modifyList(pars_out, list(psill=s)) )
  # LL_test = sapply(pars_test, function(p) sk_LL(p, g, X=X) )
  # pars_out[['psill']] = psill_test[ which.max(LL_test) ]

  # de-trend the data by subtracting GLS estimate
  #z_gls = 0
  #if( anyNA(X) | is.matrix(X) ) z_gls = sk_GLS(g, pars_out, X=X, out='z')
  #g[['gval']][is_obs] = g[['gval']][is_obs] - z_gls

  # plot the semi-variogram for de-trended data
  #vg_detrend = sk_sample_vg(g)
  #sk_plot_semi(vg_detrend, pars_out)
  return(result_optim)


  if(0)
  {

  # scale the data
  z_std = scale(g[['gval']])
  z_centre = attr(z_std, 'scaled:center')
  z_scale = attr(z_std, 'scaled:scale')
  g = modifyList(g, list(gval=as.vector(z_std)))

  # variance parameters must also be scaled
  is_v_fitted = nm_fitted %in% c('eps', 'psill')
  is_v_fixed = nm_fixed %in% c('eps', 'psill')
  initial[is_v_fitted] = initial[is_v_fitted] / (2*z_scale^2)
  if( any(is_v_fixed) ) p_fixed[nm_fixed][is_v_fixed] = p_fixed[nm_fixed][is_v_fixed] / (2*z_scale^2)
  pars = sk_pars_update(pars, p_fixed, iso=TRUE)

  # fit the model
  if( anyNA(X) ) X = NA
  result_optim = sk_optim(g, pars, X=X, iso=TRUE, initial=initial)
  pars_out = result_optim[['pars']]

  # check sequence of likely psill values
  psill_test = ( 2^(seq(5) - 3) ) * pars_out[['psill']]
  pars_test = lapply(psill_test, function(s) modifyList(pars_out, list(psill=s)) )
  LL_test = sapply(pars_test, function(p) sk_LL(p, g, X=X) )
  #print(z_scale)
  #print(data.frame(psill=psill_test, likelihood=LL_test))
  pars_out[['psill']] = psill_test[ which.max(LL_test) ]

  # de-trend the scaled data
  z_gls = sk_GLS(g, pars_out, X=X, out='z')
  if(anyNA(X)) z_gls = z_gls[1]
  g_out = modifyList(g, list(gval = z_scale * (z_std-z_gls) ) )

  # transform scaled variance parameters back to scale of input data
  v_unscaled = list(eps = 2*z_scale^2 * pars_out[['eps']], psill = 2*z_scale^2 * pars_out[['psill']])
  pars_out = modifyList(pars_out, v_unscaled)

  # plot the semi-variogram on original scale
  vg_out = sk_sample_vg(g_out)
  sk_plot_semi(vg_out, pars_out)
  return(modifyList(result_optim, list(pars=pars_out)))
  }
}



#' Fit covariance parameters to data by maximum (profile) likelihood using optim
#'
#' documentation unfinished
#'
#' @param g_obs list of form returned by `sk` (with entries 'gdim', 'gres', 'gval')
#' @param pars list of fixed kernel parameters, with NAs indicating parameters to fit
#' @param X numeric, vector, matrix, or NA, the mean or its linear predictors, passed to `sk_LL`
#' @param iso logical, indicating to constrain the y and x kernel parameters to be the same
#' @param control list, passed to `stats::optim`
#' @param quiet logical, indicating to suppress console output
#' @param ... named arguments to pass to `sk_bds`
#'
#' @return sdfsdfsd
#' @export
#'
#' @examples
#' # set up example grid and data
#' g_obs = sk(10)
#' g_obs$gval = rnorm(10^2)
#' sk_optim(g_obs, quiet=TRUE)
#'
#' # repeat with one or more parameters fixed
#' pars = sk_pars_make('gau') # NA parameters list
#' pars$psill = 1
#' sk_optim(g_obs, pars, quiet=TRUE)
#' pars$y$kp = 1
#' sk_optim(g_obs, pars, quiet=TRUE)
#'
#' # iso mode constrains x parameters to equal y parameters
#' sk_optim(g_obs, iso=T, quiet=TRUE)
#' sk_optim(g_obs, pars, iso=T, quiet=TRUE)
#'
sk_optim = function(g_obs, pars='gau', X=0, iso=FALSE, control=list(), quiet=FALSE,
                       log_scale=TRUE, method='L-BFGS-B', lower=NULL, initial=NULL, upper=NULL)
{
  # only L-BFGS-B accepts bounds, so log-scale is mandatory for other methods
  if( (method != 'L-BFGS-B' & !log_scale) )
  {
    warning('Setting log_scale=TRUE (required for methods other than L-BFGS-B)')
    log_scale = TRUE
  }

  g_obs = sk(g_obs)

  # standardize input to pars and set NAs for missing values
  pars_fix = sk_pars_make(pars)

  # extract parameter names and NA structure supplied in pars_fix
  pars_fix_vec = sk_pars_update(pars_fix, iso=iso)
  nm_fit = names(pars_fix_vec)[ is.na(pars_fix_vec) ]

  # TODO: when values are specified for all parameters in `pars`, use them as initial values
  if( length(nm_fit) == 0 )
  {
    # set initial value automatically here?
    #initial = sk_pars_update(pars_fix)
    pars_fix = sk_pars_update(pars_fix, rep(NA, length(pars_fix_vec)), iso=iso)
    pars_fix_vec = sk_pars_update(pars_fix, iso=iso)
    nm_fit = names(pars_fix_vec)
  }

  # get default initial values and bounds data frame
  pars_df = sk_bds(pars_fix, g_obs)[nm_fit,]

  # then overwrite with any user supplied settings
  if( !is.null(lower) ) pars_df[['lower']] = lower
  if( !is.null(initial) ) pars_df[['initial']] = initial
  if( !is.null(upper) ) pars_df[['upper']] = upper

  # switch to log-scale if requested
  if(log_scale)
  {
    #bds_all_df = bds_all_df |> log()
    pars_df = log(pars_df)
    pars_fix_vec = log(pars_fix_vec)
  }

  # sacling constants passed to optimizer for internal use
  optim_scale = apply(pars_df, 1L, function(p) diff(range(p)))

  # TODO: check this
  # when iso=TRUE and pars_fix contains fixed y kernel parameter(s) they must be copied to x
  pars_fix = sk_pars_update(pars_fix, pars_fix_vec, iso=iso)

  # evaluate objective at initial values and bounds
  bds_nLL = apply(pars_df, 2L, function(p) {
    sk_nLL(p, g_obs, pars_fix, X, iso, quiet=TRUE, log_scale) })

  # 1d optimization
  if( nrow(pars_df) == 1 )
  {
    tol = control[['tol']]
    optimize_out = optimize(f = sk_nLL,
                            interval = pars_df[c('lower', 'upper')],
                            tol = ifelse(is.null(tol), .Machine$double.eps^0.25, tol),
                            g_obs = g_obs,
                            pars_fix = pars_fix,
                            X = X,
                            iso = iso,
                            log_scale = log_scale,
                            quiet = quiet)

    optim_out = list(message='', par=optimize_out[['minimum']], value=optimize_out[['objective']])

  } else {

    # n-d optimization
    # set default control parameters and run the optimizer
    control = modifyList(list(maxit=1e3, parscale=optim_scale), control)
    optim_out = stats::optim(par = pars_df[['initial']],
                             lower = pars_df[['lower']],
                             upper = pars_df[['upper']],
                             f = sk_nLL,
                             method = method,
                             g_obs = g_obs,
                             pars_fix = pars_fix,
                             X = X,
                             iso = iso,
                             quiet = quiet,
                             log_scale = log_scale,
                             control = control) |> suppressWarnings()
  }

  # unpack optimizer output
  obj_val = optim_out[['value']]
  pars_fitted_v = optim_out[['par']]
  if(!quiet) cat(paste('\n', optim_out[['message']]))

  # revert to initial value if optimize/optim result was not an improvement
  if(bds_nLL[['initial']] < obj_val) pars_fitted_v = pars_df[['initial']]
  # TODO: further checking if upper/lower bounds were also better

  # reshape as list and transform back if fitting on log-scale
  pars_fitted = sk_pars_update(pars_fix, pars_fitted_v, iso=iso, na_omit=TRUE)
  pars_df['fitted'] = pars_fitted_v
  if(log_scale)
  {
    pars_df = exp(pars_df)
    pars_fitted = sk_pars_update(pars_fitted, exp(sk_pars_update(pars_fitted)))
  }

  # return parameters list and data frame in a list
  df_order = c('lower', 'initial', 'fitted', 'upper')
  return(list(pars=pars_fitted, df=pars_df[df_order], obj=obj_val))
}











