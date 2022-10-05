# bk_var.R
# Dean Koch, 2022
# Functions for matrix algebra with correlation and covariance matrices

#' Stationary 1D correlation kernels
#'
#' Computes stationary correlation function values for the n (non-negative) 1-dimensional
#' distances in `d`. Parameter list entry `pars$kp` supplies the kernel parameter(s).
#'
#' `pars$k` must be one of the following kernel names:
#'
#' * 'exp': exponential (special case of 'gex' with shape p=1)
#' * 'gau': gaussian/stable (special case of 'gex' with shape p=2)
#' * 'sph': spherical (AKA stable/Gaussian for p=2)
#'
#' * 'gex': gamma-exponential (with shape p)
#' * 'mat': Whittle-Matern (Handcock and Wallis parameterization, with shape kap)
#'
#' where the first three kernels have only a range parameters, and the last two have both a
#' range and shape parameter.
#'
#' For the 1-parameter kernels, `pars$kp` is the range parameter value ('rho'); For the
#' 2-parameter kernels, `pars$kp` is a vector whose first element is 'rho', and second
#' element is the shape parameter ('p' or 'kap'). The names in `pars$kp` are ignored and
#' only the order matters - the range parameter always comes first.
#'
#' Note that this function will not accept parameter lists `pars` of the form returned by
#' `bk_pars(...)` etc, as these include a pair of 1d kernels (however the sub-lists
#' `pars$y` and `pars$x` are accepted).
#'
#' @param pars list with elements 'k', the kernel name, and 'kp' the parameter vector
#' @param d numeric vector of length n, the distances to evaluate
#'
#' @return length-n vector or a list of parameters and bounds (see details)
#' @export
#'
#' @examples
#'
#' # define test distances, grid, and example kernel
#' n_test = 100
#' d_test = seq(n_test)-1
#' g_example = bk(n_test)
#' pars = bk_pars(g_example, c('mat', 'gau'))
#' pars_x = pars[['x']]
#'
#' # compute and plot the x component of the correlogram function
#' corr_x_example = bk_corr(pars_x, d=d_test)
#' plot(d_test, corr_x_example, pch=NA)
#' lines(d_test, corr_x_example)
#'
#' ## show how this function gets used to build more complicated objects
#'
#' # get the other component correlation, take product
#' pars_y = pars[['y']]
#' corr_y_example = bk_corr(pars_y, d=d_test)
#' corr_example = corr_y_example * corr_x_example
#'
#' # variogram
#' variogram_example = bk_vario_fun(pars, d=list(y=d_test, x=d_test))
#' variogram_compare = 2 * pars$eps + pars$psill * (1 - corr_example)
#' max(abs( variogram_example - variogram_compare ))
#'
#' # Toeplitz component matrices built entirely from these correlation vectors
#' variance_matrix_example = bk_var(g_example, pars)
#' max(abs( variance_matrix_example[['y']][,1L] - corr_y_example ))
#' max(abs( variance_matrix_example[['x']][,1L] - corr_x_example ))
#'
#'
bk_corr = function(pars, d=NA)
{
  # handle invalid pars
  if( !all( c('k', 'kp') %in% names(pars) ) ) stop('pars must be list with elements "k" and "kp"')
  k_name = pars[['k']]

  # handle five known values for kernel name
  nm_expected = c('exp', 'gau', 'sph', 'gxp', 'mat')
  msg_unrecognized = paste('unrecognized kernel name:', k_name)
  if( !(k_name %in% nm_expected) ) stop(msg_unrecognized)

  # exponential
  if(k_name == 'exp')
  {
    # new parameter list for recursive call
    pars = list(k='gxp', kp=c(rho=pars[['kp']][1], p=1L))
    return( bk_corr(pars, abs(d)) )
  }

  # gaussian/stable
  if(k_name == 'gau')
  {
    # new parameter list for recursive call
    pars = list(k='gxp', kp=c(rho=pars[['kp']][1], p=2L))
    return( bk_corr(pars, d) )
  }

  # spherical
  if(k_name == 'sph')
  {
    # assign parameter and process truncation distance
    ds = d / pars[['kp']][1]
    cvals = rep(0, length(d))
    idx.nz = ds < 1

    # evaluate on non-truncated distances and return
    cvals[idx.nz] = 1 - ( (3/2) * ds[idx.nz] ) + ( (1/2) * ( ds[idx.nz]^3 ) )
    return( cvals )
  }

  # gamma-exponential
  if(k_name == 'gxp')
  {
    # return suggested bounds and initial value if requested
    if( length( pars[['kp']] ) != 2 ) stop(paste('pars$kp must be a vector of form c(rho, p)'))

    # assign parameters and evaluate
    return( exp( -( ( d / pars[['kp']][1] )^pars[['kp']][2] ) ) )
  }

  # Whittle-Matern
  if(k_name == 'mat')
  {
    # assign parameters and compute scaling constant
    kap = pars[['kp']][2]
    rho = pars[['kp']][1] / ( 2 * sqrt(kap) )
    sc = ( (2^( 1 - kap ) ) / gamma(kap) )

    # some preprocessing addressing numerical precision issues
    cvals = rep(0, length(d))
    ds = d / rho
    idx.z = ds == 0
    cvals[idx.z] = 1
    bk = besselK(ds, kap)
    idx.big = bk == 0

    # evaluate on well-behaved inputs and return
    idx.eval = !idx.big & !idx.z
    cvals[idx.eval] = sc * (ds[idx.eval]^kap) * bk[idx.eval]
    return(cvals)
  }

  # if we got this far, input `pars$k` didn't match anything
  stop(paste('kernel name', pars$k, 'not recognized'))
}


#' Construct 1D stationary correlation matrices for regularly spaced data
#'
#' The i,jth value of the returned correlation matrix is the marginal correlation between
#' the ith and jth points in a regularly spaced sequence of `n` (1D) points, given the
#' correlation model with parameters defined in list `pars`.
#'
#' This matrix is symmetric and Toeplitz as a result of the assumption of stationarity
#' of the random field and regularity of the grid.
#'
#' The distance between adjacent points is specified by `gres`. Subsets of
#' the correlation matrix can be requested by specifying `i` and/or `j` (default
#' behaviour is to include all).
#'
#' Like `bk_corr`, this function is for computing 1D components of a 2D process.
#' The product of two matrices returned by `bk_corr_mat` is the correlation
#' matrix for a spatially separable process (see examples).
#'
#' @param pars list of kernel parameters 'k' and 'kp' (see `bk_corr`)
#' @param n positive integer, the number of points on the 1D line
#' @param gres positive numeric, the distance between adjacent grid lines
#' @param i vector, a subset of `seq(n)` indicating rows to return
#' @param j vector, a subset of `seq(n)` indicating columns to return
#'
#' @return the n x n correlation matrix, or its subset as specified in `i`, `j`
#' @export
#'
#' @examples
#'
#' # define test distances, grid, and example kernel
#' n_test = 10
#' g_example = bk(n_test)
#' pars = bk_pars(g_example, c('mat', 'gau'))
#'
#' # compute the correlation matrices and their kronecker product
#' cx = bk_corr_mat(pars[['x']], n=n_test)
#' cy = bk_corr_mat(pars[['y']], n=n_test)
#' cxy = kronecker(cx, cy)
#'
#' # bk_var returns these two matrices in a list...
#' max(abs( bk_var(g_example, pars)[['y']] - cy ))
#' max(abs( bk_var(g_example, pars)[['x']] - cx ))
#'
#' # ... or it can compute the full covariance matrix for model pars
#' var_matrix = bk_var(g_example, pars, sep=FALSE)
#' var_matrix_compare = (pars$psill*cxy) + diag(pars$eps, n_test^2)
#' max(abs( var_matrix - var_matrix_compare ))
#'
#' # extract a subgrid without computing the whole thing
#' cx_sub = bk_corr_mat(pars_x, n=n_test, i=2:4, j=2:4)
#' cx_sub - cx[2:4, 2:4]
#'
#' # gres scales distances. Increasing gres causes correlations to decrease
#' cx_long = bk_corr_mat(pars_x, n=n_test, gres=2*g_example$gres)
#' cx_long < cx
#'
bk_corr_mat = function(pars, n, gres=1, i=seq(n), j=seq(n))
{
  # compute the set of distances over which we need to evaluate kernel
  du = gres * ( seq(n) - 1 )

  # compute kernel values for these distances
  dcorr = bk_corr(pars, du)

  # build large vector to shift through in building the Toeplitz output matrix
  bigvec = c(dcorr[n:2], dcorr)

  # build and return the matrix
  return( sapply(j, function(x) bigvec[ (n-x) + i ]) )
}


#' Generate a covariance matrix or its factorization
#'
#' Computes the covariance matrix V (or one of its factorizations) for the non-NA points
#' in grid `g_obs`, given the model parameters list `pars`
#'
#' By default the output matrix is V. Alternatively, if `X` is supplied, the function
#' returns the quadratic form X^T V^{-1} X.
#'
#' When `fac_method=='eigen'` the function instead returns the eigen-decomposition of the
#' output matrix, and when `fac_method=='chol'` its lower triangular Cholesky factor is
#' returned. Supplying this factorization in argument `fac` in a subsequent call with `X`
#' can speed up calculations. `fac` is ignored when `X` is not supplied.
#'
#' `scaled=TRUE` returns the matrix scaled by the reciprocal of the partial sill,
#' `1/pars$psill`, before factorization. This is the form expected by functions
#' `bk_var_mult` and `bk_LL` in argument `fac`.
#'
#' If all grid points are observed, then the output V becomes separable. Setting `sep=TRUE`
#' in this case causes the function to returns the x and y component correlation matrices (or
#' their factorizations, as requested in `fac_method`) separately, in a list. `scaled` has no
#' effect in this output mode. Note also that `sep` has no effect when `X` is supplied, as
#' the quadratic form is not generally separable (regardless of separability in V^{-1}).
#'
#' Missing data are identified by looking for NAs in the data vector `g_obs$gval`. If all
#' are NA (or if 'gval' is missing from `g_obs`), the function behaves as though all grid
#' points are observed. For multi-layer input, NAs are instead determined from
#' `g_obs$idx_grid` and 'gval' is ignored (see `?bk`).
#'
#'
#' Note: when `pars$eps>0`, the 'eigen' factorization method will be more robust than 'chol'
#' in handling numerical precision issues with poorly conditioned covariance matrices.
#'
#' @param g_obs list of form returned by `bk` (with entries 'gdim', 'gres', 'gval')
#' @param pars list of form returned by `bk_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param scaled logical, whether to scale by `1/pars$psill`
#' @param fac_method character, the factorization to return, one of 'none', 'chol', 'eigen'
#' @param X numeric matrix, the `X` in `t(X) %*% V %*% X` (default is identity)
#' @param fac matrix or list of matrices, the variance factorization (only used with X)
#' @param sep logical, indicating to return correlation components instead of full covariance matrix
#'
#' @return either matrix `V`, or X^T V^{-1} X, or a factorization ('chol' or 'eigen')
#' @export
#'
#' @examples
#' # define example grid with NAs and example predictors matrix
#' gdim = c(12, 13)
#' n = prod(gdim)
#' n_obs = floor(n/3)
#' idx_obs = sort(sample.int(n, n_obs))
#' g = g_obs = bk(gdim)
#' g_obs$gval[idx_obs] = rnorm(n_obs)
#'
#' # example kernel
#' psill = 0.3
#' pars = bk_pars(g_obs) |> modifyList(list(psill=psill))
#'
#' # plot the covariance matrix for observed data, its cholesky factor and eigen-decomposition
#' V_obs = bk_var(g_obs, pars)
#' V_obs_chol = bk_var(g_obs, pars, fac_method='chol')
#' V_obs_eigen = bk_var(g_obs, pars, fac_method='eigen')
#' bk_plot(V_obs)
#' bk_plot(V_obs_chol)
#' bk_plot(V_obs_eigen$vectors)
#'
#' # case when there are no NAs (or no data at all)
#' g_nodata = modifyList(g_obs, list(gval=NULL))
#'
#' # get the full covariance matrix with sep=FALSE (default)
#' V_full = bk_var(g_nodata, pars)
#' max(abs( V_obs - V_full[idx_obs, idx_obs] ))
#'
#' # get 1d correlation matrices with sep=TRUE...
#' corr_components = bk_var(g_nodata, pars, sep=TRUE)
#' str(corr_components)
#'
#' # ... these are related to the full covariance matrix by psill and eps
#' corr_mat = kronecker(corr_components[['x']], corr_components[['y']])
#' V_full_compare = pars$psill * corr_mat + diag(pars$eps, n)
#' max(abs(V_full - V_full_compare))
#'
#' # ... their factorizations can be returned as (nested) lists
#' str(bk_var(g_nodata, pars, fac_method='chol', sep=TRUE))
#' str(bk_var(g_nodata, pars, fac_method='eigen', sep=TRUE))
#'
#' # compare to the full covariance matrix factorizations (default sep=FALSE)
#' str(bk_var(g_nodata, pars, fac_method='chol'))
#' str(bk_var(g_nodata, pars, fac_method='eigen'))
#'
#' # test quadratic form with X
#' nX = 3
#' X_all = cbind(1, matrix(rnorm(nX * n), ncol=nX))
#' cprod_all = crossprod(X_all, chol2inv(chol(V_full))) %*% X_all
#' abs(max(bk_var(g, pars, X=X_all) - cprod_all ))
#'
#' # test products with inverse of quadratic form with X
#' mult_test = rnorm(nX+1)
#' cprod_all_inv = chol2inv(chol(cprod_all))
#' cprod_all_inv_chol = bk_var(g, pars, X=X_all, scaled=TRUE, fac_method='eigen')
#' bk_var_mult(mult_test, pars, fac=cprod_all_inv_chol) - cprod_all_inv %*% mult_test
#'
#' # repeat with missing data
#' X_obs = X_all[idx_obs,]
#' cprod_obs = crossprod(X_obs, chol2inv(chol(V_obs))) %*% X_obs
#' abs(max(bk_var(g_obs, pars, X=X_obs) - cprod_obs ))
#' cprod_obs_inv = chol2inv(chol(cprod_obs))
#' cprod_obs_inv_chol = bk_var(g_obs, pars, X=X_obs, scaled=T, fac_method='eigen')
#' bk_var_mult(mult_test, pars, fac=cprod_obs_inv_chol) - cprod_obs_inv %*% mult_test
#'
#' # `scaled` indicates to divide matrix by psill
#' print( pars[['eps']]/pars[['psill']] )
#' diag(bk_var(g_obs, pars, scaled=TRUE)) # diagonal elements equal to 1 + eps/psill
#' ( bk_var(g_obs, pars) - psill * bk_var(g_obs, pars, scaled=TRUE) ) |> abs() |> max()
#' ( bk_var(g_obs, pars, X=X_obs, scaled=TRUE) - ( cprod_obs/psill ) ) |> abs() |> max()
#'
#' # in cholesky factor this produces a scaling by square root of psill
#' max(abs( V_obs_chol - sqrt(psill) * bk_var(g_obs, pars, fac_method='chol', scaled=TRUE) ))
#'
#' # and in the eigendecomposition, a scaling of the eigenvalues
#' vals_scaled = bk_var(g_obs, pars, fac_method='eigen', scaled=TRUE)$values
#' max(abs( bk_var(g_obs, pars, fac_method='eigen')$values - psill*vals_scaled ))
#'
bk_var = function(g_obs, pars=NULL, scaled=FALSE, fac_method='none', X=NULL, fac=NULL, sep=FALSE)
{
  # default Gaussian kernel
  if(is.null(pars)) pars = bk_pars(g_obs, 'gau')

  # check for unknown fac_method
  nm_fac_method = c('none', 'eigen', 'chol')
  msg_fac_method = paste('fac_method must be one of: NA,', paste(nm_fac_method, collapse=', '))
  if( !(fac_method %in% nm_fac_method) ) stop(msg_fac_method)

  # bk converts 1-layer matrix input to vector
  if(is.matrix(g_obs[['gval']])) g_obs = bk(g_obs)

  # copy the indexing vector in the multi-layer case
  if(is.matrix(g_obs[['gval']]))
  {
    # only the NA structure matters here
    g_obs = modifyList(g_obs, list(gval=g_obs[['idx_grid']], idx_grid=NULL))
  }

  # identify NA grid-points and flag separable case (no missing data, or all missing)
  n = prod(g_obs[['gdim']])
  is_obs = !is.na(g_obs[['gval']])
  if( !any(is_obs) | all(is_obs) ) g_obs[['gval']] = NULL
  if( is.null(g_obs[['gval']]) ) is_obs = rep(TRUE, n)
  n_obs = sum(is_obs)
  is_sep = n == n_obs

  # predictor matrix case
  if( !is.null(X) )
  {
    # call without X to get eigen-decomposition of variance
    if( is.null(fac) ) fac = bk_var(g_obs, pars, scaled=TRUE, fac_method='eigen', sep=is_sep)

    # check for invalid input
    msg_class = 'mu must be a matrix predictor columns'
    msg_mismatch = 'nrow(X) must equal the number of non-NA points in g_obs$gval'
    if( !is.matrix(X) ) stop(msg_class)
    if( nrow(X) != sum(is_obs) ) stop(msg_mismatch)

    # quadratic form of whitened data matrix (not the inverse)
    p_scale = ifelse(scaled, pars[['psill']], 1)
    X_quad = bk_var_mult(X, pars, fac=fac, quad=TRUE) / p_scale

    # return requested decomposition
    if(fac_method == 'none') return(X_quad)
    if(fac_method == 'chol' ) return(t(chol(X_quad)))
    if(fac_method == 'eigen' ) return(eigen(X_quad, symmetric=TRUE))
  }

  # unpack grid config and covariance parameters
  gres = g_obs[['gres']]
  gdim = g_obs[['gdim']]
  eps = ifelse(scaled, pars[['eps']]/pars[['psill']], pars[['eps']])
  psill = ifelse(scaled, 1, pars[['psill']])

  # complete data case
  if(is_sep)
  {
    # compute the full component correlation matrices
    cy = bk_corr_mat(pars[['y']], gdim[['y']], gres[['y']])
    cx = bk_corr_mat(pars[['x']], gdim[['x']], gres[['x']])

    # return these (or their factorizations) separately in a list...
    if(sep)
    {
      if(fac_method == 'none') return( list(y=cy, x=cx) )
      if(fac_method == 'chol') return( list( y=t(chol(cy)), x=t(chol(cx)) ) )
      if(fac_method == 'eigen') return( list( y=eigen(cy), x=eigen(cx) ) )

      # ... or construct full covariance matrix from their kronecker product
    } else {

      V = diag(eps, n_obs) + ( psill * kronecker(cx, cy) )
      if(fac_method == 'none') return(V)
      if(fac_method == 'chol') return( t(chol(V)) )
      if(fac_method == 'eigen') return( eigen(V) )
    }
  }

  # incomplete data case

  # build mapping from points to rows of component matrices
  yx_idx = which(is_obs) |> bk_vec2mat(gdim['y'], out='list')

  # draw a selection of rows/columns from component correlation matrices
  cy_obs = bk_corr_mat(pars[['y']], gdim[['y']], gres[['y']], i=yx_idx[['i']], j=yx_idx[['i']])
  cx_obs = bk_corr_mat(pars[['x']], gdim[['x']], gres[['x']], i=yx_idx[['j']], j=yx_idx[['j']])

  # Hadamard product produces the correlation matrix for observed grid points
  if(fac_method == 'none') return( ( psill * cy_obs * cx_obs ) + diag( rep(eps, n_obs) ) )

  # compute Cholesky factor (lower triangular) of correlation matrix
  if(fac_method == 'chol') return( t( chol( psill * cy_obs * cx_obs + eps * diag(rep(1, n_obs)) ) ) )

  # fac_method='eigen': compute eigen-decomposition of correlation matrix
  eigen_result = eigen(cy_obs*cx_obs, symmetric=TRUE)
  eigen_result[['values']] = ( psill * eigen_result[['values']] ) + eps
  return(eigen_result)
}

#' Multiply a vector by a power of the covariance matrix
#'
#' Computes `W %*% z`, where `z` is the vector of non-NA data in `g_obs`,
#' and `W` is the `p`th power of `V`, the covariance matrix for `z`. By default,
#' `p=-1`, so the function computes products with the inverse covariance matrix.
#'
#' Alternatively, `out='quad'` computes the quadratic form `t(z) %*% W %*% z`.
#'
#' `fac_method` specifies the covariance matrix factorization to use: either 'chol'
#' (Cholesky factorization), which only supports `p=-1`; or 'eigen' (eigen-decomposition),
#' which supports any integer or numeric power for `p`. By default, the 'eigen' method
#' is used, unless a Cholesky decomposition (matrix) is passed in `fac`.
#'
#' The 'eigen' method is required when `g_obs` has complete data (ie no NA values). Note
#' that the structure of any supplied `fac` overrules the `fac_method` argument, so if your
#' `g_obs` is complete and you supply the Cholesky decomposition, the function will throw
#' an error.
#'
#' As factorization is the slow part of the computations, it can be pre-computed
#' using `bk_var(..., scaled=TRUE)` and passed to `bk_var_mult` in argument `fac`.
#' This must be the factorization of the covariance matrix after dividing by the partial sill
#' (see `?bk_var`); either the lower Cholesky factor (eg the transposed output of
#' `base::chol`), or a list of eigen-vectors and eigen-values (eg. the output of `base::eigen`).
#' In the separable case, the eigen-decomposition is done on each of the x and y components
#' separately, and `fac` should be a list with elements 'x' and 'y', each one a list of
#' eigen-vectors and eigen-values.
#'
#' When a factorization is supplied, all entries in `pars`, except for `psill`, are ignored,
#' as they are baked into the factorization already. `g_obs` can in this case be a numeric vector
#' or matrix, containing one or more layers of observed data (with NAs omitted).
#'
#'
#' @param g_obs list of form returned by `bk` or numeric vector or matrix of non-NA data
#' @param pars list of form returned by `bk_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param fac factorization of scaled covariance matrix of z (V divided by psill)
#' @param quad logical, if TRUE the function returns the quadratic form `t(z) %*% V_inv %*% z`
#' @param p numeric, the matrix power of V^p to multiply (ignored when `method=='chol'`)
#'
#' @return numeric matrix
#' @export
#'
#' @examples
#' # relative error comparing output x to reference y
#' rel_err = \(x, y) ifelse(y == 0, 0, abs( (x - y) / y ) )
#'
#' # define example grid and data
#' gdim = c(10, 15)
#' n = prod(gdim)
#' z_all = rnorm(n)
#' g_obs = modifyList(bk(gdim), list(gval = z_all))
#'
#' # define covariance parameters
#' pars = modifyList(bk_pars(g_obs, 'gau'), list(psill=2, eps=0.5))
#'
#' # COMPLETE CASE
#'
#' # compute the full covariance matrix
#' V = bk_var(g_obs, pars, sep=FALSE)
#' V_inv = chol2inv(chol(V))
#' out_reference = V_inv %*% z_all
#' out_reference_quad = t(z_all) %*% out_reference
#' max( rel_err(bk_var_mult(g_obs, pars), out_reference) )
#' rel_err(bk_var_mult(g_obs, pars, quad=TRUE), out_reference_quad)
#'
#' # pre-computed factorization on separable components of correlation matrix
#' fac_corr = bk_var(modifyList(g_obs, list(gval=NULL)), pars, fac_method='eigen', sep=TRUE)
#' max( rel_err(bk_var_mult(g_obs, pars, fac=fac_corr), out_reference) )
#' rel_err(bk_var_mult(g_obs, pars, fac=fac_corr, quad=TRUE), out_reference_quad)
#'
#' # matrix powers
#' out_reference = V %*% z_all
#' max( rel_err(bk_var_mult(g_obs, pars, fac_method='eigen', p=1), out_reference) )
#' rel_err(bk_var_mult(g_obs, pars, fac_method='eigen', p=1, quad=T), t(z_all) %*% out_reference)
#'
#' # INCOMPLETE CASE
#'
#' n_sample = floor(n/10)
#' idx_sampled = sample.int(n, n_sample) |> sort()
#' z_obs = rep(NA, n)
#' z_obs[idx_sampled] = z_all[idx_sampled]
#' g_obs = modifyList(g_obs, list(gval = z_obs))
#' V = bk_var(g_obs, pars)
#' bk_plot(V)
#'
#' # correctness check (eigen used by default)
#' z = matrix(z_obs[idx_sampled], ncol=1)
#' V_inv = chol2inv(chol(V))
#' out_reference = (V_inv %*% z)
#' out_reference_quad = t(z) %*% out_reference
#' max(rel_err(bk_var_mult(g_obs, pars), out_reference))
#' rel_err(bk_var_mult(g_obs, pars, quad=TRUE), out_reference_quad)
#'
#' # check non-default Cholesky method
#' max( rel_err(bk_var_mult(g_obs, pars, fac_method='chol'), out_reference) )
#' rel_err(bk_var_mult(g_obs, pars, quad=TRUE, fac_method='chol'), out_reference_quad)
#'
#' # supply data as a vector instead of list by pre-computing factorization
#' fac_chol = bk_var(g_obs, pars, scaled=TRUE, fac_method='chol')
#' fac_eigen = bk_var(g_obs, pars, scaled=TRUE, fac_method='eigen')
#' max(rel_err(bk_var_mult(z, pars, fac=fac_chol), out_reference))
#' max(rel_err(bk_var_mult(g_obs, pars, fac=fac_eigen), out_reference))
#' rel_err(bk_var_mult(z, pars, fac=fac_chol, quad=TRUE), out_reference_quad)
#' rel_err(bk_var_mult(g_obs, pars, fac=fac_eigen, quad=TRUE), out_reference_quad)
#'
#' # matrix powers in eigen mode
#' out_reference = V %*% z
#' max(rel_err(bk_var_mult(g_obs, pars, p=1), out_reference))
#' rel_err(bk_var_mult(g_obs, pars, p=1, quad=TRUE), t(z) %*% out_reference)
#' max(rel_err(bk_var_mult(g_obs, pars, p=2), V %*% out_reference))
#'
#' # multiply g_obs twice by a square root of V
#' g_obs_sqrt = g_obs
#' g_obs_sqrt$gval[!is.na(g_obs$gval)] = bk_var_mult(g_obs, pars, p=1/2)
#' max( rel_err(bk_var_mult(g_obs_sqrt, pars, p=1/2), out_reference) )
#'
bk_var_mult = function(g_obs, pars, fac_method='eigen', fac=NULL, quad=FALSE, p=-1)
{
  # check for invalid factorization method name
  if( !(fac_method %in% c('chol', 'eigen')) ) stop('fac_method unrecognized')

  # reshape input as matrix
  if( is.list(g_obs) )
  {
    # unpack list g_obs as needed
    g_obs = bk(g_obs)
    n_layer = ifelse(is.matrix(g_obs[['idx_grid']]), ncol(g_obs[['idx_grid']]), 1L)

    # when there are missing data, omit from copy z
    is_obs = !is.na(g_obs[['gval']])
    z = matrix(g_obs[['gval']][is_obs], ncol=n_layer)

    # complete and empty cases trigger separability option below
    is_sep = all(is_obs) | !any(is_obs)

  } else {

    # coerce g_obs to matrix
    z = matrix(g_obs, ncol=ifelse(is.vector(g_obs), 1, ncol(g_obs)))

    # check if the supplied factorization is a Kronecker product
    if( is.null(fac) ) stop('factorization fac must be supplied if g_obs is not a list')
    is_sep = is.list(fac) & all(c('y', 'x') %in% names(fac))
  }

  # separability property and eigen-decomposition used in no-missing case
  if(is_sep)
  {
    # factorize the correlation matrix components - eigen method is forced in this case
    g_bare = g_obs[c('gres', 'gdim')]
    if( is.null(fac) ) fac = bk_var(g_bare, pars=pars, fac_method='eigen', sep=TRUE)

    # check for problems with supplied factorization type
    msg_fac = 'expected separable eigen-decomposition in fac (but got Cholesky)'
    if( !is.list(fac) ) stop(msg_fac)
    if( !all(c('y', 'x') %in% names(fac) ) ) stop('named entries "x" and "y" not found in fac')
    if( !all(sapply(fac, function(f) is.list(f)) ) ) stop(msg_fac)

    # eigenvalues of full correlation and covariance matrices
    ny = length(fac[['y']][['values']])
    nx = length(fac[['x']][['values']])
    ev_corr = kronecker(fac[['x']][['values']], fac[['y']][['values']])
    ev_p = ( pars[['eps']] + pars[['psill']] * as.numeric(ev_corr) )^p
    if( any( is.infinite(ev_p) ) ) stop('ill-conditioned covariance matrix (0 eigenvalue)')
    if( any( is.nan(ev_p) ) ) stop('ill-conditioned covariance matrix (eigenvalue < 0)')
    if( nrow(z) != (ny*nx) ) stop('factorization was inconsistent with dimensions of g_obs')

    # left multiply data vector by square root of inverse correlation matrix
    grid_evec = apply(z, 2, \(x) crossprod(fac[['y']][['vectors']], matrix(x, ny)), simplify=FALSE ) |>
      lapply(\(x) tcrossprod(x, t(fac[['x']][['vectors']]) ) )

    # quadratic form is the scaled vector multiplied by its transpose
    if(quad) return( crossprod(sqrt(ev_p) * sapply(grid_evec, as.numeric)) )

    # scale by inverse covariance eigen-values and transform back
    grid_prod = grid_evec |>
      lapply(\(x) tcrossprod(fac[['y']][['vectors']], t(ev_p * x))) |>
      sapply(\(x) tcrossprod(x, fac[['x']][['vectors']]) )

    return( grid_prod )
  }

  # code below handles non-separable case
  if( is.null(z) ) stop('data vector not found g_obs')

  # set default factorization method and switch to eigen when needed, with a warning
  if( is.null(fac_method) ) fac_method = ifelse(p==-1, 'chol', 'eigen')
  if( (p != -1) & (fac_method == 'chol') )
  {
    warning('switching to fac_method="eigen" (p=-1 required for Cholesky)')
    fac_method = 'eigen'
  }

  # factorize the variance matrix using the specified method
  if( is.null(fac) ) fac = bk_var(g_obs, pars=pars, scaled=TRUE, fac_method=fac_method)

  # detect supplied factorization type
  fac_method = ifelse(is.matrix(fac), 'chol', 'eigen')

  # computation via Cholesky factor (p=-1)
  if( fac_method == 'chol' )
  {
    # consistency check
    if( nrow(fac) != nrow(z) ) stop('factorization was inconsistent with dimensions of g_obs')

    # solve for (t(C))(V_inv)(z) whose inner product is the quadratic form
    g_chol = forwardsolve(fac, z/pars[['psill']] )
    if(quad) return( pars[['psill']] * crossprod(g_chol) )
    grid_prod = backsolve(t(fac), g_chol)
  }

  # computation via eigen-decomposition
  if( fac_method == 'eigen' )
  {
    # eigenvalues of square root of the requested power
    ev_sqrt = ( pars[['psill']] * fac[['values']] )^(p/2)

    # consistency check
    if( length(ev_sqrt) != nrow(z) ) stop('factorization was inconsistent with dimensions of g_obs')
    if( any( is.infinite(ev_sqrt) ) ) stop('ill-conditioned covariance matrix (0 eigenvalue)')
    if( any( is.nan(ev_sqrt) ) ) stop('ill-conditioned covariance matrix (eigenvalue < 0)')
    grid_evec = ev_sqrt * crossprod(fac[['vectors']], z)
    if(quad) return( crossprod(grid_evec, grid_evec) )

    # left-multiply by the remaining terms (transpose of the original transform)
    grid_prod = tcrossprod(fac[['vectors']], t( ev_sqrt * grid_evec ))
  }

  return( grid_prod )
}


#' Efficiently compute yzx for symmetric Toeplitz matrices y and x
#'
#' Computes the product `y %*% z` or `y %*% z %*% x` for symmetric Toeplitz matrices
#' `y` and `x` and any numeric matrix `z`.
#'
#' Argument(s) `y` (and `x`) can be vector(s) supplying the first row of the matrix.
#' By default, `z` is the identity matrix, so for matrix `y`, bk_toep_mult(`y`) returns
#' its argument, and for vector `y`, it returns the Toeplitz matrix generated by `y`, the
#' same as `base::toeplitz(y)`.
#'
#' Fast Fourier transforms are used to reduce the memory footprint of computations,
#' The first row(s) of `y` (and `x`) are embedded in a zero-padded vector representing a
#' circulant matrix, whose action on the zero-padded version of `z` is equivalent to
#' element-wise product in Fourier space. This allows the desired matrix product to be
#' computed without explicitly creating matrices `y` or `x` in memory.
#'
#' The function is optimized for grid data `z` that are sparse (many zeros). Before
#' computing any transformations it first scans for and removes columns and rows of
#' z which are all zero, replacing them afterwards.
#'
#' To avoid unnecessarily copying large sparse matrices, `z` can be the vector of
#' non-zero matrix entries only, where `gdim` specifies the full matrix dimensions and
#' `idx_obs` the indices of the non-zero entries.
#'
#' @param y numeric matrix or vector, the symmetric Toeplitz matrix y or its first row
#' @param z numeric matrix or vector with dimensionality conforming with y (and x)
#' @param x numeric matrix or vector, the symmetric Toeplitz matrix x or its first row
#' @param idx_obs integer vector, indices of the observed grid points
#'
#' @return numeric matrix, the product of yzx or yz (if x is NULL)
#' @export
#'
#' @examples
#' # define example matrix from 1D exponential variogram
#' n = 10
#' y = exp(1-seq(n))
#' y_mat = bk_toep_mult(y)
#' max( abs(y_mat - stats::toeplitz(y))  )
#'
#' # multiply by random matrix and compare with default matrix multiply
#' z = matrix(rnorm(n^2), n)
#' result_default = y_mat %*% z
#' max( abs( result_default - bk_toep_mult(y_mat, z) ) )
#'
#' # save memory by passing only the first row of the Toeplitz matrix
#' max( abs( result_default - bk_toep_mult(y, z) ) )
#'
#' # sparsify z and repeat
#' idx_sparse = sample.int(n^2, n^2 - n)
#' z[idx_sparse] = 0
#' result_default = y_mat %*% z
#' max( abs( result_default - bk_toep_mult(y, z) ) )
#'
#' # right-multiply with another kernel
#' x = exp( 2 *( 1-seq(n) ) )
#' x_mat = bk_toep_mult(x)
#' result_default = result_default %*% x_mat
#' max( abs( result_default - bk_toep_mult(y, z, x) ) )
#'
#' # z can also be supplied as vector of nonzero grid values
#' idx_obs = which(z != 0)
#' gdim = c(y=n, x=n)
#' max( abs( result_default - bk_toep_mult(y, z=z[idx_obs], x, idx_obs, gdim) ) )
#'
bk_toep_mult = function(y, z=NULL, x=NULL, idx_obs=NULL, gdim=NULL)
{
  # tolerance in Toeplitz check
  ttol = .Machine[['double.eps']]

  # copy first row of y and validate input
  if( is.matrix(y) )
  {
    y_first = as.numeric(y[1,])
    if( diff(dim(y)) != 0 ) stop('input matrix was not symmetric')
    if( max( abs(y - toeplitz(y_first)) ) > ttol ) stop('input matrix was not Toeplitz')
    y = y_first
  }

  # find amount of zero padding needed to get next highest composite dimension
  n_y = length(y)
  n_pad = 2*( stats::nextn(n_y) - n_y )
  z_pad = rep(0, n_y + n_pad)

  # normalization constant for inverse fft
  z_idx = seq(n_y)
  norm_constant = 2*n_y + n_pad

  # set default for z (the identity)
  if( is.null(z) ) z = diag(1, n_y)

  # reconstruct full grid with zeros for missing data
  if( !is.null(idx_obs) )
  {
    if( is.null(gdim) ) stop('gdim must be supplied with idx_obs')
    z_nz = z
    z = matrix(0, nrow=gdim['y'], ncol=gdim['x'])
    z[idx_obs] = z_nz
  }

  # omit zero columns from computation
  if( !is.matrix(z) ) z = matrix(z, n_y)
  is_col_empty = colSums(abs(z)) == 0
  n_compute = ncol(z) - sum(is_col_empty)

  # initialize zero-padded z matrix and copy input data
  z_pad = matrix(numeric(1), norm_constant, n_compute)
  z_pad[z_idx,] = z[, !is_col_empty]

  # discrete FFT of circulant matrix containing y
  fy = stats::fft( c(y, rep(0, 1 + n_pad), y[n_y:2]) )

  # add zero padding, transform by column, multiply, then transform back
  z[, !is_col_empty] = Re( stats::mvfft(fy * stats::mvfft(z_pad), inverse=TRUE)[z_idx,] )
  if( is.null(x) ) return(z/norm_constant)

  # do right-multiplication by transposing output z and passing back to bk_toep_mult
  return( t( bk_toep_mult(x, t(z/norm_constant)) ) )
}

