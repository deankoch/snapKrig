# sk_estim.R
# Dean Koch, 2022
# Functions for parameter inference and spatial prediction

#' Generalized least squares (GLS) with Kronecker covariances for sk grids
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
#' pars = utils::modifyList(sk_pars(g_empty, 'gau'), list(psill=2))
#'
#' # generate spatial noise
#' g_noise = sk_sim(g_empty, pars)
#' plot(g_noise)
#'
#' # generate more spatial noise to use as covariates
#' n_betas = 3
#' betas = stats::rnorm(n_betas, 0, 10)
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
#' # verify GLS results manually
#' X_all_obs = cbind(1, X_obs)
#' V = sk_var(g_miss, pars)
#' z = g_miss[!is.na(g_miss)]
#' X_trans = t(X_all_obs) %*% solve(V)
#' betas_compare = solve( X_trans %*% X_all_obs ) %*% X_trans %*% z
#' betas_compare - betas_est
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
    is_obs = if(is.null(g[['is_obs']])) !is.na(g[['idx_grid']]) else g[['is_obs']]
    n_layer = ncol(g[['gval']])
    z = g[['gval']]

  } else {

    # copy non-NA data as a 1-column matrix
    is_obs = if(is.null(g[['is_obs']])) !is.na(g[['gval']]) else g[['is_obs']]
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

    # unpack sk grid X as matrix
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


#' Compute kriging predictor (or variance) for an sk grid
#'
#' Evaluates the kriging prediction equations to find the expected value (mean) of
#' the spatial process for `g` at all grid points, including unobserved ones.
#'
#' This predicts a noiseless version of the random process from which grid `g` was
#' sampled, conditional on the observed data, and possibly a set of covariates. It is
#' optimal in the sense of minimizing mean squared prediction error under the covariance
#' model specified by `pars`, and assuming the predictions are a linear combination of
#' the observed data.
#'
#' The estimation method is determined by `X`. Set this to `0` and supply a de-trended
#' `g` to do simple kriging. Set it to `NA` to estimate a spatially uniform mean (ordinary
#' kriging). Or pass covariates to `X`, either as multi-layer sk grid or matrix, to do
#' universal kriging. See `sk_GLS` for details on specifying `X` in this case.
#'
#' Set `what='v'` to return the point-wise kriging variance. This usually takes much
#' longer to evaluate than the prediction, but the computer memory demands are similar.
#' A progress bar will be printed to console in this case unless `quiet=TRUE`.
#'
#'
#' Technical notes
#'
#' All calculations returned by `sk_cmean` are exact. Our implementation is based on the
#' variance decomposition suggested in section 3.4 (p. 153-155) of Cressie (1993), and uses a
#' loop over eigen-vectors (for observed locations) to compute variance iteratively.
#'
#' In technical terms, `sk_cmean` estimates the mean of the signal process behind the data.
#' The nugget effect (`eps`) is therefore added to the diagonals of the covariance matrix for
#' the observed points, but NOT to the corresponding entries of the cross-covariance matrix.
#' This has the effect of smoothing (de-noising) predictions at observed locations, which means
#' `sk_cmean` is not an exact interpolator (except in the limit `eps` -> `0`). Rather it makes
#' a correction to the observed data to make it consistent with the surrounding signal.
#'
#' This is a good thing - real spatial datasets are almost always noisy, and we are typically
#' interested in the signal, not some distorted version of it. For more on this see section 3
#' of Cressie (1993), and in particular the discussion in 3.2.1 on the nugget effect.
#'
#' The covariance factorization `fac` can be pre-computed using `sk_var` with arguments
#' `scaled=TRUE` (and, if computing variance, `fac_method='eigen'`). This will speed up
#' subsequent calls where only the observed data values have changed (same covariance structure
#' `pars`, and same `NA` structure in the data). The kriging variance does not change
#' in this case and only needs to be computed once.
#'
#' reference: "Statistics for Spatial Data" by Noel Cressie (1993)
#'
#' @param g an sk grid or list accepted by `sk` (with entries 'gdim', 'gres', 'gval')
#' @param pars list of form returned by `sk_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param X sk grid, numeric, vector, matrix, or NA: the mean, or its linear predictors
#' @param what character, what to compute: one of 'p' (predictor), 'v' (variance), or 'm' (more)
#' @param out character, the return object, either 's' (sk grid) or 'v' (vector)
#' @param fac (optional) pre-computed factorization of covariance matrix scaled by partial sill
#' @param quiet logical indicating to suppress console output
#' @param fac_method character, either 'chol' or 'eigen'
#'
#' @return numeric matrix, the predicted values (or their variance)
#'
#' @export
#' @seealso sk sk_pars
#' @family estimators
#' @family variance-related functions
#'
#' @examples
#'
#' ## set up very simple example problem
#'
#' # make example grid and covariance parameters
#' g_all = sk_sim(10)
#' pars = sk_pars(g_all)
#' plot(g_all)
#'
#' # remove most of the data
#' n = length(g_all)
#' p_miss = 0.90
#' is_miss = seq(n) %in% sample.int(n, round(p_miss*n))
#' is_obs = !is_miss
#' g_miss = g_all
#' g_miss[is_miss] = NA
#' plot(g_miss)
#'
#'
#' ## simple kriging
#'
#' # estimate the missing data conditional on what's left over
#' g_simple = sk_cmean(g_miss, pars, X=0)
#' plot(g_simple)
#'
#' # variance of the estimator
#' g_simple_v = sk_cmean(g_miss, pars, X=0, what='v', quiet=TRUE)
#' plot(g_simple_v)
#'
#' # get the same results with pre-computed variance
#' var_pc = sk_var(g_miss, pars, scaled=TRUE, fac_method='eigen')
#' g_simple_v_compare = sk_cmean(g_miss, pars, X=0, what='v', fac=var_pc, quiet=TRUE)
#' max(abs(g_simple_v_compare - g_simple_v))
#'
#' ## ordinary kriging
#'
#' # estimate spatially uniform mean - true value is 0
#' sk_GLS(g_miss, pars, out='b')
#'
#' # ordinary kriging automatically adjusts for the trend
#' g_ok = sk_cmean(g_miss, pars, X=NA)
#'
#' # additional uncertainty in estimation means variance goes up a bit
#' g_ok_v = sk_cmean(g_miss, pars, X=NA, what='v', quiet=TRUE)
#' range(g_ok_v - g_simple_v)
#'
#'
#' ## universal kriging
#'
#' # introduce some covariates
#' n_betas = 3
#' betas = stats::rnorm(n_betas, 0, 10)
#' g_X = sk_sim(g_all, pars, n_layer=n_betas-1L)
#' g_lm_all = g_all + as.vector(cbind(1, g_X[]) %*% betas)
#' g_lm_miss = g_lm_all
#' g_lm_miss[is_miss] = NA
#'
#' # prediction
#' g_uk = sk_cmean(g_lm_miss, pars, g_X)
#' g_uk_v = sk_cmean(g_lm_miss, pars, g_X, what='v', quiet=TRUE)
#'
#'
#' ## repeat with special subgrid case (faster!)
#'
#' # make g_all a subgrid of a larger example
#' g_super = sk_rescale(g_all, down=2)
#'
#' # re-generate the covariates for the larger extent
#' g_X_super = sk_sim(g_super, pars, n_layer=n_betas-1L)
#' g_lm_super = g_super + as.vector(cbind(1, g_X_super[]) %*% betas)
#'
#' # prediction
#' g_super_uk = sk_cmean(g_lm_super, pars, g_X_super)
#' g_super_uk_v = sk_cmean(g_lm_super, pars, g_X_super, what='v', quiet=TRUE)
#'
#'
#' ## verification
#'
#' # get observed variance and its inverse
#' V_full = sk_var(g_all, pars)
#' is_obs = !is.na(g_miss)
#' Vo = V_full[is_obs, is_obs]
#' Vo_inv = solve(Vo)
#'
#' # get cross covariance
#' is_diag = as.logical(diag(nrow=length(g_all))[is_obs,])
#' Vc = V_full[is_obs,]
#'
#' # nugget adjustment to diagonals (corrects the output at observed locations)
#' Vc[is_diag] = Vc[is_diag] - pars[['eps']]
#'
#' # get covariates matrix and append intercept column
#' X = g_X[]
#' X_all = cbind(1, X)
#' X_obs = X_all[is_obs,]
#'
#' # simple kriging the hard way
#' z = g_miss[is_obs]
#' z_trans = Vo_inv %*% z
#' z_pred_simple = c(t(Vc) %*% z_trans)
#' z_var_simple = pars[['psill']] - diag( t(Vc) %*% Vo_inv %*% Vc )
#'
#' # ordinary kriging the hard way
#' x_trans = ( 1 - colSums( Vo_inv %*% Vc ) )
#' m_ok = x_trans / sum(Vo_inv)
#' z_pred_ok = c( (t(Vc) + m_ok) %*% z_trans )
#' z_var_ok = z_var_simple + diag( x_trans %*% t(m_ok) )
#'
#' # universal kriging the hard way
#' z_lm = g_lm_miss[is_obs]
#' z_lm_trans = Vo_inv %*% z_lm
#' x_lm_trans = X_all - t( t(X_obs) %*% Vo_inv %*% Vc )
#' m_uk = x_lm_trans %*% solve(t(X_obs) %*% Vo_inv %*% X_obs)
#' z_pred_uk = c( (t(Vc) + t(X_obs %*% t(m_uk)) ) %*% z_lm_trans )
#' z_var_uk = z_var_simple + diag( x_lm_trans %*% t(m_uk) )
#'
#' # check that these results agree with sk_cmean
#' max(abs(z_pred_simple - g_simple), na.rm=TRUE)
#' max(abs(z_var_simple - g_simple_v))
#' max(abs(z_pred_ok - g_ok), na.rm=TRUE)
#' max(abs(z_var_ok - g_ok_v))
#' max(abs(z_pred_uk - g_uk), na.rm=TRUE)
#' max(abs(z_var_uk - g_uk_v))
#'
#'
#' ## repeat verification for sub-grid case
#'
#' # rebuild matrices
#' V_full = sk_var(g_super_uk, pars)
#' is_obs = !is.na(g_super)
#' Vo = V_full[is_obs, is_obs]
#' Vo_inv = solve(Vo)
#' is_diag = as.logical(diag(nrow=length(g_super))[is_obs,])
#' Vc = V_full[is_obs,]
#' Vc[is_diag] = Vc[is_diag] - pars[['eps']]
#' X = g_X_super[]
#' X_all = cbind(1, X)
#' X_obs = X_all[is_obs,]
#' z = g_miss[is_obs]
#'
#' # universal kriging the hard way
#' z_var_simple = pars[['psill']] - diag( t(Vc) %*% Vo_inv %*% Vc )
#' z_lm = g_lm_super[is_obs]
#' z_lm_trans = Vo_inv %*% z_lm
#' x_lm_trans = X_all - t( t(X_obs) %*% Vo_inv %*% Vc )
#' m_uk = x_lm_trans %*% solve(t(X_obs) %*% Vo_inv %*% X_obs)
#' z_pred_uk = c( (t(Vc) + t(X_obs %*% t(m_uk)) ) %*% z_lm_trans )
#' z_var_uk = z_var_simple + diag( x_lm_trans %*% t(m_uk) )
#'
#' # verification
#' max(abs(z_pred_uk - g_super_uk), na.rm=TRUE)
#' max(abs(z_var_uk - g_super_uk_v))
#'
sk_cmean = function(g, pars, X=NA, what='p', out='s', fac_method='chol', fac=NULL, quiet=FALSE)
{
  ## validate inputs

  # determine output type - sk grid or vector, prediction or variance)
  is_sk = startsWith(out, 's')
  is_var = startsWith(what, 'v')

  # variance mode forces eigen-decomposition
  if( startsWith(what, 'v') ) fac_method = 'eigen'

  # check for expected objects in list in g
  nm_expect = c('gdim', 'gres', 'gval')
  msg_expect = paste('g must be a list with entries', paste(nm_expect, collapse=', '))
  if( !is.list(g) | !all( nm_expect %in% names(g) ) ) stop(msg_expect)

  # copy grid info but not values
  g_out = g[c('gdim', 'gres', 'crs', 'gyx')]
  gdim_out = g[['gdim']]
  gres_out = g[['gres']]

  # TODO add multi-layer support

  # identify observed data points and copy their index
  is_obs = is_obs_src = if(is.null(g[['is_obs']])) !is.na(g[['gval']]) else g[['is_obs']]
  idx_obs = which(is_obs)
  n_out = prod(gdim_out)

  # copy non-NA data
  z_obs = g[['gval']][is_obs]
  if( is.null(z_obs) ) stop('No non-NA values found in g')
  n_obs = length(z_obs)

  # extract covariates matrix from sk grid input
  if(inherits(X, 'sk'))  X = matrix(X[], nrow=length(X))

  ## variance calculations

  # Check if estimating a deterministic trend, initialize the trend vector
  use_GLS = anyNA(X) | is.matrix(X)
  if( !is.matrix(X) ) mu_GLS = X

  # construct first rows of the symmetric Toeplitz correlation matrices for y and x
  cy = sqrt(pars[['psill']]) * sk_corr_mat(pars[['y']], gdim_out[['y']], gres_out[['y']], i=1)
  cx = sqrt(pars[['psill']]) * sk_corr_mat(pars[['x']], gdim_out[['x']], gres_out[['x']], i=1)

  # check for sub-grid structure and substitute simpler equivalent problem if possible
  sg = sk_sub_find(is_obs, gdim_out)
  is_subgrid = !is.null(sg)
  if(is_subgrid)
  {
    # console message indicating sub-grid size
    msg_dim = paste(sg[['gdim']], collapse=' x ')
    if(!quiet) cat(paste(msg_dim, 'complete sub-grid detected\n'))

    # overwrite g with smaller sub-grid layout (original dimensions etc are already copied)
    fac_method = 'eigen'
    is_sep = TRUE
    g = list(gval=z_obs, gdim=sg[['gdim']], gres=gres_out*sg[['res_scale']])

    # original observed index stored in is_obs_src
    is_obs = rep(TRUE, n_obs)
  }

  # compute factorization (scaled=TRUE means partial sill is factored out)
  if( is.null(fac) ) fac = sk_var(g, pars, scaled=TRUE, fac_method=fac_method, sep=is_subgrid)

  ## kriging predictor

  # de-noise observed data by left-multiplying inverse covariance
  z_tilde = sk_var_mult(z_obs, pars, fac=fac)

  # left-multiply by cross-covariance to get simple kriging predictor
  z_p_tilde = as.numeric(sk_toep_mult(cy, z_tilde, cx, idx_obs, gdim_out))
  if( !( use_GLS | is_var) )
  {
    # finished with simple kriging predictor - return as vector
    z_out = as.numeric(X) + z_p_tilde
    if(!is_sk) return(z_out)

    # or return in an sk grid object
    return(sk(utils::modifyList(g_out, list(gval=z_out))))
  }

  # compute GLS coefficients and resulting adjustments to predictor
  if(use_GLS)
  {
    # copy covariates
    if( anyNA(X) )
    {
      # any NAs in covariates matrix leads to ordinary kriging
      X = NULL
      X_obs = NA

    } else {

      # copy only subset at observed points
      X_obs = matrix(X[is_obs_src,], n_obs)
    }

    # find betas, predictor, and the data matrix with an intercept column
    gls_result = sk_GLS(g, pars, X=X_obs, fac=fac, fac_method=fac_method, out='a')
    fac_X = gls_result[['fac_X']]

    # compute bias adjustment due to estimation of linear predictor
    X_p_tilde = apply(sk_var_mult(gls_result[['x']], pars, fac=fac),
                      MARGIN=2,
                      function(x) sk_toep_mult(cy, x, cx, idx_obs, gdim_out))

    # adjusted covariates used here and in first variance loop below
    X_adj = cbind(rep(1, n_out), X) - X_p_tilde

    # estimate trend
    mu_GLS = as.numeric( tcrossprod(X_adj, t(gls_result[['b']])) )
  }

  # finished with universal and ordinary kriging predictor
  if( !is_var )
  {
    # as vector
    z_out = mu_GLS + z_p_tilde
    if(!is_sk) return(z_out)

    # or return in an sk grid object
    return(sk(utils::modifyList(g_out, list(gval=z_out))))
  }

  # verify that fac is an eigen-decomposition
  msg_fac = 'Cholesky factorization (fac) not supported in what="v" mode'
  if( !is.list(fac) ) stop(msg_fac)
  if(is_subgrid)
  {
    # check that the correct factorization was supplied (component-wise, in a list)
    if( !all( c('y', 'x') %in% names(fac) ) ) stop('supplied factorization was not separable')
    if( !all( sapply(fac[c('y', 'x')], is.list) ) ) stop(msg_fac)
  }

  # initialize variance contribution from GLS (0 in simple kriging case)
  v_gls = numeric(n_out)
  if(use_GLS)
  {
    # loop over eigen-values of X^T V^-1 X
    for(idx in seq_along(gls_result[['fac_X']][['values']]))
    {
      # copy eigen-values of inverse scaled by psill
      ev = 1 / ( pars[['psill']] * gls_result[['fac_X']][['values']][idx] )

      # multiply by adjusted covariates matrix and add square to running total in v_gls
      v_gls[] = v_gls[] + ev * tcrossprod(gls_result[['fac_X']][['vectors']][,idx], X_adj)^2
    }
  }

  # prepare a more efficient method for below when the observed points form a sub-grid
  idx_ev = seq(n_obs)
  if(is_subgrid)
  {
    # compute eigenvalues for observed covariance matrix inverse
    ev_corr = kronecker(fac[['x']][['values']], fac[['y']][['values']])
    ev = 1 / ( (pars[['psill']] * ev_corr) + pars[['eps']] )

    # directly build and multiply the relatively small component covariance matrices
    c_cross_y = sk_corr_mat(pars[['y']], gdim_out[['y']], gres_out[['y']], j=sg[['ij']][['y']])
    c_cross_x = sk_corr_mat(pars[['x']], gdim_out[['x']], gres_out[['x']], j=sg[['ij']][['x']])
    add_y2 = ( c_cross_y %*% fac[['y']][['vectors']] )
    add_x2 = ( c_cross_x %*% fac[['x']][['vectors']] )
  }

  # large loop over eigen-values of covariance matrix builds result iteratively
  v_spat = numeric(n_out)
  if(!quiet) pb = utils::txtProgressBar(max=n_obs, style=3)
  for(idx in seq(n_obs))
  {
    # separable case first
    if( is_subgrid )
    {
      # column indices in component correlation matrices corresponding to eigen-value idx
      idx_yx = sk_vec2mat(idx, sg[['gdim']], out='list')

      # Kronecker product of component columns to get large column vector
      add_yx = ( pars[['psill']] * kronecker(add_x2[, idx_yx[['j']]], add_y2[, idx_yx[['i']]]) )^2
      v_add = ev[idx] * add_yx

    } else {

      # eigen-values of inverse scaled by psill
      ev = 1 / ( pars[['psill']] * fac[['values']][idx] )

      # slow multiplication of eigen-vector with cross covariance
      v_add = ev * sk_toep_mult(cy, fac[['vectors']][, idx], cx, idx_obs, gdim_out)^2
    }

    # add result to running total
    v_spat = v_spat + as.vector(v_add)
    if(!quiet) utils::setTxtProgressBar(pb, idx)
  }
  if(!quiet) close(pb)

  # finished with kriging variance
  if(is_var)
  {
    # notice nugget effect not added to noiseless prediction
    z_out = pars[['psill']] + as.numeric(v_gls) - v_spat
    if(!is_sk) return(z_out)

    # return the sk grid
    return(sk(utils::modifyList(g_out, list(gval=z_out))))
  }

  stop('argument "what" not recognized')
}




#' Fit a covariance model to an sk grid by maximum likelihood
#'
#' This uses `stats::optim` to minimize the log-likelihood function for a grid dataset
#' `g` over the space of unknown parameters for the covariance function specified in `pars`.
#' If only one parameter is unknown, the function instead uses `stats::optimize`.
#'
#' `NA` entries in `pars` are treated as unknown parameters and fitted by the
#' function, whereas non-`NA` entries are treated as fixed parameters (and not fitted).
#' If none of the parameters in `pars` are `NA`, the function copies everything as initial
#' values, then treats all parameters as unknown. `pars` can also be a character vector
#' defining a pair of correlograms (see `?sk_pars`) in which case all covariance parameters
#' are treated as unknown.
#'
#' Bounds and initial values are set automatically using `sk_bds`, unless they are otherwise
#' specified in arguments `lower`, `initial`, `upper`. These should be vectors containing only
#' the unknown parameters, *ie.* they must exclude fixed parameters. Call
#' `sk_update_pars(pars, iso=iso)` to get a template with the expected names and order.
#'
#' All parameters in the covariance models supported by `snapKrig` are strictly positive.
#' Optimization is (by default) done on the parameter log-scale, and users can select a
#' non-constrained `method` if they wish (`?stats::optim`). As the default method 'L-BFGS-B'
#' is the only one that accepts bounds (`lower`, `initial`, `upper` are otherwise ignored)
#' `method` is ignored when `log_scale=FALSE`.
#'
#' Note that the 'gxp' and 'mat' correlograms behave strangely with very small or very large
#' shape parameters, so for them we recommended 'L-BFGS-B' only.
#'
#' When there is only one unknown parameter, the function uses `stats::optimize` instead of
#' `stats::optim`. In this case all entries of `control` with the exception of 'tol' are
#' ignored, as are bounds and initial values, and arguments to `method`.
#'
#' As a sanity check `n_max` sets a maximum for the number of observed grid points. This
#' is to avoid accidental calls with very large datasets that would cause R to hang or crash.
#' Set `n_max=Inf` (with caution) to bypass this check. Similarly the maximum number of
#' iterations is set to `1e3` but this can be changed by manually setting 'maxit' in
#' `control`.
#'
#' @param g an sk grid (or list with entries 'gdim', 'gres', 'gval')
#' @param pars covariance parameter list, with `NA`s indicating parameters to fit
#' @param X numeric (or NA), matrix, or sk grid of linear predictors, passed to `sk_LL`
#' @param iso logical, indicating to constrain the y and x kernel parameters to be the same
#' @param n_max integer, the maximum number of observations allowed
#' @param quiet logical, indicating to suppress console output
#' @param lower numeric vector, lower bounds for parameters
#' @param initial numeric vector, initial values for parameters
#' @param upper numeric vector, upper bounds for parameters
#' @param log_scale logical, indicating to log-transform parameters for optimization
#' @param method character, passed to `stats::optim` (default is 'L-BFGS-B')
#' @param control list, passed to `stats::optim`
#'
#' @return A parameter list in the form returned by `sk_pars` containing both fixed and
#' fitted parameters. The data-frame of bounds and initial values is also included in the
#' attribute 'bds'
#'
#' @export
#' @family parameter managers
#' @seealso sk sk_LL sk_nLL stats::optim stats::optimize
#'
#' @examples
#'
#' # define a grid
#' gdim = c(50, 53)
#' g_empty = sk(gdim)
#' pars_src = sk_pars(g_empty)
#' pars_src = utils::modifyList(pars_src, list(eps=runif(1, 0, 1e1), psill=runif(1, 0, 1e2)))
#' pars_src[['y']][['kp']] = pars_src[['x']][['kp']] = runif(1, 1, 50)
#'
#' # generate example data and fit to it
#' g_obs = sk_sim(g_empty, pars_src)
#' sk_plot(g_obs)
#' fit_result = sk_fit(g_obs, pars='gau')
#'
#' # compare estimates with true values
#' rbind(true=sk_pars_update(pars_src), fitted=sk_pars_update(fit_result))
#'
#' # extract bounds data frame
#' attr(fit_result, 'bds')
#'
#' # check a sequence of other psill values
#' pars_out = fit_result
#' psill_test = ( 2^(seq(5) - 3) ) * pars_out[['psill']]
#' LL_test = sapply(psill_test, function(s) sk_LL(utils::modifyList(pars_out, list(psill=s)), g_obs) )
#' plot(psill_test, LL_test)
#' lines(psill_test, LL_test)
#' print(data.frame(psill=psill_test, likelihood=LL_test))
#'
#' # repeat with most data missing
#' n = prod(gdim)
#' n_obs = 200
#' g_obs = sk_sim(g_empty, pars_src)
#' idx_miss = sample.int(length(g_empty), length(g_empty) - n_obs)
#' g_miss = g_obs
#' g_miss[idx_miss] = NA
#' sk_plot(g_miss)
#'
#' # fit and compare
#' fit_result = sk_fit(g_miss, pars='gau')
#' rbind(true=sk_pars_update(pars_src), fitted=sk_pars_update(fit_result))
#'
sk_fit = function(g, pars=NULL, X=NA, iso=TRUE, n_max=1e3, quiet=FALSE,
                  lower=NULL, initial=NULL, upper=NULL,
                  log_scale=TRUE, method='L-BFGS-B', control=list())
{
  # only L-BFGS-B accepts bounds, so log-scale is mandatory for other methods
  if( (method != 'L-BFGS-B' & !log_scale) )
  {
    warning('Setting log_scale=TRUE (required for methods other than L-BFGS-B)')
    log_scale = TRUE
  }

  # convert to sk object and check for sub-grid
  g_obs = sk(g)
  is_sg = !is.null(sk_sub_find(g_obs))

  # unpack grid input and count observations
  n_all = length(g_obs)
  is_obs = !is.na(g_obs)
  n_obs = sum(is_obs)
  if(n_obs == 0) stop('no data found in g_obs')

  # problem size sanity check
  if( !is_sg & (n_obs > n_max) ) stop('number of observed points exceeded n_max')

  # parse character string arguments to pars
  pars = ifelse(is.null(pars), 'gau', pars)
  if( !is.list(pars) ) pars = sk_pars(g_obs, pars)

  ## set up bounds and initial values

  # standardize input to pars and set NAs for missing values
  pars_fix = sk_pars_make(pars)

  # extract fixed parameter values and names of unknown parameters
  pars_fix_vec = sk_pars_update(pars_fix, iso=iso)
  nm_fit = names(pars_fix_vec)[ is.na(pars_fix_vec) ]

  # if all parameters are supplied we fix none of them (and fit all)
  if( length(nm_fit) == 0 )
  {
    # use supplied parameters as initial values
    if(is.null(initial)) initial = pars_fix_vec

    # copy NAs to all entries of fixed parameters list and update names
    pars_fix = sk_pars_update(pars_fix, rep(NA, length(pars_fix_vec)), iso=iso)
    pars_fix_vec = sk_pars_update(pars_fix, iso=iso)
    nm_fit = names(pars_fix_vec)
  }

  # get default initial values and bounds data frame
  pars_df = sk_bds(pars_fix, g_obs)[nm_fit,]
  nm_expect = rownames(pars_df)

  # collect any user supplied bounds and replace defaults
  bds_arg = list(lower=lower, initial=initial, upper=upper)
  for(arg_nm in names(bds_arg))
  {
    # skip when bounds argument not supplied
    arg = bds_arg[[arg_nm]]
    if( !is.null(arg) )
    {
      # check for valid names then replace only the named arguments
      p_nm = names(arg)
      is_nm_valid = all(p_nm %in% nm_expect)
      if( !is_nm_valid | is.null(p_nm) ) names(arg) = nm_expect[seq_along(arg)]
      pars_df[names(arg), arg_nm] = arg
    }
  }

  # switch to log-scale if requested
  if(log_scale)
  {
    # transform bounds and fixed parameters (if any)
    pars_df = log(pars_df)
    pars_fix_vec = log(pars_fix_vec)
  }

  # scaling constants passed to optimizer for internal use
  optim_scale = apply(pars_df, 1L, function(p) diff(range(p)))

  # when iso=TRUE and pars_fix contains fixed y kernel parameter(s) they must be copied to x
  pars_fix = sk_pars_update(pars_fix, pars_fix_vec, iso=iso)

  # evaluate objective at initial values and bounds
  bds_nLL = apply(pars_df, 2L, function(p) sk_nLL(p, g_obs, pars_fix, X, iso, TRUE, log_scale))

  # 1d optimization case
  if( nrow(pars_df) == 1 )
  {
    # copy user supplied tolerance or set default
    tol = ifelse(is.null(control[['tol']]), .Machine$double.eps^0.25, control[['tol']])
    optimize_out = stats::optimize(f = sk_nLL,
                                   interval = pars_df[c('lower', 'upper')],
                                   tol = tol,
                                   g_obs = g_obs,
                                   pars_fix = pars_fix,
                                   X = X,
                                   iso = iso,
                                   log_scale = log_scale,
                                   quiet = quiet)

    # reshape the output list to be like output of optim
    optim_out = list(message='', par=optimize_out[['minimum']], value=optimize_out[['objective']])

  } else {

    # n-d optimization case

    # set default control parameters and run the optimizer
    control = utils::modifyList(list(maxit=1e3, parscale=optim_scale), control)
    optim_out = suppressWarnings(stats::optim(par = pars_df[['initial']],
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
                                              control = control) )
  }

  # unpack optimizer output
  obj_val = optim_out[['value']]
  pars_fitted_v = optim_out[['par']]
  if(!quiet) cat(paste('\n', optim_out[['message']]))

  # revert to initial value if optimize/optim result was not an improvement
  if(bds_nLL[['initial']] < obj_val)
  {
    pars_fitted_v = pars_df[['initial']]
    warning('initial values had best likelihood score)')
  }

  # reshape as list and copy results to bounds data frame
  pars_fitted = sk_pars_update(pars_fix, pars_fitted_v, iso=iso, na_omit=TRUE)
  pars_df['fitted'] = pars_fitted_v

  # back transform parameters fitted on log-scale
  if(log_scale)
  {
    pars_df = exp(pars_df)
    pars_fitted = sk_pars_update(pars_fitted, exp(sk_pars_update(pars_fitted)))
  }

  # add bounds data frame as attribute of output parameter list
  df_order = c('lower', 'initial', 'fitted', 'upper')
  attr(pars_fitted, 'bds') = pars_df[df_order]
  return(pars_fitted)
}



