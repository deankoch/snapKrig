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
#' `sk_cmean` is not an exact interpolator (unless `eps=0`). Rather it makes a correction to
#' the observed data to make it consistent with the surrounding signal.
#'
#' This is a good thing - real spatial datasets are almost always noisy, and we are usually
#' interested in the signal, and not some version of it that is distorted by error. For more
#' on this see section 3 of Cressie (1993), and in particular the discussion in 3.2.1 on the
#' nugget effect.
#'
#' The covariance factorization `fac` can be pre-computed using `sk_var` with arguments
#' `scaled=TRUE` and (if computing variance) `fac_method='eigen'`. This will speed up repeated
#' calls where only the observed data values change (same covariance structure `pars`, and same
#' `NA` structure in the data). Note that the kriging variance does not change in this case and
#' only needs to be computed once.
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
#'
#' @return numeric matrix, the predicted values (or their variance)
#' @export
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
#' betas = rnorm(n_betas, 0, 10)
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
#' ## repeat verification for subgrid case
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
  g_out = g[c('gdim', 'gres', 'crs')]
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
    return(sk(modifyList(g_out, list(gval=z_out))))
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
    return(sk(modifyList(g_out, list(gval=z_out))))
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
    return(sk(modifyList(g_out, list(gval=z_out))))
  }

  stop('argument "what" not recognized')
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











