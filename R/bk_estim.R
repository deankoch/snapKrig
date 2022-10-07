# bk_estim.R
# Dean Koch, 2022
# Functions for parameter inference and spatial prediction

#' Generalized least squares (GLS) estimator
#'
#' Computes coefficients b of the linear predictor E(Z) = Xb using the GLS equation
#'
#' The GLS equation is: b = ( X^T V^{-1} X )^{-1} X^T V^{-1} z,
#'
#' where V is the covariance matrix for data z, and X is a matrix of covariates.
#' V is generated from the covariance model specified in `pars` on the grid `g_obs`,
#' where the non-NA values of `g_obs$gval` form the vector z. Operations with V^{-1}
#' are computed using the factorization `fac`, or else as specified in `fac_method`.
#'
#' When `out='z'`, the function returns the product Xb instead of b.
#'
#' Argument `X` should provide matrix X without the intercept column (a vector of 1's).
#' DO NOT include an intercept column in argument `X` or you will get collinearity errors
#' (the function appends it without checking if its there already). The columns of `X` must
#' be independent, and its rows should match the entries of `g_obs$gval`, or its non-NA data
#' values, in order. Use `X=NA` to specify an intercept-only model; ie to fit a spatially
#' constant mean. This replaces X in the GLS equation by a vector of 1's.
#'
#' `g_obs$gval` can be a matrix whose columns are multiple repetitions (layers) of the
#' same spatial process (see `bk`), in which case the covariates in `X` are recycled
#' for each layer. Layers are assumed mutually independent and the GLS equation is evaluated
#' using the corresponding block-diagonal V. Note that this is equivalent to (but faster
#' than) calling `bk_GLS` separately on each layer with the same `X` and averaging the
#' estimated b's.
#'
#' @param g_obs list of form returned by `bk` (with entries 'gdim', 'gres', 'gval')
#' @param pars list of form returned by `bk_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param X matrix or NA, the linear predictors (in columns) excluding intercept
#' @param fac_method character, factorization method: 'eigen' (default) or 'chol' (see `bk_var`)
#' @param fac matrix or list, (optional) pre-computed covariance matrix factorization
#' @param out character, either 'b' (coefficients) or 'z' (linear predictor)
#'
#' @return numeric vector, either b (default) or Xb (if `out='z'`)
#' @export
#'
#' @examples
#' # set up example grid, and covariance parameters
#' gdim = c(45, 31)
#' n = prod(gdim)
#' g = bk(gdim)
#' pars = modifyList(bk_pars(g, 'gau'), list(psill=2))
#'
#' # generate spatial noise
#' z = bk_sim(g, pars)
#' bk_plot(modifyList(g, list(gval=z)))
#'
#' # generate some covariates and data
#' n_betas = 3
#' betas = rnorm(n_betas, 0, 10)
#' X_all = cbind(1, matrix(rnorm(n*(n_betas-1)), n))
#' lm_actual = as.vector(X_all %*% betas)
#' g_obs = modifyList(g, list(gval=z+lm_actual))
#'
#' # exclude intercept column in calls to bk_GLS
#' X_pass = X_all[,-1]
#'
#' # find the GLS coefficients
#' betas_est = bk_GLS(g_obs, pars, X_pass)
#' print(betas_est)
#' print(betas)
#'
#' # compute trend as product of betas with matrix X_all, or by setting out='z'
#' lm_est = X_all %*% betas_est
#' max( abs( bk_GLS(g_obs, pars, X_pass, out='z') - lm_est ) )
#'
#' # repeat with pre-computed eigen factorization
#' fac_eigen = bk_var(g_obs, pars, fac_method='eigen', sep=TRUE)
#' betas_est_compare_eigen = bk_GLS(g_obs, pars, X_pass, fac=fac_eigen)
#' max( abs( betas_est_compare_eigen - betas_est ) )
#'
#' # missing data example
#' n_obs = 10
#' idx_rem = sort(sample.int(n, n-n_obs))
#' g_miss = g_obs
#' g_miss$gval[idx_rem] = NA
#' bk_plot(g_miss)
#' betas_est = bk_GLS(g_miss, pars, X_pass)
#' print(betas_est)
#' print(betas)
#'
#' # set X to NA to estimate the a spatially constant trend (the adjusted mean)
#' bk_GLS(g_miss, pars, X=NA)
#' mean(g_miss$gval, na.rm=TRUE)
#'
#' # generate some extra layers
#' z_extra = lapply(seq(9), function(x) bk_sim(g, pars))
#' z_multi = lm_actual + do.call(cbind, c(list(z), z_extra))
#'
#' # multi-layer example with missing data
#' is_obs = !is.na(g_miss$gval)
#' map_sparse = match(seq(n), which(is_obs))
#' g_sparse = modifyList(g_miss, list(gval=z_multi[is_obs,], idx_grid=map_sparse))
#' betas_sparse = bk_GLS(g_obs=g_sparse, pars, X=X_pass)
#' print(betas_sparse)
#' print(betas)
#'
#' bk_GLS(g_sparse, pars, NA)
#' mean(g_sparse$gval, na.rm=TRUE)
#'
bk_GLS = function(g_obs, pars, X=NA, fac=NULL, fac_method='eigen', out='b')
{
  # multi-layer support
  n_layer = 1L
  is_multi = !is.null(g_obs[['idx_grid']])
  if(is_multi)
  {
    # identify non-NA points
    is_obs = !is.na(g_obs[['idx_grid']])

    # reorder the non-NA data matrix to grid-vectorized order
    reorder_z = g_obs[['idx_grid']][is_obs]
    n_layer = ncol(g_obs[['gval']])

    # matrix(.., ncol) avoids R simplifying to vector in 1 column case
    z = matrix(g_obs[['gval']][reorder_z,], ncol=n_layer)

  } else {

    # copy non-NA data as a 1-column matrix
    is_obs = as.vector(!is.na(g_obs[['gval']]))
    if( sum(is_obs) < 2 ) stop('Not enough non-NA values in g_obs')
    z = matrix(g_obs[['gval']][is_obs], ncol=n_layer)
  }

  # set default factorization method
  is_sep = all(is_obs)
  if( is_sep & (fac_method=='chol') ) stop('eigen method is required for complete grids')

  # compute variance factorization (scaled=TRUE -> partial sill is factored out)
  if( is.null(fac) ) fac = bk_var(g_obs, pars, scaled=TRUE, fac_method=fac_method, sep=is_sep)

  # build matrix of covariate values from intercept column and (optionally) X
  n = length(is_obs)
  n_obs = sum(is_obs)
  if( anyNA(X) )
  {
    # if X is a matrix with NAs, discard it with a warning
    if( is.matrix(X) ) warning('X contained NA value(s). Setting X=NA')
    X = X_obs = matrix(1L, nrow=n_obs)

  } else {

    # check for invalid input to X and append intercept column
    if( !is.matrix(X) ) stop('X must be a matrix of covariate values')
    X = X_obs = cbind(1L, X)

    # check for incorrect length in X
    msg_expected = paste('expected', n, 'or', n_obs, 'but got', nrow(X))
    if( !( nrow(X) %in% c(n, n_obs) ) ) stop(paste('incorrect number of rows in X:', msg_expected))

    # reorder X to match z
    if( nrow(X) > n_obs ) X_obs = matrix(X[is_obs,], ncol=ncol(X))
    if( is_multi ) X_obs = matrix(X_obs[reorder_z,], ncol=ncol(X))
  }

  # find the factorization of quadratic form with X (scaling by V inverse)
  fac_X = bk_var(g_obs, pars, X=X_obs, scaled=TRUE, fac=fac, fac_method='eigen')

  # compute GLS coefficients using whitened observation data
  z_trans = bk_var_mult(z, pars, fac=fac)
  betas_gls = bk_var_mult(t(crossprod(z_trans, X_obs)), pars, fac=fac_X)
  if(is_multi) betas_gls = rowMeans(betas_gls)

  # return betas by default
  if(startsWith(out, 'b')) return(as.numeric(betas_gls))

  # or a list with everything, or else return linear predictor
  z_gls = as.numeric( tcrossprod(X, t(betas_gls)) )
  if(startsWith(out, 'a')) return(list(z=z_gls, b=as.numeric(betas_gls), x=X_obs, fac_X=fac_X))
  return(z_gls)
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
#' g_empty = bk(gdim)
#' pars_src = bk_pars(g_empty)
#' pars_src = modifyList(pars_src, list(eps=runif(1, 0, 1e1), psill=runif(1, 0, 1e2)))
#' pars_src[['y']][['kp']] = pars_src[['x']][['kp']] = runif(1, 1, 50)
#'
#' # generate example data and fit to it
#' gval = bk_sim(g_empty, pars_src, quiet=F)
#' g_obs = modifyList(g_empty, list(gval=gval))
#' bk_plot(g_obs)
#' fit_result = bk_fit(g_obs, pars='gau')
#'
#' fit_result$pars |> bk_pars_update()
#' pars_src |> bk_pars_update()
#'
#' # check sequence of other psill values
#' pars_out = fit_result$pars
#' psill_test = ( 2^(seq(5) - 3) ) * pars_out[['psill']]
#' LL_test = sapply(psill_test, function(s) bk_LL(modifyList(pars_out, list(psill=s)), g_obs) )
#' plot(psill_test, LL_test)
#' lines(psill_test, LL_test)
#' print(data.frame(psill=psill_test, likelihood=LL_test))
#'
#' # repeat with most data missing
#' n = prod(gdim)
#' n_obs = 200
#' gval = bk_sim(g_empty, pars_src, quiet=TRUE)
#' g_obs = modifyList(g_empty, list(gval=gval))
#' idx_obs = sample.int(prod(gdim), n_obs)
#' g_miss = modifyList(g_obs, list(gval=rep(NA, n)))
#' g_miss$gval[idx_obs] = g_obs$gval[idx_obs]
#' bk_plot(g_miss)
#'
#'
bk_fit = function(g_obs, pars=NULL, X=NA, iso=TRUE, initial=NULL, quiet=FALSE,
                     lower=NULL, upper=NULL, n_max=1e3)
{
  # unpack vectorized grid as list
  g_obs = bk(g_obs)
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
  #   sub_result = bk_sub_find(g_obs=g)
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
  if( !is.list(pars) ) pars = bk_pars(g_obs, ifelse(is.null(pars), 'gau', pars))
  p_fixed = bk_pars_update(pars, iso=iso)
  is_fitted = is.na(p_fixed)
  if( !any(is_fitted) ) is_fitted[] = TRUE

  # set initial value defaults
  nm_fitted = names(is_fitted)[is_fitted]
  nm_fixed = names(is_fitted)[!is_fitted]
  if( is.null(initial) ) initial = bk_bds(pars, g_obs)[nm_fitted, 'initial']
  pars = bk_pars_update(pars, p_fixed, iso=iso)

  # fit the model
  #v = var(g[['gval']], na.rm=TRUE)
  result_optim = bk_optim(g_obs, pars, X=X, iso=iso, quiet=quiet,
                             lower=lower, initial=initial, upper=upper)
  pars_out = result_optim[['pars']]

  # # check sequence of likely psill substitutions
  # psill_test = ( 2^(seq(5) - 3) ) * pars_out[['psill']]
  # pars_test = lapply(psill_test, function(s) modifyList(pars_out, list(psill=s)) )
  # LL_test = sapply(pars_test, function(p) bk_LL(p, g, X=X) )
  # pars_out[['psill']] = psill_test[ which.max(LL_test) ]

  # de-trend the data by subtracting GLS estimate
  #z_gls = 0
  #if( anyNA(X) | is.matrix(X) ) z_gls = bk_GLS(g, pars_out, X=X, out='z')
  #g[['gval']][is_obs] = g[['gval']][is_obs] - z_gls

  # plot the semi-variogram for de-trended data
  #vg_detrend = bk_sample_vg(g)
  #bk_plot_semi(vg_detrend, pars_out)
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
  pars = bk_pars_update(pars, p_fixed, iso=TRUE)

  # fit the model
  if( anyNA(X) ) X = NA
  result_optim = bk_optim(g, pars, X=X, iso=TRUE, initial=initial)
  pars_out = result_optim[['pars']]

  # check sequence of likely psill values
  psill_test = ( 2^(seq(5) - 3) ) * pars_out[['psill']]
  pars_test = lapply(psill_test, function(s) modifyList(pars_out, list(psill=s)) )
  LL_test = sapply(pars_test, function(p) bk_LL(p, g, X=X) )
  #print(z_scale)
  #print(data.frame(psill=psill_test, likelihood=LL_test))
  pars_out[['psill']] = psill_test[ which.max(LL_test) ]

  # de-trend the scaled data
  z_gls = bk_GLS(g, pars_out, X=X, out='z')
  if(anyNA(X)) z_gls = z_gls[1]
  g_out = modifyList(g, list(gval = z_scale * (z_std-z_gls) ) )

  # transform scaled variance parameters back to scale of input data
  v_unscaled = list(eps = 2*z_scale^2 * pars_out[['eps']], psill = 2*z_scale^2 * pars_out[['psill']])
  pars_out = modifyList(pars_out, v_unscaled)

  # plot the semi-variogram on original scale
  vg_out = bk_sample_vg(g_out)
  bk_plot_semi(vg_out, pars_out)
  return(modifyList(result_optim, list(pars=pars_out)))
  }
}



#' Fit covariance parameters to data by maximum (profile) likelihood using optim
#'
#' documentation unfinished
#'
#' @param g_obs list of form returned by `bk` (with entries 'gdim', 'gres', 'gval')
#' @param pars list of fixed kernel parameters, with NAs indicating parameters to fit
#' @param X numeric, vector, matrix, or NA, the mean or its linear predictors, passed to `bk_LL`
#' @param iso logical, indicating to constrain the y and x kernel parameters to be the same
#' @param control list, passed to `stats::optim`
#' @param quiet logical, indicating to suppress console output
#' @param ... named arguments to pass to `bk_bds`
#'
#' @return sdfsdfsd
#' @export
#'
#' @examples
#' # set up example grid and data
#' g_obs = bk(10)
#' g_obs$gval = rnorm(10^2)
#' bk_optim(g_obs, quiet=TRUE)
#'
#' # repeat with one or more parameters fixed
#' pars = bk_pars_make('gau') # NA parameters list
#' pars$psill = 1
#' bk_optim(g_obs, pars, quiet=TRUE)
#' pars$y$kp = 1
#' bk_optim(g_obs, pars, quiet=TRUE)
#'
#' # iso mode constrains x parameters to equal y parameters
#' bk_optim(g_obs, iso=T, quiet=TRUE)
#' bk_optim(g_obs, pars, iso=T, quiet=TRUE)
#'
bk_optim = function(g_obs, pars='gau', X=0, iso=FALSE, control=list(), quiet=FALSE,
                       log_scale=TRUE, method='L-BFGS-B', lower=NULL, initial=NULL, upper=NULL)
{
  # only L-BFGS-B accepts bounds, so log-scale is mandatory for other methods
  if( (method != 'L-BFGS-B' & !log_scale) )
  {
    warning('Setting log_scale=TRUE (required for methods other than L-BFGS-B)')
    log_scale = TRUE
  }

  g_obs = bk(g_obs)

  # standardize input to pars and set NAs for missing values
  pars_fix = bk_pars_make(pars)

  # extract parameter names and NA structure supplied in pars_fix
  pars_fix_vec = bk_pars_update(pars_fix, iso=iso)
  nm_fit = names(pars_fix_vec)[ is.na(pars_fix_vec) ]

  # TODO: when values are specified for all parameters in `pars`, use them as initial values
  if( length(nm_fit) == 0 )
  {
    # set initial value automatically here?
    #initial = bk_pars_update(pars_fix)
    pars_fix = bk_pars_update(pars_fix, rep(NA, length(pars_fix_vec)), iso=iso)
    pars_fix_vec = bk_pars_update(pars_fix, iso=iso)
    nm_fit = names(pars_fix_vec)
  }

  # get default initial values and bounds data frame
  pars_df = bk_bds(pars_fix, g_obs)[nm_fit,]

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
  pars_fix = bk_pars_update(pars_fix, pars_fix_vec, iso=iso)

  # evaluate objective at initial values and bounds
  bds_nLL = apply(pars_df, 2L, function(p) {
    bk_nLL(p, g_obs, pars_fix, X, iso, quiet=TRUE, log_scale) })

  # 1d optimization
  if( nrow(pars_df) == 1 )
  {
    tol = control[['tol']]
    optimize_out = optimize(f = bk_nLL,
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
                             f = bk_nLL,
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
  pars_fitted = bk_pars_update(pars_fix, pars_fitted_v, iso=iso, na_omit=TRUE)
  pars_df['fitted'] = pars_fitted_v
  if(log_scale)
  {
    pars_df = exp(pars_df)
    pars_fitted = bk_pars_update(pars_fitted, exp(bk_pars_update(pars_fitted)))
  }

  # return parameters list and data frame in a list
  df_order = c('lower', 'initial', 'fitted', 'upper')
  return(list(pars=pars_fitted, df=pars_df[df_order], obj=obj_val))
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
#' The covariance factorization `fac` can be pre-computed using `bk_var(..., scaled=TRUE)`
#' to speed up repeated calls where only the observed data values change (ie same covariance
#' structure `pars`, and same NA structure in the data). Note that the kriging variance does
#' not change in this case and only needs to be computed once.
#'
#' @param g_obs list of form returned by `bk` (with entries 'gdim', 'gres', 'gval')
#' @param pars list of form returned by `bk_pars` (with entries 'y', 'x', 'eps', 'psill')
#' @param X numeric, vector, matrix, or NA: the mean, or its linear predictors
#' @param out character, the return value, one of 'predictor', 'variance', or 'm'
#' @param fac (optional) pre-computed factorization of covariance matrix scaled by partial sill
#' @param quiet logical indicating to suppress console output
#'
#' @return numeric matrix, the predicted values (or their variance)
#' @export
#'
#' @examples
#' # make example grid and data
#' n = 25^2
#' n_obs = 10
#' g_obs = bk(sqrt(n))
#' idx_obs = sample.int(n, n_obs)
#' g_obs$gval[idx_obs] = rnorm(n_obs)
#' pars = bk_pars('gau', g_obs)
#' g_pred = bk_cmean(g_obs, pars)
#' g_var = bk_cmean(g_obs, pars, makev=TRUE, quiet=TRUE)
#' #g_obs |> bk_plot()
#' #g_obs |> modifyList(list(gval=g_pred)) |> bk_plot()
#' #g_obs |> modifyList(list(gval=g_var)) |> bk_plot()
#'
bk_cmean = function(g_obs, pars, X=NA, fac=NULL, out='p', fac_method='chol', quiet=FALSE)
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
  cy = cy_obs = sqrt(pars[['psill']]) * bk_corr_mat(pars[['y']], gdim[['y']], gres[['y']], i=1)
  cx = cx_obs = sqrt(pars[['psill']]) * bk_corr_mat(pars[['x']], gdim[['x']], gres[['x']], i=1)

  # set up factorization when it's not supplied. Variance mode forces eigen-decomposition
  if( startsWith(out, 'v') ) fac_method = 'eigen'

  # check for completeness and separability
  is_sep = n == n_obs

  # check for sub-grid structure and substitute simpler equivalent problem if possible
  sg = bk_sub_find(is_obs, g_obs[['gdim']])
  is_sg = !is.null(sg)
  if( is_sg )
  {
    # extract sub-grid layout and find separable covariance eigen-decomposition
    fac_method = 'eigen'
    g_obs = list(gval=z, gdim=sg[['gdim']], gres=g_obs[['gres']] * sg[['res_scale']])
    if( is.null(fac) ) fac = bk_var(g_obs, pars, scaled=TRUE, fac_method=fac_method, sep=TRUE)
    cy_obs = cy[ sg[['ij']][['y']] ]
    cx_obs = cx[ sg[['ij']][['x']] ]

    # g_obs should have no missing (NA) data points now
    is_obs = rep(TRUE, n_obs)
    # (is_obs_src still contains original observed index)

  } else {

    # compute factorization (scaled=TRUE means partial sill is factored out)
    if( is.null(fac) ) fac = bk_var(g_obs, pars, scaled=TRUE, fac_method=fac_method, sep=is_sep)
  }

  # transform the observed data by left-multiplying with inverse covariance
  z_tilde = bk_var_mult(g_obs[['gval']][is_obs], pars, fac=fac)

  # left-multiply by cross-covariance to get simple kriging predictor
  z_p_tilde = bk_toep_mult(cy, z_tilde, cx, idx_obs, gdim) |> as.numeric()
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
    gls = bk_GLS(g_obs, pars, X=X_obs, fac=fac, fac_method=fac_method, out='a')
    fac_X = gls[['fac_X']]

    # compute bias adjustment due to estimation of linear predictor
    X_p_tilde = gls[['x']] |>
      bk_var_mult(pars, fac=fac) |>
      apply(2, \(x) bk_toep_mult(cy, x, cx, idx_obs, gdim))

    # uncomment to get exact interpolator (and discontinuities at observations)
    #X_p_tilde[idx_obs,] = gls[['x']]

    # compute trend and 'm'
    X_adj = cbind(rep(1, n), X) - X_p_tilde
    mu_GLS = tcrossprod(X_adj, t(gls[['b']])) |> as.numeric()
    if( startsWith(out, 'm') ) return( as.numeric(bk_var_mult( t(X_adj), pars, fac=fac_X)) )
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
    c_cross_y = bk_corr_mat(pars[['y']], gdim[['y']], gres[['y']], j=sg[['ij']][['y']])
    c_cross_x = bk_corr_mat(pars[['x']], gdim[['x']], gres[['x']], j=sg[['ij']][['x']])
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
      v_add = ev * bk_toep_mult(cy, fac[['vectors']][, idx], cx, idx_obs, gdim)^2

    } else {

      # column indices in component correlation matrices corresponding to eigen-value idx
      idx_yx = bk_vec2mat(idx, sg[['gdim']], out='list')

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











