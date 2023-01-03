# sk_pars.R
# Dean Koch, 2022
# Functions for managing covariance model parameters
#

#' Initialize Kronecker covariance function parameters for a sk grid
#'
#' Returns a parameter list defining the Kronecker covariance model for `g`, with
#' default initial values assigned to parameters based on the grid dimensions and sample
#' variance of observed data.
#'
#' Swap `fill='initial'` with 'lower' and 'upper' to get default lower and upper bounds.
#'
#' @param g list, a sk grid list (or any other object accepted by `sk`)
#' @param pars character or list defining kernels accepted by `sk_pars_make`
#' @param fill character, either 'initial', 'lower' or 'upper'
#'
#' @return a list defining the Kronecker covariance parameters
#'
#' @export
#' @family parameter managers
#' @seealso sk sk_corr
#'
#' @examples
#' sk_pars(g=10)
#' sk_pars(c(10,15))
#' sk_pars(c(10,15), 'mat')
#' sk_pars(c(10,15), 'mat', 'upper')
sk_pars = function(g, pars='gau', fill='initial')
{
  # check for valid input and set names/NAs wherever they are missing
  pars = sk_pars_make(pars)
  g = sk(g)

  # get default values in a data frame
  bds_p = sk_bds(pars, g)

  # return in list form
  return(sk_pars_update(pars, bds_p[[fill]]))
}


#' Set default parameter covariance parameter bounds for a Kronecker covariance model
#'
#' Returns a data-frame of initial values and upper/lower bounds on covariance
#' parameters for the Kronecker covariance model defined by the correlation function
#' names in `pars`.
#'
#' Range parameters (`y.rho` and `x.rho`) are bounded by the shortest and longest
#' inter-point distances along the corresponding dimension (y or x). This is
#' computed by taking the element-wise product of dimensions and resolution, ie
#' `g$gres * dim(g)`. Ranges are initialized to the geometric mean of the upper
#' and lower bounds.
#'
#' Variance bounds are centered around `var_obs`, which by default is set to the sample
#' variance of the data in `g`. `eps` (measurement variance) and `psill` (partial
#' sill) are both initialized to one half of `var_obs`, bounded above by `var_obs`
#' times `var_mult`, and bounded below by a small positive number (`1e-6`). Note that
#' while `eps=0` produces valid models in theory, in practice `eps>0` is often
#' necessary for numerical stability.
#'
#' Shape parameter bounds are hard-coded, and are set conservatively to avoid problems
#' with numerical precision in functions like `exp` and `gamma` when evaluating very
#' large or small distances.
#'
#' @param pars list or character vector of 1-2 kernel names (see `sk_pars`)
#' @param g a sk grid (or any object accepted by `sk`)
#' @param var_obs positive numeric, the sample variance of data `g$gval`
#' @param var_mult numeric > 1, constant to multiply by `var_obs` to get upper bounds
#'
#' @return a data frame of initial values and lower/upper bounds for the parameters in `pars`
#'
#' @export
#' @keywords internal
#' @family parameter managers
#' @seealso sk
#'
#' @examples
#' gdim = c(10, 15)
#' g = sk(gdim)
#' g[] = stats::rnorm(length(g))
#' sk_bds('mat', g)
#'
#' # same result by passing in observed variance
#' sk_bds('mat', g, stats::var(g[]))
#'
#' # a less conservative bound for variance (only eps and psill affected)
#' sk_bds('mat', g, var_mult=1)
sk_bds = function(pars, g, var_obs=NULL, var_mult=2)
{
  # set up hard-coded shape parameter bounds
  kap_bds = list(gxp=list(lower=0.5, initial=1, upper=2),
                 mat=list(lower=1, initial=10, upper=30))

  # check for valid input and set names/NAs wherever they are missing
  pars = sk_pars_make(pars)
  p_vec = sk_pars_update(pars)

  # locate variance and kernel parameter indices in p_vec
  len_yx = sapply(pars[c('y', 'x')], function(p) length(p[['kp']]))
  idx_yx = list(y = 2 + seq(len_yx['y']), x = 2 + len_yx[['y']] + seq(len_yx['x']))
  idx_rho = sapply(idx_yx, function(x) x[1])
  idx_kap = sapply(idx_yx, function(x) x[-1])

  # validate g
  g = sk(g)

  # compute sample variance if not supplied, or set to default 1 when there is no data
  if( is.null(var_obs) ) var_obs = ifelse(is.null(g[['gval']]), 1, stats::var(g[], na.rm=TRUE))
  if( is.na(var_obs) ) var_obs = 1

  # set up bounds for variance components
  bds = as.list(p_vec)
  bds[['psill']] = list(lower=1e-6, initial=var_obs/2, upper=var_obs*var_mult)
  bds[['eps']] = list(lower=1e-6, initial=var_obs/2, upper=var_obs*var_mult)

  # set up bounds for kernel ranges - initial is geometric mean of upper and lower
  bds_kp = Map(function(r,d) list(lower=r, initial=NA, upper=r*d), r=g[['gres']], d=dim(g))
  bds[idx_rho] = lapply(bds_kp, function(p) {
    utils::modifyList(p, list(initial = sqrt(p[['lower']]*p[['upper']])) ) })

  # set up bounds for kernel shape parameters
  has_kap = sapply(idx_kap, function(x) length(x) == 1 )
  nm_kap = c('y', 'x')[has_kap]
  if( any(has_kap) ) for( nm in nm_kap ) bds[[ idx_kap[[nm]] ]] = kap_bds[[ pars[[nm]][['k']] ]]

  # build requested subset of output data frame
  return(do.call(rbind, lapply(bds, data.frame)))
}


#' Build a parameter list defining the 2d spatial Kronecker covariance model
#'
#' Constructs a nested list naming the expected covariance parameters for the supplied
#' correlation function names in `pars`. If `pars` is already such a list, the function
#' checks that parameter names and lengths are valid and returns a copy containing only the
#' expected parameters, properly named, with `NA`s assigned to any missing parameters.
#'
#' `pars` should be a correlation function name accepted by `sk_kp` ('exp', 'gau','sph', or
#' 'gex', 'mat'), or a vector or list of two of them in the order y, x.
#'
#' @param pars character vector of kernel name(s) or list of parameters (see details)
#'
#' @return parameter list containing sub-lists 'y', 'x', and scalars 'psill' and 'eps'
#'
#' @export
#' @keywords internal
#' @family parameter managers
#' @seealso sk_corr
#'
#' @examples
#' # pass a correlation function name to get 2d version with NAs for all parameters
#' sk_pars_make('mat')
#'
#' # pass a vector or list of correlation function names in order y, x (or else specify)
#' sk_pars_make(c('gau', 'mat'))
#' sk_pars_make(list(x='gau', y='mat'))
#'
#' # if the the x definition is missing it is copied from y, and vice versa
#' sk_pars_make(list(k='exp', kp=c(rho=1)))
#'
#' # when unnamed, kernel range and shape parameters are assigned expected names
#' sk_pars_make(list(psill=1, x=list(k='mat', kp=c(1,2))))
#'
#' # incorrectly named parameters raise errors
#' \donttest{
#' sk_pars_make(list(psill=1, x=list(k='exp', kp=c(foo=1))))
#' }
#'
#' # complete parameter definition lists are returned unchanged
#' k_list = list(k='exp', kp=c(rho=1))
#' pars = list(y=k_list, x=k_list, eps=0, psill=1)
#' identical(pars, sk_pars_make(pars))
#'
sk_pars_make = function(pars='gau')
{
  # expected spatial kernel components (as nested list)
  nm_yx = c('y', 'x')
  nm_yx_nested = c('k', 'kp')

  # expected variance parameters
  nm_var = c('eps', 'psill')

  # check for y, x spatial components
  nm_pars = names(pars)
  has_yx = (nm_yx %in% nm_pars) |> stats::setNames(nm_yx)

  # convert character input to kernel definition list
  is_char = pars |> sapply(is.character)
  if( all(is_char) )
  {
    yx_pars = pars[which(is_char)]
    if( all(nm_yx %in% names(pars)) ) yx_pars = pars[nm_yx]
    pars = yx_pars |> lapply(\(nm) list(k=nm)) |> stats::setNames( nm_yx[seq_along(pars)] )
  }

  # check for y, x spatial components
  nm_pars = names(pars)
  has_yx = (nm_yx %in% nm_pars) |> stats::setNames(nm_yx)

  # handle missing y, x entries
  if( sum(has_yx) == 0 )
  {
    # look for generic kernel definition (k, kp) in top level of list
    has_k = (nm_yx_nested %in% nm_pars) |> stats::setNames(nm_yx_nested)
    if( !has_k['k'] ) stop('kernel name not found in pars$k, pars$y or pars$x')

    # recursive call with this kernel definition in sub-list 'y'
    pars_new = utils::modifyList(pars, list(k=NULL, kp=NULL, y=pars[ nm_yx_nested[has_k] ]))
    return( sk_pars_make(pars_new) )
  }

  # if only one spatial kernel is found, copy it to get identical y, x kernels
  if( sum(has_yx) == 1 )
  {
    kernel_clone = pars[nm_yx[has_yx]] |> stats::setNames(nm_yx[!has_yx])
    pars = pars |> utils::modifyList(kernel_clone)
  }

  # checks for required spatial kernel parameters in sub-lists
  pars[nm_yx] = pars[nm_yx] |> lapply(\(p) {

    # check that each dimension has a kernel name defined
    has_k = (nm_yx_nested %in% names(p)) |> stats::setNames(nm_yx_nested)
    if( !has_k['k'] ) stop('kernel name k not found in sub-list')
    msg_kernel = paste(p[['k']], 'kernel')

    # set NA kernel parameters when they are not found
    kp_expect = sk_kp(p[['k']])
    nm_expect = names(kp_expect)
    if( !has_k['kp'] ) p[['kp']] = kp_expect

    # check for mismatch in length
    n_expect = length(kp_expect)
    n_got = length(p[['kp']])
    msg_len = paste(msg_kernel, 'expected', n_expect, 'parameter(s) but got', n_got)
    if(n_got != n_expect) stop(msg_len)

    # set default names if input is unnamed
    if( is.null(names(p[['kp']])) ) names(p[['kp']]) = nm_expect
    nm_input = names(p[['kp']])

    # warn of incorrect kernel name(s)
    msg_expect = paste0(' "', paste(nm_expect, collapse='", "'), '"')
    msg_got = paste0('"', paste(nm_input, collapse='", "'), '"')
    msg_nm = paste(msg_kernel, 'expected parameter(s)', msg_expect, '(in order) but got', msg_got)

    is_ok = (nm_expect == nm_input) | endsWith(nm_input, nm_expect) | startsWith(nm_input, nm_expect)
    if(!all(is_ok)) warning(msg_nm)
    return(p)

    })

  # set other missing parameters to NA
  for(nm in nm_var) { if( !(nm %in% nm_pars) ) pars[[nm]] = NA_real_ }
  return(pars[c(nm_yx, nm_var)])
}


#' Convert covariance parameter list to/from vectorized form
#'
#' Converts parameter list `pars` to vector `p` and back again. A convenience function for
#' numerical optimization where objective functions accept numeric vectors only.
#'
#' When `p` is not supplied, the function un-lists the numeric values of `pars`, returning
#' them as a vector. When `na_omit=TRUE`, only the non-NA values are returned; and when
#' `iso=TRUE` the x kernel parameters are also omitted from the output vector.
#'
#' When `p` is supplied, the function copies its values to the corresponding entries in
#' `pars`, returning the result as a list. In this case, when `na_omit=TRUE`, only the
#' NA entries of `pars` are filled; and when `iso=TRUE`, parameters for the y kernel
#' are recycled to fill entries for the x kernel.
#'
#' The order returned (and expected in `p`, if supplied) is:
#'
#' 'eps', 'psill', 'y.rho', ('y.kap'), 'x.rho', ('x.kap')
#'
#' where 'y.kap' and 'x.kap' are omitted in kernels without shape parameters (see
#' `sk_corr`). Only the order matters here, as names are ignored in `p`.
#'
#' With `na_omit=FALSE`, `p` should have length 3-6, the same as the vector returned by
#' `sk_pars_update(pars)`; NAs in `p` are copied over in this case, effectively inverting
#' the vectorization.
#'
#' With `na_omit=TRUE`, values in `p` are copied only to the NA entries of `pars`; ie the
#' length of `p` should equal the number of NA entries in `pars` (less any redundant x
#' kernel parameters when `iso=TRUE`). Note that this does NOT invert the vectorization
#' `p=sk_pars_update(pars, na_omit=TRUE)`, which returns non-NA entries of `pars`.
#'
#' @param pars list of kernel parameters (in form understood by `sk_pars_make`)
#' @param p numeric vector (optional) kernel parameters to update in `pars`
#' @param iso logical, indicating to treat y and x kernel parameters as equal (see details)
#' @param eps_scaled logical, indicates to drop partial sill from input/output
#' @param na_omit logical, toggles handling of `NA`s (see details)
#'
#' @return numeric vector of parameters, or, if `p` is supplied, the updated `pars` list
#'
#' @export
#' @keywords internal
#' @family parameter managers
#'
#' @examples
#'# initialize a parameter list and pass to sk_pars_update
#' k = c('gau', 'mat')
#' pars_empty = sk_pars_make(k)
#' sk_pars_update(pars_empty)
#'
#' # pars can be a character vector passed directly to sk_pars_update
#' sk_pars_update(k)
#'
#' # single kernel definitions are duplicated
#' sk_pars_update('mat')
#'
#' # example parameter vector to illustrate ordering
#' p_update = seq(5)
#'
#' # (inverse) pass a modified vector back to the function to update pars
#' pars_filled = sk_pars_update(pars_empty, p_update)
#' print(pars_filled)
#'
#' # p_update is unchanged by round trip
#' sk_pars_update(pars_filled)
#'
#' # NAs in pars can be filtered out with na_omit=TRUE
#' pars_filled$eps = NA
#' pars_filled$x$kp[1] = NA
#' sk_pars_update(pars_filled)
#' sk_pars_update(pars_filled, na_omit=TRUE)
#'
#' # when updating parameters, NAs in pars identify parameters to receive the new values
#' p_update = abs(stats::rnorm(2))
#' sk_pars_update(pars_filled, p_update, na_omit=TRUE)
#'
#' # when na_omit=FALSE, all parameters in pars are updated
#' p_update = stats::rnorm(5)
#' sk_pars_update(pars_filled, p_update)
#'
#' # iso=TRUE is for when x kernel parameters are assigned values from the y kernel
#' pars_empty = sk_pars_make('mat')
#' p_update = seq(length(sk_pars_update(pars_empty, iso=TRUE)))
#' pars_filled = sk_pars_update(pars_empty, p_update, iso=TRUE)
#' print(pars_filled)
#'
#' # notice NA shape parameter in x ignored while NA in eps is copied
#' pars_filled$eps = NA
#' pars_filled$x$kp[2] = NA
#' sk_pars_update(pars_filled, iso=TRUE)
#'
#' # update calls with iso=TRUE should omit the x kernel parameters (p is shorter)
#' sk_pars_update(pars_filled, seq(4), iso=TRUE)
#'
#' # if eps_scaled=TRUE, psill is ignored (useful for when eps is scaled by psill)
#' sk_pars_update(pars_filled)
#' sk_pars_update(pars_filled, eps_scaled=TRUE)
#'
#' # notice value of psill not updated
#' sk_pars_update(pars_filled, -seq(5), eps_scaled=TRUE)
#'
#' # compare/combine with other modes
#' pars_filled$y$kp[2] = NA
#' sk_pars_update(pars_filled)
#' sk_pars_update(pars_filled, -seq(2), iso=TRUE, na_omit=TRUE)
#' sk_pars_update(pars_filled, -seq(4), iso=TRUE)
#' sk_pars_update(pars_filled, -seq(3), iso=TRUE, eps_scaled=TRUE)
#' sk_pars_update(pars_filled, -seq(3), na_omit=TRUE)
#'
sk_pars_update = function(pars, p=NULL, iso=FALSE, eps_scaled=FALSE, na_omit=FALSE)
{
  # expected list entries
  nm_yx = c('y', 'x')
  nm_var = 'eps'
  if(!eps_scaled) nm_var = c(nm_var, 'psill')

  # check for valid input and set names/NAs wherever they are missing
  pars = sk_pars_make(pars)
  msg_iso = 'kernel names for y and x must be the same when iso=TRUE'
  if( iso & ( pars[['y']][['k']] !=  pars[['x']][['k']] ) ) stop(msg_iso)

  # extract spatial parameters in order y, x (iso mode extracts only y)
  nm_extract = nm_yx[ seq(as.integer(!iso) + 1) ]
  kp = lapply(pars[nm_extract], \(x) x[['kp']])
  names(kp) = nm_extract

  # convert pars to named numeric
  p_old = c(do.call(c, pars[nm_var]), unlist(kp))
  p_old_miss = is.na(p_old)

  # when there are no NAs, or if not omitting NAs, all entries of pars are to be replaced below
  if( !any(p_old_miss) ) p_old_miss[] = TRUE
  p_out = p_old
  if( na_omit ) { p_out = p_old[!p_old_miss] } else { p_old_miss[] = TRUE }

  # return from vectorization mode
  if( is.null(p) ) return( p_out )

  # check for wrong length input
  len_expect = sum(p_old_miss)
  len_in = length(p)
  if( len_expect != len_in ) stop(paste('expected', len_expect, 'parameter(s) but got', len_in))

  # copy from input p to full vector of replacement values
  p_old[which(p_old_miss)] = p

  # copy variance components, then remove from parameter vector
  for(nm in nm_var) pars[[nm]] = unname(p_old[nm])
  p_old = p_old[-seq_along(nm_var)]

  # copy y, x kernel parameters to list
  for(nm_dim in nm_extract)
  {
    # get the parameter names expected in pars and p
    nm = names( sk_kp(pars[[nm_dim]][['k']]) )
    idx_p = seq_along( kp[[nm_dim]] )
    nm_p = names(p_old)[idx_p]
    pars[[nm_dim]][['kp']] = stats::setNames(p_old[nm_p], nm)
    p_old = p_old[-idx_p]
  }

  # clone y parameters if requested before returning from list mode
  if(iso) pars[['x']][['kp']] = pars[['y']][['kp']]
  return(pars)
}


#' Return named vector of Kronecker covariance parameters initialized to NA
#'
#' Convenience function for looking up the number of parameters for a given
#' Kronecker covariance model, and their names. Returns a vector of `NA`s that
#' can be used as a placeholder in kernel definition lists.
#'
#' @param k character, the kernel name, one of 'exp', 'gau', 'sph', 'gxp', 'mat'
#'
#' @return named vector of NAs, a placeholder for kernel parameters
#'
#' @export
#' @keywords internal
#' @family parameter managers
#'
#' @examples
#'
#' # there are only two possible return values
#' sk_kp('gau')
#' sk_kp('mat')
sk_kp = function(k)
{
  # two parameter structures are supported
  kp_range = c(rho=NA_real_)
  kp_range_shape = c(rho=NA_real_, kap=NA_real_)

  # known kernel names
  nm_range = c('exp', 'gau', 'sph')
  nm_range_shape = c('gxp', 'mat')
  msg_unknown = paste('kernel name', k, 'not recognized')
  if( !(k %in% c(nm_range, nm_range_shape)) ) stop(msg_unknown)

  if(k %in% nm_range) return(kp_range)
  if(k %in% nm_range_shape) return(kp_range_shape)
}

#' Extract Kronecker covariance  parameters as plot-friendly strings
#'
#' Generate strings describing the model and parameter values in a Kronecker covariance
#' parameter list `pars`. These are used to fill out titles and axis labels in calls to
#' `sk_plot_pars`.
#'
#' If `pars` is a parameter list (of the form returned by `sk_pars`), the
#' function returns a list of strings: a kernel family string ('k'), dimension-wise kernel
#' parameters (sub-list 'kp', with entries 'y' and 'x'), and a title containing the
#' kernel family, nugget effect, and partial sill.
#'
#' When `pars` is a character string, the function returns it unchanged. When `pars`
#' is a vector of two character strings, it concatenates them with separator " x ".
#' When `pars` is a list of kernel parameters for a single dimension, it returns a named
#' list of two character strings: the kernel name (named entry 'k') and the parameter(s)
#' (named entry 'kp'), parenthesized, in "name = value" format where "value" is rounded
#' to `nsig` significant digits).
#'
#' @param pars character, character vector, or list (see details)
#' @param nsig number of significant figures to print
#'
#' @return a character string or list of them
#'
#' @export
#' @keywords internal
#' @family parameter managers
#' @seealso sk_plot_pars
#'
#' @examples
#' kname = 'mat'
#' sk_to_string(kname)
#' sk_to_string(rep(kname, 2))
#' sk_to_string(list(k=kname))
#'
#' gdim = c(10, 15)
#' pars = sk_pars(gdim, kname)
#' sk_to_string(pars)
#'
sk_to_string = function(pars, nsig=3)
{
  # handle character input (kernel names)
  if( is.character(pars) )
  {
    # convert to expected list format for recursive call and handle unexpected input
    if( length(pars) == 1 ) return( pars )
    if( length(pars) == 2 ) return( paste(sapply(pars, sk_to_string), collapse=' x ') )
    stop('invalid argument to pars')
  }

  # unpack nugget and sill parameters
  eps = pars[['eps']]
  psill = pars[['psill']]
  kp = ''

  # single dimension case
  if( 'k' %in% names(pars) )
  {
    k = sk_to_string(pars[['k']])
    if( 'kp' %in% names(pars) )
    {
      # get kernel parameter names, collapse name value pairs and parenthesize result
      kp_nm = pars[['k']] |> sk_kp() |> names()
      kp_all = mapply(\(nm, p) paste(nm, '=', format(p, digits=nsig)), nm=kp_nm, p=pars[['kp']])
      kp = paste0('(', paste(kp_all, collapse=', '), ')')
    }
    return( list(k=k, kp=kp) )
  }

  # 2-dimensional case: recursive call
  kp_list = lapply(pars[c('y', 'x')], sk_to_string)
  k = paste(sapply(kp_list, \(xy) xy[['k']]), collapse=' x ')
  kp = sapply(kp_list, \(xy) xy[['kp']])

  # make a title expression
  eps_val = ifelse( is.null(eps), '', paste('nugget =', format(eps, digits=nsig)))
  psill_val = ifelse( is.null(psill), '', paste('partial sill =', format(psill, digits=nsig)))
  main_extra = paste0('(', paste(c(eps_val, psill_val), collapse=', '), ')')
  main = ifelse(main_extra == '(, )', k, paste(k, 'kernel', main_extra))

  # return as named vector
  return( list(k=k, kp=kp, main=main) )
}
