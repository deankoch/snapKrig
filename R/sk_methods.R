#' sk_methods.R
#' Dean Koch, 2022
#' S3 methods for sk grid list objects
#'
#' Replace a sk list element (double-bracket assign)
#'
#' Replaces entries in the sk list object. This does no validation. If it did, then
#' `sk_validate` would have an infinite recursion problem (it uses `[[<-`). Users should
#' pass the results to `sk_validate` afterwards unless they know what they're doing.
#'
#' @param x a sk object
#' @param value the replacement object
#'
#' @return a "sk" object
#' @export
#'
#' @examples
#' # sk list elements are interrelated - for example gres must match spacing in gyx
#' g = sk_validate(list(gval=stats::rnorm(10^2), gdim=10, gres=0.5))
#' g[['gres']] = 2 * g[['gres']]
#' g[['gyx']] = lapply(g[['gyx']], function(x) 2*x)
#' sk_validate(g)
#'
`[[.sk<-` = function(x, value) {

  # construct the object directly
  sk_make(NextMethod())
}

#'
#' Extract a sk list element (single-bracket access)
#'
#' Copies the specified list element or grid point value
#'
#' Behavior depends on the class of i. For character vectors this extracts the named list
#' entries of x. For numeric, it accesses the vectorized grid data values. For multi-layer
#' objects, a layer can be specified in j.
#'
#' the default `NULL` for `i` and `j` is treated as numeric, and is shorthand for all
#' indices. For example if `x` has a single-layer `x[]` returns all grid data in a vector.
#' If `x` is multi-layer `x[,1]` all grid data from the first layer, and `x[]` returns all
#' layers, as a matrix.
#'
#' @param x a sk object
#' @param i column-vectorized index
#' @param j index of layer (only for multi-layer x)
#' @param drop ignored
#' @param ... ignored
#'
#' @return a list, vector, or matrix (see description)
#' @export
#'
#' @examples
#' # define a sk list and extract two of its elements
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' g[c('gdim', 'gres')]
#'
#' # display all the grid data as a vector or a matrix
#' g[]
#' matrix(g[], dim(g))
#'
#' # extract a particular grid point or a subset
#' g[1]
#' g[seq(5)]
#'
`[.sk` = function(x, i=NULL, j=NULL, drop=FALSE, ...) {

  # identify multi-layer rasters
  is_multi = is.matrix(x[['gval']])
  if( !is.null(j) & !is_multi ) stop('incorrect number of dimensions')

  # character indices specify sub-lists
  if(is.character(i)) { return(NextMethod()) } else {

    # calls with [] should return all data so we specify all indices
    if(is.null(i)) i = seq_along(x)

    # convert logical sk grids to indices
    if( inherits(i, 'sk') ) i = as.logical(i[])

    # handle empty gval
    if( is.null(x[['gval']]) ) {

      if( is.logical(i) ) return( rep(NA_real_, sum(i)) )
      return( rep(NA_real_, length(i)) )
    }

    # numeric indices are for the vectorized grid values
    if( !is_multi )
    {
      # non-sparse representation is direct column-vectorization
      return(x[['gval']][i])

    } else {

      # default selects all layers
      if(is.null(j)) j = seq(ncol(x[['gval']]))

      # vectorized indices mapped to rows via idx_grid in sparse representation
      return( matrix(x[['gval']][ x[['idx_grid']][i], j], ncol=length(j)) )
    }
  }
}


#'
#' Single-bracket assign
#'
#' Behavior depends on the class of i. For character vectors, this assigns to
#' the named list entries of x (as usual). For numeric indices, it assigns
#' vectorized grid data values. For multi-layer objects, specify the layer in j
#' and supply a matrix for replacement
#'
#' @param x an sk object
#' @param i column-vectorized index
#' @param j index of layer (only for multi-layer x)
#' @param value the replacement values
#'
#' @return the "sk" object with the specified subset replaced by value
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' print(g)
#' g[1] = NA
#' print(g)
#'
`[<-.sk` = function(x, i=NULL, j=NULL, value) {

  # check if we are replacing everything
  replace_all = is.null(i) & is.null(j)

  # identify multi-layer rasters
  is_multi = is.matrix(x[['gval']])
  if( !is.null(j) & !is_multi ) stop('only one index needed for single-layer grids')

  # calls with [] assign all entries/rows
  if(is.null(i)) i = seq_along(x)

  # character indices specify sub-lists (handled by default method for lists)
  if(is.character(i)) { NextMethod() } else {

    # convert logical sk grids to indices
    if( inherits(i, 'sk') ) i = as.logical(i[])

    # numeric indices are for the vectorized grid values
    if( !is_multi )
    {
      # initialize gval vector if it's not already there
      if( is.null(x[['gval']]) ) x[['gval']] = rep(NA_real_, length(x))

      # overwriting the whole object prevents unwanted coercion
      if(replace_all) { x[['gval']] = value } else { x[['gval']][i] = as.vector(value) }

    } else {

      # default selects all layers
      if(is.null(j)) j = seq(ncol(x[['gval']]))
      if(replace_all) { x[['gval']] = value } else {

        # extract all, modify the requested columns, copy back
        gval_mod = x[]
        gval_mod[i, j] = as.matrix(value)
        x[['gval']] = gval_mod
      }

      # idx_grid will be re-computed in sk_make (as needed)
      x[['idx_grid']] = NULL
    }

    # check validity, update is_na, etc then return
    return(sk(x))
  }
}


#'
#' Math group generics
#'
#' All except the `cumsum` family are supported
#'
#' @param x a sk object
#' @param ... further arguments to Math
#'
#' @return the "sk" object with data values transformed accordingly
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' summary(g)
#' summary(abs(g))
#' summary(exp(g))
#'
Math.sk = function(x, ...)
{
  # this should include everything except generics requiring an ordering in the vector
  generic_ok = c('abs', 'sign', 'sqrt', 'floor', 'ceiling', 'trunc', 'round', 'signif',
                 'exp', 'log', 'expm1', 'log1p', 'cos', 'sin', 'tan', 'cospi', 'sinpi',
                 'tanpi', 'acos', 'asin', 'atan', 'lgamma', 'gamma', 'digamma', 'trigamma')

  # dispatch to next method when generic makes sense for a grid dataset
  if(.Generic %in% generic_ok)
  {
    # replace the data values and return the validated sk object
    return( sk(utils::modifyList(x, list(gval=get(.Generic)(x[], ...)))) )

  } else { stop(paste(.Generic, 'not defined for "sk" objects')) }
}

#'
#' Operations group generics
#'
#' Applies the operation point-wise to grid data values.
#'
#' The function extracts the grid data `x[]` from all sk class arguments `x`, prior to
#' calling the default method. Before returning, the result is copied back to the grid object
#' of the second argument (or the first, if the second is not of class "sk").
#'
#' Note that the compatibility of the two arguments is not checked beyond matching
#' dimension (with vectors recycled as needed). This means for example you can do
#' operations on two grids representing different areas, so long as they have the
#' same `gdim`.
#'
#' @param e1 a "sk" object, vector or matrix
#' @param e2 a "sk" object, vector or matrix
#'
#' @return a "sk" object (a copy of `e1` or `e2`) with modified data values
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#'
#' # verify a trig identity using Ops and Math
#' summary( cos(g)^2 + sin(g)^2 )
#'
#' # create a logical grid indicating points satisfying a condition
#' plot(g < 0)
#' all( !(g > 0) == (g[] < 0) )
#'
#' # test negation
#' all( (-g)[] == -(g[]) )
#'
Ops.sk = function(e1, e2)
{
  # extract grid data as vector/matrix, saving a copy of sk object
  if(methods::is(e1, 'sk'))
  {
    g_out = e1
    e1 = e1[]
  }

  # handle missing e2 in calls like `!g` and `-g`
  if(nargs() == 1L)
  {
    # single argument case
    if(.Generic == '!')
    {
      # dispatch to default method then replace the data values in output sk object
      return( sk(utils::modifyList(g_out, list(gval=get(.Generic)(e1)))) )
    }

    # first argument is implied to be zero (as in `-e1`)
    e2 = e1
    e1 = 0

  } else {

    # same as with first argument
    if(methods::is(e2, 'sk'))
    {
      # overwrites copy of the first sk object
      g_out = e2
      e2 = e2[]
    }
  }

  # dispatch to default method then replace the data values in output sk object
  return( sk(utils::modifyList(g_out, list(gval=get(.Generic)(e1, e2)))) )
}

#'
#' Summary group generics
#'
#' Computes a summary statistic of the grid data values.
#'
#' @param ... further arguments to Summary
#' @param na.rm logical, whether to remove `NA`s before doing the computation
#'
#' @return a "sk" object with data values
#' @export
#'
#' @examples
#'
#' # make test example and verify common summary stats
#' g = sk(gval=1+abs(stats::rnorm(4^2)), gdim=4, gres=0.5)
#' sum(g) == sum(g[])
#' min(g) == min(g[])
#' max(g) == max(g[])
#'
#' # list input can include scalar, vector (but first argument must be sk class)
#' sum(g, -g)
#' sum(g, -g, 0)
#'
#' # calculate product in two ways
#' prod(g) - exp(sum(log(g)))
#'
#' # check logical conditions
#' any(g>0)
#' all(g>0)
#'
#' # count the number of grid points satisfying a condition
#' sum(g < 2)
#'
Summary.sk = function(..., na.rm=FALSE)
{
  # extract all grid data values and omit NA's if requested
  vals = do.call(c, lapply(list(...), function(x) x[]))
  if(na.rm) vals = vals[ !is.na(vals) ]
  if( length(vals) == 0 ) stop('all values are NA')

  # dispatch to default method then replace the data values in output sk object
  .Primitive(.Generic)(vals)
}


#'
#' Auto-printing
#'
#' Prints dimensions and indicates if the grid has values assigned
#'
#' This prints "(not validated)" if the sk object has no `is_na` entry,
#' to remind users to run `sk_validate`.
#'
#' @param x a sk object
#' @param ... ignored
#'
#' @return nothing
#' @export
#'
#' @examples
#' sk_make(list(gdim=10, gres=0.5))
#' sk_validate(sk_make(list(gdim=10, gres=0.5)))
print.sk = function(x, ...)
{
  # string with grid dimensions
  gdim_msg = paste(paste(x[['gdim']], collapse=' x '))

  # message about completeness
  if(is.null(x[['is_obs']])) {complete_msg = 'not validated\n'} else {

    n_miss = sum(!x[['is_obs']])
    n = prod(x[['gdim']])
    complete_msg = 'incomplete\n'
    if( n_miss == n ) complete_msg = 'empty\n'
    if( n_miss == 0 ) complete_msg = 'complete\n'
  }

  # check for matrix values
  if(is.matrix(x[['gval']]))
  {
    # indicate matrix values with layer
    n_layer = ncol(x[['gval']])
    layer_msg = paste0(n_layer, ' layer', ifelse(n_layer==1, '', 's'), '\n')
    complete_msg = paste0(complete_msg, layer_msg)
  }

  # for non-numeric types, print the type instead of the range
  if( !is.numeric(x[['gval']]) & !is.null(x[['gval']]) )
  {
    complete_msg = paste0(complete_msg, '(', typeof(x[['gval']]), ' data)\n')
  }

  # print messages
  cat(paste(gdim_msg, complete_msg))


}

#'
#' Grid summary
#'
#' Prints detailed information about a grid
#'
#' All dimensional information (`gdim`, `gres`, `gyx`) is printed in the order y, x
#'
#' @param object an sk object
#' @param ... ignored
#'
#' @return nothing
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' summary(g)
#' g[1] = NA
#' summary(g)
#'
summary.sk = function(object, ...)
{
  # validate input and check for CRS
  object = sk_validate(object)
  has_crs = !is.null(object[['crs']])

  # check grid size and number of NAs
  n = length(object)
  n_miss = sum(!object[['is_obs']])
  n_obs = n - n_miss

  # check range
  range_msg = NULL
  if( n_obs > 0 )
  {
    # for non-numeric types, print the type instead of the range
    if( !is.numeric(object[['gval']]) )
    {
      range_msg = paste0('(', typeof(object[['gval']]), ' data)')

    } else {

      range_obs = range(object, na.rm=TRUE)
      n_sig = min(max(3, round(-log(diff(range_obs), base=10))), 12)
      range_string = paste(format(range_obs, digits=n_sig, trim=TRUE), collapse=', ')
      range_msg = paste0('range [', range_string, ']')
    }
  }

  # report sample size and number observed in incomplete case
  n_msg = paste(n, 'points')
  if( !(n_miss %in% c(0, n)) ) n_msg = paste0(n_msg, '\n', n_obs, ' observed')

  # report completeness
  crs_msg = ifelse(has_crs, ' geo-referenced ', ' ')
  title_msg = paste0('incomplete', crs_msg, 'sk grid')
  if( n_miss == n ) title_msg = paste0('empty', crs_msg, 'sk grid')
  if( n_miss == 0 ) title_msg = paste0('complete', crs_msg, 'sk grid')

  # check for matrix values
  if(is.matrix(object[['gval']]))
  {
    # indicate matrix values with layer
    n_layer = ncol(object[['gval']])
    layer_msg = paste0(n_layer, ' layer', ifelse(n_layer==1, '', 's'))
    title_msg = paste0(title_msg, '\n', layer_msg)
  }

  # print top level info
  hrule_msg = '..............................\n'
  cat(paste(c(title_msg, n_msg, range_msg, hrule_msg), collapse='\n'))
  #cat('\n')

  # strings for grid dimensions and resolution
  gdim_msg = paste('dimensions :', paste(object[['gdim']], collapse=' x '))
  gres_msg = paste('resolution :', paste(object[['gres']], collapse=' x '))

  # strings for extent
  ext_list = apply(sk_coords(object, corner=TRUE, quiet=TRUE), 2, range, simplify=FALSE)
  ext_ranges = sapply(ext_list, function(e) paste0('[', paste(e, collapse=', '), ']'))
  ext_msg = paste('    extent :', paste(ext_ranges, collapse=' x '))

  # print messages
  cat(paste(c(gdim_msg, gres_msg, ext_msg), collapse='\n'))
  cat('\n')
}

#'
#' Heatmap plots
#'
#' A wrapper for sk_plot
#'
#' @param x a sk object
#' @param ... other arguments passed to sk_plot
#'
#' @return nothing
#' @export
#' @seealso sk_plot
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' plot(g)
plot.sk = function(x, ...) sk_plot(x, ...)


#' The number of grid-points
#'
#' Returns the total number of points in the grid, which is the product of the
#' number of y and x grid lines (`gdim`).
#'
#' @param x a sk object
#'
#' @return integer
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' length(g)
length.sk = function(x) as.integer(prod(x[['gdim']]))

#' Grid dimensions
#'
#' Returns `gdim`, the number of y and x grid lines, in that order.
#'
#' @param x a sk object
#'
#' @return integer vector
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' dim(g)
dim.sk = function(x) x[['gdim']]

#' Indices of grid points with missing data (NAs)
#'
#' Returns a logical vector indicating which grid points have NA values assigned
#'
#' @param x a sk object
#'
#' @return a logical vector the same length as x
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' g[c(1,3)] = NA
#' is.na(g)
is.na.sk = function(x)
{
  # handle empty grids
  n = length(x)
  if( is.null(x[['gval']]) ) return( !logical(n) )

  # return the pre-computed NAs index or tell the user to validate
  if(is.null(x[['is_obs']])) stop('invalid sk object. Pass it to sk_validate first')
  return(!x[['is_obs']])
}

#' Check for presence of grid points with missing data (NAs)
#'
#' Returns a logical indicating if any of the grid points are NA
#'
#' @param x a sk object
#' @param recursive ignored
#'
#' @return logical
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' anyNA(g)
#' g[1] = NA
#' anyNA(g)
#'
anyNA.sk = function(x, recursive)
{
  # handle empty grids
  if( is.null(x[['gval']]) ) return(TRUE)

  # check the pre-computed NAs index or tell the user to validate
  if(is.null(x[['is_obs']])) stop('invalid sk object. Pass it to sk_validate first')
  return(!all(x[['is_obs']]))
}

#' Calculate the mean value in a grid
#'
#' This calculates the mean over all layers (if any)
#'
#' @param x a sk object
#' @param ... further arguments to default mean method
#'
#' @return numeric
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' mean(g) == mean(g[])
#'
mean.sk = function(x, ...)
{
  # handle empty grids
  if( is.null(x[['gval']]) ) return(NA)

  # call the default method on the vector (or matrix)
  mean(x[], ...)
}


#' Convert grid data to vector of specified mode
#'
#' Returns a vector of the specified mode, representing the vectorized grid data.
#' For multi-layer `x`, the first layer is returned.
#'
#' For single layer `x`, and with default `mode='any'`, this is the same as `x[]`
#'
#' @param x a sk object
#' @param mode passed to as.vector
#'
#' @return a vector of the specified mode
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' as.vector(g)
as.vector.sk = function(x, mode='any')
{
  # extract data as vector or matrix and pass to default method
  as.vector(x[], mode='any')
}

#' Coerce grid values to logical
#'
#' @param x a sk object
#' @param ... further arguments to as.logical
#'
#' @return a "sk" object with logical data values
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=sample(c(0,1), 4^2, replace=TRUE), gdim=4, gres=0.5))
#' g[]
#' as.logical(g)[]
#'
#' # "range" for logical is reported as integer
#' summary(as.logical(g))
#'
as.logical.sk = function(x, ...)
{
  # replace the data values and return the sk object
  utils::modifyList(x, list(gval=as.logical(x[])))
}

#' Coerce grid values to numeric (double type)
#'
#' This also adds support for as.numeric
#'
#' @param x a sk object
#' @param ... further arguments to as.double
#'
#' @return an "sk" object with numeric data values
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=sample(c(FALSE, TRUE), 4^2, replace=TRUE), gdim=4, gres=0.5))
#' g[]
#' as.numeric(g)[]
#'
as.double.sk = function(x, ...)
{
  # replace the data values and return the sk object
  utils::modifyList(x, list(gval=as.double(x[], ...)))
}


#' Coerce grid values to integer
#'
#' @param x a sk object
#' @param ... further arguments to as.integer
#'
#' @return an "sk" object with integer data values
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' g[]
#' as.integer(g)[]
as.integer.sk = function(x, ...)
{
  # replace the data values and return the sk object
  utils::modifyList(x, list(gval=as.integer(x[], ...)))
}

#' convert to matrix
#'
#' Returns a matrix representation of the grid data. This is shorthand for
#' extracting the data using `x[]` (single layer) or `x[,j]` (multi-layer),
#' then passing the result to `matrix` along with `dim(x)`.
#'
#' @param x a sk object
#' @param layer integer, for multi-layer grids, the layer number to return
#' @param rownames.force ignored
#' @param ... further arguments to as.matrix
#'
#' @return the grid data as a matrix
#' @export
#'
#' @examples
#' g = sk_validate(list(gval=stats::rnorm(4^2), gdim=4, gres=0.5))
#' plot(g)
#' as.matrix(g)
as.matrix.sk = function(x, rownames.force=NA, layer=1, ...)
{
  if( is.null(x[['idx_grid']]) ) { return(matrix(x[], dim(x))) } else {

    # for multi-layer objects we have to select a particular column
    return(matrix(x[, layer], dim(x)))
  }
}

