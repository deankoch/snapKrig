#' bk_methods.R
#' Dean Koch, 2022
#' methods for bk objects (grid lists)
#'
#' ## INDEXING METHODS
#'
#' Replace a bk list element (double-bracket assign)
#'
#' Replaces entries in the bk list object. This does no validation (since `bk_validate`
#' uses `[[<-`, so there would be an infinite recursion problem) so users are advised
#' to pass the results to `bk_validate` afterwards unless they know what they're doing.
#'
#' @param x a bk object
#' @param value the replacement object
#'
#' @return a "bk" object
#' @export
#'
#' @examples
#' # bk list elements are interrelated - for example gres must match spacing in gyx
#' g = bk_validate(list(gval=rnorm(10^2), gdim=10, gres=0.5))
#' g[['gres']] = 2 * g[['gres']]
#' g[['gyx']] = lapply(g[['gyx']], function(x) 2*x)
#' bk_validate(g)
#'
`[[.bk<-` = function(x, value) {

  # construct the object directly
  bk_make(NextMethod())
}

#'
#' Extract a bk list element (single-bracket access)
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
#' @param x a bk object
#' @param i column-vectorized index
#' @param j index of layer (only for multi-layer x)
#' @param ... ignored
#'
#' @return a list, vector, or matrix (see description)
#' @export
#'
#' @examples
#' # define a bk list and extract two of its elements
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
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
`[.bk` = function(x, i=NULL, j=NULL, drop=FALSE, ...) {

  # identify multi-layer rasters
  is_multi = is.matrix(x[['gval']])
  if( !is.null(j) & !is_multi ) stop('incorrect number of dimensions')

  # character indices specify sub-lists
  if(is.character(i)) { return(NextMethod()) } else {

    # calls with [] return all data so we specify all indices
    if(is.null(i)) i = seq_along(x)

    # handle empty gval
    if( is.null(x[['gval']]) ) return( rep(NA_real_, length(i)) )

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
#' @param x a bk object
#' @param value the replacement values
#' @param i column-vectorized index
#' @param j index of layer (only for multi-layer x)
#' @param ... ignored
#'
#' @return the "bk" object with the specified subset replaced by value
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' print(g)
#' g[1] = NA
#' print(g)
#'
`[<-.bk` = function(x, i=NULL, j=NULL, value, ...) {

  # check if we are replacing everything
  replace_all = is.null(i) & is.null(j)

  # identify multi-layer rasters
  is_multi = is.matrix(x[['gval']])
  if( !is.null(j) & !is_multi ) stop('only one index needed for single-layer grids')

  # calls with [] assign all entries/rows
  if(is.null(i)) i = seq_along(x)

  # character indices specify sub-lists (handled by default method for lists)
  if(is.character(i)) { NextMethod() } else {

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

      # idx_grid will be re-computed in bk_make (as needed)
      x[['idx_grid']] = NULL
    }

    # check validity, update n_missing, etc then return
    return(bk_validate(bk_make(x)))
  }
}

#'
#' ## GROUP GENERICS
#'
#' Math group generics
#'
#' All except the `cumsum` family are supported
#'
#' @param x a bk object
#' @param ... further arguments to Math
#'
#' @return the "bk" object with data values transformed accordingly
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' summary(g)
#' summary(abs(g))
#' summary(exp(g))
#'
Math.bk = function(x, ...)
{
  # this should include everything except generics requiring an ordering in the vector
  generic_ok = c('abs', 'sign', 'sqrt', 'floor', 'ceiling', 'trunc', 'round', 'signif',
                 'exp', 'log', 'expm1', 'log1p', 'cos', 'sin', 'tan', 'cospi', 'sinpi',
                 'tanpi', 'acos', 'asin', 'atan', 'lgamma', 'gamma', 'digamma', 'trigamma')

  # dispatch to next method when generic makes sense for a grid dataset
  if(.Generic %in% generic_ok)
  {
    # replace the data values and return the bk object
    return( modifyList(x, list(gval=get(.Generic)(x[], ...))) )

  } else { stop(paste(.Generic, 'not defined for "bk" objects')) }
}

#'
#' Operations group generics
#'
#' Applies the operation point-wise to grid data values.
#'
#' The function extracts the grid data `x[]` from all bk class arguments `x`, prior to
#' calling the default method. Before returning, the result is copied back to the grid object
#' of the second argument (or the first, if the second is not of class "bk").
#'
#' Note that the compatibility of the two arguments is not checked beyond matching
#' dimension (with vectors recycled as needed). This means for example you can do
#' operations on two grids representing different areas, so long as they have the
#' same `gdim`.
#'
#' @param e1 a "bk" object, vector or matrix
#' @param e2 a "bk" object, vector or matrix
#'
#' @return a "bk" object (a copy of `e1` or `e2`) with modified data values
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#'
#' # verify a trig identity using Ops and Math
#' summary( cos(g)^2 + sin(g)^2 )
#'
#' # create a logical grid indicating points satisfying a condition
#' plot(g < 0)
#'
#' # test negation
#' all( (-g)[] == -(g[]) )
#'
Ops.bk = function(e1, e2)
{
  # extract grid data as vector/matrix, saving a copy of bk object
  if(is(e1, 'bk'))
  {
    g_out = e1
    e1 = e1[]
  }

  # handle missing e2 in calls like `-g` and `+g`
  if(nargs() == 1L)
  {
    # produces 0 +/- e2
    e2 = e1
    e1 = 0

  } else {

    # same as with first argument
    if(is(e2, 'bk'))
    {
      # overwrites copy of the first bk object
      g_out = e2
      e2 = e2[]
    }
  }

  # dispatch to default method then replace the data values in output bk object
  return( modifyList(g_out, list(gval=get(.Generic)(e1, e2))) )
}

#'
#' Summary group generics
#'
#' Computes a summary statistic of the grid data values.
#'
#' @param ... further arguments to Summary
#' @param na.rm logical, whether to remove `NA`s before doing the computation
#'
#' @return a "bk" object with data values
#' @export
#'
#' @examples
#'
#' # make test example and verify common summary stats
#' g = bk_validate(list(gval=1+abs(rnorm(4^2)), gdim=4, gres=0.5))
#' sum(g) == sum(g[])
#' min(g) == min(g[])
#' max(g) == max(g[])
#'
#' # list input can include scalar, vector (but first argument must be bk class)
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
Summary.bk = function(..., na.rm=FALSE)
{
  # extract all grid data values and omit NA's if requested
  vals = do.call(c, lapply(list(...), function(x) x[]))
  if(na.rm) vals = vals[ !is.na(vals) ]
  if( length(vals) == 0 ) stop('all values are NA')

  # dispatch to default method then replace the data values in output bk object
  .Primitive(.Generic)(vals)
}


#'
#' ## OTHER METHODS
#'
#' Auto-printing
#'
#' Prints dimensions and indicates if the grid has values assigned
#'
#' This prints "(not validated)" if the bk object has no `n_missing` entry.
#' This is to remind users to run `bk_validate`.
#'
#' @param x a bk object
#' @param ... ignored
#'
#' @return nothing
#' @export
#'
#' @examples
#' bk_make(list(gdim=10, gres=0.5))
#' bk_validate(bk_make(list(gdim=10, gres=0.5)))
print.bk = function(x, ...)
{
  # string with grid dimensions
  gdim_msg = paste(paste(x[['gdim']], collapse=' x '))

  # message about completeness
  n_miss = x[['n_missing']]
  if(is.null(n_miss)) {complete_msg = 'not validated\n'} else {

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

  # for non-numeric types, print the class instead of the range
  if( !is.numeric(x[['gval']]) & !is.null(x[['gval']]) )
  {
    complete_msg = paste0(complete_msg, '(', class(x[['gval']])[1], ' data)\n')
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
#' @param x a bk object
#' @param ... ignored
#'
#' @return nothing
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' summary(g)
#' g[1] = NA
#' summary(g)
#'
summary.bk = function(x, ...)
{
  # validate input and check for CRS
  x = bk_validate(x)
  has_crs = !is.null(x[['crs']])

  # check grid size and number of NAs
  n = length(x)
  n_miss = x[['n_missing']]
  n_obs = n - n_miss

  # check range
  range_msg = NULL
  if( n_obs > 0 )
  {
    # for non-numeric types, print the class instead of the range
    if( !is.numeric(x[['gval']]) )
    {
      range_msg = paste0('(', class(x[['gval']])[1], ' data)')

    } else {

      range_obs = range(x, na.rm=TRUE)
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
  title_msg = paste0('incomplete', crs_msg, 'bk grid')
  if( n_miss == n ) title_msg = paste0('empty', crs_msg, 'bk grid')
  if( n_miss == 0 ) title_msg = paste0('complete', crs_msg, 'bk grid')

  # check for matrix values
  if(is.matrix(x[['gval']]))
  {
    # indicate matrix values with layer
    n_layer = ncol(x[['gval']])
    layer_msg = paste0(n_layer, ' layer', ifelse(n_layer==1, '', 's'))
    title_msg = paste0(title_msg, '\n', layer_msg)
  }

  # print top level info
  hrule_msg = '..............................\n'
  cat(paste(c(title_msg, n_msg, range_msg, hrule_msg), collapse='\n'))
  #cat('\n')

  # strings for grid dimensions and resolution
  gdim_msg = paste('dimensions :', paste(x[['gdim']], collapse=' x '))
  gres_msg = paste('resolution :', paste(x[['gres']], collapse=' x '))

  # strings for extent
  ext_list = apply(bk_coords(x, corner=TRUE, quiet=TRUE), 2, range, simplify=FALSE)
  ext_ranges = sapply(ext_list, function(e) paste0('[', paste(e, collapse=', '), ']'))
  ext_msg = paste('    extent :', paste(ext_ranges, collapse=' x '))

  # print messages
  cat(paste(c(gdim_msg, gres_msg, ext_msg), collapse='\n'))
  cat('\n')
}

#'
#' Heatmap plots
#'
#' A wrapper for bk_plot
#'
#' @param x a bk object
#' @param ... other arguments passed to bk_plot
#'
#' @return nothing
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' plot(g)
plot.bk = function(x, ...) bk_plot(x, ...)


#' The number of grid-points
#'
#' Returns the total number of points in the grid, which is the product of the
#' number of y and x grid lines (`gdim`).
#'
#' @param x a bk object
#'
#' @return integer
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' length(g)
length.bk = function(x) as.integer(prod(x[['gdim']]))

#' Grid dimensions
#'
#' Returns `gdim`, the number of y and x grid lines, in that order.
#'
#' @param x a bk object
#'
#' @return integer vector
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' dim(g)
dim.bk = function(x) x[['gdim']]

#' Indices of grid points with missing data (NAs)
#'
#' Returns a logical vector indicating which grid points have NA values assigned
#'
#' @param x a bk object
#'
#' @return a logical vector the same length as x
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' g[c(1,3)] = NA
#' is.na(g)
is.na.bk = function(x)
{
  # handle empty grids
  n = length(x)
  if( is.null(x[['gval']]) ) return( !logical(n) )

  # identify multi-layer rasters
  is_multi = is.matrix(x[['gval']])
  if( is_multi ) { return(is.na(x[['idx_grid']])) } else { return(is.na(x[])) }
}

#' Check for presence of grid points with missing data (NAs)
#'
#' Returns a logical indicating if any of the grid points are NA
#'
#' @param x a bk object
#'
#' @return logical
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' anyNA(g)
#' g[1] = NA
#' anyNA(g)
#'
anyNA.bk = function(x)
{
  # handle empty grids
  if( is.null(x[['gval']]) ) return(TRUE)

  # identify multi-layer rasters
  is_multi = is.matrix(x[['gval']])

  # call the default method on the vectors
  if( is_multi ) {  return( anyNA(x[['idx_grid']]) ) } else {
   return( anyNA(x[]) )
  }
}

#' Calculate the mean value in a grid
#'
#' This calculates the mean over all layers (if any)
#'
#' @param x a bk object
#' @param ... further arguments to default mean method
#'
#' @return numeric
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' mean(g) == mean(g[])
#'
mean.bk = function(x, ...)
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
#' @param x a bk object
#' @param mode passed to as.vector
#'
#' @return a vector of the specified mode
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' as.vector(g)
as.vector.bk = function(x, mode='any')
{
  # extract data as vector or matrix and pass to default method
  as.vector(x[], mode='any')
}

#' Coerce grid values to logical
#'
#' @param x a bk object
#' @param ... further arguments to as.logical
#'
#' @return a "bk" object with logical data values
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=sample(c(0,1), 4^2, replace=TRUE), gdim=4, gres=0.5))
#' g[]
#' as.logical(g)[]
#'
#' # "range" for logical is reported as integer
#' summary(as.logical(g))
#'
as.logical.bk = function(x, ...)
{
  # replace the data values and return the bk object
  modifyList(x, list(gval=as.logical(x[])))
}

#' Coerce grid values to numeric (double type)
#'
#' This also adds support for as.numeric
#'
#' @param x a bk object
#' @param mode passed to as.double
#'
#' @return a "bk" object with numeric data values
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=sample(c(F,T), 4^2, replace=TRUE), gdim=4, gres=0.5))
#' g[]
#' as.numeric(g)[]
#'
as.double.bk = function(x, ...)
{
  # replace the data values and return the bk object
  modifyList(x, list(gval=as.double(x[], ...)))
}


#' Coerce grid values to integer
#'
#' @param x a bk object
#' @param mode passed to as.vector
#'
#' @return a "bk" object with integer data values
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' g[]
#' as.integer(g)[]
as.integer.bk = function(x, ...)
{
  # replace the data values and return the bk object
  modifyList(x, list(gval=as.integer(x[], ...)))
}

#' convert to matrix
#'
#' Returns a matrix representation of the grid data. This is shorthand for
#' extracting the data using `x[]` (single layer) or `x[,j]` (multi-layer),
#' then passing the result to `matrix` along with `dim(x)`.
#'
#' @param x a bk object
#' @param layer integer, for multi-layer grids, the layer number to return
#' @param rownames.force ignored
#' @param ... further arguments to as.matrix
#'
#' @return the grid data as a matrix
#' @export
#'
#' @examples
#' g = bk_validate(list(gval=rnorm(4^2), gdim=4, gres=0.5))
#' plot(g)
#' as.matrix(g)
as.matrix.bk = function(x, rownames.force=NA, layer=1, ...)
{
  if( is.null(x[['idx_grid']]) ) { return(matrix(x[], dim(x))) } else {

    # for multi-layer objects we have to select a particular column
    return(matrix(x[, layer], dim(x)))
  }
}

