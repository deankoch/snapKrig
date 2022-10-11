# sk_index.R
# Dean Koch, 2022
# Helper functions for indexing grid points

#' Up or down-scale a sk grid by an integer factor
#'
#' Changes the resolution of a sk grid by a factor of `up` or `down`. For down-scaling, this
#' introduces `NA`s at unobserved grid points (and does no interpolation).
#'
#' Users should specify a sk grid `g` to re-scale and an integer scaling factor; either `up`
#' or `down` (and not both). This effects the scaling of resolution (`g[['gres']]`) by `up`
#' or `1/down`.
#'
#' `up` (or `down`) should be a vector of two positive integers, the desired re-scaling
#' factors in the y and x dimensions, in that order, or a single value to be used for both.
#'
#' When `up` is supplied, a lower resolution grid is returned comprising every `up`th grid
#' line of `g` along each dimension. All other grid lines, and any data values lying on them,
#' are ignored. `up` should be no greater than `dim(g) - 1`. Note that if `up` does not
#' evenly divide this number, the bounding box will shrink slightly.
#'
#' When `down` is supplied, the function returns a higher resolution grid (say `g_fine`) with
#' the same bounding box as `g`. Along each dimension, every `down`th grid line of `g_fine`
#' coincides with a grid line of `g`. Any non-NA values found in `g[]` are copied to `g_fine`,
#' and `g` can be recovered from `g_fine` with `sk_rescale(g_fine, up=down)`.
#'
#' @param g a sk grid or any grid object accepted by `sk`
#' @param up integer > 0, or vector of two, the up-scaling factors
#' @param down integer > 0, or vector of two, the down-scaling factors
#'
#' @return a sk grid of the requested resolution
#'
#' @export
#' @seealso sk sk_cmean
#' @family indexing functions
#' @family sk constructors
#'
#' @examples
#'
#' # example data
#' gdim = c(50, 53)
#' g = sk(gdim)
#' pars = modifyList(sk_pars(g), list(eps=1e-6))
#' g[] = sk_sim(g, pars)
#' plot(g)
#'
#' # upscale
#' plot(sk_rescale(g, up=1)) # does nothing
#' plot(sk_rescale(g, up=2))
#'
#' # downscale
#' sk_plot(sk_rescale(g, down=1)) # does nothing
#' sk_plot(sk_rescale(g, down=2))
#'
#' # length-2 vectors to rescale differently in x and y directions
#' plot(sk_rescale(g, up=c(2,3)))
#' plot(sk_rescale(g, down=c(2,3)))
#'
#' # invert a down-scaling
#' g_compare = sk_rescale(sk_rescale(g, down=c(5,3)), up=c(5,3))
#' all.equal(g, g_compare)
#'
#' # multi-layer example with missing data
#' n_pt = prod(gdim)
#' n_layer = 3
#'
#' # generate some more data and omit 50% of it
#' gval_multi = sk_sim(sk(list(gdim=gdim, gval=matrix(NA, n_pt, n_layer))), pars)
#' idx_miss = sample.int(n_pt, round(0.5*n_pt))
#' gval_multi[idx_miss,] = NA
#'
#' # plot third layer, then down-scaled and up-scaled versions
#' g_sim_multi = sk(gdim=gdim, gval=gval_multi)
#' sk_plot(g_sim_multi, layer=3)
#' sk_plot(sk_rescale(g=g_sim_multi, down=2), layer=3)
#' sk_plot(sk_rescale(g=g_sim_multi, up=2), layer=3)
#'
sk_rescale = function(g, up=NULL, down=NULL)
{
  # user has to pick one or the other
  is_up = !is.null(up)
  is_down = !is.null(down)
  if(is_up & is_down) stop('both up and down were specified')
  if( !(is_up | is_down) ) stop('either up or down must be specified')

  # unpack the grid object as list
  g = sk(g)
  gdim = dim(g)
  is_empty = is.null(g[['gval']])

  # multi-layer support
  if( is.matrix(g[['gval']]) )
  {
    # make simple grid with same dimensions and mapping vector in place of data
    g_first = sk(dim(g), vals=FALSE)
    g_first[] = g[['idx_grid']]

    # re-scale mapping vector to get new grid and mapping from the old
    g_result = sk_rescale(g_first, up=up, down=down)
    is_obs_first = !is.na(g_result)
    idx_keep = g_result[is_obs_first]

    # copy gval, omitting rows no longer mapped to grid (applies only to up-scaling)
    g_result[is_obs_first] = seq_along(idx_keep)
    g_result[['idx_grid']] = g_result[]
    g_result[['gval']] = g[['gval']][idx_keep,]

    # validate before returning
    return(sk_validate(g_result))
  }

  # up-scaling
  if(is_up)
  {
    # check for invalid up-scaling arguments
    msg_max = paste(gdim, collapse=', ')
    up = as.integer(up)
    if( any(up < 1) ) stop('upscaling factors cannot be less than 1')
    if( !any(up < gdim) ) stop( paste('upscaling factors must be less than', msg_max) )

    # set up dimensions of sub-grid
    ij_values = Map(function(d, r) seq(1, d, r), d=gdim, r=up)
    names(ij_values) = c('i', 'j')

    # build and return the sk sub-grid
    return(sk(gdim = stats::setNames(sapply(ij_values, length), c('y', 'x')),
              gyx = Map(function(yx, idx) yx[idx], yx=g[['gyx']], idx=ij_values),
              gval = g[ sk_sub_idx(gdim, ij_values, idx=TRUE) ]))
  }

  # check for invalid down-scaling arguments
  down = as.integer(down)
  if( any(down < 1) ) stop('downscaling factors cannot be less than 1')

  # set up dimensions of super-grid and copy CRS as needed
  gdim_new = gdim + (down - 1L) * (gdim - 1L)
  gres_new = g[['gres']] / down
  g_new = sk(gdim=gdim_new, gres=gres_new, vals=!is_empty)
  g_new[['gyx']] = Map(function(old, new) new + old[1], old=g[['gyx']], new=g_new[['gyx']])
  g_new[['crs']] = g[['crs']]

  # copy data to the sub-grid of the output and return
  ij_sub = stats::setNames(Map(function(d, b) seq(1, d, b), d=gdim_new, b=down), c('i', 'j'))
  idx_sub = sk_sub_idx(gdim_new, ij=ij_sub, idx=TRUE)
  if( !is_empty ) g_new[idx_sub] = g[]
  return(g_new)
}

#' Column-vectorization indices
#'
#' Maps matrix indices i, j to a single vectorized index, k
#'
#' Column vectorization (as in `base::as.vector`) builds a length(mn) vector by stacking
#' the columns of an m X n matrix, with the leftmost column appearing first in the vector
#' and the rightmost column last. Matrix element i,j gets mapped to element k = i + m * (j-1)
#' in the vector. This function returns that index.
#'
#' `ij` can be a matrix or a list of length-n vectors 'i' and 'j', or a vector
#' representing a single point at the given row (i) and column (j) number (in that order).
#'
#' `gdim` should either be an integer number of rows in the matrix, or a vector of the form
#' `c(ni, nj)` (the return value of `dim` for example) in which case its first element is used.
#'
#' @param ij n x 2 matrix, the row and column indices
#' @param gdim integer (or vector with first element equal to) the number of rows in the matrix
#' @param simplified, if FALSE, the function returns an n x 1 matrix
#'
#' @return integer vector, the vectorized `ij` indices
#'
#' @export
#' @keywords internal
#' @family indexing functions
#'
#' @examples
#' # define matrix dimensions and look up a specific index
#' gdim = c(4, 5)
#' ij = c(i=3, j=2)
#' sk_mat2vec(ij, gdim)
#'
#' # display all matrix indices in column-vectorized order
#' gyx = expand.grid(i=seq(gdim[1]), j=seq(gdim[2]))
#' result = sk_mat2vec(gyx, gdim)
#' data.frame(k=result, gyx)
#'
sk_mat2vec = function(ij, gdim, simplified=TRUE)
{
  # handle vector input to gdim
  if( length(gdim) > 1 ) gdim = gdim[1]

  # coerce list to matrix
  if( is.list(ij) )
  {
    if( !all( diff(sapply(ij, length)) == 0 ) ) stop('elements of list ij must have equal length')
    ij = do.call(cbind, ij)
  }

  # handle vector input (single point)
  if( is.vector(ij) ) ij = matrix(ij, nrow=1)

  # coerce input to matrix
  ij = as.matrix(ij)

  # check for invalid input
  if( any(ij[,1] > gdim) ) stop('ij contains "i" indices exceeding gdim')

  # return the vectorized index
  idx = ( ij %*% c(1, gdim) ) - gdim
  if( !simplified ) return( idx )
  return( as.vector(idx) )
}


#' Invert column-vectorization indices
#'
#' Inverts the function `sk_mat2vec`, returning matrix row and column numbers i, j,
#' given the column-vectorized index `k` and matrix dimensions `gdim`.
#'
#' Output indices are returned in a matrix with columns 'i', 'j' and rows in same
#' order as the input `k`. When `out='list'` list of vectors 'i' and 'j' (with entries
#' in the same order) is returned instead.
#'
#' The entries of `k` can be any permutation with replacement from `seq(prod(gdim))`
#'
#' @param k a vector of positive integers, the vector indices to look up
#' @param gdim integer (or vector with first element equal to) the number of rows in the matrix
#' @param out either 'matrix' or 'list'
#'
#' @return a two column matrix of integers (row and column numbers) with `length(k)` rows
#'
#' @export
#' @keywords internal
#' @family indexing functions
#'
#' @examples
#'
#' # show how elements are ordered in `base::matrix`
#' gdim = c(5, 6)
#' matrix_order = matrix(1:prod(gdim), gdim)
#' print(matrix_order)
#'
#' # identify the row and column numbers for specific entry, or several
#' sk_vec2mat(2, gdim)
#' sk_vec2mat(c(2, 10, 5), gdim)
#' sk_vec2mat(c(2, 10, 5), gdim, out='list')
#'
sk_vec2mat = function(k, gdim, out='matrix')
{
  # handle vector input to ni
  if( length(gdim) > 1 ) gdim = gdim[1]

  # compute column and row numbers
  cnum = as.integer( ceiling( k / gdim ) )
  rnum = as.integer( k - ( gdim * (cnum - 1) ) )

  # return as matrix
  if(out == 'matrix') return( cbind(i=rnum, j=cnum) )
  if(out == 'list') return( list(i=rnum, j=cnum) )
}


#' Find column-vectorized index of a sub-grid
#'
#' Returns a logical vector indicating all grid points lying on the specified sub-grid.
#' A grid point is `TRUE` only if both its i and j grid lines are found in `ij`.
#'
#' `ij` should be a list containing integer vectors named 'i', 'j', enumerating the
#' i and j grid lines of the desired sub-grid. If 'i' (or 'j') is missing, the function
#' automatically specifies all rows (or columns).
#'
#' If `idx=TRUE`, the function computes the vectorized index of the sub-grid points
#' with respect to the full grid `gdim` (see `sk_mat2vec`). Letting `i = c(i1, i2, ...im)`
#' and `j = c(j1, j2, ...in)` the function orders the sub-grid points as follows:
#'
#' (i1, j1), (i2, j1), ... (im, j1), (i1, j2), (i2, j2), ..., (i1, j3), ... (in, jm)
#'
#' This is the column-major vectorized order for the sub-grid (with y descending and
#' x ascending), provided the input grid line numbers in `ij` are in ascending order.
#' When `nosort=FALSE`, the function orders the input grid lines automatically.
#'
#' @param gdim integer vector, the number rows and columns in the full grid (in that order)
#' @param ij list containing vectors "i" and "j", the sub-grid row and column numbers
#' @param nosort logical, skips sorting the input vectors in `ij`
#' @param idx logical, indicates to return indices (default TRUE) versus logical vector
#'
#' @return integer or logical vector
#'
#' @export
#' @keywords internal
#' @family indexing functions
#'
#' @examples
#'
#' # example grid and a particular grid point
#' gdim = c(i=10, j=13)
#' ij_list = list(i=6, j=3)
#'
#' # sk_sub_idx returns a logical vector indexing the point (or the index itself)
#' is_pt = sk_sub_idx(gdim, ij_list)
#' idx_sub = sk_sub_idx(gdim, ij_list, idx=TRUE)
#' sk_plot(is_pt, gdim, col_grid='white', zlab='index', breaks=c('other', idx_sub))
#'
#' # equivalent call when ij_list is a single point
#' sk_mat2vec(ij_list, gdim) == idx_sub
#'
#' # if i or j are omitted from ij, the function returns the full row or column
#' is_col2 = sk_sub_idx(gdim, ij_list['i'])
#' is_row3 = sk_sub_idx(gdim, ij_list['j'])
#' sk_plot(is_col2, gdim, col_grid='white', breaks=c('other', paste('row', ij_list['i'])))
#' sk_plot(is_row3, gdim, col_grid='white', breaks=c('other', paste('col', ij_list['j'])))
#'
#' # indices in column-vectorized order
#' sk_sub_idx(gdim, list(i=2), idx=TRUE)
#' sk_sub_idx(gdim, list(j=3), idx=TRUE)
#' sk_sub_idx(gdim, idx=TRUE) # call without arguments returns all indices
#'
#' # bigger sub-grid example
#' origin_sg = c(5, 2) # assign i,j of top left point
#' gdim_sg = c(3, 4) # sub-grid dimensions (make sure this fits in gdim!)
#' ij_list = stats::setNames(Map(\(d, o) o + seq(d) - 1, gdim_sg, origin_sg), c('i', 'j'))
#' is_sg = sk_sub_idx(gdim, ij=ij_list)
#' sk_plot(is_sg, gdim, col_grid='white', zlab='sub-grid')
#'
#' # plot the index values: column major vectorization with y descending, x ascending
#' idx_sg = sk_sub_idx(gdim, ij=ij_list, idx=TRUE)
#' vec_order = rep(NA, prod(gdim))
#' vec_order[is_sg] = as.character(idx_sg)
#' sk_plot(vec_order, gdim, col_grid='black', zlab='vector idx')
#'
#' # example with j indices supplied in reverse (descending) order
#' ij_list_xflip = modifyList(ij_list, list(j=rev(ij_list[['j']])))
#'
#' # ordering in the vectors ij$i and ij$j doesn't matter if `nosort=FALSE` or `idx=FALSE`
#' identical(is_sg, sk_sub_idx(gdim, ij=ij_list, nosort=TRUE))
#' all.equal(which(is_sg), sk_sub_idx(gdim, ij=ij_list_xflip, idx=TRUE))
#'
#' # when `nosort=TRUE` and `idx=TRUE` we get the same indices but in a different order
#' idx_sg_xflip = sk_sub_idx(gdim, ij=ij_list_xflip, idx=TRUE, nosort=TRUE)
#' all.equal(sort(idx_sg), sort(idx_sg_xflip))
#' all.equal(idx_sg, idx_sg_xflip)
#' vec_order[is_sg] = as.character(idx_sg_xflip)
#' sk_plot(vec_order, gdim, col_grid='black', zlab='vector index')
#'
sk_sub_idx = function(gdim, ij=NULL, idx=FALSE, nosort=FALSE)
{
  # check input and set expected order
  ij_nm = c(y='i', x='j')
  if( !all( c('y', 'x') %in% names(gdim) ) ) gdim = stats::setNames(gdim, c('y', 'x'))
  if( is.null(ij) ) ij = stats::setNames(lapply(gdim, seq), ij_nm)
  if( is.null(names(ij)) & (length(ij) == 2) ) names(ij) = ij_nm
  if( sum( ij_nm %in% names(ij) ) != length(ij) ) stop('unexpected names in ij')

  # set default i to select all rows
  if( is.null(ij[['i']]) ) { ij[['i']] = seq( gdim['y'] ) } else {

    # otherwise coerce to integer and sort as needed
    ij[['i']] = as.integer(ij[['i']])
    if( !nosort ) ij[['i']] = sort(ij[['i']])
  }

  # set default j to select all columns
  if( is.null(ij[['j']]) ) { ij[['j']] = seq( gdim['x'] ) } else {

    # otherwise coerce to integer and sort as needed
    ij[['j']] = as.integer(ij[['j']])
    if( !nosort ) ij[['j']] = sort(ij[['j']])
  }

  # count desired sub-grid dimensions and compute the indexing vector
  n_ij = sapply(ij, length)
  idx_out = rep(ij[['i']], n_ij['j']) + rep(gdim['y'] * (ij[['j']] - 1L), each=n_ij['i'])
  if(idx) return(as.integer(unname(idx_out)))

  # compute the logical vector
  is_out = logical(prod(gdim))
  is_out[idx_out] = TRUE
  return(is_out)
}

#' Return a sub-grid of a sk grid object
#'
#' Creates a "sk" object containing only the grid-lines specified in `idx_keep`. Alternatively,
#' grid lines to remove can be specified in `idx_rem`.
#'
#' One of `idx_keep` or `idx_rem` (but not both) can be specified, and the grid line numbers
#' (not intercepts) should be supplied in ascending order in list entries named "i" and "j".
#'
#' If `idx_rem` is specified, `mirror=TRUE` will cause the selection in `idx_rem` to be
#' reflected about the central grid line (useful for specifying outer grid lines). `mirror`
#' is ignored if `idx_keep` is specified instead.
#'
#' Default `idx=FALSE` causes the function to return the sub-grid as a sk grid object.
#' If `idx=TRUE`, the function instead returns a list containing `idx_keep` and `idx_rem` as
#' specified above.
#'
#' If neither `idx_keep` nor `idx_rem` is supplied, the function removes outer grid lines
#' iteratively (selecting the one with highest proportion of `NA`s), attempting to find a
#' complete sub-grid (having no `NA`s) somewhere in the interior. This heuristic is designed
#' for rasters with few `NA`s, all located around the perimeter.
#'
#' @param g sk grid or any grid-like object accepted by `sk`
#' @param ij_keep list of grid line numbers ("i" and "j") forming regular sub-grid
#' @param ij_rem list of grid line numbers ("i" and "j") whose exclusion forms regular sub-grid
#' @param idx logical, if TRUE the function returns a list containing `ij_keep` and `ij_rem`
#' @param mirror logical, whether to mirror the selection in `ij_rem` (see details)
#'
#' @export
#' @keywords internal
#' @family sk constructors
#'
#' @examples
#'
#' # make an example grid
#' g = sk(c(50, 100))
#' g[] = apply(expand.grid(g[['gyx']]), 1, \(z) cos( 2*sum(z^2) ) )
#' plot(g)
#'
#' # subset by specifying grid lines to keep
#' ij_keep = list(i=seq(1, 50, by=2), j=seq(1, 50, by=2))
#' g_keep = sk_sub(g, ij_keep)
#' plot(g_keep)
#'
#' # get the indices kept and removed
#' idx = sk_sub(g, ij_keep, idx=TRUE)
#'
#' # equivalent call specifying grid lines to omit
#' g_rem = sk_sub(g, ij_rem=idx[['rem']])
#' identical(g_rem, g_keep)
#'
#' # remove some data around the edges of the grid
#' idx = sk_sub(g, ij_rem=list(i=seq(10), j=seq(10)), mirror=TRUE, idx=TRUE)
#' idx_y_pts = sk_sub_idx(dim(g), idx[['rem']]['i'], idx=TRUE)
#' idx_x_pts = sk_sub_idx(dim(g), idx[['rem']]['j'], idx=TRUE)
#' idx_pts = c(idx_y_pts, idx_x_pts)
#' idx_na = sort(sample(idx_pts, 0.6*length(idx_pts)))
#' g[idx_na] = NA
#' plot(g)
#'
#' # identify the interior sub-grid that is complete
#' g_sub = sk_sub(g)
#' print(g_sub)
#' plot(g_sub)
#'
#' # verify it is as large as expected
#' ( dim(g) - dim(g_sub) ) == sapply(idx[['rem']], length)
#'
sk_sub = function(g, ij_keep=NULL, ij_rem=NULL, idx=FALSE, mirror=FALSE)
{
  # expected names of vectors and lists in ij
  ij_nm = c('i', 'j')

  # open grid as snapKrig list object and copy original dimensions
  g = sk(g)
  gdim = dim(g)

  # sort out the supplied arguments, preferring ij_keep when both supplied
  has_keep = !is.null(ij_keep)
  has_rem = !is.null(ij_rem)
  if(has_keep & has_rem)
  {
    warning('ignoring ij_rem')
    ij_rem = NULL
    has_rem = FALSE
  }

  # construct ij_rem if it wasn't supplied but ij_keep was
  ij_keep = ij_keep[ij_nm]
  if(has_keep & !has_rem)
  {
    #ij_keep[[i]]
    ij_rem = Map(function(d, i) seq(d)[!(seq(d) %in% i)], d=gdim, i=ij_keep)
    names(ij_rem) = ij_nm
  }

  # handle case where user supplied either ij_keep or ij_rem
  if( !is.null(ij_rem) )
  {
    # mirror the grid lines if requested (only when idx_keep not specified)
    if(mirror & has_rem) ij_rem = Map(function(d, i) c(i, d - i + 1L), d=gdim, i=ij_rem)
    ij_ok = Map(function(d, i) !any( (i < 1) | (i > d) ), d=gdim, i=ij_rem)

    # a validity check
    if( any(!unlist(ij_ok)) ) stop('ij_rem contained indices less than 1 or greater than gdim')

    # update list of grid lines to keep
    ij_keep = Map(function(d, i) seq(d)[!(seq(d) %in% i)], d=gdim, i=ij_rem)
    gdim_new = sapply(ij_keep, length)

    # another validity check
    if( any( gdim_new < 1 ) ) stop('request resulted in an empty sub-grid')

    # assign names lost in the `Map` calls and sort the removed grid lines
    names(ij_keep) = ij_nm
    names(ij_rem) = ij_nm
    ij_rem = lapply(ij_rem, sort)

    # identify points to include in the sub-grid
    is_sub = sk_sub_idx(gdim, ij_keep, nosort=TRUE, idx=FALSE)

    # verify that the result will be a regular sub-grid
    sub_result = sk_sub_find(is_sub, gdim)
    if(is.null(sub_result)) stop('request resulted in an irregular sub-grid')

    # return either the grid line indices or the sub-grid itself as a sk object
    if(idx) return( list(rem=ij_rem, keep=ij_keep) )
    gyx_new = Map(function(gl, i) gl[i], gl=g[['gyx']], i=sub_result[['ij']])
    return(sk(gdim=gdim_new, gyx=gyx_new, gval=g[is_sub]))
  }

  # indicator for NA's as a matrix representation of the grid data
  is_NA_mat = matrix(is.na(g), gdim)
  has_NA = any(is_NA_mat)

  # initialize lists of grid line numbers that we may examine
  i_seq = seq(gdim[['y']])
  j_seq = seq(gdim[['x']])

  # initialize indexing variables for the grid lines to count (start with outer ones)
  ij_count = lapply(gdim, function(d) c(1L, d))
  names(ij_count) = ij_nm

  # loop while there are still NAs in the sub-grid
  while(has_NA)
  {
    # count NAs around on outer grid lines
    ni_NA = apply(is_NA_mat[ij_count[['i']], j_seq], 1, sum)
    nj_NA = apply(is_NA_mat[i_seq, ij_count[['j']]], 2, sum)

    # score each grid line and identify the worst
    i_score = 1 - ( ni_NA/length(j_seq) )
    j_score = 1 - ( nj_NA/length(i_seq) )
    dim_worst = ij_nm[ which.min(c(min(i_score), min(j_score))) ]
    side_worst = which.min( list(i=i_score, j=j_score)[[ dim_worst ]] )

    # remove the worst grid line
    inc = ifelse(side_worst==2, -1, 1)
    ij_count[[dim_worst]][side_worst] = inc + ij_count[[dim_worst]][side_worst]

    # copy the new grid line vectors then check again for NAs
    i_seq = seq(ij_count[['i']][1], ij_count[['i']][2])
    j_seq = seq(ij_count[['j']][1], ij_count[['j']][2])
    has_NA = any(is_NA_mat[i_seq, j_seq])
  }

  # copy results to ij_keep then generate corresponding ij_rem
  ij_keep = list(i=i_seq, j=j_seq)
  ij_rem = Map(function(d, i) seq(d)[ !(seq(d) %in% i) ], d=gdim, i=ij_keep)
  names(ij_rem) = ij_nm
  if(idx) return( list(remove=ij_rem, keep=ij_keep) )

  # recursive call to make the sk object
  return(sk_sub(g, ij_rem=ij_rem))
}

#' Find complete regular sub-grids in a sk grid object
#'
#' If a sk grid `g` has missing values (`NA`s) but the set of non-`NA` points form a
#' complete (regular) sub-grid, this function finds its grid lines, resolution, and
#' dimensions. If no eligible sub-grids are found, the function returns `NULL`.
#'
#' A sub-grid is only eligible if it contains ALL of the non-`NA` points in `g` and none
#' of the `NA`s. For example if a single point missing from the sub-grid, or a single non-`NA`
#' point lies outside the sub-grid, the function will fail to detect any sub-grids and return
#' `NULL`. If no points are `NA`, the function returns indices for the full grid.
#'
#' The returned list contains the following named elements:
#'
#'  * `ij` the grid line numbers of the sub-grid with respect to `g`
#'  * `res_scale` the resolution scaling factor (relative increase in grid line spacing of `g`)
#'  * `gdim` the number of y and x grid lines in the sub-grid
#'
#' As in `sk`, each of these is given in the 'y', 'x' order.
#'
#' Users can also pass the logical vector returned by `!is.na(g)` instead of `g`, in which
#' case argument `gdim` must also be specified. This can be much faster with large grids.
#'
#' @param g logical vector, sk grid, or any grid object accepted by `sk`
#' @param gdim integer vector, the grid dimensions (in order 'y', 'x')
#'
#' @return `NULL` or list of information about the location and spacing of the sub-grid
#' within `g` (see details)
#'
#' @export
#' @keywords internal
#' @family indexing functions
#'
#' @examples
#'
#' # define a grid and example data
#' gdim = c(50, 53)
#' g = sk(gdim)
#' g[] = sk_sim(g, modifyList(sk_pars(g), list(eps=1e-12)))
#' plot(g)
#'
#' # define a super-grid containing the original data and make sure we can find it
#' g_big = sk_rescale(g, down=3)
#' plot(g_big)
#' print(sk_sub_find(g_big))
#'
#' # define a smaller sub-grid at random
#' spacing = sapply(floor(gdim/10), function(x) 1 + sample.int(x, 1))
#' gdim_sg = sapply(floor( (gdim - 1) / spacing), function(x) sample.int(x, 1))
#' ij_first = sapply(gdim - ( spacing * gdim_sg ), function(x) sample.int(x, 1))
#'
#' # find index of sub-grid lines and vectorized index of points
#' ij_sg = Map(\(idx, r, n) seq(idx, by=r, length.out=n), idx=ij_first, r=spacing, n=gdim_sg)
#' names(ij_sg) = c('i', 'j')
#' is_sg = sk_sub_idx(gdim, ij_sg, idx=FALSE)
#'
#' # assign values to the sub-grid points
#' g_sub = sk(gdim)
#' g_sub[is_sg] = g[is_sg]
#' plot(g_sub, zlab='sub-grid')
#'
#' # call the function and check for expected results
#' sub_result = sk_sub_find(g_sub)
#' all.equal(unname(sub_result[['gdim']]), gdim_sg)
#' all.equal(unname(sub_result[['ij']]), unname(ij_sg))
#'
#' # sub grids with side length 1 have no spacing defined along that dimension
#' spacing[gdim_sg==1] = NA
#'
#' # check consistency in spacing
#' all.equal(unname(sub_result[['res_scale']]), spacing)
#'
#' # can also call on the vector and supply gdim separately
#' identical(sub_result, sk_sub_find(!is.na(g_sub), dim(g_sub)))
#'
sk_sub_find = function(g, gdim=NULL)
{
  # expected order for dimensional info
  nm_dim = c('y', 'x')

  # open various grid objects with sk
  if(!is.logical(g))
  {
    # overwrite g with logical NAs indicator, copying gdim first
    g = sk(g)
    gdim = dim(g)
    g = !is.na(g)

  } else {

    # validity checks
    if( anyNA(g) ) stop('logical vector g cannot have NAs')
    if( is.null(gdim) ) stop('full grid dimensions gdim must be supplied when g is a vector')
    msg_len = paste('expected logical vector g to have length', prod(gdim), 'but got', length(g))
    if( prod(gdim) != length(g)) stop(msg_len)

    # set names for user-supplied gdim
    gdim = stats::setNames(as.integer(gdim), c('y', 'x'))
  }

  # need this to get indices of first, second, and last elements in sub-grid
  idx_obs = which(g)
  n_obs = sum(g)

  # find the dimensions of the smallest sub-grid enclosing all observed points
  ij_bbox = sk_vec2mat(c(idx_obs[1], idx_obs[n_obs]), gdim)
  gdim_bbox = apply(ij_bbox, 2, function(x) diff(x) + 1L)
  if( !all(gdim_bbox > 0) ) return(NULL)

  # compute number of rows in sub-grid and do first existence check
  skip_i = diff(idx_obs[1:2]) - 1L
  ni = as.integer( ifelse(skip_i > -1L, ( gdim_bbox[['i']] + skip_i ) / (1L + skip_i), NA) )
  if( ni == 1 ) skip_i = NA
  if( (ni == 0) | ( ni %% 1L != 0 ) ) return(NULL)

  # compute number of columns in sub-grid and do second existence check
  nj = n_obs / ni
  if( nj %% 1L != 0 ) return(NULL)
  nj = as.integer(nj)
  skip_j = as.integer( ifelse(nj > 1, (gdim_bbox[['j']] - nj) / (nj - 1L), NA) )

  # sub-grid resolution scaling factor and dimensions
  res_ij = 1L + c(skip_i, skip_j)
  res_scale = 1L / res_ij

  # sub-grid line indices
  ij_obs = Map(\(idx, r, n) seq(idx, by=r, length.out=n), idx=ij_bbox[1,], r=res_ij, n=c(ni, nj))
  ij_na = sapply(ij_obs, anyNA)
  ij_obs[ij_na] = as.list(ij_bbox[1,])[ ij_na ]

  # final existence check
  idx_sub = sk_sub_idx(gdim, ij_obs)
  if( !all( g[idx_sub] ) ) return(NULL)

  # return sub-grid info in a list
  gdim_result = sapply(ij_obs, length)
  names(gdim_result) = nm_dim
  names(ij_obs) = nm_dim
  names(res_ij) = nm_dim
  return(list(ij=ij_obs, res_scale=res_ij, gdim=gdim_result))
}



