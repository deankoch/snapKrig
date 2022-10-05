# bk_index.R
# Dean Koch, 2022
# Convenience functions for indexing grid points


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
#' @export
#'
#' @examples
#' # define matrix dimensions and look up a specific index
#' gdim = c(4, 5)
#' ij = c(i=3, j=2)
#' bk_mat2vec(ij, gdim)
#'
#' # display all matrix indices in column-vectorized order
#' gyx = expand.grid(i=seq(gdim[1]), j=seq(gdim[2]))
#' result = bk_mat2vec(gyx, gdim)
#' data.frame(k=result, gyx)
#'
bk_mat2vec = function(ij, gdim, simplified=TRUE)
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
#' Inverts the function `bk_mat2vec`, returning matrix row and column numbers i, j,
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
#' @export
#'
#' @examples
#'
#' # show how elements are ordered in `base::matrix`
#' gdim = c(5, 6)
#' matrix_order = matrix(1:prod(gdim), gdim)
#' print(matrix_order)
#'
#' # identify the row and column numbers for specific entry, or several
#' bk_vec2mat(2, gdim)
#' bk_vec2mat(c(2, 10, 5), gdim)
#' bk_vec2mat(c(2, 10, 5), gdim, out='list')
#'
bk_vec2mat = function(k, gdim, out='matrix')
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
#' with respect to the full grid `gdim` (see `bk_mat2vec`). Letting `i = c(i1, i2, ...im)`
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
#' @export
#'
#' @examples
#'
#' # example grid and a particular grid point
#' gdim = c(i=10, j=13)
#' ij_list = list(i=6, j=3)
#'
#' # bk_sub_idx returns a logical vector indexing the point (or the index itself)
#' is_pt = bk_sub_idx(gdim, ij_list)
#' idx_sub = bk_sub_idx(gdim, ij_list, idx=TRUE)
#' bk_plot(is_pt, gdim, col_grid='white', zlab='index', breaks=c('other', idx_sub))
#'
#' # equivalent call when ij_list is a single point
#' bk_mat2vec(ij_list, gdim) == idx_sub
#'
#' # if i or j are omitted from ij, the function returns the full row or column
#' is_col2 = bk_sub_idx(gdim, ij_list['i'])
#' is_row3 = bk_sub_idx(gdim, ij_list['j'])
#' bk_plot(is_col2, gdim, col_grid='white', breaks=c('other', paste('row', ij_list['i'])))
#' bk_plot(is_row3, gdim, col_grid='white', breaks=c('other', paste('col', ij_list['j'])))
#'
#' # indices in column-vectorized order
#' bk_sub_idx(gdim, list(i=2), idx=TRUE)
#' bk_sub_idx(gdim, list(j=3), idx=TRUE)
#' bk_sub_idx(gdim, idx=TRUE) # call without arguments returns all indices
#'
#' # bigger sub-grid example
#' origin_sg = c(5, 2) # assign i,j of top left point
#' gdim_sg = c(3, 4) # sub-grid dimensions (make sure this fits in gdim!)
#' ij_list = stats::setNames(Map(\(d, o) o + seq(d) - 1, gdim_sg, origin_sg), c('i', 'j'))
#' is_sg = bk_sub_idx(gdim, ij=ij_list)
#' bk_plot(is_sg, gdim, col_grid='white', zlab='sub-grid')
#'
#' # plot the index values: column major vectorization with y descending, x ascending
#' idx_sg = bk_sub_idx(gdim, ij=ij_list, idx=TRUE)
#' vec_order = rep(NA, prod(gdim))
#' vec_order[is_sg] = as.character(idx_sg)
#' bk_plot(vec_order, gdim, col_grid='black', zlab='vector idx')
#'
#' # example with j indices supplied in reverse (descending) order
#' ij_list_xflip = modifyList(ij_list, list(j=rev(ij_list[['j']])))
#'
#' # ordering in the vectors ij$i and ij$j doesn't matter if `nosort=FALSE` or `idx=FALSE`
#' identical(is_sg, bk_sub_idx(gdim, ij=ij_list, nosort=TRUE))
#' all.equal(which(is_sg), bk_sub_idx(gdim, ij=ij_list_xflip, idx=TRUE))
#'
#' # when `nosort=TRUE` and `idx=TRUE` we get the same indices but in a different order
#' idx_sg_xflip = bk_sub_idx(gdim, ij=ij_list_xflip, idx=TRUE, nosort=TRUE)
#' all.equal(sort(idx_sg), sort(idx_sg_xflip))
#' all.equal(idx_sg, idx_sg_xflip)
#' vec_order[is_sg] = as.character(idx_sg_xflip)
#' bk_plot(vec_order, gdim, col_grid='black', zlab='vector index')
#'
bk_sub_idx = function(gdim, ij=NULL, idx=FALSE, nosort=FALSE)
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

#' Return a sub-grid of a blitzKrig grid object
#'
#' Removes grid lines as specified in `ij` and returns the result as a blitzKrig grid
#' list; If `ij` is empty, the function attempts to find a set of grid lines that forms
#' a complete regular sub-grid of `g_obs`.
#'
#' `ij` can be a list of integer vectors, named `i` and `j`, specifying grid line numbers to
#' remove, in ascending order. `ij` can also be a nested list, with the grid lines
#' to remove in element `remove`. Alternatively, users can specify the grid lines to NOT
#' remove in element `keep`, and the function will remove the others. If both `keep` and
#' `remove` are specified, the function ignores `keep` and uses `remove`.
#'
#' Default `idx=FALSE` causes the function to return the sub-grid as a blitzKrig grid object.
#' If `idx=TRUE`, the function instead returns a list containing `keep` and `remove` as
#' specified above. This can be passed back to bk_sub when repeatedly cropping multiple
#' grids in the same way.
#'
#' `mirror=TRUE` indicates to mirror the supplied selection of grid lines in `ij`
#' about the central grid line. For example `ij_rem=list(i=1)` removes both the first
#' (left-most) column and the last (right-most).
#'
#' When `ij` is an empty `list()`, the function sets `ij` automatically using a
#' greedy algorithm that iteratively checks all (four) outer grid lines and removes
#' the one with the highest proportion of NAs. Ties are broken in clockwise order, starting
#' from the left vertical. The algorithm stops when the remaining sub-grid is complete.
#'
#'
#'
bk_sub = function(g_obs, ij=list(), mirror=FALSE, idx=FALSE)
{
  # expected names of vectors and lists in ij
  ij_nm = c('i', 'j')

  # open grid as blitzKrig list object
  g_obs = bk(g_obs)
  gdim_old = stats::setNames(g_obs[['gdim']], ij_nm)

  # compute grid lines to remove from keep
  if( 'keep' %in% names(ij) ) ij_rem = Map(function(g, i) seq(g)[!(seq(g) %in% i)],
                                           g = gdim_old, i = ij[['keep']])

  # copy remove list or set it to NA when not supplied
  if( any(ij_nm %in% names(ij) ) ) ij = list(remove=ij)
  if( 'remove' %in% names(ij) ) ij_rem = ij[['remove']]
  if( !any(c('remove', 'keep') %in% names(ij)) ) ij_rem = NA

  # ij_rem is specified
  if( !anyNA(ij_rem) )
  {
    # validate and fill in missing ij_rem arguments
    ij_supplied = stats::setNames(ij_nm %in% names(ij_rem), nm=ij_nm)
    msg_ij = 'ij_rem must be a list with named element(s) i and/or j'
    if( !is.list(ij_rem) | !any(ij_supplied) ) stop(msg_ij)
    if( !ij_supplied['i'] ) ij_rem[['i']] = integer()
    if( !ij_supplied['j'] ) ij_rem[['j']] = integer()

    # mirror the grid lines if requested, and another validity check
    if(mirror) ij_rem = Map(function(g, i) c(i, g - i + 1L), g=gdim_old, i=ij_rem)
    ij_ok = Map(function(g, i) !any( (i < 1) | (i > g) ), g=gdim_old, i=ij_rem)
    if( any(!unlist(ij_ok)) ) stop('ij_rem contained indices less than 1 or greater than gdim')

    # identify the grid lines to keep
    ij_keep = Map(function(i, g) seq(g)[!(seq(g) %in% i)], g=gdim_old, i=ij_rem)
    gdim_new = sapply(ij_keep, length)
    if( any( gdim_new < 1 ) ) stop('the supplied ij_rem resulted in an empty sub-grid')
    ij_keep = stats::setNames(ij_keep, ij_nm)
    if(idx) return( list(remove=ij_rem, keep=ij_keep) )

    # modify and return the list object
    g_obs[['gval']] = g_obs[['gval']][ bk_sub_idx(gdim_old, ij_keep, nosort=TRUE) ]
    g_obs[['gyx']] = Map(function(g, i) g[i], g=g_obs[['gyx']], i=ij_keep)
    g_obs[['gdim']] = gdim_new
    return(bk(g_obs))
  }

  is_NA_mat = matrix(is.na(g_obs[['gval']]), gdim_old)

  # initialize indexing variables
  has_NA = any(is_NA_mat)
  ij_count = stats::setNames(lapply(gdim_old, function(g) c(1L, g)), ij_nm)
  i_seq = seq(gdim_old[[1]])
  j_seq = seq(gdim_old[[2]])

  # loop while there are NAs in the sub-grid
  while( has_NA )
  {
    # check edge grid lines for missing data
    ni_NA = apply(is_NA_mat[ij_count[['i']], j_seq], 1, sum)
    nj_NA = apply(is_NA_mat[i_seq ,ij_count[['j']]], 2, sum)

    # score each grid line and identify the worst
    i_score = 1 - ( ni_NA/length(j_seq) )
    j_score = 1 - ( nj_NA/length(i_seq) )
    dim_worst = ij_nm[ which.min(c(min(i_score), min(j_score))) ]
    side_worst = which.min( list(i=i_score, j=j_score)[[ dim_worst ]] )

    # remove the worst grid line from sub-grid and check again for NAs
    inc = ifelse(side_worst==2, -1, 1)
    ij_removed = ij_count[[dim_worst]][side_worst]
    ij_count[[dim_worst]][side_worst] = inc + ij_removed

    #cat(paste('\nremoved grid line', ij_removed, 'from dimension', dim_worst))

    # copy the new grid line vectors
    i_seq = seq(ij_count[['i']][1], ij_count[['i']][2])
    j_seq = seq(ij_count[['j']][1], ij_count[['j']][2])
    has_NA = any(is_NA_mat[i_seq, j_seq])
  }

  ij_keep = list(i=i_seq, j=j_seq)
  ij_rem = Map(function(g, i) seq(g)[ !(seq(g) %in% i) ], g=gdim_old, i=ij_keep)
  if(idx) return( list(remove=ij_rem, keep=ij_keep) )
  return(bk_sub(g_obs, ij=list(remove=ij_rem)))
}



#' Up or down-scale a grid
#'
#' Changes the resolution of a grid by a factor of `up` or `down`.
#'
#' Users should specify a grid `g` to re-scale and an integer scaling factor; either `up`
#' or `down`. This effects the scaling of resolution (`g$gres`) by `up` or `1/down`.
#'
#' `up` (or `down`) should be a vector of two positive integers, supplying the re-scaling
#' factors in the y and x dimensions in that order, or a single value to be used for both.
#'
#' When `up` is supplied, a lower resolution grid is returned comprising every `up`th grid
#' line of `g` along each dimension. All other grid lines, and any data values lying on them,
#' are ignored. `up` should be no greater than `g$gdim - 1`. Note that if `up` does not
#' evenly divide this number, the bounding box will shrink slightly.
#'
#' When `down` is supplied, the function returns a higher resolution grid (`g_fine`) with
#' the same bounding box as `g`. Along each dimension, every `down`th grid line of `g_fine`
#' coincides with a grid line of `g`. Any values found in `g$gval` are copied to `g_fine`,
#' and un-mapped grid lines in `g_fine` are initialized to `NA`. Recover `g` from `g_fine`
#' with `bk_rescale(g_fine, up=down)`.
#'
#' @param g any object accepted or returned by `bk`
#' @param up integer > 0, or vector of two, the up-scaling factor(s)
#' @param down integer > 0, or vector of two, the down-scaling factor(s)
#'
#' @return a grid list of the form returned by `bk`
#' @export
#'
#' @examples
#'
#' # example data
#' gdim = c(50, 53)
#' g = bk(gdim)
#' pars = modifyList(bk_pars(g), list(eps=1e-6))
#' gval = bk_sim(g, pars)
#' g_obs = modifyList(g, list(gval=gval))
#' bk_plot(g_obs)
#'
#' # upscale
#' bk_plot(bk_rescale(g=g_obs, up=1)) # does nothing
#' bk_plot(bk_rescale(g=g_obs, up=2))
#'
#' # downscale
#' bk_plot(bk_rescale(g=g_obs, down=1)) # does nothing
#' bk_plot(bk_rescale(g=g_obs, down=2))
#'
#' # length-2 vectors to rescale differently in x and y directions
#' bk_plot(bk_rescale(g=g_obs, up=c(2,3)))
#' bk_plot(bk_rescale(g=g_obs, down=c(2,3)))
#'
#' # invert a down-scaling
#' g_obs_compare = bk_rescale(bk_rescale(g=g_obs, down=c(5,3)), up=c(5,3))
#' identical(g_obs, g_obs_compare)
#'
#' # multi-layer example with missing data
#' n_pt = prod(gdim)
#' n_layer = 3
#'
#' # generate some data and omit 50% of it
#' gval_multi = bk_sim(bk(list(gdim=gdim, gval=matrix(NA, n_pt, n_layer))), pars)
#' idx_miss = sample.int(n_pt, round(0.5*n_pt))
#' gval_multi[idx_miss,] = NA
#'
#' # plot third layer, then down-scaled and up-scaled versions
#' g_sim_multi = modifyList(g, list(gval=gval_multi))
#' bk_plot(g_sim_multi, layer=3)
#' bk_plot(bk_rescale(g=g_sim_multi, down=2), layer=3)
#' bk_plot(bk_rescale(g=g_sim_multi, up=2), layer=3)
#'
bk_rescale = function(g, up=NULL, down=NULL)
{
  # user has to pick one or the other
  is_up = !is.null(up)
  is_down = !is.null(down)
  if(is_up & is_down) stop('both up and down were specified')
  if( !(is_up | is_down) ) stop('either up or down must be specified')

  # unpack the grid object as list
  g = bk(g)
  gdim = g[['gdim']]

  # multi-layer support
  if( !is.null(g[['idx_grid']]) )
  {
    # re-scale mapping vector to get new grid and mapping from the old
    g_first = modifyList(g, list(gval=g[['idx_grid']], idx_grid=NULL))
    g_result = bk_rescale(g_first, up=up, down=down)
    is_obs_first = !is.na(g_result[['gval']])

    # copy gval, omit rows no longer mapped to grid (applies only to up-scaling)
    idx_keep = g_result[['gval']][is_obs_first]
    g_result[['gval']] = g[['gval']][idx_keep,]

    # compute and copy the new indexing vector
    g_result[['idx_grid']] = match(seq(prod(g_result[['gdim']])), which(is_obs_first))
    return(g_result)
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
    ij_values = stats::setNames(Map(function(d, r) seq(1, d, r), d=gdim, r=up), c('i', 'j'))
    gdim_new = stats::setNames(sapply(ij_values, length), c('y', 'x'))

    # overwrite data vector with the subset then update other fields in g
    idx_sub = bk_sub_idx(gdim, ij_values, idx=TRUE)
    g[['gval']] = g[['gval']][idx_sub]
    g[['gres']] = g[['gres']] * up
    g[['gdim']] = gdim_new
    g[['gyx']] = Map(function(yx, idx) yx[idx], yx=g[['gyx']], idx=ij_values)
    return(g)
  }

  # check for invalid down-scaling arguments
  down = as.integer(down)
  if( any(down < 1) ) stop('downscaling factors cannot be less than 1')

  # set up dimensions of super-grid
  gdim_new = gdim + (down - 1L) * (gdim - 1L)
  gres_new = g[['gres']] / down
  g_new = bk(list(gdim=gdim_new, gres=gres_new), vals=!is.null(g[['gval']]))
  g_new[['gyx']] = Map(function(old, new) new + old[1], old=g[['gyx']], new=g_new[['gyx']])
  g_new[['crs']] = g[['crs']]

  # copy data to the sub-grid of the output and return
  ij_sub = stats::setNames(Map(function(d, b) seq(1, d, b), d=gdim_new, b=down), c('i', 'j'))
  idx_sub = bk_sub_idx(gdim_new, ij=ij_sub, idx=TRUE)
  if( !is.null(g[['gval']]) ) g_new[['gval']][idx_sub] = g[['gval']]
  return(g_new)
}


#' Check vectorized grid data for non-NA points that form a complete sub-grid
#'
#' If a gridded data-set `g_obs` has missing values (NAs), but the set of non-NA points
#' form a complete sub-grid, this function finds its grid lines, resolution scaling factor,
#' and dimensions. If no eligible sub-grids are found, the function returns NULL.
#'
#' A sub-grid is only eligible if it contains all of the non-NA points in `g_obs` and none
#' of the NAs; eg if a single point missing from the sub-grid, or a single non-NA point lies
#' outside the sub-grid, the function will fail to detect the sub-grid and return NULL. If no
#' points are NA, the function returns indices for the full grid.
#'
#' In the special case that `g_obs` is a logical vector, it should indicate the the non-NA
#' locations in a grid with dimensions `gdim`. Otherwise, grid dimensions are extracted
#' from `g_obs`, overriding any argument to `gdim`.
#'
#' @param g_obs logical vector, or any other object accepted by `bk`
#' @param gdim integer vector, the grid dimensions (ny, nx)
#' @param g_out logical, indicates to return a grid list object
#'
#' @return NULL or list of information about the location and spacing of the sub-grid
#' within `g` (see details)
#' @export
#'
#' @examples
#'
#' # define a grid and example data
#' gdim = c(50, 53)
#' g_bare = bk(gdim)
#' gval = bk_sim(g_bare, modifyList(bk_pars(g), list(eps=1e-12)))
#' g_obs = modifyList(g_bare, list(gval=gval))
#' bk_plot(g_obs)
#'
#' # define a super-grid containing the original data and make sure we can find it
#' g_obs_big = bk_rescale(g_obs, down=3)
#' bk_plot(g_obs_big)
#' str(bk_sub_find(g_obs_big))
#'
#' # define a smaller sub-grid at random
#' spacing = sapply(floor(gdim/10), function(x) 1 + sample.int(x, 1))
#' gdim_sg = sapply(floor( (gdim - 1) / spacing), function(x) sample.int(x, 1))
#' ij_first = sapply(gdim - ( spacing * gdim_sg ), function(x) sample.int(x, 1))
#'
#' # find index of sub-grid lines and vectorized index of points
#' ij_sg = Map(\(idx, r, n) seq(idx, by=r, length.out=n), idx=ij_first, r=spacing, n=gdim_sg)
#' is_sg = bk_sub_idx(gdim, ij_sg, idx=FALSE)
#'
#' # assign values to the sub-grid points
#' g_obs_sub = g_bare
#' g_obs_sub$gval[is_sg] = gval[is_sg]
#' bk_plot(g_obs)
#' bk_plot(g_obs_sub, zlab='sub-grid')
#'
#' # call the function and check for expected results
#' subgrid_result = bk_sub_find(g_obs_sub)
#' all.equal(unname(subgrid_result$gdim), gdim_sg)
#' all.equal(unname(subgrid_result$ij), ij_sg)
#'
#' # sub grids with side length 1 have no spacing defined along that dimension
#' spacing[gdim_sg==1] = NA
#' all.equal(unname(subgrid_result$res_scale), spacing)
#'
#' # or call on the vector and supply gdim separately
#' identical(subgrid_result, bk_sub_find(g_obs_sub$gval, g_obs_sub$gdim))
#' identical(subgrid_result, bk_sub_find(!is.na(g_obs_sub$gval), g_obs_sub$gdim))
#'
bk_sub_find = function(g_obs, gdim=NULL)
{
  # handle vector input
  if( is.logical(g_obs) & is.vector(g_obs) )
  {
    # logical vectors interpreted as indicating non-NAs
    if( anyNA(g_obs) ) stop('g_obs vector of logical class cannot have NAs')
    gdim = stats::setNames(as.integer(gdim), c('y', 'x'))

  } else {

    # open as blitzKrig list object
    if( is.vector(g_obs) & !is.list(g_obs) ) g_obs = list(gval=g_obs, gdim=gdim)
    g_result = bk(g_obs)

    # process only the first column of multi-layer input
    if( is.matrix(g_result[['gval']]) ) g_result[['gval']] = as.vector(g_obs[['gval']][,1])
    g_obs = !is.na( g_result[['gval']] )
    gdim = g_result[['gdim']]

    # handle sparse indexing
    if( !is.null(g_result[['idx_obs']]) )
    {
      # decompress and replace NAs with FALSE
      g_obs = g_obs[ g_result[['idx_obs']] ]
      g_obs[is.na(g_obs)] = FALSE
    }
  }

  # checks for valid arguments
  n_obs = sum(g_obs)
  if(n_obs < 2) return(NULL)
  if( is.null(gdim) ) stop('full grid dimensions gdim must be supplied when g_obs is a vector')
  msg_len = paste('Expected', prod(gdim), 'but got', length(g_obs))
  if( prod(gdim) != length(g_obs)) stop(paste('input g_obs was the wrong length.', msg_len))

  # need this to get indices of first, second, and last elements in sub-grid
  idx_obs = which(g_obs)

  # find the dimensions of the smallest sub-grid enclosing all observed points
  ij_bbox = bk_vec2mat(c(idx_obs[1], idx_obs[n_obs]), gdim)
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
  skip_j = as.integer( ifelse(nj > 1, (gdim_bbox[['j']] - nj) / (nj - 1), NA) )

  # sub-grid resolution scaling factor and dimensions
  nm_dim = c('y', 'x')
  res_ij = setNames(1L + c(skip_i, skip_j), nm_dim)
  res_scale = 1L / res_ij
  gdim_sub = setNames(as.integer(c(ni, nj)), nm_dim)

  # sub-grid line indices
  ij_obs = Map(\(idx, r, n) seq(idx, by=r, length.out=n), idx=ij_bbox[1,], r=res_ij, n=gdim_sub)
  ij_na = sapply(ij_obs, anyNA)
  ij_obs[ij_na] = as.list(ij_bbox[1,])[ ij_na ]

  # final existence check
  idx_sub = bk_sub_idx(gdim, ij_obs)
  if( !all( g_obs[idx_sub] ) ) return(NULL)

  # return sub-grid info in a list
  gdim_result = setNames(sapply(ij_obs, length), c('y', 'x'))
  sub_result = list(ij=setNames(ij_obs, nm_dim), res_scale=res_ij, gdim=gdim_result)
  return(sub_result)
}


#' Return coordinates of a grid of points in column-vectorized order
#'
#' Expands a set of y and x grid line numbers in the column-vectorized order returned
#' by `bk`.
#'
#' This is similar to `base::expand.grid` but with the first dimension (y) descending
#' instead of ascending.
#'
#' `out='sf'` returns an `sf` simple features object containing points in the same order,
#' with data (if any) copied from `g$gval` into column 'gval'. Note that `prod(g$gdim)`
#' points are created, which can be slow for large grids.
#'
#' @param g any object accepted by `bk`
#' @param out character indicating return value type, either 'list', 'matrix', or 'sf'
#' @param corners logical, indicates to only return the corner points
#'
#' @return a matrix, list, or sf POINT collection, in column vectorized order
#' @export
#'
#' @examples
#' gdim = c(5,3)
#' g_example = bk(list(gdim=gdim, gres=c(0.5, 0.7), gval=seq(prod(gdim))))
#' bk_coords(g_example)
#' bk_coords(g_example, out='list')
#'
#' # corner points
#' bk_coords(g_example, corner=TRUE)
#' bk_coords(g_example, corner=TRUE, out='list')
#'
#' # sf output type
#' if( requireNamespace('sf') ) {
#'
#' # make the points
#' sf_coords = bk_coords(g_example, out='sf')
#'
#' # data are copied to variable 'gval'
#' plot(sf_coords, pch=16)
#'
#' }
#'
bk_coords = function(g, out='matrix', corner=FALSE, quiet=FALSE)
{
  # unpack input and slice multi-layer input
  g = bk(g)
  if( !is.null(g[['idx_grid']]) )
  {
    # keep only the first layer
    g[['gval']] = as.vector(g[['gval']][g[['idx_grid']], 1L])
    g[['idx_grid']] = NULL
  }

  # take subset of corner points if requested
  if(corner)
  {
    g[['gdim']] = c(y=2L, x=2L)
    g[['gyx']] = lapply(g[['gyx']], range)
  }

  # sort grid lines
  g[['gyx']][['y']] = sort(g[['gyx']][['y']], decreasing=TRUE)
  g[['gyx']][['x']] = sort(g[['gyx']][['x']], decreasing=FALSE)

  # compute the coordinates
  ij = bk_vec2mat(seq(prod(g[['gdim']])), g[['gdim']], out='list')
  out_list = stats::setNames(Map(function(gl, idx) gl[idx], g[['gyx']], idx=ij), c('y', 'x'))

  # return as list
  if( out == 'list' ) return(out_list)

  # return as matrix
  out_mat = do.call(cbind, out_list)
  if( out == 'matrix' ) return( out_mat )

  # sf return mode
  if( !startsWith(out, 'sf') ) stop('Argument `out` must be either "list", "matrix", or "sf"')
  sf_loaded = requireNamespace('sf', quietly=TRUE)
  if( !sf_loaded ) stop('sf package not loaded. Try library(sf)')
  if( !quiet ) cat(paste('processing', prod(g[['gdim']]), 'grid points...\n'))

  # create the points object
  if( is.null(g[['crs']]) ) g[['crs']] = ''
  sf_out = sf::st_as_sf(as.data.frame(out_mat), coords=c('x', 'y'), crs=g[['crs']])

  # copy any data and return
  if( !is.null(g[['gval']]) & !corner ) sf_out['gval'] = g[['gval']]
  return(sf_out[2:1])
}

