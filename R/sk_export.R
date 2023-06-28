# sk_export.R
# Dean Koch, October 2022
# Functions for exporting grids

#' Convert "sk" grid to SpatRaster
#'
#' @param g any object accepted or returned by `sk`
#' @param template character or RasterLayer/SpatRaster to set output type
#'
#' Converts a vector or matrix to a SpatRaster or RasterLayer. Multi-layer outputs are
#' supported for terra but not raster.
#'
#' @return a RasterLayer or SpatRaster containing the data from `g` (or a sub-grid)
#'
#' @family exporting functions
#' @seealso sk
#' @export
#'
#' @examples
#' if( requireNamespace('raster') ) {
#'
#' # open example file as RasterLayer
#' r_path = system.file('ex/logo.tif', package='terra')
#' r = raster::raster(r_path, band=1)
#' g = sk(r)
#'
#' # convert back to RasterLayer and compare
#' r_from_g = sk_export(g, 'raster')
#' print(r)
#' print(r_from_g)
#'
#' # NOTE: layer name, band number, and various other metadata are lost
#' all.equal(r_from_g, r)
#'
#' }
#'
#' # same with terra
#' if( requireNamespace('terra') ) {
#'
#' # convert all layers
#' r = terra::rast(r_path)
#' g = sk(r)
#' r_from_g = sk_export(g)
#'
#' # NOTE: various metadata are lost
#' all.equal(r_from_g, r)
#'
#' }
#'
sk_export = function(g, template='terra')
{
  # stop with an error message if raster/terra package is unavailable
  pkg_check = stats::setNames(nm=c('raster', 'terra'))
  pkg_msg = paste(pkg_check, collapse=' or ')
  msg_dep = paste(pkg_msg, 'package must be loaded first. Try `install.packages(terra)`')
  is_loaded = sapply(pkg_check, function(pkg) requireNamespace(pkg, quietly=TRUE))
  if( !any(is_loaded) ) stop(msg_dep)

  # load the input as snapKrig list and set empty CRS if not supplied
  g = sk(g)
  is_multi = is.matrix(g[['gval']])
  n_layer = ifelse(is_multi, ncol(g[['gval']]), 1L)
  if( is.null(g[['crs']]) ) g[['crs']] = ''

  # extract grid cell boundaries as defined in raster/terra
  yx_bbox = Map(\(g, s) range(g) + (c(-1,1) * s/2), g=g[['gyx']], s=g[['gres']])

  # handle grid objects as templates
  if( !is.character(template) )
  {
    template = class(g)[1]

    # set template class name
    if( inherits(g, 'SpatRaster') ) template = 'terra'
    if( inherits(g, c('RasterLayer', 'RasterStack')) ) template = 'raster'
  }

  # terra is preferred when available
  if( template == 'terra' )
  {
    # initialize raster
    r_ext = terra::ext(do.call(c, rev(yx_bbox)))
    r_out = terra::rast(extent=r_ext, resolution=rev(g[['gres']]), crs=g[['crs']], nlyr=n_layer)

    # reorder values then write to cells
    idx_reorder = t(matrix(seq_along(g), dim(g)))
    idx_write = seq_along(idx_reorder)
    if(!is_multi) return( terra::init(r_out, as.matrix(g)[idx_reorder]) )
    for(i in seq(n_layer)) terra::set.values(r_out,
                                             cells = idx_write,
                                             values = g[c(idx_reorder),i],
                                             layer = i)
    return(r_out)
  }

  # revert to raster if requested
  if( template == 'raster' )
  {
    # initialize raster
    r_ext = raster::extent(do.call(c, rev(yx_bbox)))
    r_out = raster::raster(ext=r_ext, resolution=rev(g[['gres']]), crs=g[['crs']])

    # warn about lack of support for RasterStack output
    if(n_layer > 1) warning('only the first layer was exported')

    # write values to cells
    if( !is.null(g[['gval']]) ) r_out = raster::setValues(r_out, as.matrix(g))
    return(r_out)
  }

  # error if we didn't get an expected class
  stop(paste('unrecognized template class', template))
}


#' Snap a set of points to a "sk" grid
#'
#' Maps the input points in `from` to the closest grid points in the lattice of which
#' `g` is a sub-grid. In cases of duplicate mappings, the function returns the first
#' matches only.
#'
#' `from` can be a geometry collection from packages `sf` or `sp`, or a matrix or list
#' of y and x coordinates. When `from` is a matrix, its first two columns should be the
#' y and x coordinates (in that order), and the (optional) third column should be the
#' data. When `from` is a list, the function expects (two or three) vectors of equal
#' length, ordered as above.
#'
#' When `from` is a geometry collection with a coordinates reference system (CRS) string,
#' points are first transformed to the CRS of `g`. If one or both of `from` and `g` are
#' missing a CRS definition, the function assumes the same one is shared in both.
#'
#' `g` can be a raster geometry object (such as SpatRaster), in which case the function
#' behaves like `terra::rasterize`, or an sk grid object. It can also be a matrix (supplying
#' dimensions) or a list containing either `gdim` or`gres`, from which an appropriately
#' spaced set of grid lines is derived, centered under the bounding box of the points.
#' If `g` is not supplied, it is automatically set to equal `nrow(from)`, so that there
#' there is one grid line along each dimension for each input point.
#'
#' `crop_from` and `crop_g` control the extent of the output grid. If both are `FALSE`
#' (the default) the function returns the smallest regular grid containing both `g`
#' and the snapped `from` points. If `crop_from=TRUE` and `crop_g=FALSE` the output
#' grid will match `g` exactly. If `crop_from=FALSE` and `crop_g=TRUE` the output
#' grid will include all snapped points, and possibly omit some or all of `g`. And if
#' both are `TRUE`, the output grid encloses the intersection of the points with the
#' bounding box of `g`.
#'
#' @param from matrix, data frame, or points object from `sp` or `sf`, the source points
#' @param g any object accepted or returned by `sk`, the destination grid
#' @param crop_from logical, indicating to omit points not overlying `g`.
#' @param crop_g logical, indicating to trim `g` to the extent of `from`.
#' @param quiet logical, suppresses console output
#'
#' @return sk object, a grid containing the snapped points. These are assigned
#' the corresponding data value in `from`, or if  `from` has no data, an integer mapping
#' to the points in `from`. Un-mapped grid points are set to NA.
#' @export
#' @family sk constructors
#' @seealso sk sk_coords
#'
#' @examples
#'
#' # functions to scale arbitrary inverval to (1, 2,... 100) and make color palettes
#' num_to_cent = function(x) 1L + floor(99*( x-min(x) ) / diff(range(x)))
#' my_pal = function(x) grDevices::hcl.colors(x, 'Spectral', rev=TRUE)
#' my_col = function(x) my_pal(1e2)[ num_to_cent(x) ]
#'
#' # create a grid object
#' gdim = c(40, 30)
#' g = sk(gdim=gdim, gres=1.1)
#'
#' # randomly position points within bounding box of g
#' n_pts = 10
#' from = lapply(g$gyx, function(yx) runif(n_pts, min(yx), max(yx)) )
#'
#' # translate away from g (no overlap is required)
#' from[['y']] = from[['y']] + 5
#' from[['x']] = from[['x']] + 15
#'
#' # add example data values and plot
#' from[['z']] = stats::rnorm(length(from[['y']]))
#' plot(g, reset=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col(from[['z']]))
#' graphics::points(from[c('x', 'y')])
#'
#' # snap only the points overlying the input grid
#' g_snap = sk_snap(from, g, crop_from=TRUE)
#' plot(g_snap, col_grid='black', reset=FALSE, leg=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col(from[['z']]))
#' graphics::points(from[c('x', 'y')])
#'
#' # snap all points to grid extension (default settings)
#' g_snap = sk_snap(from, g, crop_from=FALSE, crop_g=FALSE)
#' plot(g_snap, col_grid='black', reset=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col(from[['z']]))
#' graphics::points(from[c('x', 'y')])
#'
#' # find smallest subgrid enclosing all snapped grid points
#' g_snap = sk_snap(from, g, crop_g=TRUE)
#' plot(g_snap, col_grid='black', reset=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col(from[['z']]))
#' graphics::points(from[c('x', 'y')])
#'
#' # create a new grid of different resolution enclosing all input points
#' g_snap = sk_snap(from, g=list(gres=c(0.5, 0.5)))
#' plot(g_snap, reset=FALSE, col_grid='black')
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col(from[['z']]))
#' graphics::points(from[c('x', 'y')])
#'
#' if( requireNamespace('sf') ) {
#'
#' # a different example, snapping mis-aligned subgrid
#' g_pts = sk(list(gdim=c(15, 8), gres=1.7), vals=FALSE)
#' g_pts[['gyx']][['y']] = g_pts[['gyx']][['y']] + 5
#' g_pts[['gyx']][['x']] = g_pts[['gyx']][['x']] + 5
#' from = sk_coords(g_pts, out='list')
#'
#' # convert to sf
#' eg_sfc = sf::st_geometry(sk_coords(g_pts, out='sf'))
#' plot(g, reset=FALSE)
#' plot(eg_sfc, add=TRUE)
#'
#' # generate example data and plot
#' eg_sf = sf::st_sf(data.frame(z=stats::rnorm(length(g_pts))), geometry=eg_sfc)
#' plot(g, reset=FALSE)
#' plot(eg_sf, pch=16, add=TRUE, pal=my_pal)
#' plot(eg_sfc, add=TRUE)
#'
#' # snap points
#' g_snap = sk_snap(from=eg_sf, g)
#' plot(g_snap, reset=FALSE, col_grid='black')
#' plot(eg_sf, pch=16, add=TRUE, pal=my_pal)
#' plot(eg_sfc, add=TRUE)
#'
#' # snapping points without data produces the mapping (non-NA values index "from")
#' g_snap = sk_snap(from=eg_sfc, g)
#' plot(g_snap, ij=TRUE, reset=FALSE, col_grid='black')
#' plot(eg_sfc, add=TRUE)
#'
#' # with crop_g=TRUE)
#' g_snap = sk_snap(from=eg_sfc, g, crop_g=TRUE)
#' plot(g_snap, reset=FALSE, col_grid='black')
#' plot(eg_sfc, add=TRUE)
#'
#' # test with sp class
#' eg_sp = as(eg_sf,'Spatial')
#' g_snap = sk_snap(from=eg_sp, g)
#' plot(g_snap, reset=FALSE, col_grid='black')
#' plot(eg_sf, pch=16, add=TRUE, pal=my_pal)
#' plot(eg_sfc, add=TRUE)
#'
#' }
#'
sk_snap = function(from, g=nrow(from), crop_from=FALSE, crop_g=FALSE, quiet=FALSE)
{
  # set up more appropriate grid lines later when input g doesn't specify them
  gres_input = NULL
  if( is.matrix(g) ) g = list(gdim=dim(g))
  auto_gyx = (!is.list(g) & is.vector(g)) | ( is.list(g) & !('gyx' %in% names(g)) )

  # set up grid dimensions later when only gres specified (set dummy dimensions now)
  auto_gdim = ifelse(auto_gyx & is.list(g), ('gres' %in% names(g)) & !('gdim' %in% names(g)), FALSE)
  if(auto_gdim & is.list(g) ) g = utils::modifyList(g, list(gdim=c(y=2L, x=2L)))

  # unpack g as snapKrig grid list but don't copy data values
  g = sk(g, vals=FALSE)
  to_crs = g[['crs']]
  vals = NULL

  # expected point object classes and names
  nm_yx = c('y', 'x')
  is_sf = inherits(from, c('sf','sfc', 'sfg'))
  is_sp = inherits(from, c('SpatialPoints', 'SpatialPointsDataFrame'))

  # get coordinate(s) from sf objects as matrix/vector
  if(is_sf)
  {
    # transform to destination coordinate system, if supplied (otherwise assume same as source)
    from_crs = sf::st_crs(from)[['wkt']]
    is_geo = !is.null(to_crs) & !is.na(from_crs)
    if( is_geo ) { from = sf::st_transform(from, to_crs) } else { to_crs = from_crs }

    # copy values if possible
    if( length(names(from)) > 1 ) vals = unlist(sf::st_drop_geometry(from)[,1])

    # reverse column order to get y, x coordinates as matrix
    from = apply(sf::st_coordinates(from)[,2:1], 2L, identity)
  }

  # the same for sp objects
  if(is_sp)
  {
    from_crs = sp::wkt(from)
    is_geo = !is.null(to_crs) & !is.null(from_crs)
    if( is_geo ) { from = sp::spTransform(from, to_crs) } else { to_crs = from_crs }

    # grab first column if there are any
    if( length(names(from)) > 0 ) vals = unlist(from[[names(from)]])
    from = apply(sp::coordinates(from)[,2:1], 2L, identity)
  }

  # convert (length-2) vector, matrix, and dataframe to list
  if( is.data.frame(from) ) from = as.matrix(from)
  if( is.list(from) ) from = do.call(cbind, from)
  if( is.vector(from) ) from = as.matrix(from, nrow=1)
  if( is.matrix(from) ) from = apply(from, 2, identity, simplify=FALSE)

  # crop input points if requested
  if( crop_from )
  {
    omit_y = ( from[[1]] < min(g[['gyx']][['y']]) ) | ( from[[1]] > max(g[['gyx']][['y']]) )
    omit_x = ( from[[2]] < min(g[['gyx']][['x']]) ) | ( from[[2]] > max(g[['gyx']][['x']]) )
    omit_from = omit_y | omit_x
    from = lapply(from, function(x) x[!omit_from])
    vals = vals[!omit_from]
  }

  # copy data and find bounding box
  if( length(from) > 2 ) vals = from[[3]]
  from_yx = stats::setNames(from[1:2], nm_yx)
  from_n = length(from_yx[['x']])
  from_bds = lapply(list(min, max), function(f) sapply(from_yx, f))

  # automatically set grid lines to enclose all input points
  if(auto_gyx)
  {
    # automatically set grid dimensions based on requested resolution
    if(auto_gdim)
    {
      # given gres, compute number of grid lines neeeded to span the bounding box
      g[['gdim']] = 1L + sapply(from_yx, function(yx) diff(range(yx))) %/% g[['gres']]

      # find offsets that center the grid lines under point bounding box
      auto_pad = ( sapply(from_yx, function(yx) diff(range(yx))) %% g[['gres']] ) / 2
      gyx = Map(function(yx, p, r, n) seq(yx + p, by=r, length.out=n),
                yx=from_bds[[1]], p=auto_pad, r=g[['gres']], n=g[['gdim']])

      # make grid object
      g = sk(list(gdim=g[['gdim']], gyx=gyx, gres=g[['gres']]))

    } else {

      # place grid lines to coincide with bounding box edge points, then recompute gres
      gyx = Map(function(yx, n) seq(min(yx), max(yx), length.out=n), yx=from_yx, n=g[['gdim']])
      gres = sapply(gyx, function(yx) diff(yx[1:2]))
      g = sk(list(gdim=g[['gdim']], gyx=gyx, gres=gres))

    }
  }

  # in this case the output configuration should match g exactly
  if(crop_from & !crop_g)
  {
    # copy grid info from g
    to_yx = g[['gyx']]
    gdim_out = g[['gdim']]

  } else {

    # construct the grid line locations based on from

    # find bounding box of destination grid template and reshape to list of min and max
    g_bbox = lapply(sk_coords(g, out='list', corner=TRUE, quiet=TRUE), range)
    g_bds = lapply(list(min, max), function(f) sapply(g_bbox, f))

    # find the offsets between the two bounding boxes
    to_pad = Map(function(a, b) abs((a - b) %% g[['gres']]), a=from_bds, b=g_bds)
    to_min = from_bds[[1]] - to_pad[[1]] + as.integer(to_pad[[1]] > (g[['gres']]/2)) * g[['gres']]
    to_max = from_bds[[2]] - to_pad[[2]] + as.integer(to_pad[[2]] > (g[['gres']]/2)) * g[['gres']]

    # crop the grid to the extent of g
    if( !crop_g )
    {
      # if not cropping, extend these grid lines to include all of g_bbox
      to_min = pmin(to_min, sapply(g_bbox, min))
      to_max = pmax(to_max, sapply(g_bbox, max))
    }

    # compute new grid line locations
    to_yx = Map(function(a, b, r) seq(a, b, by=r), a=to_min, b=to_max, r=g[['gres']])
    gdim_out = sapply(to_yx, length)
  }


  # find cross-distance matrices for point coordinates and grid lines
  d_yx_all = Map(function(a, b) outer(a, b, '-')^2, a=from_yx, b=to_yx)

  # find minimum distance mappings in each dimension
  ij_min = stats::setNames(lapply(d_yx_all, function(d) max.col(-d, ties.method='f')), c('i', 'j'))

  # print maximum snapping distance
  if(!quiet)
  {
    # find component squared distances then print the maximum
    dy = mapply(function(i, j) d_yx_all[['y']][i, j], i=seq_along(ij_min[['i']]), j=ij_min[['i']])
    dx = mapply(function(i, j) d_yx_all[['x']][i, j], i=seq_along(ij_min[['j']]), j=ij_min[['j']])
    cat(paste('maximum snapping distance:', max(sqrt(dx+dy)), '\n'))
  }

  # reflect indexing of vertical axis for snapKrig
  ij_min[['i']] = gdim_out['y'] + 1L - ij_min[['i']]
  to_idx = sk_mat2vec(ij_min, gdim_out)

  # handle multiple points mapping to a single grid-point
  is_dupe = duplicated(to_idx)
  if( !quiet & any(is_dupe) ) cat( paste('omitting', sum(is_dupe), 'duplicate mapping(s)\n') )

  # match magic to get NAs at unassigned grid points
  from_idx = match(seq(prod(gdim_out)), to_idx)
  #if( is.null(vals) ) { gval = from_idx } else { gval = vals[from_idx] }
  return(sk(gyx = to_yx,
            crs = to_crs,
            gval = if(is.null(vals)) from_idx else vals[from_idx]))

  # if( is.null(vals) ) { g_out[['gval']] = from_idx } else { g_out[['gval']] = vals[from_idx] }
  # return(sk_validate(g_out))
}


#' Return coordinates of a grid of points in column-vectorized order
#'
#' Expands a set of y and x grid line numbers in the column-vectorized order returned
#' by `sk`. This is similar to `base::expand.grid` but with the first dimension (y)
#' descending instead of ascending.
#'
#' `out='sf'` returns an `sf` simple features object containing points in the same order,
#' with data (if any) copied from `g[['gval']]` into column 'gval'. Note that `length(g)`
#' points are created, which can be slow for large grids.
#'
#' If `na_omit` is `TRUE` the function omits points with `NA` data (in `gval`) and only
#' returns the coordinates for observations. This argument is ignored when `corners=TRUE`
#' (which always returns the four corner points) or when the grid contains no observations
#' (all points returned).
#'
#' @param g any object accepted by `sk`
#' @param out character indicating return value type, either 'list', 'matrix', or 'sf'
#' @param na_omit logical, indicates to return only locations of observed points
#' @param corner logical, indicates to return only the four corner points
#' @param quiet logical, suppresses console output
#'
#' @return a matrix, list, or sf point collection in column vectorized order
#' @export
#' @family exporting functions
#' @seealso sk sk_snap base::expand.grid
#'
#' @examples
#' gdim = c(5,3)
#' g_complete = sk(gdim=gdim, gres=c(0.5, 0.7), gval=seq(prod(gdim)))
#' sk_coords(g_complete)
#' sk_coords(g_complete, out='list')
#'
#' # missing data example
#' idx_obs =  sort(sample.int(length(g_complete), 5))
#' g = sk(gdim=gdim, gres=c(0.5, 0.7))
#' g[idx_obs] = g_complete[idx_obs]
#' all.equal(sk_coords(g, na_omit=TRUE), sk_coords(g_complete)[idx_obs,])
#'
#' # corner points
#' sk_coords(g, corner=TRUE)
#' sk_coords(g, corner=TRUE, out='list')
#'
#' # repeat with multi-layer example
#' g_multi = sk(utils::modifyList(g, list(gval = cbind(g[], 2*g[]))))
#' all.equal(sk_coords(g_multi, na_omit=TRUE), sk_coords(g_complete)[idx_obs,])
#' sk_coords(g_multi, corner=TRUE)
#'
#' # sf output type
#' if( requireNamespace('sf') ) {
#'
#' # gather all points but don't copy data
#' sf_coords_all = sk_coords(sk(g, vals=FALSE), out='sf')
#'
#' # extract non-NA data
#' sf_coords = sk_coords(g, out='sf', na_omit=TRUE)
#'
#' # plot everything together
#' plot(g, reset=FALSE)
#' plot(sf_coords_all, add=TRUE)
#' plot(sf_coords, pch=16, cex=2, add=TRUE)
#'
#' }
#'
sk_coords = function(g, out='matrix', corner=FALSE, na_omit=FALSE, quiet=FALSE)
{
  # unpack input and check for empty grids
  g = sk(g)
  is_empty = !any(g[['is_obs']])
  if(is_empty) na_omit = FALSE

  # take subset of corner points if requested
  if(corner)
  {
    g[['gdim']] = c(y=2L, x=2L)
    g[['gyx']] = lapply(g[['gyx']], range)
    g[['gres']] = lapply(g[['gyx']], diff)
    na_omit = FALSE
  }

  # create vectorized index of points to return (if all locations are NA, return all)
  idx_get = if(na_omit) { which(!is.na(g)) } else { seq_along(g) }
  if( length(idx_get) == 0 )
  {
    idx_get = seq_along(g)
    is_empty = TRUE
  }

  # copy data values at these locations
  values_out = NULL
  if( !is_empty )
  {
   # for vector gval, get subset of values using normal extract function
    if( !is.matrix(g[['gval']]) ) { values_out = g[idx_get] } else {

      # grab first layer of matrix gval directly (it is already subsetted)
      values_out = as.vector(g[['gval']][, 1L])
    }
  }

  # sort grid lines positions so they match column vectorization
  gyx = Map(function(gl, dec) sort(gl, decreasing=dec), gl=g[['gyx']], dec=c(y=TRUE, x=FALSE))

  # console message
  if( !quiet )
  {
    n_get = length(idx_get)
    msg_get = ifelse(length(g) == n_get, '', paste(n_get, 'of'))
    cat(paste('processing', paste0(msg_get, length(g)), 'grid points...\n'))
  }

  # compute the coordinates
  ij = sk_vec2mat(idx_get, dim(g), out='list')
  coords_out = Map(function(gl, idx) gl[idx], gl=gyx, idx=ij)

  # return as list
  if( out == 'list' ) return(coords_out)

  # return as matrix
  coords_out_mat = do.call(cbind, coords_out)
  if( out == 'matrix' ) return(coords_out_mat)

  # handle invalid `out` arguments
  if( !startsWith(out, 'sf') ) stop('Argument `out` must be either "list", "matrix", or "sf"')
  sf_loaded = requireNamespace('sf', quietly=TRUE)
  if( !sf_loaded ) stop('sf package not loaded. Try library(sf)')

  # create the sf points object
  if( is.null(g[['crs']]) ) g[['crs']] = ''
  sf_out = sf::st_as_sf(as.data.frame(coords_out_mat), coords=c('x', 'y'), crs=g[['crs']])

  # copy any data to it and return - making sure geometry column is second
  if( !is.null(values_out) )
  {
    sf_out['gval'] = values_out
    return(sf_out[2:1])
  }

  return(sf_out)
}


