# bk_grid.R
# Dean Koch, 2022
# Functions for building grids out of various inputs


#' Make a blitzKrig grid list object
#'
#' Constructs a list representing a 2-dimensional regular spatial grid of data
#'
#' This function constructs blitzKrig grid list objects, accepting 'RasterLayer' and 'RasterStack'
#' inputs from the `raster` package, 'SpatRaster' objects from `terra`, as well as any
#' non-complex matrix, or a list containing the vectorization of one. Empty grids can be
#' initialized by specifying dimensions (and or setting `vals=FALSE`)
#'
#' The function returns a list with the following 3-6 elements:
#'
#' * gval: the data (if any) in column-major order with y descending, x ascending
#' * crs: character, the WKT string (if available) describing coordinate reference system
#' * idx_grid: integer vector, mapping rows of matrix `gval` to points on the grid
#' * gyx: a list containing the coordinates of the y and x grid lines in vectors `y` and `x`
#' * gres: the (numeric) y and x distances between grid lines
#' * gdim: the (integer) number of y and x grid lines
#'
#' The last three elements are required to define a valid `blitzKrig` grid list object. Note that
#' regular grids must have equally spaced grid lines in `gyx`.
#'
#' The input `g` can itself be a list containing some/all of these elements (including at least
#' one of 'gdim' or 'gyx'), and the function will fill in missing entries wherever possible:
#' If 'gval' is missing in the single-layer case, the function sets NAs; If 'res' is missing,
#' it is computed from the first two grid lines in 'gyx'. If 'gyx' is missing, it is assigned
#' the sequence `1:n` (scaled by 'res', if available) for each `n` in 'gdim'; and if 'gdim'
#' is missing, it is set to equal the number of grid lines specified in (each vector of) 'gyx'.
#'
#' `gval` can be a vector, for single-layer grids, or a matrix whose columns are a set of grid
#' layers. In the matrix case, `gval` stores the observed data only, with NAs indicating by the
#' mapping `idx_grid`. This mapping is assumed to be the same in all layers, but is only computed
#' for the first layer. If a point is missing from one layer, it should be missing from all layers.
#'
#' `idx_grid` is a length `prod(gdim)` integer vector with NAs for unobserved points, and
#' otherwise the row number (in `gval`) of the observed point. These non-NA values must
#' comprise `seq(nrow(gval))` (ie all rows must be mapped), but they can be in any order.
#' If `gval` is a matrix but `idx_grid` is missing, it is computed automatically (from the
#' first column); and if `idx_grid` is supplied, but `gval` is a vector, it coerced to a 1-column
#' matrix.
#'
#' Scalar inputs to 'gdim', 'gres' are duplicated for both dimensions, and for convenience
#' 'gdim' can be specified directly in `g` to initialize a simple grid; For example the call
#' `bk_grid(list(gdim=c(5,5)))` can be simplified to `bk_grid(list(gdim=5))` or
#' `bk_grid(5)`.
#'
#' @param g raster, matrix, numeric vector, or list (see details)
#' @param vals logical indicating to include the data vector 'gval' in return list
#'
#' @return named list containing at least 'gyx', 'gres', and 'gdim' (see details)
#' @export
#'
#' @examples
#'
#' # simple grid construction from dimensions
#' gdim = c(12, 10)
#' g = bk_grid(g=list(gdim=gdim))
#' str(g)
#' str(bk_grid(gdim, vals=FALSE))
#'
#' # pass result to bk_grid and get the same thing back
#' identical(g, bk_grid(g))
#'
#' # supply grid lines instead and get the same result
#' all.equal(g, bk_grid(g=list(gyx=lapply(gdim, function(x) seq(x)-1L))) )
#'
#' # display coordinates and grid line indices
#' bk_plot(g)
#' bk_plot(g, ij=TRUE)
#'
#' # gres argument is ignored if a non-conforming gyx is supplied
#' gres_new = c(3, 4)
#' bk_plot(bk_grid(g=list(gyx=lapply(gdim, seq), gres=gres_new)))
#' bk_plot(bk_grid(g=list(gdim=gdim, gres=gres_new)))
#'
#' # shorthand for square grids
#' all.equal(bk_grid(2), bk_grid(g=c(2,2)))
#'
#' # example with random data
#' gdim = c(25, 25)
#' yx = as.list(expand.grid(lapply(gdim, seq)))
#' eg_vec = as.numeric( yx[[2]] %% yx[[1]] )
#' eg_mat = matrix(eg_vec, gdim)
#' g = bk_grid(eg_mat)
#' bk_plot(g, ij=T, zlab='j mod i')
#'
#' # y is in descending order
#' bk_plot(g, xlab='x = j', ylab='y = 26 - i', zlab='j mod i')
#'
#' # data vectors should be in R's default matrix vectorization order
#' all.equal(eg_vec, as.vector(eg_mat))
#' all.equal(g, bk_grid(list(gdim=gdim, gval=as.vector(eg_mat))))
#'
#' # multi-layer example from matrix
#' n_pt = prod(gdim)
#' n_layer = 3
#' mat_multi = matrix(rnorm(n_pt*n_layer), n_pt, n_layer)
#' g_multi = bk_grid(list(gdim=gdim, gval=mat_multi))
#' str(g_multi)
#'
#' # repeat with missing data (note all columns must have consistent NA structure)
#' mat_multi[sample.int(n_pt, 0.5*n_pt),] = NA
#' g_multi_miss = bk_grid(list(gdim=gdim, gval=mat_multi))
#' str(g_multi_miss)
#'
#' # only observed data points are stored, idx_grid maps them to the full grid vector
#' max(abs( g_multi$gval - g_multi_miss$gval[g_multi_miss$idx_grid,] ), na.rm=TRUE)
#'
#' # vals=FALSE drops multi-layer information
#' bk_grid(g=list(gdim=gdim, gval=mat_multi), vals=FALSE)
#'
#' if( requireNamespace('raster') ) {
#'
#' # open example file as RasterLayer
#' r_path = system.file('external/rlogo.grd', package='raster')
#' r = raster::raster(r_path)
#'
#' # convert to blitzKrig list (notice only first layer was loaded by raster)
#' g = bk_grid(r)
#' str(g)
#' bk_plot(g)
#'
#' # open a RasterStack - gval becomes a matrix with layers in columns
#' r_multi = raster::stack(r_path)
#' g_multi = bk_grid(r_multi)
#' str(g_multi)
#' bk_plot(g_multi, layer=1)
#' bk_plot(g_multi, layer=2)
#' bk_plot(g_multi, layer=3)
#'
#' # repeat with terra
#' if( requireNamespace('terra') ) {
#'
#' # open example file as SpatRaster (all layers loaded by default)
#' r_multi = terra::rast(r_path)
#' g_multi = bk_grid(r_multi)
#' str(g_multi)
#' bk_plot(g_multi, layer=1)
#' bk_plot(g_multi, layer=2)
#' bk_plot(g_multi, layer=3)
#'
#' # open first layer only
#' g = bk_grid(r[[1]])
#' str(g)
#' bk_plot(g)
#'
#' }
#' }
bk_grid = function(g, vals=TRUE)
{
  # names for dimensional components and required entries of g, and a preferred order
  nm_dim = c('y', 'x')
  nm_g = c('gyx', 'gres', 'gdim')
  nm_order = c('gval', 'idx_grid', 'crs', nm_g)

  # handle raster objects
  is_terra = any(c('SpatRaster') %in% class(g))
  is_raster = any(c('RasterLayer', 'RasterStack') %in% class(g))
  if( is_terra | is_raster )
  {
    # terra and raster use another ordering for vectorization
    gdim = dim(g)[1:2] # order y, x
    if(vals) idx_reorder = matrix(seq(prod(gdim)), gdim, byrow=TRUE)

    # package-specific calls
    gval = NULL
    if(is_terra)
    {
      # terra class
      gcrs = terra::crs(g)
      gyx = list(y=terra::yFromRow(g, seq(gdim[1])), x=terra::xFromCol(g, seq(gdim[2])))
      gres = terra::res(g)[2:1] # order dy, dx
      if(vals)
      {
        # single layer gval is a vector
        n_layer = terra::nlyr(g)
        if(n_layer == 1) { gval = terra::values(g)[idx_reorder]  } else {

          # multi-layer gval is a matrix
          gval = terra::values(g)[idx_reorder,]
        }
      }

    } else {

      # raster class
      gcrs = raster::wkt(g)
      gyx = list(y=raster::yFromRow(g, seq(gdim[1])), x=raster::xFromCol(g, seq(gdim[2])))
      gres = raster::res(g)[2:1] # order dy, dx
      if(vals)
      {
        # single layer gval is a vector
        n_layer = raster::nlayers(g)
        if(n_layer == 1) { gval = raster::getValues(g)[idx_reorder] } else {

          # multi-layer gval is a matrix
          gval = raster::getValues(g)[idx_reorder,]
        }
      }
    }

    # sort both sets of grid lines into ascending order
    gyx = lapply(gyx, sort)

    # recursive call to validate and set names
    bk_grid(list(gval=gval, crs=gcrs, gyx=gyx, gres=gres, gdim=gdim), vals=vals)
    return( bk_grid(list(gval=gval, crs=gcrs, gyx=gyx, gres=gres, gdim=gdim), vals=vals) )
  }

  # handle matrix objects
  if( is.matrix(g) )
  {
    # set unit resolution by default
    g_out = list(gres=c(y=1, x=1), gdim=dim(g))
    g_out[['gyx']] = lapply(g_out[['gdim']], function(d) as.numeric(seq(d)-1L))

    # build named list and return
    g_out = stats::setNames(lapply(g_out[nm_g], function(r) stats::setNames(r, nm_dim)), nm_g)
    if(!vals) return(g_out)
    return( c(list(gval=as.vector(g)), g_out) )
  }

  # handle list objects
  if( is.list(g) )
  {
    # check for values
    is_gval = names(g) == 'gval'
    is_indexed = !is.null(g[['idx_grid']])
    if( is_indexed & !any(is_gval) ) stop('g$idx_grid supplied without g$gval')
    if( any(is_gval) )
    {
      # check for multi-layer input
      is_multi = is.matrix(g[['gval']])

      # create missing indexing vector when it is expected
      if( is_multi & !is_indexed )
      {
        # identify observed data in first layer and build an indexing vector from it
        is_obs_first = !is.na(g[['gval']][,1L])
        g[['idx_grid']] = match(seq(nrow(g[['gval']])), which(is_obs_first))

        # handle no observed data case
        if(length(g[['idx_grid']]) == 0) g[['idx_grid']] = rep(NA, nrow(g[['gval']]))

        # keep only the non-NA rows
        g[['gval']] = g[['gval']][is_obs_first,]
        is_indexed = TRUE
      }

      # vector gval overwritten with 1-column matrix when it is indexed
      if( !is_multi & is_indexed ) g[['gval']] = matrix(g[['gval']], ncol=1L)

      # if idx_grid is supplied, gval should have no NAs
      gval_ok = ifelse(is_indexed, !anyNA(g[['gval']]), TRUE)
      if(!gval_ok) stop('inconsistent pattern of NAs among layers')
    }

    # validate or compute gdim
    if( is.null(g[['gdim']]) )
    {
      # require gyx when gdim missing
      if( is.null(g[['gyx']]) ) stop('both gdim and gyx missing from input grid g')
      g[['gdim']] = sapply(g[['gyx']], length)

    } else {

      # coerce to integer and duplicate scalar input
      g[['gdim']] = as.integer(g[['gdim']])
      if( length(g[['gdim']]) == 1 ) g[['gdim']] = rep(g[['gdim']], 2)
    }

    # when grid resolution is missing calculate it from gyx or set up default
    if( is.null(g[['gres']]) )
    {
      # compute from gyx where available (or set unit default)
      g[['gres']] = c(1, 1)
      if( !is.null(g[['gyx']]) ) g[['gres']] = as.numeric(sapply(g[['gyx']], function(r) diff(r)[1]))

    } else {

      # duplicate scalar input
      if( length(g[['gres']]) == 1 ) g[['gres']] = rep(g[['gres']], 2)
    }

    # set up grid line positions if they're missing
    if( is.null(g[['gyx']]) )
    {
      g[['gyx']] = Map(function(d, r) as.numeric(r*(seq(d)-1L)),
                       d=g[['gdim']],
                       r=g[['gres']])
    } else {

      # TODO: figure out a better threshold or produce a more informative warning
      # check that the grid is regular
      g_spacing = lapply(g[['gyx']], function(r) diff(r) )
      is_regular = sapply(g_spacing, function(s) max(abs(s[1] - s)) < .Machine[['double.eps']] )
      #if( any(!is_regular) ) warning('grid line positions gyx were not regular')
    }

    # check number of dimensions
    if( !all(sapply(g[nm_g], length) == 2) ) stop('gdim, gres, and gyx must each have length 2')

    # initialize with NAs if no data supplied
    if( vals & !any(is_gval) ) g[['gval']] = rep(NA_real_, prod(g[['gdim']]))

    # set dimension names and order of output
    g[nm_g] = lapply(g[nm_g], function(r) stats::setNames(r, nm_dim))
    nm_known = names(g)[ match(nm_order, names(g)) ]
    nm_unknown = names(g)[ !(names(g) %in% nm_order) ]
    nm_out = c(nm_unknown, nm_known[!is.na(nm_known)])

    # finished if not returning grid values
    if(!vals)
    {
      nm_out = nm_out[ !(nm_out %in% c('idx_grid', 'gval')) ]
      return(g[nm_out])
    }

    # check that index (if any) has correct length
    index_ok = ifelse(is_indexed, prod(g[['gdim']]) == length(g[['idx_grid']]), TRUE)
    if(!index_ok) stop('indexing vector idx_grid must have length prod(gdim)')

    # compare number of grid points found in data matrix/vector, and total expected
    n_got = ifelse(is_indexed, nrow(g[['gval']]), length(g[['gval']]))
    n_expect = ifelse(is_indexed, sum(!is.na(g[['idx_grid']])), prod(g[['gdim']]))
    msg_gval = ifelse(is_indexed, 'number of rows in gval', 'length of gval')
    msg_gdim = ifelse(is_indexed, 'number of non-NAs in idx_grid', 'prod(gdim)')
    if(n_got != n_expect) stop(paste('mismatch between', msg_gdim, 'and', msg_gval))

    return(g[nm_out])
  }

  # handle numeric vectors
  if( is.vector(g) )
  {
    if( length(g) > 2 ) stop('numeric vector g must be of length 1 or 2')
    if( length(g) == 1 ) g = stats::setNames(rep(g, 2), nm_dim)
    return( bk_grid(g=list(gdim=as.integer(g)), vals) )
  }

  # unrecognized objects are returned unchanged with a warning
  warning('input g was not recognized')
  return(NA)
}


#' Convert column-vectorized grid to SpatRaster
#'
#' @param g any object accepted or returned by `bk_grid`
#' @param template character or RasterLayer/SpatRaster to set output type
#'
#' Converts a column-vectorized vector or matrix to a SpatRaster, or if terra is
#' unavailable, a RasterLayer.
#'
#' @return a RasterLayer or SpatRaster containing the data from `g` (or a sub-grid)
#' @export
#'
#' @examples
#'
#' if( requireNamespace('raster') ) {
#'
#' # open example file as RasterLayer
#' r_path = system.file('external/rlogo.grd', package='raster')
#' r = raster::raster(r_path, band=1)
#' g = bk_grid(r)
#'
#' # convert back to RasterLayer and compare
#' r_from_g = bk_export(g, 'raster')
#' print(r_from_g)
#' print(r)
#'
#' # layer name, band number, and various other metadata are lost
#' all.equal(r_from_g, r)
#'
#' # same with terra
#' if( requireNamespace('terra') ) {
#'
#' # convert all layers
#' r = terra::rast(r_path)
#' g = bk_grid(r)
#' r_from_g = bk_export(g)
#'
#' # various metadata are lost
#' all.equal(r_from_g, r)
#'
#' }
#' }
#'
bk_export = function(g, template='terra')
{
  # stop with an error message if raster/terra package is unavailable
  pkg_check = stats::setNames(nm=c('raster', 'terra'))
  pkg_msg = paste(pkg_check, collapse=' or ')
  msg_dep = paste(pkg_msg, 'package must be loaded first. Try `install.packages(terra)`')
  is_loaded = sapply(pkg_check, function(pkg) requireNamespace(pkg, quietly=TRUE))
  if( !any(is_loaded) ) stop(msg_dep)

  # load the input as blitzKrig list
  g = bk_grid(g)
  g[['crs']] = ifelse(is.null(g[['crs']]), '', g[['crs']])
  n_layer = ifelse(is.null(g[['idx_grid']]), sum(!is.null(g[['gval']])), ncol(g[['gval']]))

  # extract grid cell boundaries as defined in raster/terra
  yx_bbox = Map(\(g, s) range(g) + (c(-1,1) * s/2), g=g[['gyx']], s=g[['gres']])

  # handle grid objects as templates
  if( !is.character(template) )
  {
    # check if template is a sub-grid or super-grid of `g`?
    # g_template = bk_grid(template)
    # TODO: implement crop?

    # set template class name
    template = class(g)
    if( 'SpatRaster' %in% template ) template = 'terra'
    if( any(c('RasterLayer', 'RasterStack') %in% template ) ) template = 'raster'
    template = paste(template, collapse=', ')
  }

  # terra is preferred when available
  if( template == 'terra' )
  {
    g_ext = terra::ext(do.call(c, rev(yx_bbox)))
    r_out = terra::rast(extent=g_ext, resolution=rev(g[['gres']]), crs=g[['crs']], nlyr=n_layer)
    if( n_layer == 1 ) { r_out = terra::setValues(r_out, matrix(g[['gval']], g[['gdim']])) } else {

      # pad with NAs to recover full grid
      gval = g[['gval']][ g[['idx_grid']], ]

      # multi-layer assignments require a different vectorization ordering!
      gval_list = apply(gval, 2, function(x) t(matrix(x, g[['gdim']])), simplify=FALSE)
      r_out = terra::setValues(r_out, do.call(cbind, gval_list))
    }

    return(r_out)
  }

  # check for unknown class
  if( template == 'raster' )
  {
    # attempt to use raster if terra unavailable
    g_ext = raster::extent(do.call(c, rev(yx_bbox)))
    r_out = raster::raster(ext=g_ext, resolution=rev(g[['gres']]), crs=g[['crs']])
    if( !is.null(g[['gval']]) )
    {
      r_out = raster::setValues(r_out, matrix(g[['gval']], g[['gdim']]))


    }
    return(r_out)
  }

  # error if we didn't get an expected class
  stop(paste('unrecognized template class', template))
}


#' Snap a set of points to a grid
#'
#' Maps the input points in `from` to the closest grid points in the extension of `g`
#' covering the bounding box of `from` (ie. the lattice of which `g` is a sub-grid).
#' In cases of duplicate mappings, the function returns the first matches only.
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
#' behaves like `terra::rasterize`. It can also be a matrix (supplying dimensions) or a
#' list containing either `gdim` or`gres`, from which an appropriately spaced set of grid
#' lines is derived, centered under the bounding box of the points.
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
#' @param g any grid object accepted or returned by `bk_grid`, the destination grid
#' @param crop_from logical, indicating to omit points not overlying `g`.
#' @param crop_g logical, indicating to trim `g` to the extent of `from`.
#'
#' @return list of form returned by `bk_grid`, defining a grid containing the snapped
#' points. These are assigned the corresponding data value in `from`, or if  `from` has no
#' data, an integer mapping to the points in `from`. Un-mapped grid points are set to NA.
#' @export
#'
#' @examples
#'
#' # functions to scale arbitrary inverval to (1, 2,... 100) and make color palettes
#' num_to_cent = function(x) 1L + floor(99*( x-min(x) ) / diff(range(x)))
#' my_pal = function(x) hcl.colors(x, 'Spectral', rev=T)
#' my_col = function(x) my_pal(1e2)[ num_to_cent(x) ]
#'
#' # create a grid object
#' gdim = c(40, 30)
#' g = bk_grid(list(gdim=gdim, gres=1.1))
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
#' from[['z']] = rnorm(length(from[['y']]))
#' bk_plot(g, reset=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col(from[['z']]))
#' graphics::points(from[c('x', 'y')])
#'
#' # snap only the points overlying the input grid
#' g_snap = bk_snap(from, g, crop_from=TRUE)
#' bk_plot(g_snap, col_grid='black', reset=FALSE, leg=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col(from[['z']]))
#' graphics::points(from[c('x', 'y')])
#'
#' # snap points with  and plot (default settings)
#' g_snap = bk_snap(from, g, crop_from=FALSE, crop_g=FALSE)
#' bk_plot(g_snap, col_grid='black', reset=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col(from[['z']]))
#' graphics::points(from[c('x', 'y')])
#'
#' # find smallest subgrid enclosing all snapped grid points
#' g_snap = bk_snap(from, g, crop_g=TRUE)
#' bk_plot(g_snap, col_grid='black', reset=FALSE)
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col(from[['z']]))
#' graphics::points(from[c('x', 'y')])
#'
#' # create a new grid of different resolution enclosing all input points
#' g_snap = bk_snap(from, g=list(gres=c(0.5, 0.5)))
#' bk_plot(g_snap, reset=FALSE, col_grid='black')
#' graphics::points(from[c('x', 'y')], pch=16, col=my_col(from[['z']]))
#' graphics::points(from[c('x', 'y')])
#'
#' if( requireNamespace('sf') ) {
#'
#' # a different example, snapping mis-aligned subgrid
#' g_pts = bk_grid(list(gdim=c(15, 8), gres=1.7), vals=FALSE)
#' g_pts[['gyx']][['y']] = g_pts[['gyx']][['y']] + 5
#' g_pts[['gyx']][['x']] = g_pts[['gyx']][['x']] + 5
#' n_pts = prod(g_pts$gdim)
#' from = bk_coords(g_pts, out='list')
#'
#' # convert to sf
#' eg_sfc = sf::st_geometry(bk_coords(g_pts, out='sf'))
#' bk_plot(g, reset=FALSE)
#' plot(eg_sfc, add=TRUE)
#'
#' # generate example data and plot
#' eg_sf = sf::st_sf(data.frame(z=rnorm(n_pts)), geometry=eg_sfc)
#' bk_plot(g, reset=FALSE)
#' plot(eg_sf, pch=16, add=TRUE, pal=my_pal)
#' plot(eg_sfc, add=TRUE)
#'
#' # snap points
#' g_snap = bk_snap(from=eg_sf, g)
#' bk_plot(g_snap, reset=FALSE, col_grid='black')
#' plot(eg_sf, pch=16, add=TRUE, pal=my_pal)
#' plot(eg_sfc, add=TRUE)
#'
#' # snapping points without data produces the mapping (non-NA values index "from")
#' g_snap = bk_snap(from=eg_sfc, g)
#' bk_plot(g_snap, ij=TRUE, reset=FALSE, col_grid='black')
#' plot(eg_sfc, add=TRUE)
#'
#' # with crop_g=TRUE)
#' g_snap = bk_snap(from=eg_sfc, g, crop_g=TRUE)
#' bk_plot(g_snap, reset=FALSE, col_grid='black')
#' plot(eg_sfc, add=TRUE)
#'
#' # test with sp class
#' eg_sp = as(eg_sf,'Spatial')
#' g_snap = bk_snap(from=eg_sp, g)
#' bk_plot(g_snap, reset=FALSE, col_grid='black')
#' plot(eg_sf, pch=16, add=TRUE, pal=my_pal)
#' plot(eg_sfc, add=TRUE)
#'
#' }
#'
bk_snap = function(from, g=NULL, crop_from=FALSE, crop_g=FALSE, quiet=FALSE)
{
  # set up more appropriate grid lines later when input g doesn't specify them
  gres_input = NULL
  if( is.matrix(g) ) g = list(gdim=dim(g))
  auto_gyx = (!is.list(g) & is.vector(g)) | ( is.list(g) & !('gyx' %in% names(g)) )

  # set up grid dimensions later when only gres specified (set dummy dimensions now)
  auto_gdim = ifelse(auto_gyx & is.list(g), ('gres' %in% names(g)) & !('gdim' %in% names(g)), FALSE)
  if(auto_gdim & is.list(g) ) g = modifyList(g, list(gdim=c(y=2L, x=2L)))

  # unpack g as blitzKrig grid list but don't copy data values
  g = bk_grid(g, vals=FALSE)
  to_crs = g[['crs']]
  vals = NULL

  # expected point object classes and names
  nm_yx = c('y', 'x')
  is_sf = any( c('sf','sfc', 'sfg') %in% class(from) )
  is_sp = any( c('SpatialPoints', 'SpatialPointsDataFrame') %in% class(from) )
  #to_crs = NULL

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
      g = bk_grid(list(gdim=g[['gdim']], gyx=gyx, gres=g[['gres']]))

    } else {

      # place grid lines to coincide with bounding box edge points, then recompute gres
      gyx = Map(function(yx, n) seq(min(yx), max(yx), length.out=n), yx=from_yx, n=g[['gdim']])
      gres = sapply(gyx, function(yx) diff(yx[1:2]))
      g = bk_grid(list(gdim=g[['gdim']], gyx=gyx, gres=gres))

    }
  }

  # find bounding box of destination grid template and reshape to list of min and max
  g_bbox = lapply(bk_coords(g, out='list', corner=TRUE), range)
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

  # compute new grid line locations and initialize the output grid list object
  to_yx = Map(function(a, b, r) seq(a, b, by=r), a=to_min, b=to_max, r=g[['gres']])
  g_out = bk_grid(list(gyx=to_yx), vals=FALSE)
  if( !is.null(to_crs) ) g_out[['crs']] = to_crs

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

  # reflect indexing of vertical axis for blitzKrig
  ij_min[['i']] = g_out[['gdim']]['y'] + 1L - ij_min[['i']]
  to_idx = bk_mat2vec(ij_min, g_out[['gdim']])

  # handle multiple points mapping to a single grid-point
  is_dupe = duplicated(to_idx)
  if( !quiet & any(is_dupe) ) warning( paste('omitting', sum(is_dupe), 'duplicate mapping(s)') )

  # match magic to get NAs at unassigned grid points
  to_all = seq( prod(g_out[['gdim']]) )
  from_idx = match(to_all, to_idx)
  if( is.null(vals) ) { gval = from_idx } else { gval = vals[from_idx] }
  return( c(list(gval=gval), g_out) )
}

