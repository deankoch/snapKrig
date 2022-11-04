# sk_plot.R
# Dean Koch, 2022
# graphics


#' Plot grid data
#'
#' Plots a matrix or raster as a heatmap with an optional color bar legend. This is a wrapper
#' for `graphics::image` similar to `terra::plot` but with tighter margins to increase the
#' image area on screen; and with a different layout for heat-maps of matrices, which are
#' plotted with i and j axes (row and column numbers) replacing y and x.
#'
#' This method is called automatically when a "sk" object is passed to `plot`.
#'
#' `g` can be a vector, in which case `gdim` supplies the y, x dimensions (ie number
#' of rows, number of columns), in that order. `g` can also can be a matrix, raster or any
#' other object understood by `sk` (in which case `gdim` can be omitted).
#'
#' The data in `g` can be of numeric, integer, logical, or factor class. The numeric class is
#' plotted with a continuous color bar legend, while the others get a categorical legend.
#'
#' Category names (ie tick labels) for the legend can be supplied in argument `breaks`, with
#' one character string for each unique non-NA value in `g`, as integer, in ascending order.
#' If the data are continuous, `breaks` can either be the desired number of bins for coloring,
#' or a vector of break points delineating bins (passed to `graphics::image`). Note that a set
#' of `n` bins has `n+1` break points.
#'
#' `pal` should be one of the palette names returned by `grDevices::hcl.pals`, or else a
#' vector of color names with length one fewer than the number of break points.
#'
#' If the data are all `NA`, the function omits the heatmap and legend, and draws grid lines
#' instead. `col_grid` can be specified to enable/disable grid lines more generally. These
#' lines can sometimes appear misaligned due to anti-aliasing. If this is a problem, try
#' switching to a different graphics device back-end (eg. in Windows Rstudio, try changing from
#' the default to Cairo with `options(RStudioGD.backend = 'cairo')`).
#'
#' The function sets the graphical parameters 'oma' (to `c(0,0,0,0)`) and 'mar' (to values
#' between 1.1 and 5.1, depending on whether titles, axes, and legend are plotted), then calls
#' `graphics::image`, which sets other graphical parameters such as 'usr'. By default all
#' graphical parameters are reset to their original values at the end of the function call.
#' `reset=FALSE` prevents this, so that additional elements can be added to the plot later
#' (such as by calling `sf::st_plot(.., add=TRUE)` or `graphics:lines`).
#'
#' @param g vector, matrix, or any object understood by `sk`
#' @param gdim numeric vector, (optional) grid dimensions of the data when `g` is a vector
#' @param ... plotting parameters (see details)
#'
#' @return The function returns a vector of suggested plot height and width values in units of
#'  inches wich minimize the unused margin space. For example to save a trim version of your
#'  plot as png, call `sk_plot` first to get the suggested height and width, say `y` and `x`,
#'  then pass the result to `png(filename, height=m*y, width=m*x, pointsize=m*12, ...)`,
#'  where `m` is any positive scaling factor.
#'
#' @section Plotting parameters:
#'
#' The following style parameters are optional:
#'
#' \describe{
#'
#' \item{adj, leg_just}{numeric in [0,1]: respectively, the horizontal justification
#' of the title and vertical justification of the color bar legend (default 0.5 for both) }
#'
#' \item{asp}{ numeric or NA: the aspect ratio parameter passed to graphics::image (default 1) }
#'
#' \item{axes, leg}{logical: respectively indicates to draw axes (y and x, or i and j),
#' the color bar legend (default TRUE)}
#'
#' \item{axes_w, leg_w, lab_w, main_w, oma_w}{numeric: respectively, the number of lines of
#' margin to reserve for axes (default 2), legend (default 5), axis labels (default 2), main
#' title (default 2), and outer white-space (default 0.5)}
#'
#' \item{breaks}{numeric (vector) or character vector: the color break points (see details)}
#'
#' \item{cex, cex.main, cex.x, cex.y, cex.z}{numeric: controls the size of text elements in
#' the plot (default 1), those for title, x/y labels and ticks, and legend title and ticks
#' all inherit the value assigned to cex (unless otherwise specified).}
#'
#' \item{col_box, col_grid}{character: respectively, the colors to use for drawing a box
#' around the image border and for drawing grid cell boundaries (NA to omit)}
#'
#' \item{col_invert, col_rev}{logical: respectively, inverts (default FALSE), and reverses
#' the color scale (default TRUE}
#'
#' \item{ij}{logical: enables/disables matrix style plot with j axis annotations on top
#' (default TRUE for vector and matrix input, otherwise FALSE)}
#'
#' \item{layer}{integer: the layer (column) to plot (default 1)}
#'
#' \item{lwd_axis, lwd_ticks, lwd_grid}{numeric: respectively, line widths for the axis
#' lines, ticks, and grid lines (default 1)}
#'
#' \item{main, zlab, ylab, xlab}{character: respectively, a title to put on top in bold,
#' a legend title to put over the color bar, and axis titles for dimensions y and x.
#' Setting to '' omits both the label and its margin space}
#'
#' \item{minimal}{logical: produces a stripped down plot showing only the grid data. This
#' omits all annotation (unless otherwise specified by `axes` and/or `leg`) and removes
#' all margin space (unless otherwise specified by `leg_w` and/or `mar_w`)}
#'
#' \item{pal}{character (vector): one of `graphics::hcl.pals` (default 'Spectral') or a
#'  vector of colors}
#'
#' \item{pxmax}{integer: the maximum number of pixels to draw along each dimension
#' (default 2000). If either x or y dimension exceeds this limit, the grid is up-scaled
#' before plotting}
#'
#' \item{reset}{logical: indicates to restore original graphical parameters after plot is
#'  finished (default TRUE)}
#'
#' \item{zlim}{numeric vector: range in the data to plot (ignored for discrete plots)}
#'
#' \item{x_ontop}{logical: toggles the placement of the horizontal dimension axis on
#' top vs bottom }
#'
#' }
#'
#' @export
#' @seealso sk
#' @family plotting functions
#'
#' @examples
#' # example grid
#' gdim = c(50, 100)
#' n = prod(gdim)
#' g = sk(gdim)
#'
#' # plot the grid layout as raster then as matrix
#' plot(g)
#' plot(g, ij=TRUE)
#'
#' # example data: cosine of squared distance to top left corner
#' z = apply(expand.grid(g$gyx), 1, \(z) cos( 2*sum(z^2) ) )
#' g_example = modifyList(g, list(gval=z))
#' plot(g_example)
#'
#' # plot as matrix (changes default palette)
#' plot(g_example, ij=T)
#'
#' # alignment
#' plot(g_example, ij=T, main='Centered title and legend by default')
#' plot(g_example, ij=T, main='adj: left-right justification of title', adj=0)
#' plot(g_example, ij=T, main='leg_just: top-bottom justification of color bar', leg_just=0)
#'
#' # set the palette - see hcl.pals() for valid names
#' pal = 'Zissou 1'
#' plot(g_example, pal=pal, main=pal)
#' plot(g_example, pal=pal, main=pal, col_invert=TRUE)
#' plot(g_example, pal=pal, main=pal, col_invert=TRUE, col_rev=TRUE)
#'
#' # example data: cosine of distance to top left corner
#' g[] = apply(expand.grid(g$gyx), 1, \(z) cos( sqrt(sum(z^2))/50 ) )
#' plot(g)
#'
#' # specify the layer for multi-layer objects (default is first layer)
#' g_multi = sk(list(gdim=gdim, gval=cbind(z, z^2)))
#' plot(g_multi)
#' plot(g_multi, layer=2)
#'
#' # reduce number of color breaks or specify a factor for discrete value plots
#' plot(g, breaks=50)
#' plot(g, breaks=3)
#' g[] = cut(g[], breaks=3, dig.lab=1)
#' plot(g)
#'
#' # pass color bar labels for discrete plots in breaks (in order low to high)
#' plot(g, breaks=c('a', 'b', 'c'), zlab='group')
#'
#' # select some "observed" points and make a covariance matrix
#' idx_obs = match(seq(n), sort(sample.int(prod(gdim), 1e2)))
#' g[] = idx_obs
#' plot(g)
#' v = sk_var(g)
#'
#' # matrix display mode is automatic when first argument is a matrix or vector
#' sk_plot(v, zlab=expression(V[ij]))
#' sk_plot(c(v), dim(v), zlab=expression(V[ij]))
#'
#' # or pass the matrix to `sk` first to turn it into a sk grid object
#' g_v = sk(v)
#' plot(g_v, zlab=expression(V[ij]))
#'
#' # minimal versions
#' plot(g_v, minimal=T)
#' plot(g_v, minimal=T, leg=T)
#' plot(g_v, minimal=T, col_grid='white', leg=TRUE)
#'
#' # logical matrix plots are gray-scale by default
#' plot(g_v > 1e-2, main='logical matrix')
#'
#' # logical, integer and factor class matrices get a discrete color-bar
#' interval = 1e-2 # try also 1e-3 to see behaviour with large number of bins
#' v_discrete = cut(v, seq(0, ceiling(max(v)), by=interval), dig.lab=2)
#' g_v[] = cut(v, seq(0, ceiling(max(v)), by=interval), dig.lab=2)
#' plot(g_v)
#'
#' # labels are preserved for character matrices
#' z_char = rep(c('foo', 'bar'), n/2)
#' z_char[sample.int(n, n/2)] = NA
#' sk_plot(z_char, gdim)
#'
#' # multi-pane plot
#' g_sim = sk_sim(c(100, 50), n_layer=3)
#' split.screen(c(1,3))
#' screen(1)
#' plot(g_sim, main='layer 1', layer=1, minimal=T, col_box='black')
#' screen(2)
#' plot(g_sim, main='layer 2', layer=2, minimal=T, col_box='black')
#' screen(3)
#' plot(g_sim, main='layer 3', layer=3, minimal=T, col_box='black')
#'
sk_plot = function(g, gdim=NULL, ...)
{
  ###################################################################################
  ## case by case defaults

  # determine if we are in ij and/or minimal mode, and whether legend is requested
  ij = ifelse( is.null( list(...)[['ij']] ), is.matrix(g), list(...)[['ij']])
  minimal = ifelse( is.null( list(...)[['minimal']] ), FALSE, list(...)[['minimal']])
  leg = ifelse( is.null( list(...)[['leg']] ), !minimal, list(...)[['leg']])

  # default plot styling
  ylab = 'y'
  xlab = 'x'
  zlab = 'z(x,y)'
  pal = 'Spectral'
  col_box = 'black'
  axes = TRUE

  # adjusted defaults for matrix plot (ij mode)
  if(ij)
  {
    ylab = 'row i'
    xlab = 'column j'
    zlab = expression(z[ij])
    pal = 'Inferno'
  }

  # adjusted defaults for minimal mode
  if(minimal)
  {
    xlab = ''
    ylab = ''
    col_box = NA
    leg = FALSE
    axes = FALSE
  }

  ###################################################################################
  ## set more defaults and overwrite with any user-defined plot arguments

  # check for title and axis labels
  main = ifelse( is.null( list(...)[['main']] ), '', list(...)[['main']])
  ylab = ifelse( is.null( list(...)[['ylab']] ), ylab, list(...)[['ylab']])
  xlab = ifelse( is.null( list(...)[['xlab']] ), xlab, list(...)[['xlab']])
  has_main = nchar(main) > 0
  has_ylab = nchar(ylab) > 0
  has_xlab = nchar(xlab) > 0

  # leave these NULL for now if not supplied by user
  zlim = list(...)[['zlim']]
  breaks = list(...)[['breaks']]

  # lots more defaults, some of which depend on the above
  asp = ifelse( is.null( list(...)[['asp']] ), 1, list(...)[['asp']])
  axes = ifelse( is.null( list(...)[['axes']] ), axes, list(...)[['axes']])
  axes_w = ifelse( is.null( list(...)[['axes_w']] ), ifelse(axes, 2, 0), list(...)[['axes_w']])
  cex = ifelse( is.null( list(...)[['cex']] ), 1, list(...)[['cex']])
  cex.main = ifelse( is.null( list(...)[['cex.main']] ), cex, list(...)[['cex.main']])
  cex.y = ifelse( is.null( list(...)[['cex.y']] ), cex, list(...)[['cex.y']])
  cex.x = ifelse( is.null( list(...)[['cex.x']] ), cex, list(...)[['cex.x']])
  cex.z = ifelse( is.null( list(...)[['cex.z']] ), cex, list(...)[['cex.z']])
  lwd_axis = ifelse( is.null( list(...)[['lwd_axis']] ), 1, list(...)[['lwd_axis']])
  lwd_ticks = ifelse( is.null( list(...)[['lwd_ticks']] ), 1, list(...)[['lwd_ticks']])
  lwd_grid = ifelse( is.null( list(...)[['lwd_grid']] ), 1, list(...)[['lwd_grid']])
  col_box = ifelse( is.null( list(...)[['col_box']] ), col_box, list(...)[['col_box']])
  col_grid = ifelse( is.null( list(...)[['col_grid']] ), NA, list(...)[['col_grid']])
  col_grid_y = ifelse( is.null( list(...)[['col_grid_y']] ), col_grid, list(...)[['col_grid_y']])
  col_grid_x = ifelse( is.null( list(...)[['col_grid_x']] ), col_grid, list(...)[['col_grid_x']])
  col_invert = ifelse( is.null( list(...)[['col_invert']] ), FALSE, list(...)[['col_invert']])
  col_rev = ifelse( is.null( list(...)[['col_rev']] ), TRUE, list(...)[['col_rev']])
  col_w = ifelse( is.null( list(...)[['col_w']] ), 0.7, list(...)[['col_w']])
  pal = ifelse( is.null( list(...)[['pal']] ), pal, list(...)[['pal']])
  layer = ifelse( is.null( list(...)[['layer']] ), 1L, list(...)[['layer']])
  leg_just = ifelse( is.null( list(...)[['leg_just']] ), 0.5, list(...)[['leg_just']])
  oma_w = ifelse( is.null( list(...)[['oma_w']] ), 0.5, list(...)[['oma_w']])
  lab_w = ifelse( is.null( list(...)[['lab_w']] ), 2, list(...)[['lab_w']])
  leg_w = ifelse( is.null( list(...)[['leg_w']] ), ifelse(leg, 5, 0), list(...)[['leg_w']])
  main_w = ifelse( is.null( list(...)[['main_w']] ), ifelse(has_main, 2, 0), list(...)[['main_w']])
  zlab = ifelse( is.null( list(...)[['zlab']] ), zlab, list(...)[['zlab']])
  adj = ifelse( is.null( list(...)[['adj']] ), 0.5, list(...)[['adj']])
  x_ontop = ifelse( is.null( list(...)[['x_ontop']] ), ij, list(...)[['x_ontop']])
  reset = ifelse( is.null( list(...)[['reset']] ), FALSE, list(...)[['reset']])
  px_max = ifelse( is.null( list(...)[['px_max']] ), 2e3, list(...)[['px_max']])


  # additional settings that depend on the above: where to put the axess
  axis_s = c(y=2L, x=ifelse(x_ontop, 3L, 1L))


  ###################################################################################
  ## process input data

  # convert vector input to matrix
  if( is.vector(g) & !is.list(g) )
  {
    # convert vectors to matrix
    if( is.null(gdim) ) stop('grid dimensions gdim must be supplied when g is a vector')
    if( length(g) != prod(gdim) ) stop('length of grid vector g was not equal to prod(gdim)')
    g = matrix(g, gdim)
  }

  # convert matrix and raster to sk
  g = sk(g)

  # slice off one layer of multi-layer input
  is_multi = is.matrix(g[['gval']])
  if( is_multi )
  {
    # overwrite with the specified layer
    g[['gval']] = as.vector(g[, layer])
    g[['idx_grid']] = NULL
  }

  # upscale as needed then overwrite new grid dimensions
  up_fac = ceiling( dim(g)/px_max )
  if( any( up_fac > 1 ) ) g = sk_rescale(g, up=up_fac)
  gdim = dim(g)

  # compute row-first ordering expected by graphics/grDevices
  idx_image = matrix(seq_along(g), gdim)[gdim['y']:1, ]

  # copy the vectorized data
  z = g[]
  z_is_na = is.na(g)
  na_image = all(z_is_na)

  # handle logical, character, integer types
  if(!na_image)
  {
    # coerce logical to integer
    bool_image = is.logical(z)
    if( bool_image )
    {
      # these annoying edge cases require different breaks
      all_true = all(z[!z_is_na])
      all_false = !any(z[!z_is_na])
      z = as.integer(z)
      if( is.null(breaks) )
      {
        breaks = c('false', 'true')
        if(all_true) breaks = 'true'
        if(all_false) breaks = 'false'
      }
    }

    # coerce factor and character to integer, preserving names in `breaks`
    if( is.character(z) ) z = as.factor(z)
    if( is.factor(z) )
    {
      # this can be slow for large vectors
      bins = sort(unique(z))
      z = match(z, bins)
      if( is.null(breaks) ) breaks = as.character(bins)
    }
  }

  # construct matrix representation expected by graphics/grDevices
  z_image = matrix(z[idx_image], rev(gdim), byrow=TRUE)


  ###################################################################################
  ## call graphics::image

  # find limits from data if they are not all NA
  if( !na_image ) { if( is.null(zlim) ) zlim = range(z_image, na.rm=TRUE) } else {

    # turn off legend and show grid lines unless user specifies otherwise
    zlim = c(0,1)
    if( is.null( list(...)[['leg']] ) ) leg = FALSE
    if( is.null( list(...)[['col_grid']] ) ) col_grid = 'black'
  }

  # handle plots with a single unique non-NA value
  if( diff(zlim) == 0 )
  {
    # assign the (single) break point label to the non-NA value
    if( is.null(breaks) ) breaks = as.character(zlim[1])

    # convert matrix to integer (1's indicating the non-NA value)
    z_image = matrix(as.integer(z_image == zlim[1]), rev(gdim))
    zlim = c(1, 1)
  }

  # integer data case (discrete)
  discrete_image = !na_image & is.integer(z_image)
  if( discrete_image )
  {
    # use a better default palette for small numbers of categories
    if( is.null( list(...)[['pal']] ) ) pal = 'Zissou'

    # discretize z limits and find unique levels for color scale
    zlim = c(ceiling(zlim[1]), floor(zlim[2]))
    z_lvl = as.numeric(zlim[1]:zlim[2])

    # add padding so extremes are included
    zlim = zlim + c(-0.5, 0.5)

    # character input to `breaks` interpreted as labels (else use integers, and `breaks` ignored)
    z_label = as.character(z_lvl)
    if( is.character(breaks) ) z_label = breaks
    breaks = length(z_lvl)

    # set gray-scale for unary/binary cases
    if( breaks == 1 ) pal = '#1A1A1A'
    if( breaks == 2 & !any(z_is_na) ) pal = c('#E5E5E5', '#1A1A1A')

  } else {

    # set default set number of breaks high enough to look continuous
    if( is.null(breaks) ) breaks = 1e3
    z_lvl = pretty(zlim)
    z_label = as.character(z_lvl)

  }

  # set up a color palette
  pal_known = any(startsWith(grDevices::hcl.pals(), pal[1]))
  if( length(breaks) == 1 ) breaks = seq(zlim[1], zlim[2], length.out=1L+breaks)
  if( pal_known ) pal = grDevices::hcl.colors(length(breaks)-1L, pal[1], rev=col_rev)
  if( col_invert ) pal = rgb( ( 255 - t(col2rgb(pal)) ) / 255 )

  # define margins with space for color bar, axes, titles
  mar_new = c(bottom = ( axes_w * (!x_ontop) ) + ( lab_w * has_xlab ),
              left = axes_w + ( lab_w * has_ylab ),
              top = ( ( axes_w + ( lab_w * has_xlab) ) * x_ontop ) + ( main_w * has_main ),
              right = leg_w)

  # set new graphical parameters but keep a backup
  if(reset) par_existing = par(no.readonly=TRUE)
  par(mar=mar_new, oma=rep(oma_w, 4L))

  # draw the raster plot without axes or annotations
  graphics::image(g[['gyx']], z=z_image, axes=FALSE, ann=FALSE, xpd=NA,
                   asp=asp, useRaster=TRUE, zlim=zlim, col=pal, breaks=breaks)


  ###################################################################################
  ## draw annotations

  # compute bounding box of image
  plot_bbox = list(min=sapply(g[['gyx']], min), max=sapply(g[['gyx']], max))
  half_pixel = sapply(g[['gyx']], function(yx) diff(yx[1:2])) / 2

  # position and draw axes and tick marks
  if(axes)
  {
    # axes are snapped to raster image edge
    axis_xy = rev(plot_bbox[['min']] - half_pixel)

    # find row and column indices for tick mark locations, then y/x coordinates
    tick_pos = tick_lab = lapply(g[['gyx']], function(yx) pretty(range(yx)))
    if(ij)
    {
      # force labeling of first and last row/column in matrices and set new default labels
      tick_ij = lapply(gdim, function(d) unique(round(c(1, seq(1,d, length.out=5), d))))
      tick_lab = modifyList(tick_ij, list( y = gdim['y'] + 1 - tick_ij[['y']] ))
      tick_pos = Map(\(yx, idx) yx[idx], yx=g[['gyx']], idx=tick_ij)

      # reverse ordering of y axis for i labels, j labels go on top
      axis_xy['y'] = plot_bbox[['max']]['y'] + half_pixel['y']

      # different default aesthetics for matrix plot
      draw_box = TRUE
      lwd_axis = 0
    }

    # draw axes and tick marks (loop is over y, x)
    Map(\(s, a, lab, yx, cex) graphics::axis(s, at=a, labels=lab, pos=yx,
                                        lwd=lwd_axis, lwd.ticks=1, cex.axis=cex),
        s = axis_s,
        a = stats::setNames(tick_pos, c('y', 'x')),
        lab = tick_lab,
        yx = axis_xy,
        cex = c(cex.y, cex.x))
  }

  # title and axis label line positions
  main_nline = par('mar')[3] - 0.7 * main_w
  y_nline = par('mar')[2L] - 0.7 * lab_w
  x_nline = par('mar')[ axis_s['x'] ] - ( main_w * (ij & x_ontop) ) - 0.7 * lab_w

  # print any titles
  graphics::mtext(main, side=3, line=main_nline, font=2, xpd=NA, adj=adj, cex=cex.main)
  graphics::mtext(ylab, side=axis_s['y'], line=y_nline, xpd=NA, cex=cex.y)
  graphics::mtext(xlab, side=axis_s['x'], line=x_nline, xpd=NA, cex=cex.x)

  # finding bounding box for the image
  bmin = plot_bbox[['min']] - half_pixel
  bmax = plot_bbox[['max']] + half_pixel

  # find line height for horizontal text, in same units as user coordinates
  y_inch = diff(graphics::grconvertY(0:1, 'inches', 'user'))
  y_line = y_inch * par('cin')[2] * par('cex') * par('lheight')

  # same for vertical text
  x_inch = diff(graphics::grconvertX(0:1, 'inches', 'user'))
  x_line = x_inch * par('cin')[2] * par('cex') * par('lheight')

  # draw x grid lines
  if( !is.na(col_grid_x) )
  {
    gl_x = g[['gyx']][['x']] - half_pixel['x']
    graphics::segments(x0=gl_x, y0=bmin['y'], y1=bmax['y'], col=col_grid_x, lwd=lwd_grid)
  }

  # draw y grid lines
  if( !is.na(col_grid_y) )
  {
    gl_y = g[['gyx']][['y']] - half_pixel['y']
    graphics::segments(y0=gl_y, x0=bmin['x'], x1=bmax['x'], col=col_grid_y, lwd=lwd_grid)
  }

  # draw a frame around the image
  if( !is.na(col_box) ) graphics::rect(bmin['x'], bmin['y'], bmax['x'], bmax['y'], border=col_box)

  # add a color bar legend
  # TODO: add option to draw this invisibly
  if(leg)
  {
    # set x coordinates with width equal to height of one line of text
    bar_x = bmax['x'] + ( col_w * x_line * c(1,2) )

    # set y coordinates, attempting to match breaks to line height
    n_bin = length(breaks) - 1L
    n_pad = 2L
    n_y_line_max = floor( ( bmax['y'] - bmin['y'] ) / y_line ) - 2L * n_pad
    n_y_line = pmax(pmin(n_bin, n_y_line_max), 2)
    y_pad = ( bmax['y'] - bmin['y'] - y_line * n_y_line )
    bar_y = c(bmin['y'], bmax['y']) + y_pad * c(leg_just, -(1-leg_just))

    # map scale colors to y position
    bar_z = stats::approxfun(range(breaks), bar_y)(breaks)

    # draw the color bar strips then a black border around the whole thing
    graphics::rect(bar_x[1], head(bar_z, -1), bar_x[2], bar_z[-1], col=pal, border=NA, xpd=NA)
    graphics::rect(bar_x[1], bar_y[1], bar_x[2], bar_y[2], xpd=NA)

    # draw annotations as another axis
    bar_ticks_y = stats::approxfun(zlim, bar_y)(z_lvl)
    graphics::axis(4L,
                   lwd = 0,
                   at = bar_ticks_y,
                   labels = z_label,
                   pos = bar_x[2],
                   lwd.ticks = lwd_ticks,
                   las = 1L,
                   cex.axis = cex.z)

    # draw the legend title one half-line above the color bar
    graphics::text(bar_x[1], bar_y[2] + y_line / 2, zlab, xpd=NA, adj=c(0,0), cex=cex.z)
  }

  # find minimal y, x lengths required for image plot plus margins in inches
  y_graphics = ( ( bmax['y'] - bmin['y'] ) + sum( mar_new[c(1,3)] ) * y_line ) / y_inch
  x_graphics = ( ( bmax['x'] - bmin['x'] ) + sum( mar_new[1+c(1,3)] ) * x_line ) / x_inch

  ###################################################################################
  ## tidy up

  # restore old margin settings before the user's next plot call
  if(reset) par(par_existing)
  return(invisible(c(y=y_graphics, x=x_graphics)))
}


#' Plot the covariance structure of a snapKrig model
#'
#' Visualization of the footprint of a covariance kernel as a heatmap of approximate
#' size `dim(g)`, where each grid cell is colored according to its covariance with
#' the central grid point.
#'
#' If `g` is not supplied, the function sets a default with dimensions 100 x 100.
#' A default resolution is computed such that the maximum nugget-free covariance along
#' the outer edge of the plot is 5% of `pars$psill`.
#'
#' When `simple=FALSE` (the default), covariance parameters are printed in the
#' title and axis labels with values rounded to 3 decimal places. This can be
#' customized by passing arguments 'main', 'ylab', 'xlab' (and any others
#' accepted `sk_plot` apart from `gdim`).
#'
#' @param pars list of the form returned by `sk_pars` with entries 'y', 'x', ('eps', 'psill')
#' @param g any object understood by `sk`
#' @param ... additional arguments passed to `sk_plot`
#'
#' @return the same as `sk_plot`
#'
#' @export
#' @seealso sk sk_pars
#' @family plotting functions
#'
#' @examples
#' gdim = c(100, 100)
#' g = sk(gdim)
#' pars = sk_pars(g, 'mat')
#'
#' # plot with default grid
#' sk_plot_pars(pars)
#'
#' # plot with a predefined grid
#' sk_plot_pars(pars, g)
#'
#' # zoom in/out by passing a grid object with suitably modified resolution
#' gres = g[['gres']]
#' sk_plot_pars(pars, sk(gdim=gdim, gres=0.5*gres))
#' sk_plot_pars(pars, sk(gdim=gdim, gres=2*gres))
#'
#' # change plot style settings (all parameters of sk_plot accepted)
#' sk_plot_pars(pars, simple=TRUE)
#' sk_plot_pars(pars, minimal=TRUE)
#'
sk_plot_pars = function(pars, g=NULL, simple=FALSE, ...)
{
  # if a grid is not supplied, the function sets a default based on range parameters
  if( !is.null(g) ) { g = sk(g) } else {

    # approximate the range at which correlation is around 0.05
    c_target = sqrt(0.05) * pars[['psill']]
    d_max = 100 * max(sapply(pars[c('y', 'x')], function(p) p[['kp']]['rho']))
    d_test = seq(0, d_max, length.out=1e3)
    dist_score = sapply(pars[c('y', 'x')], function(p) sk_corr(p, d_test) - c_target)
    d_target = d_test[which.min(rowSums(dist_score))]

    # by default plot 5 times the largest range
    if( !(d_target > 0) ) stop('invalid range parameter(s)')

    # set default plot dimensions and compute resolution
    gdim = c(100, 100)
    g = sk(gdim=gdim, gres=d_target/gdim)
  }

  # copy important parameters
  gdim = dim(g)
  gres = g[['gres']]
  psill = pars[['psill']]
  eps = pars[['eps']]

  # increment dimensions to get odd number (having a central row and column)
  gdim = gdim + 1 - (g[['gdim']] %% 2)
  ij_mid = 1L + floor(gdim/2)

  # construct the required rows of the correlation matrices
  cy = sk_corr_mat(pars[['y']], gdim[['y']], gres[['y']], i=ij_mid['y'])
  cx = sk_corr_mat(pars[['x']], gdim[['x']], gres[['x']], i=ij_mid['x'])

  # find the covariance values and create a grid list object for them
  cov_values = eps + psill * as.vector(kronecker(cx, cy))
  g_plot = sk(gdim=gdim, gres=gres, gval=cov_values)

  # replace grid line positions with displacement variable
  g_plot[['gyx']] = Map(function(gyx, ij) gyx - ij + 1L, gyx=g_plot[['gyx']], ij=ij_mid)

  # make a plot title
  titles = sk_to_string(pars)
  ylab = paste('y displacement', ifelse(simple, '', titles[['kp']][['y']]))
  xlab = paste('x displacement', ifelse(simple, '', titles[['kp']][['x']]))
  main = ifelse(simple, paste(titles[['k']]), paste(titles[['main']]))
  zlab = expression(V[ij])

  # user arguments override titles set above in call to sk_plot
  if( !is.null( list(...)[['main']] ) ) main = list(...)[['main']]
  if( !is.null( list(...)[['zlab']] ) ) ylab = list(...)[['zlab']]
  if( !is.null( list(...)[['ylab']] ) ) ylab = list(...)[['ylab']]
  if( !is.null( list(...)[['xlab']] ) ) xlab = list(...)[['xlab']]
  return( sk_plot(g_plot, main=main, zlab=zlab, ylab=ylab, xlab=xlab, ...) )
}


#' Plot a semi-variogram
#'
#' Plots a sample semi-variogram using the point pair difference data in `vg`.
#' Binned summary statistics are drawn as circles with size scaled to the sample sizes.
#' A covariance model (`pars`) is optionally drawn over the sample data as a ribbon plot.
#'
#' If `vg` is a data frame, it should contain absolute differences (numeric 'dabs'),
#' inter-point distances (numeric 'd'), and, optionally, an assignment into distance bins
#' (integer 'bin') for a sample of point pairs. If 'bin' is missing, the function calls
#' `sk_add_bins` to assign them automatically.
#'
#' Function `fun` is the statistic to use for estimating the variogram (ie twice the
#' semi-variogram) from the distance-binned absolute differences in `vg`. If `fun` is a
#' function, it must accept sub-vectors of the numeric `vg$dabs` as its only argument,
#' returning a non-negative numeric scalar. `fun` can also be set to one of the names
#' 'root_median', 'root_mean' (the default), or 'classical', as shorthand for the robust
#' fourth-root-based methods in section 2.4 of Cressie (1993), or the classical mean of
#' squares method of Matheron.
#'
#' Optional list `pars` defines a theoretical semi-variogram to draw over the sample data.
#' When `add=TRUE`, the function overlays it on an existing plot (without changing the
#' legend, plot limits etc). Anisotropic models, which may assume a range of semi-variances
#' for any given distance, are drawn as a ribbon plot.
#'
#' `add=TRUE` can only be used in combination with an earlier call to `sk_plot_semi`
#' where `reset=FALSE` (which allows the function to change R's graphical parameters)
#'
#' `vg` can be a grid object (anything understood by `sk`) rather than a
#' variogram data frame. When `add=FALSE`, the function uses it to set the distance limits
#' for an initial empty plot (the model semi-variance is then drawn if `pars` is supplied).
#'
#' @param vg data frame of sample absolute differences, with columns 'dabs', 'd' (and 'bin')
#' @param pars list of the form returned by `sk_pars` with entries 'y', 'x', 'eps', 'psill'
#' @param add logical, indicates to draw on an existing plot rather than create a new one
#' @param fun character or function, the aggregation function (see details)
#' @param ... further plotting parameters (see below)
#'
#' @section Plotting parameters:
#'
#' The following style parameters are optional:
#'
#' \describe{
#'
#' \item{alpha_bin, alpha_model, alpha_bin_b, alpha_model_b}{numeric in [0,1]: respectively,
#' the transparency of the fill color for circles and ribbons (default 0.3), and their borders
#' (default 0.5 )}
#'
#' \item{bty}{character: plot frame border type, passed to base::plot (default 'n' for no border)}
#'
#' \item{col_bin, col_model}{character: respectively, the color to use for circles (default
#' 'black') and ribbons (default 'blue')}
#'
#' \item{cex_bin}{numeric > 0: scaling factor for circles (default 1.5)}
#'
#' \item{d_max}{numeric > 0: x (distance) limit for plotting in input units}
#'
#' \item{leg_main}{character: title for the sample bin legend (default 'model')}
#'
#' \item{lwd}{numeric: line width for the model semi-variance}
#'
#' \item{main}{character: a title}
#'
#' \item{n_bin, n_test}{integer: respectively, the number of distance bins for the sample
#' (optional if `vg` has a 'bin' column, and ignored if `vg` is a grid object), and the
#' number of distances at which to evaluate the semi-variogram for model `pars` (default 5e3,
#' ignored if `pars` not supplied)}
#'
#' \item{reset}{logical: indicates to reset graphical parameters to their original values
#' when finished (default TRUE)}
#'
#' \item{unit_in, unit_out}{character: udunits2 strings specifying input distance units and
#' the desired units for distances in the plot (default for both is meters 'm')}
#'
#' \item{xlab, ylab}{character: titles for the y and x axes. The default for y is
#' 'semi-variance (gamma)', and for x 'distance (<unit_out>)'}
#'
#' }
#'
#' @return nothing
#' @export
#'
#' @examples
#'
#' # make example grid and reference covariance model
#' gdim = c(10, 15)
#' g_empty = sk(gdim)
#' pars = sk_pars(g_empty, 'mat')
#'
#' # plot a semivariance model
#' sk_plot_semi(g_empty)
#' sk_plot_semi(g_empty, pars)
#'
#' # change annotations, sharpen ribbon border
#' sk_plot_semi(g_empty, pars, main='title', xlab='x', ylab='y')
#' sk_plot_semi(g_empty, pars, alpha_model_b=1, main='example title', xlab='x', ylab='y')
#'
#' # input and output units are 'm' by default
#' sk_plot_semi(g_empty, pars, unit_out='km')
#' sk_plot_semi(g_empty, pars, unit_in='km', unit_out='km')
#'
#' # generate sample data and sample semivariogram
#' g_obs = sk_sim(g_empty, pars)
#' vg = sk_sample_vg(g_obs)
#' sk_plot_semi(vg)
#'
#' # different aggregation methods produce variety of results
#' sk_plot_semi(vg, fun='root_median')
#' sk_plot_semi(vg, fun='root_mean')
#' sk_plot_semi(vg, fun='classical') # default
#' sk_plot_semi(vg, fun=function(x) mean(x^2)) # same as classical
#'
#' # plot again with reference model and adjust distance limits, number of bins
#' sk_plot_semi(vg, pars)
#' sk_plot_semi(vg, pars, d_max=10)
#' sk_plot_semi(vg, pars, d_max=10, n_bin=1e2)
#'
#' # add dashed line for half sample variance (this tends to underestimate the sill)
#' sk_plot_semi(vg, pars)
#' sample_var = var(g_obs[['gval']], na.rm=TRUE)
#' abline(h=sample_var/2, lty='dashed')
#'
#' # initial call with reset=FALSE, then use add=TRUE to overlay the same model with a green fill
#' sk_plot_semi(vg, pars, lwd=2, reset=FALSE)
#' sk_plot_semi(vg, pars, add=TRUE, col_model='green', alpha_model_b=0)
#'
#' # overlay several models with varying nugget effect
#' pars_vary = pars
#' for(i in seq(3))
#' {
#'   pars_vary$eps = 0.8 * pars_vary$eps
#'   sk_plot_semi(vg, pars_vary, add=TRUE, alpha_model_b=0)
#' }
#' dev.off()
#'
sk_plot_semi = function(vg, pars=NULL, add=FALSE, fun='classical', ...)
{
  # flags if plotting a sample semi-variogram, adding a model function
  input_vg = is.data.frame(vg)
  draw_model = !is.null(pars)

  # assign defaults
  n_test = ifelse( is.null( list(...)[['n_test']] ), 5e3, list(...)[['n_test']])
  cex_main = ifelse( is.null( list(...)[['cex_main']] ), 1, list(...)[['cex_main']])
  cex_bin = ifelse( is.null( list(...)[['cex_bin']] ), 1.5, list(...)[['cex_bin']])
  col_bin = ifelse( is.null( list(...)[['col_bin']] ), 'black', list(...)[['col_bin']])
  col_model = ifelse( is.null( list(...)[['col_model']] ), 'blue', list(...)[['col_model']])
  alpha_bin = ifelse( is.null( list(...)[['alpha_bin']] ), 0.3, list(...)[['alpha_bin']])
  alpha_model = ifelse( is.null( list(...)[['alpha_model']] ), 0.3, list(...)[['alpha_model']])
  alpha_bin_b = ifelse( is.null( list(...)[['alpha_bin_b']] ), 0.5, list(...)[['alpha_bin_b']])
  alpha_model_b = ifelse( is.null( list(...)[['alpha_model_b']] ), 0.5, list(...)[['alpha_model_b']])
  unit_in = ifelse( is.null( list(...)[['unit_in']] ), 'm', list(...)[['unit_in']])
  unit_out = ifelse( is.null( list(...)[['unit_out']] ), 'm', list(...)[['unit_out']])
  ylab = ifelse( is.null( list(...)[['ylab']] ), 'semi-variance (gamma)', list(...)[['ylab']])
  leg = ifelse( is.null( list(...)[['leg']] ), TRUE, list(...)[['leg']])
  leg_main = ifelse( is.null( list(...)[['leg_main']] ), 'model', list(...)[['leg_main']])
  n_bin = ifelse( is.null( list(...)[['n_bin']] ), NA, list(...)[['n_bin']])
  reset = ifelse( is.null( list(...)[['reset']] ), TRUE, list(...)[['reset']])
  bty = ifelse( is.null( list(...)[['bty']] ), 'n', list(...)[['bty']])
  lwd = ifelse( is.null( list(...)[['lwd']] ), 1, list(...)[['lwd']])

  # set up the more complicated defaults
  d_max = list(...)[['d_max']]
  main_def = paste(ifelse(input_vg, 'sample', 'model'), 'semi-variogram')
  xlab_def = paste0('distance (', unit_out, ')')
  xlab = ifelse( is.null( list(...)[['xlab']] ), xlab_def, list(...)[['xlab']])
  plot_max = 0

  # handle named variogram aggregation functions
  if(is.character(fun))
  {
    # list of named functions (see section 2.4 on variogram estimation in Cressie, 1993)
    fun_lookup = list(root_mean = function(x) mean( sqrt(x) )^4 / ( 0.457 + 0.494/length(x) ),
                      root_median = function(x) stats::median( sqrt(x) )^4 / 0.457,
                      classical = function(x) mean(x^2) )

    # catch unknown names
    info_agg = paste(names(fun_lookup), collapse=', ')
    msg_agg = 'unrecognized fun name. Try one of:' |> paste(info_agg)
    if( is.null(fun) ) stop(msg_agg)
    fun = fun_lookup[[fun]]
  }

  # set up new margins but keep a backup of existing graphical parameters
  mar_new = mar_existing = par('mar')
  mar_new[4] = ifelse(!add & input_vg, 6.1, mar_existing)
  par(mar=mar_new, oma=c(0,0,0,0))

  # input is a variogram data frame
  if(input_vg)
  {
    # find distance range and take corresponding subset of data frame as needed
    d_max = ifelse(is.null(d_max), max(vg[['d']]), d_max)
    vg = vg[vg[['d']] < d_max, ]

    # calculate new distance bins only if we need them
    n_bin_new = ifelse(is.na(n_bin), 25L, n_bin)
    if( !is.na(n_bin) | is.null(vg[['bin']]) ) vg = sk_add_bins(vg, n_bin_new)

    # compute summary stats by distance bin: sample size, distance mean, semi-var stat
    s_bin = sapply(split(rep(1, nrow(vg)), vg[['bin']]), sum)
    d_bin = sapply(split(vg[['d']], vg[['bin']]), mean)
    v_bin = sapply(split(vg[['dabs']], vg[['bin']]), fun)
    plot_max = 1.1 * ( max(v_bin)/2 )

    # convert size of circles to plot and distances in output units
    size_plot = 1 + cex_bin * ( s_bin / max(s_bin) )
    d_bin_src = units::set_units(d_bin, unit_in, mode='s')
    d_bin_out = units::drop_units( units::set_units(d_bin_src, unit_out, mode='s') )

  } else {

    # get distance range from grid layout and initialize semi-variance range
    g = sk(vg)
    d_max = ifelse(is.null(d_max), sqrt( sum( ( g[['gdim']] * g[['gres']] )^2 ) ), d_max)

    # do nothing if no model was supplied in add mode
    if( add & !draw_model ) return( invisible() )
  }

  # compute model semi-variance values
  if(draw_model)
  {
    # component distances to test along each axis
    d_test = seq(0, d_max, length.out=n_test)

    # compute model semi-variance range
    vario_result = sk_vario_fun(pars, d_test)
    sv_min = vario_result[['min']]/2
    sv_max = vario_result[['max']]/2
    plot_max = max(c(plot_max, 1.1 * sv_max))
  }

  # initialize the plot
  if(!add)
  {
    # set up axis limits, leaving space for bin size legend
    ylim_out = c(0, plot_max)
    xlim_src = units::set_units(c(0, d_max) * c(1, 1 + ifelse(input_vg, 1e-1, 0)), unit_in, mode='s')
    xlim_src = units::set_units(c(0, d_max), unit_in, mode='s')
    xlim_out = units::drop_units( units::set_units(xlim_src, unit_out, mode='s') )

    # set up title
    pars_out = ifelse(draw_model, paste(sk_to_string(pars)[['k']], 'kernel'), '')
    main_def = ifelse(input_vg, paste(main_def, paste0('(n = ', nrow(vg), ')')), main_def)
    if(draw_model) main_def = paste(c(main_def, pars_out), collapse=', ')

    # initialize the plot area without points
    main_out = ifelse( is.null( list(...)[['main']] ), main_def, list(...)[['main']])
    plot(xlim_out, ylim_out, main=main_out, xlab=xlab, ylab=ylab, pch=NA, bty=bty, cex.main=cex_main)

    # find line height in same units as user coordinates
    y_inch = diff(graphics::grconvertY(0:1, 'lines', 'user'))
    x_inch = diff(graphics::grconvertX(0:1, 'lines', 'user'))
    y_line = y_inch * par('cin')[2] * par('cex') * par('lheight')
    x_line = x_inch * par('cin')[2] * par('cex') * par('lheight')

    # add points for sample semi-variance bins
    if(input_vg)
    {
      # draw and color circles for sample semi-variance aggregated by distance
      col_bin_fill = grDevices::adjustcolor(col_bin, alpha.f=alpha_bin)
      col_bin_border = grDevices::adjustcolor(col_bin, alpha.f=alpha_bin_b)
      points(d_bin_out, v_bin/2, pch=16, col=col_bin_fill, cex=size_plot)
      points(d_bin_out, v_bin/2, col=col_bin_border, cex=size_plot)

      # add legend for the circles
      if(input_vg & leg)
      {
        # compute legend ticks (bin sizes to show)
        s_bin_plot = pretty(s_bin, n=3)
        if( any(s_bin_plot == 0) ) s_bin_plot = s_bin_plot[s_bin_plot != 0]
        s_bin_ticks = length(s_bin_plot)
        if(s_bin_ticks > 3) s_bin_plot = s_bin_plot[c(1, s_bin_ticks)]
        size_leg_plot = 1 + cex_bin * ( s_bin_plot / max(s_bin) )

        # draw the legend
        leg_plot = paste(' ', s_bin_plot)
        legend(x=2*x_line + par('usr')[2], y=plot_max/2 + y_line, yjust=0, y.intersp=1,
               title='sample size', legend=leg_plot, pch=1, pt.cex=size_leg_plot,
               border=col_bin, bty='n', xpd=NA)
      }
    }
  }

  # add the theoretical variogram lines
  if(draw_model)
  {
    # convert distances to output units before drawing lines
    d_test_src = units::set_units(d_test, unit_in, mode='s')
    d_test_out = units::drop_units( units::set_units(d_test_src, unit_out, mode='s') )

    # draw upper and lower as polygon
    x_out = c(d_test_out, rev(d_test_out))
    y_out = c(sv_max, rev(sv_min))
    col_model_fill = grDevices::adjustcolor(col_model, alpha.f=alpha_model)
    col_model_border = grDevices::adjustcolor(col_model, alpha.f=alpha_model_b)
    graphics::polygon(x_out, y_out, col=col_model_fill, border=col_model_border, lwd=lwd)

    # add a legend
    if(!add & input_vg & leg) legend(x=2*x_line + par('usr')[2], y=plot_max/2 - y_line,
                               yjust=1,
                               title=leg_main,
                               fill=col_model_fill,
                               border=col_model,
                               xpd=NA, legend='', bty='n')
  }

  # restore old margin settings before the user's next plot call
  if(reset) par(mar=mar_existing)
}
