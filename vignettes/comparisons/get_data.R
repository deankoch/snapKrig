#'
#' ## Example datasets for snapKrig
#'
#' These are three earth science data examples to use in demonstration code for the `snapKrig` package.
#' They are drawn from open and easily accessible sources: two are included in the `fields` and `sp`
#' packages, and the other can be downloaded from a permanent archive at FRDR (LINK HERE) using
#' `rasterbc`.
#'
#' We have variables on the atmosphere (ozone), the biosphere (forest health), and the lithosphere
#' (soil contamination), so maybe I should add a hydrosphere example at some point.
#'
#' ## Ozone data in USA (requires `fields`)
#'
#' These point datasets are included with the `fields` package, and appear frequently in its example
#' code (see for example `?fields::SpatialProcess`). They shows maximum recorded ozone concentrations
#' in parts per billion, over 8-hour periods, on several dates in the summer of 1987, at 147 stations
#' in the US Midwest
#'
#' The function by default returns this full station dataset for a particular date (June 19th).
#' Optionally a smaller subset of 20 station summer averages centered around Chicago can be returned
#' by setting `chi=TRUE `
#'
#' The data are in geographical coordinates (latitude, longitude), so we transform to a projected UTM
#' coordinate system so that Euclidean distances can be more easily computed. We also remove NA
#' points.
get_ozone = function(chi=FALSE)
{
  # define output projection
  epsg_UTM15 = 26715
  epsg_UTM15_string = paste0('EPSG:', epsg_UTM15)

  # loads smaller Chicago example
  if(chi)
  {
    # ChicagoO3 dataset included with fields package
    if( !requireNamespace('fields') ) { stop('fields package not found') }
    ozone_df = data.frame(ozone = fields::ChicagoO3[['y']], fields::ChicagoO3[['lon.lat']])
    row.names(ozone_df) = NULL

  } else {

    # bigger ozone2 dataset included with fields package
    data(ozone2, package='fields', envir = environment())

    # pull a particular date from ozone table and form into a points geometry (omit NAs)
    ozone_date_idx = 16
    ozone_df = data.frame(ozone = ozone2[['y']][ozone_date_idx,], ozone2[['lon.lat']])
    names(ozone_df) = c('ozone', 'lon', 'lat')
  }

  # convert to sf object and transform to a projected coordinate system
  ozone_geo_sf = st_as_sf(ozone_df[!is.na(ozone_df[['ozone']]),], coords=c('lon', 'lat'), crs=4326)
  ozone_sf = st_transform(ozone_geo_sf, crs=epsg_UTM15)
  return(ozone_sf)
}

#'
#' ## Pine beetle data in Canada (requires `rasterbc` and `snapKrig`)
#'
#' This is a raster dataset on insect outbreak damage in pine forests of British Columbia (BC), Canada,
#' appearing in the rasterbc vignette for `snapKrig`. We use the `rasterbc` package to download and save
#' a local copy for the year 2006. The data are model output at 100m resolution, covering the area
#' around the town of Merritt, in central BC, showing a damage severity index, proportional to the
#' number of trees killed.
#'
#' Alternatively, set out='treed' to return the an estimate of pine tree density (%) or out='dem'
#' to the elevation (DEM), covering the same area.
#'
#' The function optionally returns an up-scaled version of the data, in which every 20th pixel is
#' included along each axis and the others dropped. Set `up` to any positive integer to change this
#' sampling interval, or set `up=NA` to get the original data at full resolution. Users may also change
#' the `SNRC` code (see `?rasterbc::ntspoly_bc`) to select a different area of BC.
#'
#' `rasterbc` needs to know where to look for files (and download missing ones) - users should replace
#' the default argument `rasterbc_storage_path` with an appropriate path on their local machine. If
#' you have not run the script before, it will attempt to download about 20MB of GeoTIFF files to that
#' location.
get_mpb = function(up=NA, SNRC='092I', rasterbc_storage_path=NA, out='mpb')
{
  # check for the packages
  if( !requireNamespace('snapKrig') ) { stop('snapKrig package not found') }
  if( !requireNamespace('rasterbc') ) { stop('rasterbc package not found') }

  # set storage directory for TIFF files (default uses same directory as the other vignette)
  rasterbc::datadir_bc(rasterbc_storage_path, quiet=TRUE)

  # DEM requests
  if(out=='dem') output_rast = rasterbc::opendata_bc(SNRC, 'dem', 'dem')

  # treed requests (nearest available year)
  if(out %in% c('treed', 'mpb')) output_rast = rasterbc::opendata_bc(SNRC, 'pine', 'vegTreed', 2001)

  # mpb requests
  if(out=='mpb')
  {
    # download other forest layers to get pine estimate
    needle_full_rast = rasterbc::opendata_bc(SNRC, 'pine', 'needle', 2001)
    pinus_full_rast = rasterbc::opendata_bc(SNRC, 'pine', 'pinusTotal', 2001)

    # estimate Pinus density as fraction of pixel area (output_rast should be treed %)
    pine_rast = (output_rast/100) * (needle_full_rast/100) * (pinus_full_rast/100)

    # download midpoint of MPB damage (fraction) estimate in 2006
    IBM_full_rast = rasterbc::opendata_bc(SNRC, 'fids', 'IBM_mid', 2006)

    # estimate pine beetle damage (fraction of area affected)
    output_rast = pine_rast * IBM_full_rast
  }

  # convert to snapKrig grid object and trim outer NA rows and columns
  output_g = sk_sub(sk(output_rast))

  # upscale as needed before returning
  if( is.na(up) ) {  return( sk_export(output_g)) } else {

    # export converts back to SpatRaster
    return( sk_export(sk_rescale(output_g, up=up)) )
  }
}

#'
#' ## Soil quality in Europe (requires `sp`)
#'
#' This is point data on heavy metal concentrations on a river floodplain in Northern Europe.
#' It is included with the `sp` and `gstat` packages and appears in their vignettes and docs.
#'
#' The helper function by default simply opens the data, assigns the CRS code, and coerces to
#' an sf POINT collection containing only the logarithm of the zinc concentration.
#'
#' We also have the option to return some additional data layers related to the river. These
#' are used in the Meuse vignette but not the comparison vignette. Activate this additional output
#' by setting `dfMaxLength=NA` (or assign a units object to replace the default 50m)
# supply dfMaxLength to get additional output
get_meuse = function(dfMaxLength=NULL)
{
  # dfMaxLength = (optional) interval for river line geometries re-sampling
  # see below in comments

  # EPSG code for the coordinate system
  epsg_meuse = 28992
  crs_meuse = sf::st_crs(epsg_meuse)[['wkt']]

  # load soil points data
  data(meuse, package='sp', envir = environment())
  meuse_soils = sf::st_as_sf(meuse, coords=c('x', 'y'), crs=epsg_meuse)

  # add 'logzinc' column
  meuse_soils[['log_zinc']] = log(meuse_soils[['zinc']])

  # this argument triggers a more complicated multi-layer output that uses river geometry
  if( is.null(dfMaxLength) ) { return(meuse_soils['log_zinc']) } else {

    # dfMaxLength sets the interval used to sample line geometries of the river
    # using Voronoi tiles. This is a fussy and not well-tested algorithm for finding the
    # centre line of a river polygon, but it seems to work well enough for the example here

    # set default tuning parameters - 50m seems to work well in Meuse example
    if( is.na(dfMaxLength) ) dfMaxLength = units::set_units(50, m)

    # open river location data
    data(meuse.riv, package='sp', envir = environment())

    # reshape the river (edge) point data as a more densely segmented polygon
    colnames(meuse.riv) = c('x', 'y')
    meuse_river_points = sf::st_as_sf(as.data.frame(meuse.riv), coords=c('x', 'y'), crs=crs_meuse)
    meuse_river_seg = sf::st_cast(sf::st_combine(meuse_river_points), 'LINESTRING')
    meuse_river_poly = sf::st_cast(sf::st_segmentize(meuse_river_seg, dfMaxLength), 'POLYGON')

    # skeletonization trick to get a single linestring at center of the river
    meuse_river_voronoi = sf::st_cast(sf::st_voronoi(meuse_river_poly, bOnlyEdges=TRUE), 'POINT')
    meuse_river_skele = sf::st_intersection(meuse_river_voronoi, meuse_river_poly)
    n_skele = length(meuse_river_skele)

    # compute distance matrix
    dmat_skele = units::drop_units(sf::st_distance(meuse_river_skele))

    # re-order to start from northernmost point
    idx_first = which.max(st_coordinates(meuse_river_skele)[,2])
    idx_reorder = c(idx_first, integer(n_skele-1L))
    for(idx_skele in seq(n_skele-1L))
    {
      # find least distance match
      idx_tocheck = seq(n_skele) != idx_first
      idx_first = which(idx_tocheck)[ which.min(dmat_skele[idx_tocheck, idx_first]) ]
      idx_reorder[1L+idx_skele] = idx_first

      # modify distance matrix so the matching point is not selected again
      dmat_skele[idx_first, ] = Inf
    }

    # connect the points to get the spine
    meuse_river = sf::st_cast(sf::st_combine(meuse_river_skele[idx_reorder]), 'LINESTRING')

    # add 'distance' to river to soils data frame
    meuse_soils[['distance']] = units::drop_units( sf::st_distance(meuse_soils, meuse_river))

    # crop the river objects to buffered bounding box of soils data
    bbox_padded = st_buffer(sf::st_as_sfc(sf::st_bbox(meuse_soils)), units::set_units(500, m))
    meuse_river_poly = sf::st_crop(meuse_river_poly, bbox_padded)
    meuse_river = sf::st_crop(meuse_river, bbox_padded)

    # return three geometry objects in a list
    return( list(soils=meuse_soils, river_poly=meuse_river_poly, river_line=meuse_river) )
  }
}
