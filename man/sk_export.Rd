% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sk_export.R
\name{sk_export}
\alias{sk_export}
\title{Convert "sk" grid to SpatRaster}
\usage{
sk_export(g, template = "terra")
}
\arguments{
\item{g}{any object accepted or returned by \code{sk}}

\item{template}{character or RasterLayer/SpatRaster to set output type

Converts a vector or matrix to a SpatRaster or RasterLayer. Multi-layer outputs are
supported for terra but not raster.}
}
\value{
a RasterLayer or SpatRaster containing the data from \code{g} (or a sub-grid)
}
\description{
Convert "sk" grid to SpatRaster
}
\examples{
if( requireNamespace('raster') ) {

# open example file as RasterLayer
r_path = system.file('ex/logo.tif', package='terra')
r = raster::raster(r_path, band=1)
g = sk(r)

# convert back to RasterLayer and compare
r_from_g = sk_export(g, 'raster')
print(r)
print(r_from_g)

# NOTE: layer name, band number, and various other metadata are lost
all.equal(r_from_g, r)

}

# same with terra
if( requireNamespace('terra') ) {

# convert all layers
r = terra::rast(r_path)
g = sk(r)
r_from_g = sk_export(g)

# NOTE: various metadata are lost
all.equal(r_from_g, r)

}

}
\seealso{
sk

Other exporting functions: 
\code{\link{sk_coords}()}
}
\concept{exporting functions}
