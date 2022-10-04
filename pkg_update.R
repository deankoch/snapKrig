# run me once in a while!

# A list of packages connected in some way to the project.
pkg_to_update = c('raster', 'terra', 'sf', 'sp', 'Rcpp',
                  'gstat', 'geoR', 'fields', 'RandomFields', 'microbenchmark', 'R.utils',
                  'devtools', 'sloop', 'ggplot2', 'gridExtra')

# installs/updates
install.packages(pkg_to_update)
