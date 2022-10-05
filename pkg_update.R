# run me once in a while!

# A list of packages connected in some way to the project.
pkg_to_update = c('raster', 'terra', 'sf', 'sp', 'Rcpp',
                  'gstat', 'geoR', 'fields', 'microbenchmark', 'R.utils',
                  'devtools', 'sloop', 'ggplot2', 'gridExtra',
                  'here', 'rmarkdown', 'smoothr', 'Rtools')

# not listed: RandomFields (not available for newer R version) and rjtools (see below)

# installs/updates
install.packages(pkg_to_update)

# rjtools requires Rtools, which is installed from a binary here:
# https://cran.r-project.org/bin/windows/Rtools/

# once this is up to date you can then run
remotes::install_github('rjournal/rjtools')
