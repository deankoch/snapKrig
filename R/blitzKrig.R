#' ## INTRODUCTION
#'
#' A "bk" object is just a list of vectors representing a grid and its data, similar to
#' (but much simpler than) the "raster" and "SpatRaster" classes in terra and raster.
#'
#' Our S3 class doesn't accomplish anything performance-wise. We use it to link blitzKrig
#' grid methods to R's many generic functions (`print`, `plot`, etc) in a sensible way, and
#' to point R's internal generic functions (like `[`, `[<-`, and `length`) to the contents
#' of `bk[['gval']]`, so that "bk" objects can behave like vectors.
#'
#' At minimum, a "bk" object contains three entries: `gdim`, `gres`, and `gyx`, defining the
#' grid dimensions, spacing, and extent. All are given in the order y, x. This is so that if
#' we view observations at grid points as *matrix* data, then `gdim` will be consistent with
#' `base::dim`.
#'
#' Optionally, geo-referenced data can be accompanied by metadata describing how the grid
#' is mapped to the globe. This goes in the `crs` entry, and is copied over automatically
#' when importing from `sf`, `terra`, and `raster`.
#'
#' * `crs`: character, the WKT representation of the CRS for the grid
#'
#' Observations at some or all of the grid points are stored in the entry `gval`
#'
#' * `gval`: numeric vector or matrix, the grid data
#'
#' In the single-layer case this is a vector with as many entries as there are grid points.
#' We use column-vectorized ordering, which stacks the columns of the data matrix (with the
#' left-most, or first column appearing first in the vector). This is the ordering we get
#' when coercing a matrix to vector with `base::as.vector`, for example.
#'
#' In the multi-layer case we have one such vector per layer, and these are stored as
#' columns of a matrix. For example, the column vectorization of the first layer is the
#' vector `gval[,1]`. **It is assumed that each layer has the same NA structure**.
#'
#' To save memory, when `gval` is a matrix, "bk" objects use a sparse representation that
#' omits NAs. This means the matrix `gval` only stores the observed data values, so it will
#' have as many rows as there are observed grid points. This requires an additional
#' indexing vector:
#'
#' * `idx_grid`: length-n numeric vector mapping rows of `gval` to grid points
#'
#' Users can supply (the shortened) `gval` matrix together with the corresponding,
#' `idx_grid`, or just pass the complete `gval` matrix (with NAs) on its own, and `bk`
#' will do the indexing and simplification for you.
#'
#' ## CREATION
#'
#' Typical usage is to pass a grid-like object to the helper function `bk`, which extracts
#' the list entries discussed above and passes them to `bk_make`, the constructor, then
#' `bk_validate` (for sanity checking, and to fill in missing entries).
#'
