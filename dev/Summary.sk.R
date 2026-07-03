#' Summary group generics
#'
#' Computes a summary statistic of the grid data values.
#'
#' @param ... further arguments to Summary
#' @param na.rm logical, whether to remove `NA`s before doing the computation
#'
#' @return a "sk" object with data values
#' @export
#'
#' @examples
#'
#' # make test example and verify common summary stats
#' g = sk(gval=1+abs(stats::rnorm(4^2)), gdim=4, gres=0.5)
#' sum(g) == sum(g[])
#' min(g) == min(g[])
#' max(g) == max(g[])
#'
#' # list input can include scalar, vector (but first argument must be sk class)
#' sum(g, -g)
#' sum(g, -g, 0)
#'
#' # calculate product in two ways
#' prod(g) - exp(sum(log(g)))
#'
#' # check logical conditions
#' any(g>0)
#' all(g>0)
#'
#' # count the number of grid points satisfying a condition
#' sum(g < 2)
#'
Summary.sk = function(..., na.rm=FALSE)
{
  # extract all grid data values and omit NA's if requested
  vals = do.call(c, lapply(list(...), function(x) x[]))
  if(na.rm) vals = vals[ !is.na(vals) ]
  if( length(vals) == 0 ) stop('all values are NA')

  # dispatch to default method then replace the data values in output sk object
  .Primitive(.Generic)(vals)
}
