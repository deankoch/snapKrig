% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sk_index.R
\name{sk_sub_find}
\alias{sk_sub_find}
\title{Find complete regular sub-grids in a sk grid object}
\usage{
sk_sub_find(g, gdim = NULL)
}
\arguments{
\item{g}{logical vector, sk grid, or any grid object accepted by \code{sk}}

\item{gdim}{integer vector, the grid dimensions (in order 'y', 'x')}
}
\value{
\code{NULL} or list of information about the location and spacing of the sub-grid
within \code{g} (see details)
}
\description{
If a sk grid \code{g} has missing values (\code{NA}s) but the set of non-\code{NA} points form a
complete (regular) sub-grid, this function finds its grid lines, resolution, and
dimensions. If no eligible sub-grids are found, the function returns \code{NULL}.
}
\details{
A sub-grid is only eligible if it contains ALL of the non-\code{NA} points in \code{g} and none
of the \code{NA}s. For example if a single point missing from the sub-grid, or a single non-\code{NA}
point lies outside the sub-grid, the function will fail to detect any sub-grids and return
\code{NULL}. If no points are \code{NA}, the function returns indices for the full grid.

The returned list contains the following named elements:
\itemize{
\item \code{ij} the grid line numbers of the sub-grid with respect to \code{g}
\item \code{res_scale} the resolution scaling factor (relative increase in grid line spacing of \code{g})
\item \code{gdim} the number of y and x grid lines in the sub-grid
}

As in \code{sk}, each of these is given in the 'y', 'x' order.

Users can also pass the logical vector returned by \code{!is.na(g)} instead of \code{g}, in which
case argument \code{gdim} must also be specified. This can be much faster with large grids.
}
\examples{

# define a grid and example data
gdim = c(50, 53)
pars = utils::modifyList(sk_pars(gdim), list(eps=1e-12))
g = sk_sim(gdim, pars)
plot(g)

# define a super-grid containing the original data and make sure we can find it
g_big = sk_rescale(g, down=3)
plot(g_big)
print(sk_sub_find(g_big))

# define a smaller sub-grid at random
spacing = sapply(floor(gdim/10), function(x) 1 + sample.int(x, 1))
gdim_sg = sapply(floor( (gdim - 1) / spacing), function(x) sample.int(x, 1))
ij_first = sapply(gdim - ( spacing * gdim_sg ), function(x) sample.int(x, 1))

# find index of sub-grid lines and vectorized index of points
ij_sg = Map(function(idx, r, n) seq(idx, by=r, length.out=n), idx=ij_first, r=spacing, n=gdim_sg)
names(ij_sg) = c('i', 'j')
is_sg = sk_sub_idx(gdim, ij_sg, idx=FALSE)

# assign values to the sub-grid points
g_sub = sk(gdim)
g_sub[is_sg] = g[is_sg]
plot(g_sub, zlab='sub-grid')

# call the function and check for expected results
sub_result = sk_sub_find(g_sub)
all.equal(unname(sub_result[['gdim']]), gdim_sg)
all.equal(unname(sub_result[['ij']]), unname(ij_sg))

# sub grids with side length 1 have no spacing defined along that dimension
spacing[gdim_sg==1] = NA

# check consistency in spacing
all.equal(unname(sub_result[['res_scale']]), spacing)

# can also call on the vector and supply gdim separately
identical(sub_result, sk_sub_find(!is.na(g_sub), dim(g_sub)))

}
\seealso{
Other indexing functions: 
\code{\link{sk_mat2vec}()},
\code{\link{sk_rescale}()},
\code{\link{sk_sub_idx}()},
\code{\link{sk_vec2mat}()}
}
\concept{indexing functions}
\keyword{internal}
