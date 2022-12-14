% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sk_pars.R
\name{sk_bds}
\alias{sk_bds}
\title{Set default parameter covariance parameter bounds for a Kronecker covariance model}
\usage{
sk_bds(pars, g, var_obs = NULL, var_mult = 2)
}
\arguments{
\item{pars}{list or character vector of 1-2 kernel names (see \code{sk_pars})}

\item{g}{a sk grid (or any object accepted by \code{sk})}

\item{var_obs}{positive numeric, the sample variance of data \code{g$gval}}

\item{var_mult}{numeric > 1, constant to multiply by \code{var_obs} to get upper bounds}
}
\value{
a data frame of initial values and lower/upper bounds for the parameters in \code{pars}
}
\description{
Returns a data-frame of initial values and upper/lower bounds on covariance
parameters for the Kronecker covariance model defined by the correlation function
names in \code{pars}.
}
\details{
Range parameters (\code{y.rho} and \code{x.rho}) are bounded by the shortest and longest
inter-point distances along the corresponding dimension (y or x). This is
computed by taking the element-wise product of dimensions and resolution, ie
\code{g$gres * dim(g)}. Ranges are initialized to the geometric mean of the upper
and lower bounds.

Variance bounds are centered around \code{var_obs}, which by default is set to the sample
variance of the data in \code{g}. \code{eps} (measurement variance) and \code{psill} (partial
sill) are both initialized to one half of \code{var_obs}, bounded above by \code{var_obs}
times \code{var_mult}, and bounded below by a small positive number (\code{1e-6}). Note that
while \code{eps=0} produces valid models in theory, in practice \code{eps>0} is often
necessary for numerical stability.

Shape parameter bounds are hard-coded, and are set conservatively to avoid problems
with numerical precision in functions like \code{exp} and \code{gamma} when evaluating very
large or small distances.
}
\examples{
gdim = c(10, 15)
g = sk(gdim)
g[] = stats::rnorm(length(g))
sk_bds('mat', g)

# same result by passing in observed variance
sk_bds('mat', g, stats::var(g[]))

# a less conservative bound for variance (only eps and psill affected)
sk_bds('mat', g, var_mult=1)
}
\seealso{
sk

Other parameter managers: 
\code{\link{sk_fit}()},
\code{\link{sk_kp}()},
\code{\link{sk_pars_make}()},
\code{\link{sk_pars_update}()},
\code{\link{sk_pars}()},
\code{\link{sk_to_string}()}
}
\concept{parameter managers}
\keyword{internal}
