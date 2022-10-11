# snapKrig

an R package for computationally fast grid-based Gaussian processes models and kriging

<img src="https://raw.githubusercontent.com/deankoch/snapKrig/master/vignettes/meuse_vignette_files/figure-gfm/ordinary_kriging-1.png" width="30%"></img>
<img src="https://raw.githubusercontent.com/deankoch/snapKrig/master/vignettes/meuse_vignette_files/figure-gfm/predictor_plot-1.png" width="30%"></img>
<img src="https://raw.githubusercontent.com/deankoch/snapKrig/master/vignettes/meuse_vignette_files/figure-gfm/variance_plot-1.png" width="30%"></img>


# Update October 2022

This is a continuation of the [`pkern`](https://github.com/deankoch/pkern) package, just with a catchier
name and some updates. Development on `pkern` has ended, with `snapKrig` picking up where it left off.

The most important update in this version was to define an S3 class, "sk" for grid list objects (which are
usually named `g` in examples). The base class for "sk" still a list, and square-bracket replace/access
calls with entry names, like `g[['gdim']]` and `g[c('gdim', 'gres')]`, still behave the same way. However,
grid data values can now be accessed/assigned directly via square-bracket indexing, such as
with `g[i]` instead of `g[["gval"]][i]`. For numeric indices, `g` now behaves like a numeric vector.

We also now have a large number of methods defined for "sk", including common generics like `summary`,
`plot`, `length`, `range`, `mean`, `sum`, and much more. Group generics like `exp`, `*`, and `>` also have
methods for `g` and will return a "sk" object directly (no more fussing with `modifyList`).

# Coming soon

A CRAN submission is coming up soon. For now you can install the package using devtools and try out
the [Meuse River vignette](https://github.com/deankoch/snapKrig/blob/master/vignettes/meuse_vignette.md).
An paper introducing `snapKrig` has also been drafted for submission to the R Journal around the same time.


# Overview

`snapKrig` provides a computationally lean implementation of a 2-dimensional spatial correlation model for
gridded data. This can be useful when working with geo-referenced data, such as in earth sciences, where 
users often need a tool for interpolation or down-scaling.

The package offers an fast and simple back-end for modeling with spatially correlated errors.
It works much faster than alternatives like `gstat` and `fields` (sometimes by several orders of magnitude),
at the expense of slightly restricting the type of model users can select, and forcing the use of gridded
datasets (possibly requiring the data to be snapped).

`snapKrig` depends only on packages included by default in R (like `graphics` and `stats`), but supports 
raster and geometry classes from `sf` and `terra`. While we recommend the newest version of R, the package
is written for compatibility with older versions of R (*eg.* we steer clear of syntax like `|>` and `\(x)`).


# Technical Features

* models anisotropic Gaussian processes on 2-dimensional regular grids for a choice of covariance kernels
* optimized computation of the likelihood function, generalized least squares, and kriging predictor/variance
* fast computations with missing data problems, and even faster in the complete data case 
* automated maximum likelihood model fitting and support for sample semi-variograms
* user friendly helper functions for raster down-scaling and point interpolation

# Technical Background

This is an R implementation of some methods I developed in [my thesis](https://doi.org/10.7939/r3-91zn-v276)
for speeding up geostatistical computations involving large covariance matrices. The central idea is to model
spatial dependence using a separable 2-dimensional covariance kernel, defined as the product of two (1-dimensional)
univariate covariance kernels. This introduces special symmetries and structure in the covariance matrix, which are
exploited in this package for fast and memory-efficient computations.

This package will accompany a paper on fast and efficient downscaling of weather grids so the focus is on a particular
application of kriging, but the underlying methods are applicable much more generally. See [[1](https://doi.org/10.7939/r3-g6qb-bq70)],
where I use product kernels to study directions of anisotropy in a nonstationary random fields, and
[[2](https://doi.org/10.1007/s11538-021-00899-z), [3](https://doi.org/10.1098/rsif.2020.0434)], where I apply it to
fit acovariance structure, and to speed up calculations of dispersal kernel convolutions.

