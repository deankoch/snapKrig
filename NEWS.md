# snapKrig (development version)

# snapKrig 0.0.3

*2026-07-03*

* add missing dependency on R >= 4.1.0 (for pipe)

* update maintainer email contact address

# snapKrig 0.0.2

*2022-05-05*

* new features

  * assignments with scalars now recycle values. Example `g[] = 1` is now equivalent to `g[] = rep(1, length(g))`

  * `sk` now accepts calls of the form `sk(g, ...)` where `g` is an existing `sk` object to modify with the named elements in `...`. For example `g = sk(g, gval=1))` is equivalent to `g[] = 1`

  * added some new example applications to the github repo
  
* bugfixes

  * `anyNA` now reports the correct result instead of its negation

  * single-bracket extract with logical index now returns `NA` vector of correct length 

  * fixed issue with multi-layer input to `sk_export` losing all but first layer

  * fixed issue with pars argument to `sk_fit` not accepting lists

  * The log-likelihood plot in the vignette now uses `reset=FALSE` so the circle is visible


## Version 0.0.1

*2022-12-19*

* new submission to CRAN
