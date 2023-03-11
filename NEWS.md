## 2022-12-19 Version 0.0.2

### bug fixes

* `anyNA` now reports the correct result instead of its negation
* The log-likelihood plot in the vignette now uses `reset=FALSE` so the circle is visible

### new

* assignments with scalars now recycle values. Example `g[] = 1` is now equivalent to `g[] = rep(1, length(g))`

* `sk` now accepts calls of the form `sk(g, ...)` where `g` is an existing `sk` object to modify with the named elements in `...`. For example `sk(g, gval=1))` is equivalent to `g[] = 1`

## 2022-12-19 Version 0.0.1

New submission to CRAN
