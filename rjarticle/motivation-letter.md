---
output: pdf_document
fontsize: 12pt
---

\thispagestyle{empty}
\today

Catherine Hurley  
Editor-in-Chief   
The R Journal  
\bigskip

Dear Professor Hurley,
\bigskip

Please consider our article "snapKrig: An R Package for Fast Geostatistics with Kronecker Covariances" for publication in the R Journal. This introduces the `snapKrig` package, newly published on CRAN, which does likelihood-based model fitting, simulation, inference, and prediction for 2-dimensional spatial data lying on a grid.

snapKrig uses unconventional methods for exact computations, making it substantially faster than standard alternative packages like `gstat`, `fields`, `geoR`, and `spatial`. This is of particular importance for kriging prediction (interpolation), a very common problem in earth sciences which can be computationally prohibitive even for moderately large raster datasets.

Our motivation in developing this package was to reduce computation time when interpolating weather station data and weather forecast grids onto a common high resolution grid. However we think `snapKrig` will be useful in a wide range of applications. We illustrate this in the manuscript by including three different real-life data examples, looking at soil contamination, ozone concentration, and forest health.

The manuscript also reports on a benchmarking experiment where we compared computation time for kriging by `gstat`, `fields`, `geoR`, and `snapKrig`. This highlights situations where snapKrig can be faster by *several orders of magnitude*, with similarly reduced memory requirements. Our submission includes a pair of well-documented R scripts for reproducing this workflow, and for repeating our model-selection. We also include a "readme.txt" with instructions for readers who wish to run these scripts.

We have worked hard to make `snapKrig` both feature-rich, user-friendly, and instructive. The package defines its own S3 class for grids, handling conversions to and from common alternatives (like `terra`'s SpatRaster), and it offers a wide variety of methods for grid objects including a custom plot method. 

We also provide a highly optimized likelihood function so that users can exploit our computational methods on related problems, like Bayesian MCMC for hierarchical models. We are currently exploring this idea in a separate paper. Our hope is that an R Journal publication might similarly inspire others to research new applications for our computational methods.


\bigskip

Regards,
    
Dean Koch  
Department of Mathematics and Statistics  
University of Alberta, Canada    
dkoch@ualberta.ca  



\bigskip

