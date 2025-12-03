
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tosccamm

<!-- badges: start -->

![Status:
Finished](https://img.shields.io/badge/status-finished-brightgreen)
![Build
Status](https://github.com/nuria-sv/toscca-mm/actions/workflows/r.yml/badge.svg)
![License](https://img.shields.io/github/license/nuria-sv/tosccamm)
<!-- badges: end -->


_tosccamm_ is the package to implement the Thresholded Ordered Sparse CCA for Multiple Measurements (TOSCCA-MM) method described in [Senar et al. (2025)](https://doi.org/10.1002/bimj.70090).

## Installation

NOT READY

You can install the development version of tosccamm like so:

``` r
devtools::install_github("nuria-sv/toscca")    # dependencies
devtools::install_github("nuria-sv/tosccamm")
```

## TOSCCA-MM

<img src="diagram/tosccamm_tikz.png" height="350" width="700" align="center">

TOSCCA-MM is a novel extension of sparse CCA that incorporates time
dynamics at the latent variable level through longitudinal models, such
as autoregressive models or linear mixed effect models. This approach
addresses the correlation of repeated measurements while drawing latent
paths, for each component. To aid interpretability and computational
efficiency, we implement an $\ell_0$ penalty to enforce fixed sparsity
levels. We estimate these trajectories fitting longitudinal models to
the low-dimensional latent variables, (i.e.: linear mixed effects
model). By leveraging the clustered structure of high-dimensional
datasets, we are able to explore the shared longitudinal latent
mechanisms. The sparse canonical weights, yield interpretable outcomes
on variable contribution to the estimated correlated trajectories.
Furthermore, modelling time in the latent space significantly reduces
computational burden.

## Example

This is a basic example over simulated data of TOSCCA-MM

``` r
# library(tosccamm)
source("C:/Users/PC/OneDrive/github/tosccamm/R/tosccam_permut.R")
source("C:/Users/PC/OneDrive/github/tosccamm/R/tosccamm_core.R")
source("C:/Users/PC/OneDrive/github/tosccamm/R/tosccamm_folds.R")
source("C:/Users/PC/OneDrive/github/tosccamm/R/toscca_helpers.R")
source("C:/Users/PC/OneDrive/github/tosccamm/R/general_functions.R")
# for plots
library(grid)
library(ggplot2)
library(gridExtra)
library(viridis)
#> Loading required package: viridisLite
```

Estimate the canonical weights and latent paths for $K$ components.

``` r
res_k = list()

X.temp = XX2
Y.temp = YY2
for (k in 1:5) {
  if(k > 1) {
    # residualise for subsequent components
    X.temp = data.frame(X.temp[,c(1,2)],toscca::residualisation(as.matrix(X.temp[,-c(1,2)]), res_k[[k-1]]$alpha, type = "basic") )
    Y.temp = data.frame(Y.temp[,c(1,2)],toscca::residualisation(as.matrix(Y.temp[,-c(1,2)]), res_k[[k-1]]$beta, type = "basic") )

    nz_a_gen = as.numeric(table(res_k[[k-1]]$alpha != 0)[2])
    nz_b_gen = as.numeric(table(res_k[[k-1]]$beta != 0)[2])
  }

  res_k[[k]] <- tosccamm(X.temp, Y.temp, folds = 2,
                                            nonzero_a = nonz_a, nonzero_b = nonz_b,
                                            model = "lme", lmeformula = " ~ 0 + poly(time,3) + (1|id)")

}
#>  Common convergence error: 0 & Iterations: 5  Common convergence error: 0 & Iterations: 5 
#> k-fold cv max. cancor 
#>             0.6347814 
#> 
#>  ........................................ 
#>  # nonzero A: 10
#>  # nonzero B: 5
#>  ........................................ 
#>  Common convergence error: 0 & Iterations: 14  Common convergence error: 0.00814 & Iterations: 21 
#> k-fold cv max. cancor 
#>             0.4786131 
#> 
#>  ........................................ 
#>  # nonzero A: 15
#>  # nonzero B: 50
#>  ........................................ 
#>  Common convergence error: 0.0432 & Iterations: 21  Common convergence error: 0.06226 & Iterations: 21 
#> k-fold cv max. cancor 
#>             0.2137972 
#> 
#>  ........................................ 
#>  # nonzero A: 45
#>  # nonzero B: 50
#>  ........................................ 
#>  Common convergence error: 0.02561 & Iterations: 21  Common convergence error: 0.0193 & Iterations: 21 
#> k-fold cv max. cancor 
#>             0.2098927 
#> 
#>  ........................................ 
#>  # nonzero A: 30
#>  # nonzero B: 16
#>  ........................................ 
#>  Common convergence error: 0.01477 & Iterations: 21  Common convergence error: 0.04615 & Iterations: 21 
#> k-fold cv max. cancor 
#>             0.2233089 
#> 
#>  ........................................ 
#>  # nonzero A: 40
#>  # nonzero B: 39
#>  ........................................
```

### Results

#### Latent paths for $k=1$ and $k=2$

<img src="man/figures/README-plots-1.png" width="100%" />

#### Canonical weights for $k=1$ and $k=2$

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

#### Latent path and canonical weights for $k=3$, noise

<img src="man/figures/README-plotNoise-1.png" width="100%" />

<img src="man/figures/README-gridPlots-1.png" width="100%" /><img src="man/figures/README-gridPlots-2.png" width="100%" /><img src="man/figures/README-gridPlots-3.png" width="100%" />

    #>  Common convergence error: 0 & Iterations: 5  Common convergence error: 0 & Iterations: 5 
    #> k-fold cv max. cancor 
    #>             0.6347814 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 10
    #>  # nonzero B: 5
    #>  ........................................ 
    #>  Common convergence error: 0 & Iterations: 8  Common convergence error: 0 & Iterations: 11 
    #> k-fold cv max. cancor 
    #>              0.431531 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 10
    #>  # nonzero B: 5
    #>  ........................................ 
    #>  Common convergence error: 0.00695 & Iterations: 21  Common convergence error: 0.09078 & Iterations: 21 
    #> k-fold cv max. cancor 
    #>            0.08403192 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 10
    #>  # nonzero B: 5
    #>  ........................................ 
    #>  Common convergence error: 1e-05 & Iterations: 21  Common convergence error: 0.1628 & Iterations: 21 
    #> k-fold cv max. cancor 
    #>             0.1082644 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 10
    #>  # nonzero B: 5
    #>  ........................................ 
    #>  Common convergence error: 0.2184 & Iterations: 21  Common convergence error: 2e-05 & Iterations: 21 
    #> k-fold cv max. cancor 
    #>            0.07707604 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 10
    #>  # nonzero B: 5
    #>  ........................................
    #> Loading required package: iterators
    #> Loading required package: parallel
    #> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    #> ℹ Please use `linewidth` instead.
    #> This warning is displayed once every 8 hours.
    #> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    #> generated.
    #> Warning: The dot-dot notation (`..count..`) was deprecated in ggplot2 3.4.0.
    #> ℹ Please use `after_stat(count)` instead.
    #> This warning is displayed once every 8 hours.
    #> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    #> generated.

<img src="man/figures/README-permutationTesting-1.png" width="100%" />

    #> Empirical p-values:
    #> 0
    #> 0
    #> 0.34
    #> 0.213
    #> 0.395
    #> NULL
