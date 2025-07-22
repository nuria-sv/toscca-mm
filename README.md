
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tosccamm

<!-- badges: start -->

<!-- badges: end -->

The goal of tosccamm is to …

## Installation

You can install the development version of tosccamm like so:

``` r
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
#>  Common convergence error: 0 & Iterations: 5  Common convergence error: 0 & Iterations: 6 
#> k-fold cv max. cancor 
#>             0.4191102 
#> 
#>  ........................................ 
#>  # nonzero A: 10
#>  # nonzero B: 5
#>  ........................................ 
#>  Common convergence error: 2e-05 & Iterations: 21  Common convergence error: 0.00034 & Iterations: 21 
#> k-fold cv max. cancor 
#>             0.3968518 
#> 
#>  ........................................ 
#>  # nonzero A: 10
#>  # nonzero B: 28
#>  ........................................ 
#>  Common convergence error: 0.0324 & Iterations: 21  Common convergence error: 0.0282 & Iterations: 21 
#> k-fold cv max. cancor 
#>             0.1791438 
#> 
#>  ........................................ 
#>  # nonzero A: 10
#>  # nonzero B: 5
#>  ........................................ 
#>  Common convergence error: 0.04798 & Iterations: 21  Common convergence error: 0.03029 & Iterations: 21 
#> k-fold cv max. cancor 
#>             0.3115394 
#> 
#>  ........................................ 
#>  # nonzero A: 5
#>  # nonzero B: 50
#>  ........................................ 
#>  Common convergence error: 0.06816 & Iterations: 21  Common convergence error: 0.06314 & Iterations: 21 
#> k-fold cv max. cancor 
#>             0.3007874 
#> 
#>  ........................................ 
#>  # nonzero A: 5
#>  # nonzero B: 16
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

    #>  Common convergence error: 0 & Iterations: 4  Common convergence error: 0 & Iterations: 6 
    #> k-fold cv max. cancor 
    #>             0.4191102 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 10
    #>  # nonzero B: 5
    #>  ........................................ 
    #>  Common convergence error: 0 & Iterations: 9  Common convergence error: 0 & Iterations: 21 
    #> k-fold cv max. cancor 
    #>             0.1683788 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 10
    #>  # nonzero B: 5
    #>  ........................................ 
    #>  Common convergence error: 0 & Iterations: 9  Common convergence error: 0 & Iterations: 16 
    #> k-fold cv max. cancor 
    #>             0.3118109 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 10
    #>  # nonzero B: 5
    #>  ........................................ 
    #>  Common convergence error: 0 & Iterations: 20  Common convergence error: 0.09223 & Iterations: 21 
    #> k-fold cv max. cancor 
    #>             0.0630195 
    #> 
    #>  ........................................ 
    #>  # nonzero A: 10
    #>  # nonzero B: 5
    #>  ........................................ 
    #>  Common convergence error: 0.01518 & Iterations: 21  Common convergence error: 0.00273 & Iterations: 21 
    #> k-fold cv max. cancor 
    #>           0.007440166 
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
    #> 0.062
    #> 0.001
    #> 0.476
    #> 0.969
    #> NULL
