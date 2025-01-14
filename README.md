
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tosccamm

<!-- badges: start -->
<!-- badges: end -->

The goal of tosccamm is to …

## Installation

You can install the development version of tosccamm like so:

``` r
devtools::install_github("nuria-sv/toscamm")
```

## TOSCCA-MM

<img src="vignettes/tosccamm_tikz.png" height="350" width="700" align="center">

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
library(tosccamm)
## basic example code
```
