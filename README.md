
<!-- README.md is generated from README.Rmd. Please edit that file -->

# INLAcomp

<!-- badges: start -->

[![R-CMD-check](https://github.com/jmartinez-minaya/INLAcomp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jmartinez-minaya/INLAcomp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of **INLAcomp** is to analyze compositional data with Logistic
Normal distribution regression with Dirichlet residuals (LNDM) using the
integrated nested Laplace approximation via the [R-INLA
package](https://www.r-inla.org/). Method is in [Martinez-Minaya and
Rue,
2024](https://link.springer.com/article/10.1007/s11222-024-10427-3).

- For a **tutorial about how to fit a simple LNDM model using R-INLA**,
  click
  [here](https://github.com/jmartinez-minaya/INLAcomp/blob/main/vignettes/Dirichlet-CoDa.Rmd)
  or
  [here](https://github.com/hrue/r-inla/commit/3577c1b030a8460ff7194893bd97e57a62a1399d)

- The Spatial LNDM simulation conducted in the paper is available
  [here](https://github.com/jmartinez-minaya/INLAcomp/tree/main/simulations).
  The generated html is availabe
  [here](https://jmartinez-minaya.github.io/en/supplementary/supplementary/INLAComp/simulations.html)

- The code for the Real example depicted in the paper is available in
  the
  [vignette](https://github.com/jmartinez-minaya/INLAcomp/blob/main/vignettes/my-vignette.Rmd).
  The generated html is available
  [here](https://jmartinez-minaya.github.io/en/supplementary/supplementary/INLAComp/my-vignette.html)

## Installation

It is not yet in CRAN, but you can install the latest bugfix release of
**INLAcomp** from [github](https://github.com/jmartinez-minaya/INLAcomp)
with:

``` r
remotes::install_github("https://github.com/jmartinez-minaya/INLAcomp")
```
