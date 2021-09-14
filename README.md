
# sGBJ

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/sGBJ)](https://CRAN.R-project.org/package=sGBJ)
[![R-CMD-check](https://github.com/lauravillain/sGBJ/workflows/R-CMD-check/badge.svg)](https://github.com/lauravillain/sGBJ/actions)
<!-- badges: end -->

The goal of `sGBJ` is to provide an extension of the Generalized Berk-Jones test for set-based inference in genomics studies with a censored time-to-event outcome. In particular, it is geared towards 

This packages implements an extension of the GBJ statistic for censored time-to-event data. It computes the sGBJ statistic and its p-value for testing the association between a gene set and a time-to-event outcome with possible adjustment on additional covariates.

The main function of the package is `sGBJ()`

The methods implemented in this package are detailed in the following article:

> Villain L, Ferté T, Thiébaut R, Hejblum BP (2021). Gene Set Analysis for time-to-event outcome with the Generalized Berk–Jones statistic. *bioRxiv* 2021.09.07.459329.  
DOI: [10.1101/2021.09.07.459329](https://www.biorxiv.org/content/10.1101/2021.09.07.459329)

## Installation

You can install the released version of sGBJ from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sGBJ")
```

Laura Villain, Thomas Ferté, Rodolphe Thiébaut & Boris Hejblum

