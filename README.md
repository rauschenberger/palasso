
<!--
[![Travis-CI Build
Status](https://travis-ci.org/rauschenberger/palasso.svg)](https://travis-ci.org/rauschenberger/palasso)
-->

[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/rauschenberger/palasso?svg=true)](https://ci.appveyor.com/project/rauschenberger/palasso)
[![R-CMD-check](https://github.com/rauschenberger/palasso/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rauschenberger/palasso/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/rauschenberger/palasso/graph/badge.svg)](https://app.codecov.io/gh/rauschenberger/palasso)

## Scope

Sparse regression with paired covariates (extending
[glmnet](https://CRAN.R-project.org/package=glmnet)).

## Installation

Install the current release from
[CRAN](https://CRAN.R-project.org/package=palasso), or the latest
development version from
[GitHub](https://github.com/rauschenberger/palasso):

``` r
install.packages("palasso")
```

``` r
#install.packages("remotes")
remotes::install_github("rauschenberger/palasso")
```

Optionally, install [ashr](https://github.com/stephens999/ashr) and
[CorShrink](https://github.com/kkdey/CorShrink) for the correlation
shrinkage:

``` r
remotes::install_github("stephens999/ashr")
remotes::install_github("kkdey/CorShrink")
```

## Reference

A Rauschenberger, I Ciocanea-Teodorescu, RX Menezes, MA Jonker, and MA
van de Wiel (2019). “Sparse classification with paired covariates”.
*Advances in Data Analysis and Classification.* In press.
[doi: 10.1007/s11634-019-00375-6](https://doi.org/10.1007/s11634-019-00375-6)

[![CRAN
version](https://www.r-pkg.org/badges/version/palasso)](https://CRAN.R-project.org/package=palasso)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/palasso)](https://CRAN.R-project.org/package=palasso)
[![Total CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/palasso)](https://CRAN.R-project.org/package=palasso)
