
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

Armin Rauschenberger
[![AR](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-6498-4801),
Iuliana Ciocănea-Teodorescu
[![ICT](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-2489-9025),
Renée X. Menezes
[![RXM](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-1378-2721),
Marianne A. Jonker
[![MAJ](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-0134-8482),
and Mark A. van de Wiel
[![MvdW](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-4780-8472)
(2020). "Sparse classification with paired covariates".
*Advances in Data Analysis and Classification* 14:571-588.
[doi: 10.1007/s11634-019-00375-6](https://doi.org/10.1007/s11634-019-00375-6)

[![CRAN
version](https://www.r-pkg.org/badges/version/palasso)](https://CRAN.R-project.org/package=palasso)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/palasso)](https://CRAN.R-project.org/package=palasso)
[![Total CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/palasso)](https://CRAN.R-project.org/package=palasso)
