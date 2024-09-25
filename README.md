
<!--
[![Travis-CI Build
Status](https://travis-ci.org/rauschenberger/palasso.svg)](https://travis-ci.org/rauschenberger/palasso)
-->

[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/rauschenberger/palasso?svg=true)](https://ci.appveyor.com/project/rauschenberger/palasso)
[![R-CMD-check](https://github.com/rauschenberger/palasso/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rauschenberger/palasso/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/rauschenberger/palasso/graph/badge.svg)](https://app.codecov.io/gh/rauschenberger/palasso)

# Sparse regression with paired covariates

The paired lasso is designed for situations where each covariate in one set forms a pair with a covariate in the other set (i.e., settings with a one-to-one correspondence between two covariate sets).

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
Marianne A. Jonker
[![MAJ](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-0134-8482),
Renée X. Menezes
[![RXM](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-1378-2721),
and Mark A. van de Wiel
[![MvdW](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-4780-8472)
(2020). "Sparse classification with paired covariates".
*Advances in Data Analysis and Classification* 14:571-588.
[doi: 10.1007/s11634-019-00375-6](https://doi.org/10.1007/s11634-019-00375-6).

[![CRAN
version](https://www.r-pkg.org/badges/version/palasso)](https://CRAN.R-project.org/package=palasso)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/palasso)](https://CRAN.R-project.org/package=palasso)
[![Total CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/palasso)](https://CRAN.R-project.org/package=palasso)

## Disclaimer

The R package `palasso` implements sparse regression with paired covariates ([Rauschenberger et al., 2020](https://doi.org/10.1007/s11634-019-00375-6)).

Copyright &copy; 2017 Armin Rauschenberger, Department of Epidemiology and Biostatistics, Amsterdam UMC, VU University Amsterdam, Amsterdam, The Netherlands

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
