
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LoveER

<!-- badges: start -->
<!-- badges: end -->

LoveER provides statistical guaranteed algorithms to

-   perform overlapping clustering of features under a structured latent
    factor model;
-   perform prediction, estimation and inference under a structured
    latent factor regression model.

## Installation

You can install the released version of LoveER from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("LoveER")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bingx1990/LoveER")
```

## Example

This is a basic example which shows you how to use LoveER. We start by
generating a data set.

The following code calls the LOVE function to perform overlapping
clustering of the columns of the matrix.

``` r
# library(LoveER)
## basic example code
# res_LOVE <- LOVE(X)
```

To predict , estimate the coefficient and provide confidence intervals
of , the following code provides an example.

``` r
# res_ER <- ER(Y, X, res_LOVE, CI = TRUE)
```
