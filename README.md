
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LoveER

<!-- badges: start -->
<!-- badges: end -->

LoveER provides algorithms with statistical guarantees to

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

This is a basic example which shows you how to use two main functions of
LoveER: LOVE and ER. We start by generating a synthetic data set.

``` r
p <- 6
n <- 50
K <- 2
A <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, 1), c(1/3, 2/3), c(1/2, -1/2))
Z <- matrix(rnorm(n * K, sd = 2), n, K)
E <- matrix(rnorm(n * p), n, p)
X <- Z %*% t(A) + E

eps <- rnorm(n)
beta <- c(1, -0.5)
Y <- Z %*% beta + eps
```

The following code calls the LOVE function to perform overlapping
clustering of the columns of the matrix.

``` r
library(LoveER)
# basic example code
res_LOVE <- LOVE(X)
#> Selecting optimal delta by using data splitting...
#> Finishing selecting optimal delta = 0.0279715 with leading constant 1 ...
#> Estimating pure rows...
#> Estimating C and Sigma_IJ...
#> Selecting the tuning parameter for estimating Omega...
#> Selecting the optimal lambda = 0.01398575 with leading constant 1 ...
#> Estimating Omega...
#> Estimating non-pure rows by Hard Thresholding ...
```

To predict , estimate the coefficient and provide confidence intervals
of , the following code provides an example.

``` r
res_ER <- ER(Y, X, res_LOVE, CI = TRUE)
```
