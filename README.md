
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LoveER

<!-- badges: start -->
<!-- badges: end -->

LoveER performs prediction of *Y*, estimation and inference of *β* with
statistical guarantees under the structured latent factor regression
model *X* = *A* *Z* + *E* and *Y* = *Z*<sup>⊤</sup>*β* + *ε*.

## Installation

<!-- the released version of LoveER from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("LoveER") -->
<!-- ``` -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("bingx1990/LoveER")
```

It requires to pre-install the [LOVE](https://github.com/bingx1990/LOVE)
package with

``` r
devtools::install_github("bingx1990/LOVE")
```

## Example

This is a basic example which shows you how to use the ER package. We
start by generating a synthetic data set.

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
clustering of the columns of the **X** matrix. It should be noted that
`pure_homo = TRUE` is required for downstream prediction of *Y* and
inference of *β*.

``` r
library(LOVE)
library(LoveER)
# overlapping clustering
res_LOVE <- LOVE(X, pure_homo = TRUE, delta = seq(0.1, 1.1 ,0.1))
```

To predict *Y*, estimate the coefficient *β* and provide confidence
intervals of *β*, the following code provides an example.

``` r
# use the output of LOVE to perform prediction and inference.
res_ER <- ER(Y, X, res_LOVE, CI = TRUE)
res_ER <- ER(Y, X, res_LOVE, beta_est = "Dantzig", CI = FALSE)
```
