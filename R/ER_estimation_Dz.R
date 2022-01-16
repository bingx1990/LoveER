################################################################################
###                 Estimate beta by Dantzig-type estimator                  ###
################################################################################
# library(linprog)


#' Estimate \eqn{\beta} via the Dantzig approach.
#'
#' @inheritParams ER
#' @param A_hat The estimated loading matrix.
#' @param C_hat The estimated covariance matrix of \eqn{Z}.
#' @param I_hat The estimated set of pure variables.
#' @param delta_opt The selected constant of \code{delta} returned by \code{\link{LOVE}}.
#'
#' @return A numeric vector of length \eqn{K}.
#' @noRd

ER_est_beta_dz <- function(Y, X, A_hat, C_hat, I_hat, delta_opt, mu = 0.5, lbd = 0.5) {
  n <- nrow(X); p <- ncol(X)
  AI <- A_hat[I_hat,]
  h <- solve(crossprod(AI), t(AI) %*% crossprod(X[ ,I_hat], Y) / n)
  as.vector(Solve_dz(C_hat, h, mu * delta_opt, lbd * delta_opt))
}


#' Solve the linear program to estimate \eqn{\beta}.
#'
#' @param C A \eqn{K} by \eqn{K} matrix.
#' @param h A vector of length \eqn{K}.
#' @param mu,lbd Two numeric constants.
#' @noRd

Solve_dz <- function(C, h, mu, lbd) {
  K <- nrow(C)
  cvec <- c(1, rep(0, 2*K))
  Amat <- -cvec
  Amat <- rbind(Amat, c(-1, rep(1, 2*K)))
  tmp_constr <- C %x% t(c(1,-1))
  Amat <- rbind(Amat, cbind(-1 * lbd, rbind(tmp_constr, -tmp_constr)))
  bvec <- c(0, 0, mu + h, mu - h)

  lpResult <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  while (length(lpResult) == 0) {
    cat("The penalty lambda =", lbd, "is too small and increased by 0.01...\n")
    lbd <- lbd + 0.01
    Amat[-(1:2), 1] <- lbd
    lpResult <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution[-1]
  }
  ind <- seq(2, 2*K, 2)
  return(lpResult[ind] - lpResult[ind + 1])
}
