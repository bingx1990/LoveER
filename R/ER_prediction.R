########################################################################
####          Estimation of the coefficients for prediction       ######
########################################################################


#' @title Estimation for prediction under Essential Regression
#'
#' @description Estimate the p-dimensional coefficient between \eqn{Y} and \eqn{X}
#'   and the \eqn{p} by \eqn{K} matrix to predict \eqn{Z}.
#'
#' @inheritParams ER
#' @param Theta_hat A \eqn{p} by \eqn{K} matrix.
#'
#' @return A list including: \itemize{
#'   \item \code{theta} The estimated \eqn{p}-dimensional coefficients of \eqn{X}.
#'   \item \code{mat_trans_to_Z} The \eqn{p} by \eqn{K} matrix used to predict \eqn{Z}.
#'   \item \code{fitted_val} The fitted values.
#'   \item \code{Z_pred} The predicted \strong{Z} matrix.
#' }

ER_prediction <- function(Y, X, Theta_hat) {
  n <- nrow(X)
  Q <- try(Theta_hat %*% solve(crossprod(X %*% Theta_hat) / n, t(Theta_hat)), silent = T)
  if (class(Q)[1] == "try-error")
    Q <- Theta_hat %*% ginv(crossprod(X %*% Theta_hat) / n) %*% t(Theta_hat)

  theta_hat <- Q %*% crossprod(X, Y) / n
  fitted_val <- X %*% theta_hat
  mat_trans_to_Z <- Q %*% Theta_hat
  Z_hat <- X %*% mat_trans_to_Z
  return(list(theta = theta_hat,
              mat_trans_to_Z = mat_trans_to_Z,
              fitted_val = fitted_val,
              Z_pred = Z_hat))
}












################################################################################
###           Functions for re-fitting by using the support of beta          ###
################################################################################


# Pred_Z_BLP <- function(X, A_hat, C_hat, est_Gamma, S_beta) {
#   est_Gamma_inv <- diag(est_Gamma ** (-1))
#   G_hat <- crossprod(A_hat, est_Gamma_inv) %*% A_hat + solve(C_hat)
#   Z_hat <- X %*% est_Gamma_inv %*% A_hat %*% ginv(G_hat)
#   Z_hat[,S_beta,drop = F]
# }

# Pred_Z_BLP_avg <- function(X, A_hat, C_hat, S_beta) {
#   G_hat <- crossprod(A_hat) + solve(C_hat)
#   Z_hat <- X %*% A_hat %*% ginv(G_hat)
#   Z_hat[,S_beta,drop = F]
# }

# Pred_avg <- function(Y, X, A_hat, S, S_beta) {
#   B_hat <- t(solve(crossprod(A_hat[S,]), t(A_hat[S, ])))[, S_beta]
#   Z_hat <- X[ ,S] %*% B_hat
#   eta_hat <- try(solve(crossprod(Z_hat), crossprod(Z_hat, Y)), silent = T)
#   if (class(eta_hat) == "try-error")
#     eta_hat <- ginv(crossprod(Z_hat)) %*% crossprod(Z_hat, Y)
#   theta_hat <- B_hat %*% eta_hat
#   return(list(pred = X[ ,S] %*% theta_hat, theta = theta_hat))
# }
