########################################################################
########                  Essential Regression                ##########
########################################################################
# source("CV.R")
# source("Utilities.R")
# source("EstNonpure.R")
# source("EstPure.R")
# source("EstOmega.R")
# source("ER_estimation_Dz.R")
# source("ER_estimation_LS.R")
# source("LOVE.R")
# source("ER_prediction.R")
# library(linprog)


#' @title Essential Regression
#'
#' @description Under the Essential
#'  Regression framework \deqn{X = AZ+E,    Y = Z' \beta + \epsilon,}
#'  perform prediction of \eqn{Y}, estimation and inference of \eqn{\beta}.
#'
#' @param Y A vector of response with length \eqn{n}.
#' @param X A \eqn{n} by \eqn{p} data matrix.
#' @param res_LOVE The returned object from \code{\link{LOVE}}.
#' @param beta_est The procedure used for estimating \eqn{\beta}. One of \{\code{NULL}, "LS", "Dantzig"\}
#' @param mu,lbd The tuning parameters used for estimating \eqn{\beta} via the Dantzig approach. The default value is 0.5.
#' @param CI Logical. TRUE if confidence intervals are constructed.
#' @param alpha_level The significance level. The default set to 0.05.
#' @param correction The approach used for addressing the multiple testing problem.
#'  Either \code{NULL} or "Bonferroni". The default value is "Bonferroni".
#
#' @return A list of objects including: \itemize{
#'   \item \code{beta} The estimated coefficients of \eqn{\beta}.
#'   \item \code{beta_CIs} The coordinate-wise confidence intervals of \eqn{\beta}.
#'   \item \code{beta_var} The variances of the \code{beta}.
#'   \item \code{coef_X} The estimated p-dimensional coefficient between \eqn{Y} and \eqn{X}.
#'   \item \code{mat_trans_to_Z} The \eqn{p} by \eqn{K} matrix used to predict \eqn{Z}.
#'   \item \code{fitted_val} The fitted values of length \eqn{n}.
#'   \item \code{Z_pred} The predicted \strong{Z} matrix.
#'   \item \code{X_center} Centers of the input \code{X}.
#'   \item \code{Y_center} Center of the input \code{Y}.
#' }


#' @examples
#' p <- 6
#' n <- 50
#' K <- 2
#' A <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, 1), c(1/3, 2/3), c(1/2, -1/2))
#' Z <- matrix(rnorm(n * K, sd = 2), n, K)
#' E <- matrix(rnorm(n * p), n, p)
#' X <- Z %*% t(A) + E
#'
#' eps <- rnorm(n)
#' beta <- c(1, -0.5)
#' Y <- Z %*% beta + eps
#'
#' res_LOVE <- LOVE(X)
#' res_ER <- ER(Y, X, res_LOVE, CI = TRUE)
#'
#' @export


ER <- function(Y, X, res_LOVE, beta_est = "LS", mu = 0.5, lbd = 0.5,
               CI = F, alpha_level = 0.05, correction = "Bonferroni") {

  n <- nrow(X); p <- ncol(X)

  # center both Y and X
  Y_center <-  mean(Y)
  Y <- Y - Y_center
  X <- scale(X, T, F)

  A_hat <- res_LOVE$A;  C_hat <- res_LOVE$C; I_hat <- res_LOVE$pureVec;
  I_hat_partition <- res_LOVE$pureInd; Gamma_hat <- res_LOVE$Gamma;
  optDelta <- res_LOVE$optDelta

  Sigma <- crossprod(X) / n

  Theta_hat <- matrix(0, nrow = p, ncol = ncol(A_hat))
  BI <- t(solve(crossprod(A_hat[I_hat,]), t(A_hat[I_hat,])))
  Theta_hat[I_hat, ] = (Sigma[I_hat, I_hat] - diag(Gamma_hat[I_hat])) %*% BI
  Theta_hat[-I_hat, ] =  Sigma[-I_hat, I_hat] %*% BI


  pred_result <- ER_prediction(Y, X, Theta_hat)


  if (is.null(beta_est) || beta_est == "Dantzig") {
    beta_hat <- beta_CIs <- beta_var <- NULL
    if (beta_est == "Dantzig") {
      beta_hat <- ER_est_beta_dz(Y, X, A_hat, C_hat, I_hat, optDelta, mu, lbd)
    }
  } else {
    # least squares estimation
    res_beta <- ER_est_beta_LS(Y, X, Theta_hat, C_hat, Gamma_hat, I_hat, I_hat_partition, BI,
                               CI = CI, alpha_level = alpha_level, correction = correction)
    beta_hat <- res_beta$beta
    beta_CIs <- res_beta$CIs
    beta_var <- res_beta$beta_var
  }
  return(list(beta = beta_hat, beta_CIs = beta_CIs, beta_var = beta_var,
              coef_X = pred_result$theta,
              mat_trans_to_Z = pred_result$mat_trans_to_Z,
              fitted_val = pred_result$fitted_val,
              Z_pred = pred_result$Z_pred,
              X_center = attr(X, "scaled:center"), Y_center = Y_center))
}





#' Prediction under Essential Regression
#'
#' Predict either the response or the latent factors under Essential Regression.
#'
#' @param fitted_ER The object returned from \code{\link{ER}}.
#' @param newX A new data matrix of \eqn{X}. If not provided, the fitted values of
#'   either \strong{Y} or \strong{Z} are returned, depending on \code{type}.
#' @param type Either "response" or "factor". For "response", predicted values of
#'   \eqn{Y} are returned. For "factor", predicted \eqn{Z} is returned.
#'
#' @export

Predict.ER <- function(fitted_ER, newX = NULL, type = "response") {
  # newdata is a n by p matrix / array
  # type is one of {"response", "factor"}

  if (is.null(newX)) {
    if (type == "response")
      return(fitted_ER$fitted_val + fitted_ER$Y_center)
    else
      return(fitted_ER$Z_pred)
  } else {
    newX <- t(t(newX) - fitted_ER$X_center)
    if (type == "response")
      return(newX %*% fitted_ER$coef_X + fitted_ER$Y_center)
    else if (type == "factor")
      return(newX %*% fitted_ER$mat_trans_to_Z)
    else
      stop("Unknown prediction type.\n")
  }
}





