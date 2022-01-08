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
#' @description Perform prediction, estimation and inference under Essential Regressions
#'
#' @param Y a vector of response with length n.
#' @param X n by p data matrix.
#' @param res_LOVE the fitted model from \code{\link{LOVE}}.
#' @param beta_est the procedure of estimating beta. One of \{"NULL", "LS", "Dantzig"\}
#' @param mu the tuning parameter used for estimating beta via the Dantzig approach. Default value is 0.5.
#' @param lbd the tuning parameter used for estimating beta via the Dantzig approach. Default value is 0.5.
#' @param CI construct confidence intervals if TRUE, otherwise FALSE.
#' @param alpha_level significance level. Default set to 0.05.
#' @param correction correction for multiple hypothesis testing. Either NULL or "Bonferroni". Default set to "Bonferroni".
#
#' @return a list of objects including: \itemize{
#'   \item \code{beta} the estimated coefficients beta
#'   \item \code{beta_CIs} the coordinate-wise confidence intervals of beta
#'   \item \code{beta_var} the variances of the estimated beta
#'   \item \code{coef_X} the coefficients between Y and X
#'   \item \code{mat_trans_to_Z} the p by K matrix used to predict Z
#'   \item \code{fitted_val} the fitted values
#'   \item \code{Z_pred} the predicted Z matrix
#'   \item \code{X_center} the centers of the input X
#'   \item \code{Y_center} the sample mean of the input Y
#' }


#' @examples
#' data(toydata)
#' output_table <- ER(X = toydata)
#' @export


ER <- function(Y, X, res_LOVE, beta_est = "NULL", mu = 0.5, lbd = 0.5,
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


  if (beta_est == "NULL" || beta_est == "Dantzig") {
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





