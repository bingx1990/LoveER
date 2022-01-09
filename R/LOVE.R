########################################################################################
##########                                                                ##############
##########              Latent-model based OVErlap clustering             ##############
##########                                                                ##############
########################################################################################

#' @title LOVE: Latent-model based OVErlapping clustering
#'
#' @description Perform overlapping (variable) clustering of a \eqn{p}-dimensional feature
#' generated from the latent factor model \deqn{X = AZ + E} with identifiability
#' conditions on \eqn{A} and \eqn{Cov(Z)}.
#'
#' @param X A \eqn{n} by \eqn{p} data matrix.
#' @param delta The grid of leading constant of \eqn{\delta}.
#' @param lbd The grid of leading constant of \eqn{\lambda}.
#' @param mu The leading constant used for thresholding the loading matrix.
#' @param merge Logical. If TRUE, take the union of all candidate pure variables;
#'   otherwise, take the intersection.
#' @param diagonal Logical. If TRUE, the covariance matrix of \eqn{Z} is diagonal; else FALSE.
#' @param est_non_pure_row String. Procedure used for estimating the non-pure rows. One of \{"HT", "ST", "Dantzig"\}.
#' @param center_flag Logical. If TRUE, center the features.
#' @param equal_var Logical. TRUE if all features have equal variance.
#' @param rep_CV The number of repetitions used for cross validation.
#' @param verbose Logical. Set FALSE to suppress printing the progress.
#
#' @return A list of objects including: \itemize{
#'   \item \code{K} The estimated number of clusters.
#'   \item \code{pureVec} The estimated set of pure variables.
#'   \item \code{pureInd} The estimated partition of pure variables.
#'   \item \code{group} The estimated clusters (indices of each cluster).
#'   \item \code{A} The estimated \eqn{p} by \eqn{K} assignment matrix.
#'   \item \code{C} The covariance matrix of \eqn{Z}.
#'   \item \code{Omega} The precision matrix of \eqn{Z}.
#'   \item \code{Gamma} The diagonal of the covariance matrix of \eqn{E}.
#'   \item \code{optDelta} The selected value from \code{delta}.
#' }

#' @examples
#' data(toydata)
#' output_table <- LOVE(X = toydata)
#' @export



LOVE <- function(X, delta = seq(0.1, 1.1 ,0.1), lbd = 0.5, mu = 0.5, merge = FALSE,
                 diagonal = FALSE, est_non_pure_row = "HT", center_flag = TRUE,
                 equal_var = TRUE, rep_CV = 50, verbose = TRUE) {

  n <- nrow(X);  p <- ncol(X)

  if (center_flag)   # centering
    X <- scale(X, TRUE, FALSE)

  if (equal_var)
    se_est <- rep(1,p)
  else
    se_est <- apply(X, 2, sd)   # estimate the standard errors of each feature

  deltaGrids <- delta * sqrt(log(max(p, n)) / n)
  if (verbose)
    cat("Selecting optimal delta by using data splitting...\n")
  optDelta <- ifelse(length(deltaGrids) > 1,
                     median(replicate(rep_CV, CV_Delta(X, deltaGrids, diagonal, se_est, merge))), deltaGrids)

  if (verbose) {
    cat("Finishing selecting optimal delta =", optDelta, "with leading constant", min(which(deltaGrids >= optDelta)),"...\n")
    cat("Estimating pure rows...\n")
  }

  Sigma <- crossprod(X) / n
  resultAI <- EstAI(Sigma, optDelta, se_est, merge)

  ### Check if there is any group with ONLY ONE pure variable
  pure_numb <- sapply(resultAI$pureSignInd, FUN = function(x) {length(c(x$pos, x$neg))})
  if (sum(pure_numb == 1) > 0) {
    cat("Changing 'merge' to 'union' and reselecting delta ... \n")
    optDelta <- ifelse(length(deltaGrids) > 1,
                       median(replicate(rep_CV, CV_Delta(X, deltaGrids, diagonal, se_est, merge = F))),
                       deltaGrids)
    resultAI <- EstAI(Sigma, optDelta, se_est, merge = F)
  }

  A_hat <- resultAI$AI
  I_hat <- resultAI$pureVec

  if (is.null(I_hat)) {
    cat("Algorithm fails due to the non-existence of any pure variable.\n")
    stop()
  }

  if (verbose)
    cat("Estimating C and Sigma_IJ...\n")

  C_hat <- EstC(Sigma, A_hat, diagonal)

  # Estimate the covariance matrix of the error
  Gamma_hat <- rep(0, p)
  Gamma_hat[I_hat] <- diag(Sigma[I_hat, I_hat]) - diag(A_hat[I_hat,] %*% C_hat %*% t(A_hat[I_hat,]))
  Gamma_hat[Gamma_hat < 0] <- 0


  if (diagonal)
    Omega <- diag(diag(C_hat) ** (-1))
  else {
    if (verbose)
      cat("Selecting the tuning parameter for estimating Omega...\n")
    lbdGrids <- lbd * optDelta
    optLbd <- ifelse(length(lbd) > 1,
                     median(replicate(rep_CV, CV_lbd(X, lbdGrids, A_hat, I_hat, diagonal))),
                     lbdGrids)
    if (verbose) {
      cat("Selecting the optimal lambda =", optLbd, "with leading constant", min(which(lbdGrids >= optLbd)),"...\n")
      cat("Estimating Omega...\n")
    }
    Omega <- estOmega(optLbd, C_hat)
  }

  if (length(I_hat) == p) ### non-overlapping case
    group <- resultAI$pureSignInd
  else {
    Y <- EstY(Sigma, A_hat, I_hat)
    threshold <- mu * optDelta * norm(Omega, "I")

    if (verbose) {
      cat("Estimating non-pure rows by", switch(est_non_pure_row, "HT" = "Hard Thresholding", "ST" = "Soft Thresholding", "Dantzig" = "Dantzig"), "...\n")
    }
    if (est_non_pure_row == "HT")
      AJ <- threshA(t(Omega %*% Y), threshold)
    else if (est_non_pure_row == "ST")
      AJ <- EstAJInv(Omega, Y, threshold)
    else if (est_non_pure_row == "Dantzig") {
      AI <- abs(A_hat[I_hat, ])
      sigma_bar_sup <- max(solve(crossprod(AI), t(AI)) %*% se_est[I_hat])
      AJ <- EstAJDant(C_hat, Y, mu * optDelta * sigma_bar_sup, sigma_bar_sup + se_est[-I_hat])
    } else
      cat("Unknown method of estimating the non-pure rows.\n")

    A_hat[-I_hat, ] <- AJ
    group <-recoverGroup(A_hat)
  }
  return(list(K = ncol(A_hat),
              pureVec = I_hat,
              pureInd = resultAI$pureSignInd,
              group = group,
              A = A_hat,
              C = C_hat,
              Omega = Omega,
              Gamma = Gamma_hat,
              optDelta = delta[min(which(deltaGrids >= optDelta))]))
}


