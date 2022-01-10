#################################################################################
#####                                                                       #####
#####      Functions related to the estimation of non-pure rows of A        #####
#####                                                                       #####
#################################################################################


#' Estimate the covariance matrix of \eqn{Z}.
#'
#' @inheritParams EstAI
#' @param AI A \eqn{p} by \eqn{K} matrix.
#' @inheritParams LOVE
#'
#' @return A \eqn{K} by \eqn{K} matrix.

EstC <- function(Sigma, AI, diagonal) {
  K <- ncol(AI)
  C <- diag(0, K, K)
  for (i in 1:K) {
    groupi <- which(AI[ ,i] != 0)
    sigmai <- as.matrix(abs(Sigma[groupi,groupi]))
    tmpEntry <- sum(sigmai) - sum(diag(sigmai))
    C[i,i] <- tmpEntry / (length(groupi) * (length(groupi) - 1))
    if (!diagonal && i < K) {
      for (j in (i+1):K) {
        groupj <- which(AI[ ,j]!=0)
        # adjust the sign for each row
        sigmaij <- AI[groupi,i] * as.matrix(Sigma[groupi, groupj])
        sigmaij <- t(AI[groupj, j] * t(sigmaij))
        C[i,j] <- C[j,i] <- sum(sigmaij) / (length(groupi) * length(groupj))
      }
    }
  }
  return(C)
}



#' Function to estimate the \eqn{K} by \eqn{|J|} submatrix of \eqn{\Sigma}.
#'
#' @inheritParams EstC
#' @param pureVec The estimated set of pure variables.
#'
#' @return A \eqn{K} by \eqn{|J|} matrix.

EstY <- function(Sigma, AI, pureVec) {
   AI_sub <- AI[pureVec,,drop = F]
   return(solve(crossprod(AI_sub), t(AI_sub) %*% Sigma[pureVec, -pureVec, drop = F]))
}



#' Estimates non-pure rows via soft-thresholding
#'
#' This function estimates the \eqn{|J|} by \eqn{K} submatrix \eqn{A_J} by using soft thresholding.
#'
#' @param Omega The estimated precision matrix of \eqn{Z}.
#' @param Y A \eqn{K} by \eqn{|J|} reponse matrix.
#' @param lbd A tuning parameter for soft-thresholding.
#'
#' @return A \eqn{|J|} by \eqn{K} matrix.

EstAJInv <- function(Omega, Y, lbd) {
  AJ <- matrix(0, ncol(Y), nrow(Y))
  for (i in 1:ncol(Y)) {
    Atilde <- Omega %*% as.matrix(Y[ ,i])
    AJ[i, ] <- LP(Atilde, lbd)
    if (sum(abs(AJ[i, ])) > 1)
      AJ[i,] <- AJ[i,] / sum(abs(AJ[i, ]))
  }
  return(AJ)
}



#' Cast the procedure as a linear program
#   beta^+ - beta^- \leq lbd + y
#   - beta^+ + beta^- \leq lbd - y
#   beta^+ \ge 0; beta^- \ge 0
#
#' @param y A vector of length \eqn{K}.
#' @inheritParams EstAJInv
#'
#' @return A vector of length \eqn{K}.

LP <- function(y, lbd) {
  K <- length(y)
  cvec <- rep(1, 2 * K)
  bvec <- c(lbd + y, lbd - y, rep(0, 2 * K))
  C <- matrix(0, K, 2 * K)
  for (i in 1:K) {
    indices <- c(i,i + K)
    C[i,indices] = c(1,-1)
  }
  Amat <- rbind(C, -C, diag(-1, nrow = 2 * K))
  LPsol <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  beta <- LPsol[1:K] - LPsol[(K + 1):(2 * K)]
  return(beta)
}


#' Estimate the non-pure rows via the Dantzig approach
#'
#' @inheritParams estOmega
#' @inheritParams EstAJInv
#' @param se_est_J The estimated standard errors of the non-pure variables.
#'
#' @return A \eqn{|J|} by \eqn{K} matrix.

EstAJDant <- function(C_hat, Y, lbd, se_est_J) {
  AJ <- matrix(0, ncol(Y), nrow(Y))
  for (i in 1:ncol(Y)) {
    AJ[i, ] <- Dantzig(C_hat, Y[,i], lbd * se_est_J[i])
    if (sum(abs(AJ[i, ])) > 1)
      AJ[i,] <- AJ[i,] / sum(abs(AJ[i, ]))
    # cat("Finishing estimating the", i, "th row...\n")
  }
  return(AJ)
}


#' The Dantzig approach of solving one non-pure row
#'
#' @inheritParams LP
#' @inheritParams EstAJDant
#'
#' @return A vector of length \eqn{K}.

Dantzig <- function(C_hat, y, lbd) {
  K <- length(y)
  cvec <- rep(1, 2 * K)
  bvec <- c(lbd + y, lbd - y, rep(0, 2 * K))
  new_C_hat <- matrix(0, K, 2 * K)
  for (i in 1:K) {
    new_C_hat[i, ] <- c(C_hat[i,], -C_hat[i,])
  }
  Amat <- rbind(new_C_hat, -new_C_hat, diag(-1, nrow = 2 * K))
  LPsol <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  beta <- LPsol[1:K] - LPsol[(K + 1):(2 * K)]
  return(beta)
}







# EstY <- function(Sigma, AI, pureVec) {
#   SigmaS <- AdjustSign(Sigma, AI)
#   SigmaJ <- matrix(SigmaS[ ,-pureVec], nrow = nrow(Sigma))
#   SigmaTJ <- matrix(0, ncol(AI), nrow(AI) - length(pureVec))
#   for (i in 1:ncol(AI)) {
#     groupi <- which(AI[ ,i] != 0)
#     SigmaiJ <- as.matrix(SigmaJ[groupi, ])
#     SigmaTJ[i, ] <- apply(SigmaiJ, 2, mean) # Average columns along the rows.
#   }
#   return(SigmaTJ)
# }


# AdjustSign <- function(Sigma, AI) {
#   SigmaS <- matrix(0, nrow(AI), nrow(AI))
#   for (i in 1:nrow(AI)) {
#     index <- which(AI[i, ] != 0)
#     if (length(index) != 0)
#       SigmaS[i, ] = sign(AI[i,index]) * Sigma[i, ]
#   }
#   return(SigmaS)
# }



