#####     Functions to select tuning parameters via cross validation


### Functions to select delta


#' @title Cross validation to select \eqn{\delta}
#'
#' @description Cross validation for choosing the tuning parameter \eqn{\delta}.
#'   For each value of \code{deltaGrids}, first split the data into two parts and
#'   calculate \eqn{I}, \eqn{A_I} and \eqn{Cov(Z)}. Then calculate the fit \eqn{A_I Cov(Z) A_I'}
#'   to find the value which minimizes the loss criterion \deqn{||\Sigma - A_I Cov(Z) A_I'||_{F-off}/(|I|(|I|-1))}.
#'
#' @inheritParams LOVE
#' @param se_est The vector of the standard deviations of \eqn{p} features.
#' @param deltaGrids A vector of numerical constants.
#'
#' @return A numeric constant. The selected optimal \eqn{\delta}.

CV_Delta <- function(X, deltaGrids, diagonal, se_est, merge) {
  n <- nrow(X); p <- ncol(X)
  sampInd <- sample(n, floor(n / 2))
  X1 <- X[sampInd, ]
  X2 <- X[-sampInd, ]
  Sigma1 <- crossprod(X1) / nrow(X1);
  diag(Sigma1) <- 0
  Sigma2 <- crossprod(X2) / nrow(X2)

  result_Ms <- FindRowMax(abs(Sigma1))
  Ms <- result_Ms$M
  arg_Ms <- result_Ms$arg_M

  loss <- c()
  for (i in 1:length(deltaGrids)) {
    resultFitted <- CalFittedSigma(Sigma1, deltaGrids[i], Ms, arg_Ms, se_est,
                                   diagonal, merge)
    fittedValue <- resultFitted$fitted
    estPureVec <- resultFitted$pureVec
    if (is.null(dim(fittedValue)) && fittedValue == -1)
      loss[i] <- Inf
    else {
      denom <- length(estPureVec) * (length(estPureVec) - 1)
      # loss[i] <- 2 * offSum(Sigma2[estPureVec, estPureVec], fittedValue, 1) / denom
      loss[i] <- 2 * offSum(Sigma2[estPureVec, estPureVec] - fittedValue, se_est[estPureVec]) / denom
    }
  }
  # cat(loss)
  return(deltaGrids[which.min(loss)])
}




#' Calculate the fitted value \eqn{A_I Cov(Z) A_I'}.
#'
#' @inheritParams LOVE
#' @inheritParams EstAI
#' @inheritParams FindPureNode
#'
#' @return A list including: \itemize{
#'    \item \code{pureVec} A vector of the indices of the estimated pure variables.
#'    \item \code{fitted} The fitted value \eqn{A_I Cov(Z) A_I'}.
#' }
#' @noRd

CalFittedSigma <- function(Sigma, delta, Ms, arg_Ms, se_est, diagonal, merge) {
  resultPureNode <- FindPureNode(abs(Sigma), delta, Ms, arg_Ms, se_est, merge)

  estPureIndices <- resultPureNode$pureInd
  # lapply(estPureIndices, function(x) cat(x, "\n"))

  if (singleton(estPureIndices))
    return(list(pureVec = NULL, fitted = -1))

  estSignPureIndices <- FindSignPureNode(estPureIndices, Sigma)
  AI <- RecoverAI(estSignPureIndices, length(se_est))
  C <- EstC(Sigma, AI, diagonal)

  if (length(estPureIndices) == 1)
    fitted <- -1
  else {
    subAI <- AI[resultPureNode$pureVec, ]
    fitted <- subAI %*% C %*% t(subAI)
  }
  return(list(pureVec = resultPureNode$pureVec, fitted = fitted))
}



### Functions to select lambda for estimating the precision matrix


#' @title Cross validation to select \eqn{\lambda}
#'
#' @description Cross-validation to select \eqn{\lambda} for estimating the precision
#'   matrix of \eqn{Z}. Split the data into two parts. Estimating \eqn{Cov(Z)} on two datasets.
#'   Then, for each value in \code{lbdGrids}, calculate \eqn{Omega} on the first dataset
#'   and calculate the loss on the second dataset. Choose the value which minimizes
#'    \deqn{<Cov(Z), \Omega> - log(det(\Omega)).}
#'
#' @inheritParams LOVE
#' @param lbdGrids A vector of numerical constants.
#' @param AI A \eqn{p} by \eqn{K} matrix.
#' @param pureVec The estimated set of pure variables.
#'
#' @return The selected \eqn{\lambda}.

CV_lbd <- function(X, lbdGrids, AI, pureVec, diagonal) {
  sampInd <- sample(nrow(X), floor(nrow(X) / 2))
  X1 <- X[sampInd, ]
  X2 <- X[-sampInd, ]
  Sigma1 <- crossprod(X1) / nrow(X1)
  Sigma2 <- crossprod(X2) / nrow(X2)
  C1 <- EstC(Sigma1, AI, diagonal)
  C2 <- EstC(Sigma2, AI, diagonal)
  loss <- c()
  for (i in 1:length(lbdGrids)) {
    Omega <- estOmega(lbdGrids[i], C1)
    det_Omega <- det(Omega)
    loss[i] <- ifelse(det_Omega <= 0, Inf, sum(Omega * C2) - log(det_Omega))
  }
  return(lbdGrids[which.min(loss)])
}
