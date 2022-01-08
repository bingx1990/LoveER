########################################################################################
##########                                                                ############## 
##########              Latent model based OVErlap clustering             ############## 
##########                                                                ##############
########################################################################################

###    Require modules "Tuning.R", "EstPure.R", "EstOmega.R", "EstNonpure.R"


LOVE <- function(X, delta = seq(0.1, 1.1 ,0.1), lbd = 0.5, mu = 0.5, merge = FALSE, 
                 diagonal = FALSE, est_non_pure_row = "HT", center_flag = TRUE, equal_var = TRUE, 
                 rep_CV = 50, verbose = TRUE) {
  # Non-overlapping/Overlappig variable clustering for p features collected in n samples
  #
  # Args: 
  #   X: n by p data matrix.
  #   delta: grid of leading constants for [delta]. 
  #   lbd: the grid of leading constants for [labmda]. Default value is 0.5.
  #   mu: the leading constant used in thresholding. Default value is 0.5.
  #   merge: if True, taking the union of all candidate pure variables; otherwise, taking the intersection.
  #   diagonal: TRUE if the cov(Z) is diagonal; else FALSE. 
  #   est_non_pure_row: procedure used for estimating the non-pure rows. One of {"HT", "ST", "Dantzig"}
  #   center_flag: center X if TRUE.
  #   eq_var: True if all features have equal variance.
  #   rep_CV: the number of repetitions used for cross validation.
  # 
  # Returns: 
  #   a list including: 
  #   K: the number of clusters.
  #   pureVec: the set of pure variables.
  #   pureInd: the partition of pure variables.
  #   group: the group structure(indices for each cluster).
  #   A: the p by K membership matrix.
  #   C: the covariance matrix of latent variables
  #   Omega: the precision of latent variables
  #   Gamma: the covariance matrix of the error

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
  return(list(K = ncol(A_hat), pureVec = I_hat, pureInd = resultAI$pureSignInd,
              group = group, A = A_hat, C = C_hat, Omega = Omega, Gamma = Gamma_hat,
              optDelta = delta[min(which(deltaGrids >= optDelta))]))
}


