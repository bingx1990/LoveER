######### This script contains some utilities function to test the performance ########

checkElement <- function(element, groupList) {
  # Check if an element is in a list. If it does exist in some group, return the group
  # index in that list and its sublist index. Otherwise, return c(0,0).
  for (i in 1:length(groupList)) {
    for (j in 1:length(groupList[[i]])) {
      if (element %in% groupList[[i]][[j]])
        return(c(i,j))
    }
  }
  return(c(0,0))
}


pureRowInd <- function(A) {
  # For given matrix A, find the row indices which correspond to pure nodes
  #
  # Args:
  #   A: p by K matrix.
  #
  # Returns:
  #   vector of pure row indices
  pureVec <- c()
  for (i in 1:ncol(A)) {
    pureVec <- c(pureVec, which(abs(A[ ,i]) == 1))
  }
  return(pureVec)
}



recoverGroup <- function(A) {
  # Recover group structure based given p by K matrix A and perform thresholding
  #
  # Args:
  #   A: estimated matrix.
  #   thresh: constant > 0.
  #
  # Returns:
  #   list of group indices with sign subpartition
  Group <- list()
  for (i in 1:ncol(A)) {
    column <- A[,i]
    posInd <- which(column > 0)
    negInd <- which(column < 0)
    Group[[i]] <- list(pos = posInd, neg = negInd)
  }
  return(Group)
}

singleton <- function(estPureIndices) {
  # Check if there exists an element of the given list has length equal to 1
  # If exists at least one, return TRUE; otherwise return FALSE
  if (length(estPureIndices) == 0)
    return(T)
  else
    ifelse(sum(sapply(estPureIndices, FUN = function(x) {length(x)}) == 1) > 0, T, F)
}


threshA <- function(A, mu, scale = FALSE) {
  # Threshold the estimated {@code A} based on the given {@code mu}. If {@code scale} is true,
  # then normalize each row of A such that the l-1 norm of each row is not larger than 1.
  scaledA <- A
  for (i in 1:nrow(A)) {
    colInd <- abs(A[i, ]) <= mu
    scaledA[i,colInd] = 0
    if (scale && sum(abs(scaledA[i, ])) > 1)
      scaledA[i, ] <- scaledA[i, ] / sum(abs(scaledA[i, ]))
  }
  return(scaledA)
}

offSum <- function(M, N, weights) {
  # Calculate the sum of squares of the upper off-diagonal elements of two matrices
  # require: M and N have the same dimensions
  tmp <- (M-N) / weights
  tmp <- t(t(tmp) / weights)
  return(sum((tmp[row(tmp) <= (col(tmp) - 1)])^2))
}




