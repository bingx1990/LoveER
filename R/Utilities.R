#####      This script contains some helper functions

#' Recover clusters based on a given \eqn{p} by \eqn{K} loading matrix.
#'
#' @param A A \eqn{p} by \eqn{K} matrix.
#'
#' @return A list of group indices with sign sub-partition.

recoverGroup <- function(A) {
  Group <- list()
  for (i in 1:ncol(A)) {
    column <- A[,i]
    posInd <- which(column > 0)
    negInd <- which(column < 0)
    Group[[i]] <- list(pos = posInd, neg = negInd)
  }
  return(Group)
}


#' Function to check if there exists any element in the given list that has length
#' equal to 1
#'
#' @param estPureIndices A list of indices of the estimated pure variables.
#'
#' @return Logical. If exists at least one, return TRUE; otherwise return FALSE

singleton <- function(estPureIndices) {
  if (length(estPureIndices) == 0)
    return(T)
  else
    ifelse(sum(sapply(estPureIndices, FUN = function(x) {length(x)}) == 1) > 0, T, F)
}


#' @title Function to hard-threshold
#'
#' @description Threshold the estimated \code{A} based on the given \code{mu}.
#'   If \code{scale} is TRUE, then normalize each row of \code{A} such that the
#'   \eqn{\ell-1} norm of each row is no larger than 1.
#'
#' @inheritParams recoverGroup
#' @param mu A numeric value.
#' @param scale Logical. Normalize the row-wise \eqn{\ell-1} norm if TRUE.
#'
#' @return A matrix with the same dimension as \code{A}.
#' @noRd

threshA <- function(A, mu, scale = FALSE) {
  scaledA <- A
  for (i in 1:nrow(A)) {
    colInd <- abs(A[i, ]) <= mu
    scaledA[i,colInd] = 0
    if (scale && sum(abs(scaledA[i, ])) > 1)
      scaledA[i, ] <- scaledA[i, ] / sum(abs(scaledA[i, ]))
  }
  return(scaledA)
}



#' Calculate the weighted sum of squares of the upper off-diagonal elements
#'
#' @param M A given symmetric matrix
#' @param weights A vector of length equal to the number of rows of \code{M}.
#'
#' @return A numeric value.

offSum <- function(M, weights) {
  tmp <- M / weights
  tmp <- t(t(tmp) / weights)
  return(sum((tmp[row(tmp) <= (col(tmp) - 1)])^2))
}




# checkElement <- function(element, groupList) {
#   # Check if an element is in a list. If it does exist in some group, return the group
#   # index in that list and its sublist index. Otherwise, return c(0,0).
#   for (i in 1:length(groupList)) {
#     for (j in 1:length(groupList[[i]])) {
#       if (element %in% groupList[[i]][[j]])
#         return(c(i,j))
#     }
#   }
#   return(c(0,0))
# }


# pureRowInd <- function(A) {
#   # For given matrix A, find the row indices which correspond to pure nodes
#   #
#   # Args:
#   #   A: p by K matrix.
#   #
#   # Returns:
#   #   vector of pure row indices
#   pureVec <- c()
#   for (i in 1:ncol(A)) {
#     pureVec <- c(pureVec, which(abs(A[ ,i]) == 1))
#   }
#   return(pureVec)
# }



