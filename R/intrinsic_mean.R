#' Intrinsic mean
#'
#' @param data list or array [\code{n,n,B}] of \code{B} symmetric positive definite matrices of dimension \code{nxn}
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot"
#' @param metric_ts metric used on the tangent space. It must be either "Frobenius" or "FrobeniusScaled"
#' @param tolerance tolerance for the computation of the intrinsic_mean
#' @param weight vector of length \code{N} to weight the locations in the computation of the intrinsic mean. If NULL
#' a vector of ones is used
#' @description Evaluate the intrinsic mean of a given set of symmetric positive definite matrices
#' @details ...
#' @examples
#' data_manifold_tot <- Manifoldgstat::fieldCov
#' Sigma <-intrinsic_mean(data_manifold_tot, metric_manifold = "Frobenius", metric_ts = "Frobenius")
#' print(Sigma)
#' @useDynLib Manifoldgstat
#' @export

intrinsic_mean = function(data, metric_manifold = "Frobenius", metric_ts = "Frobenius",
                     tolerance = 1e-6, weight= NULL){
  if( is.array(data)){
    data = alply(data,3)
  }

  if(is.null(weight)) weight = rep(1, length(data))

  result =.Call("intrinsic_mean",data, N=length(data), metric_manifold, metric_ts, tolerance, weight)

  return (result)
}
