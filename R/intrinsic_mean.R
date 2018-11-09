#' Intrinsic mean
#'
#' @param data list or array [\code{n,n,B}] of \code{B} symmetric positive definite matrices of dimension \code{nxn}
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot", "Correlation"
#' @param metric_ts metric used on the tangent space. It must be chosen among "Frobenius", "FrobeniusScaled", "Correlation"
#' @param tolerance tolerance for the computation of the intrinsic_mean
#' @param weight_intrinsic vector of length \code{B} to weight the matrices in the computation of the intrinsic mean. If NULL
#' a vector of ones is used
#' @param weight_extrinsic vector of length \code{B} to weight the matrices in the computation of the extrinsic mean. If NULL
#' weight_intrinsic are used
#' @param tolerance_map_cor tolerance to use in the maps. Required only if \code{metric_manifold== "Correlation"}
#' @return A matrix representing the intrinsic mean of the \code{data}
#' @description Evaluate the intrinsic mean of a given set of symmetric positive definite matrices
#' @examples
#' data_manifold_tot <- Manifoldgstat::fieldCov
#' Sigma <-intrinsic_mean(data_manifold_tot, metric_manifold = "Frobenius", metric_ts = "Frobenius")
#' print(Sigma)
#' @useDynLib Manifoldgstat
#' @export

intrinsic_mean = function(data, metric_manifold = "Frobenius", metric_ts = "Frobenius",
                     tolerance = 1e-6, weight_intrinsic= NULL, weight_extrinsic= weight_intrinsic, tolerance_map_cor=1e-6){
  if( is.array(data)){
    data = alply(data,3)
  }

  if(is.null(weight_intrinsic)) {
    weight_intrinsic = rep(1, length(data))
    weight_extrinsic = weight_intrinsic
  }

  result =.Call("intrinsic_mean",data, N=length(data), metric_manifold, metric_ts, tolerance, weight_intrinsic, weight_extrinsic, tolerance_map_cor)

  return (result)
}
