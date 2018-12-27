#' Distance on the manifold
#'
#' @param data1 list or array [\code{p,p,B}] of \code{B} symmetric positive definite matrices of dimension \code{p*p}
#' @param data2 a list or array [\code{p,p,B}] of \code{B} symmetric positive definite matrices of dimension \code{p*p}.
#' Or a single \code{nxn} matrix.
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot"
#' @description Compute the manifold distance between symmetric positive definite matrices
#' @details If \code{B2}=\code{B1} then the result is a vector of length \code{B1=B2} containing in position \code{i} the manifold distance beetween \code{data1[,,i]} and \code{data2[,,i]}.
#' Instead if \code{B2}=1 the result is a vector of length \code{B1} containing in position \code{i} the manifold distance between \code{data1[,,i]} and \code{data2[,,1]}
#' @examples
#' data_manifold_model <- Manifoldgstat::rCov
#' distances <-distance_manifold(data_manifold_model, diag(2), metric_manifold = "Frobenius")
#' print(distances)
#' @useDynLib Manifoldgstat
#' @export

distance_manifold = function(data1, data2, metric_manifold = "Frobenius", metric_ts = "Frobenius"){
  if(is.array(data1)){
    data1 = alply(data1,3)
  }
  if(is.array(data2)){
    if(length(dim(data2))==3) data2 = alply(data2,3)
    if(length(dim(data2))==2) data2 = list(data2)
  }

  N1=length(data1)
  N2=length(data2)
  if(N1!=N2 & N2!=1) stop("Data must be either arrays of the same dimension or an array and a single matrix")

  result =.Call("distance_manifold",data1,data2, N1,N2, metric_manifold)

  return (result)
}
