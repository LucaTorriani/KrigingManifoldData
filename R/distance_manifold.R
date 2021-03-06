#' Distance on the manifold
#'
#' @param data1 Either a list/array [\code{p,p,B1}] of \code{B1} symmetric positive definite matrices of dimension \code{p*p}, or a single \code{p*p} matrix
#' @param data2 Either a list/array [\code{p,p,B2}] of \code{B2} symmetric positive definite matrices of dimension \code{p*p}, or a single \code{p*p} matrix.
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot", "Correlation"
#' @return A vector of distances, or a double if \code{data1} and \code{data2} are single matrices.
#' @description Compute the manifold distance between symmetric positive definite matrices
#' @details If \code{B2}=\code{B1} then the result is a vector of length \code{B1=B2} containing in position \code{i} the manifold distance beetween \code{data1[,,i]} and \code{data2[,,i]}.
#' Instead if \code{B2}=1 and \code{B1}!=1 the result is a vector of length \code{B1} containing in position \code{i} the manifold distance between \code{data1[,,i]} and \code{data2[,,1]}
#' @examples
#' data_manifold_model <- Manifoldgstat::rCov
#' distances <-distance_manifold(data_manifold_model, diag(2), metric_manifold = "Frobenius")
#' print(distances)
#' @useDynLib Manifoldgstat
#' @export

distance_manifold = function(data1, data2, metric_manifold = "Frobenius"){
  if(is.array(data1)){
    if(length(dim(data1))==3) data1 = alply(data1,3)
    if(length(dim(data1))==2) data1 = list(data1)
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
