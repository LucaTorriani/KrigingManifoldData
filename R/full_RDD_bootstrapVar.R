#' Bootstrap variance for method full_RDD
#' @useDynLib Manifoldgstat
#' @export
full_RDD_bootstrapVar <- function(res_RDD_OOK, K,n, metric_manifold, dim_matrix)
{
  # INPUT :
  # res_RDD_OOK = result of the function RDD_OOK
  # K = number of centers used
  # n = number of point evaluations of each function
  # t_step = step of discretization of the function
  # metric_manifold = metric to use on the manifold
  # dim_matrix = p
  # OUTPUT :
  # bootstrap variance at the predicted locations
  res.boot = res_RDD_OOK$resBootstrap$fpred
  res.aggr = res_RDD_OOK$resAggregated
  ngrid=dim(res.boot[[1]])[1]
  B=length(res.boot)
  boot.v <- rep(0,ngrid)

  if(K==1)
    return(boot.v)
  if(K>1)
  {
    for(i in 1:ngrid)
    {
      if(!is.na(res.aggr[i,1]))
      {
        imat = matrix_to_matrixArray(map(res.boot, return_ith_row, i)%>%ldply, dim_matrix)
        iaggr = vec_to_matrix(res.aggr[i,])
        boot.v[i] = 1.0/(B-1)*sum(Manifoldgstat::distance_manifold(imat, iaggr))
      }
    }
  }
  return(Re(boot.v))

}
