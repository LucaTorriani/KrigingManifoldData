#' Bootstrap variance 
#' @useDynLib Manifoldgstat
#' @export

# res.boot = res_RDD_OOK$resBootstrap$fmean 
# res.boot = res_RDD_OOK$resBootstrap$fpred
# res.aggr = res_RDD_OOK$resAggregated
bootstrapVar <- function(res.boot, res.aggregated, K, metric_manifold)
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
  
  # res.boot = res_RDD_OOK$resBootstrap$fpred
  
  # dimgrid=dim(res.boot[[1]])[1]
  dimgrid = length(res.boot[[1]])
  B=length(res.boot)
  boot.v <- rep(0,dimgrid)
  
  if(K==1)
    return(boot.v)
  if(K>1)
  {
    for(i in 1:dimgrid)
    {
      if(!is.na(res.aggr[[i]][1,1]))
      {
        # imat = matrix_to_matrixArray(map(res.boot, return_ith_row, i)%>%ldply, dim_matrix)
        # iaggr = vec_to_matrix(res.aggr[i,])
        # boot.v[i] = 1.0/(B-1)*sum(Manifoldgstat::distance_manifold(imat, iaggr))
        boot.v[i] = 1.0/(B-1)*sum(Manifoldgstat::distance_manifold(map(res.boot,return_ith_list_element, i=i), res.aggr[[i]]))
      }
    }
  }
  # return(Re(boot.v))
  return(boot.v)
}