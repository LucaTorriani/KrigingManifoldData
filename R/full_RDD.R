#' Perform full_RDD
#' Create a GLS model
#'
#' @param data_coords \code{N*2} or \code{N*3} matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to
#' be provided in signed decimal degrees
#' @param data_val array [\code{p,p,N}] of \code{N} symmetric positive definite matrices of dimension \code{p*p}
#' @param K number of neighborhood (i.e., centers) to sample at each iteration
#' @param grid prediction grid
#' @param nk_min minimum number of observations within a neighborhood
#' @param B number of \texit{divide} iterations to perform
#' @param suppressMes {\code{TRUE}, \code{FALSE}} controls the level of interaction and warnings given
#' @param tol tolerance for the main loop of model_kriging
#' @param max_it maximum number of iterations for the main loop of model_kriging
#' @param n_h number of bins in the empirical variogram
#' @param tolerance_intrinsic tolerance for the computation of the intrinsic mean. Not needed if Sigma is provided
#' @param X matrix (N rows and unrestricted number of columns) of additional covariates for the tangent space model, possibly NULL
#' @param X_new matrix (with the same number of rows of \code{new_coords}) of additional covariates for the new locations, possibly NULL
#' @param ker.width.intrinsic parameter controlling the width of the Gaussian kernel for the computation of the local mean (if 0, no
#' kernel is used)
#' @param ker.width.vario parameter controlling the width of the Gaussian kernel for the computation of the empirical variogram (if 0,
#' no kernel is used)
#' @param graph.distance.complete \code{N*N} distance matrix (the [i,j] element is the length of the shortest path between points i and j)
#' @param data.grid.distance \code{N*dim(grid)[1]} distance matrix between locations where the datum has been observed and locations where
#' the datum has to be predicted
#' @param N_samples number of samples
#' @param p dimension of the manifold matrices
#' @param aggregation_mean "Weighted" ...
#' @param aggregation_kriging "Weighted" if the prediction must be aggregated using different weights, "Equal" to use equal weights
#' @param method.analysis "Local mean" to predict just with the mean, "Kriging" to predict via Kriging procedure
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot"
#' @param metric_ts metric used on the tangent space. It must be chosen among "Frobenius", "FrobeniusScaled", "Correlation"
#' @param model_ts type of model fitted on the tangent space. It must be chosen among "Intercept", "Coord1", "Coord2", "Additive"
#' @param vario_model type of variogram fitted. It must be chosen among "Gaussian", "Spherical", "Exponential"
#' @param distance type of distance between coordinates. It must be either "Eucldist" or "Geodist"
#' @return According to the analysis chosen:
#' \itemize{
#'   \item If \code{method.analysis} = "Local mean" it returns a list with the following fields
#'     \itemize{
#'       \item \code{resBootstrap}{...}
#'       \item \code{resAggregated}{...}
#'     }
#'   \item If \code{method.analysis} = "Kriging" it returns a list with the following fields
#'     \itemize{
#'       \item \code{resBootstrap}{...}
#'       \item \code{resAggregated}{...}
#'       \item \code{resLocalMean}{...}
#'     }
#' }
#' @useDynLib Manifoldgstat
#' @export
#'
full_RDD = function(data_coords, data_val, K, grid, nk_min=1, B=100,
                    # spdist='euclidean',
                    suppressMes=F,
                    tol=1e-12, max_it=100, n_h=15, tolerance_intrinsic =10^(-6), X=NULL, X_new=NULL, # plot=FALSE,
                    ker.width.intrinsic = 0, ker.width.vario = 1.5,
                    # mesh,
                    graph.distance.complete,
                    # assign.matrix, no.assg.grid,
                    data.grid.distance,
                    # is.observed, border.length,
                    N_samples, p,  aggregation_mean, aggregation_kriging,
                    method.analysis = 'Local mean',
                    metric_manifold, metric_ts, model_ts,
                    vario_model, distance = NULL)

{
  if(aggregation_mean != "Weighted" && aggregation_mean != "Equal") stop("aggregation_mean expected to be: Equal or Weighted")
  if(aggregation_kriging != "Weighted" && aggregation_kriging != "Equal") stop("aggregation_mean expected to be: Equal or Weighted")

  if(ker.width.intrinsic==0 & aggregation_mean=="Weighted") warning("Aggregation for the mean will use equal weights since ker.width.intrinsic=0")
  if(ker.width.vario==0 & aggregation_kriging=="Weighted") warning("Aggregation for the kriging predictors will use equal weights since ker.width.vario=0")

  # This function implements RDD-OOK procedure, by using the functions RDD_OOK_boot_man and RDD_OOK_aggr_man
  if(K==1)
  {
    resBootstrap=RDD_OOK_boot_man(data_coords=data_coords, data_val=data_val, K=K, grid=grid, nk_min=nk_min, B=1,
                                  # spdist=spdist,
                                  suppressMes=suppressMes, tol=tol, max_it = max_it,
                                  n_h=n_h, tolerance_intrinsic=tolerance_intrinsic, X=X, X_new=X_new,  # plot=plot,
                                  ker.width.intrinsic = 0,
                                  ker.width.vario = 0,
                                  # mesh=mesh,
                                  graph.distance.complete= graph.distance.complete,
                                  # assign.matrix=assign.matrix, no.assg.grid=no.assg.grid,
                                  data.grid.distance = data.grid.distance,
                                  # is.observed=is.observed, border.length= border.length,
                                  p=p, method.analysis = method.analysis,
                                  metric_manifold = metric_manifold, metric_ts = metric_ts, model_ts = model_ts,
                                  vario_model = vario_model, distance = distance)


    # resAggregated=simplify2array(resBootstrap$fpred)[,,1]
    resAggregated=resBootstrap$fpred[[1]]
    if(method.analysis == 'Kriging'){
      # resLocalMean=simplify2array(resBootstrap$fmean)[,,1]
      resLocalMean=resBootstrap$fmean[[1]]
      return(list(resBootstrap = resBootstrap,
                  resAggregated = resAggregated,
                  resLocalMean = resLocalMean))
    }
  }
  if(K>1)
  {
    resBootstrap=RDD_OOK_boot_man(data_coords=data_coords, data_val=data_val, K=K, grid=grid, nk_min=nk_min, B=B,
                                  # spdist=spdist,
                                  suppressMes=suppressMes, tol=tol, max_it = max_it,
                                  n_h=n_h, tolerance_intrinsic =tolerance_intrinsic, X=X, X_new=X_new,  # plot=plot,
                                  ker.width.intrinsic = ker.width.intrinsic,
                                  ker.width.vario = ker.width.vario,
                                  # mesh=mesh,
                                  graph.distance.complete= graph.distance.complete,
                                  # assign.matrix=assign.matrix, no.assg.grid=no.assg.grid,
                                  data.grid.distance = data.grid.distance,
                                  # is.observed=is.observed, border.length= border.length,
                                  p=p,  method.analysis = method.analysis,
                                  metric_manifold = metric_manifold, metric_ts = metric_ts, model_ts = model_ts,
                                  vario_model = vario_model, distance = distance)

    if(method.analysis == 'Local mean') {
      if (aggregation_mean== "Equal") ker.width.intrinsic=0
      resAggregated=RDD_OOK_aggr_man(fOKBV = resBootstrap$fmean, weights_intrinsic = resBootstrap$kervalues_mean,
                                     ker.width.intrinsic=  ker.width.intrinsic)
    }


    if(method.analysis == 'Kriging')
    {
      if (aggregation_kriging== "Equal") ker.width.vario = 0
      resAggregated=RDD_OOK_aggr_man(fOKBV = resBootstrap$fpred, weights_intrinsic = resBootstrap$kervalues_krig,
                                     ker.width.intrinsic=  ker.width.vario)

      if (aggregation_mean== "Equal") ker.width.intrinsic = 0
      resLocalMean=RDD_OOK_aggr_man(fOKBV = resBootstrap$fmean, weights_intrinsic = resBootstrap$kervalues_mean,
                                    ker.width.intrinsic =  ker.width.intrinsic)


      return(list(resBootstrap = resBootstrap,
                  resAggregated = resAggregated,
                  resLocalMean = resLocalMean))
    }
  }
  return(list(resBootstrap = resBootstrap,
              resAggregated = resAggregated  ))
}
