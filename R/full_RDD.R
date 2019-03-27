#' Perform full_RDD
#'
#' @param data_coords \code{N*2} or \code{N*3} matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to
#' be provided in signed decimal degrees
#' @param data_val array [\code{p,p,N}] of \code{N} symmetric positive definite matrices of dimension \code{p*p}
#' @param K number of cells the domain is subdivided in
#' @param grid prediction grid, i.e. \code{M*2} or \code{M*3} matrix of coordinates where to predict
#' @param nk_min minimum number of observations within a cell
#' @param B number of \emph{divide} iterations to perform
#' @param suppressMes \{\code{TRUE}, \code{FALSE}\} controls the level of interaction and warnings given
#' @param tol tolerance for the main loop of \code{model_kriging}
#' @param max_it maximum number of iterations for the main loop of \code{model_kriging}
#' @param n_h number of bins in the empirical variogram
#' @param tolerance_intrinsic tolerance for the computation of the intrinsic mean
#' @param X matrix (N rows and unrestricted number of columns) of additional covariates for the tangent space model, possibly NULL
#' @param X_new matrix (with the same number of rows of \code{new_coords}) of additional covariates for the new locations, possibly NULL
#' @param ker.width.intrinsic parameter controlling the width of the Gaussian kernel for the computation of the local mean (if 0, 
#' a "step kernel" is used, giving weight 1 to all the data within the cell and 0 to those outside of it)
#' @param ker.width.vario parameter controlling the width of the Gaussian kernel for the computation of the empirical variogram (if 0, 
#' a "step kernel" is used, giving weight 1 to all the data within the cell and 0 to those outside of it)
#' @param graph.distance.complete \code{N*N} distance matrix (the [i,j] element is the length of the shortest path between points i and j)
#' @param data.grid.distance \code{N*M} distance matrix between locations where the datum has been observed and locations where
#' the datum has to be predicted
#' @param aggregation_mean "Weighted" to aggregate the mean predictions using kernel-based weights, "Equal" to use equal weights
#' @param aggregation_kriging "Weighted" to aggregate the Kriging predictions using kernel-based weights, "Equal" to use equal weights
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
#'       \item \code{resBootstrap} {A list consisting of}
#'          \itemize{
#'           \item \code{fmean} {list of length \code{B}. Each field contains the prediction (at iteration \code{b}) for each new location, obtained 
#'                             as the intrinsic mean of the data within the tile it belongs to}
#'           \item \code{kervalues_mean} {Weights used for aggregating \code{fmean}}
#'          }
#'       \item \code{resAggregated} {Predictions, for each new location, obtained aggregating \code{fmean} using \code{kervalues_mean} as weights}
#'     }
#'   \item If \code{method.analysis} = "Kriging" it returns a list with the following fields
#'     \itemize{
#'       \item \code{resBootstrap} {A list consisting of}
#'          \itemize{
#'           \item \code{fmean} {list of length \code{B}. Each field contains the prediction (at iteration \code{b}) for each new location, obtained 
#'                             as the intrinsic mean of the data within the tile it belongs to}
#'           \item \code{fpred} {list of length \code{B}. Each field contains the prediction (at iteration \code{b}) for each new location, obtained 
#'                             through kriging}
#'           \item \code{kervalues_mean} {Weights used for aggregating \code{fmean}}
#'           \item \code{kervalues_krig} {Weights used for aggregating \code{fpred}}
#'           \item \code{variofit} {list of length \code{B}. Each field contains, for each datum, the parameters of the variogram fitted in the tile it belongs to}
#'          }
#'       \item \code{resAggregated} {Predictions, for each new location, obtained aggregating \code{fpred} using \code{kervalues_krig} as weights}
#'       \item \code{resLocalMean} {Predictions, for each new location, obtained aggregating \code{fmean} using \code{kervalues_mean} as weights}
#'     }
#' }
#' @details It uses a repetition of local analyses, through a \emph{divide} et \emph{impera} strategy. In the \emph{divide} step, the domain
#' is randomly decomposed in subdomains where local tangent-space models are estimated in order to predict at new locations (in each subregion is 
#' performed exactly the analysis described in the \code{model_kriging} function). This is repeated \code{B} times with different 
#' partitions of the domain. Then, in the \emph{impera} step, the results of these iterations are aggregated, by means of the intrinsic mean,
#' to provide a final prediction.
#' @description Perform kriging prediction using FullRDD procedure
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
                    # is.observed, border.length, N_samples, p,
                    aggregation_mean, aggregation_kriging,
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
                                  # is.observed=is.observed, border.length= border.length, p=p, 
                                  method.analysis = method.analysis,
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
                                  # is.observed=is.observed, border.length= border.length, p=p,
                                  method.analysis = method.analysis,
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
