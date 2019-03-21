#' Perform mixed_RDD
#' @param data_coords \code{N*2} or \code{N*3} matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to
#' be provided in signed decimal degrees
#' @param data_val array [\code{p,p,N}] of \code{N} symmetric positive definite matrices of dimension \code{p*p}
#' @param K number of neighborhood (i.e., centers) to sample at each iteration
#' @param grid prediction grid
#' @param nk_min minimum number of observations within a neighborhood
#' @param B number of \texit{divide} iterations to perform
#' @param suppressMes \{\code{TRUE}, \code{FALSE}\} controls the level of interaction and warnings given
#' @param ker.width.intrinsic parameter controlling the width of the Gaussian kernel for the computation of the local mean (if 0, no
#' kernel is used)
#' @param graph.distance.complete \code{N*N} distance matrix (the [i,j] element is the length of the shortest path between points i and j)
#' @param data.grid.distance \code{N*dim(grid)[1]} distance matrix between locations where the datum has been observed and locations where
#' @param N_samples number of samples
#' @param p dimension of the manifold matrices
#' @param aggregation_mean "Weighted" ...
#' @param metric_ts metric used on the tangent space. It must be chosen among "Frobenius", "FrobeniusScaled", "Correlation"
#' @param tol tolerance for the main loop of model_kriging
#' @param max_it maximum number of iterations for the main loop of model_kriging
#' @param n_h number of bins in the empirical variogram
#' @param tolerance_intrinsic tolerance for the computation of the intrinsic mean. Not needed if Sigma is provided
#' @param X matrix (N rows and unrestricted number of columns) of additional covariates for the tangent space model, possibly NULL
#' @param X_new matrix (with the same number of rows of \code{new_coords}) of additional covariates for the new locations, possibly NULL
#' @param create_pdf_vario boolean. If \code{TRUE} the empirical and fitted variograms are plotted in a pdf file
#' @param pdf_parameters list with the fields \code{test_nr} and \code{sample_draw}. Additional parameters to name the pdf
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot"
#' @param model_ts type of model fitted on the tangent space. It must be chosen among "Intercept", "Coord1", "Coord2", "Additive"
#' @param vario_model type of variogram fitted. It must be chosen among "Gaussian", "Spherical", "Exponential"
#' @param distance type of distance between coordinates. It must be either "Eucldist" or "Geodist"
#' @return it returns a list with the following fields
#' \itemize{
#'       \item \code{resBootstrap}{...}
#'       \item \code{resAggregated}{...}
#'       \item \code{model_pred}{...}
#' }
#' @useDynLib Manifoldgstat
#' @export
#'
mixed_RDD = function(data_coords, data_val, K, grid, nk_min=1, B=100,
                       # spdist='euclidean',
                       suppressMes=F,
                       ker.width.intrinsic = 0,
                       # mesh,
                       graph.distance.complete,
                       # assign.matrix, no.assg.grid,
                       data.grid.distance,
                       # is.observed, border.length,
                       N_samples, p,
                       aggregation_mean, metric_ts,
                       tol=1e-12, max_it=100, n_h=15, tolerance_intrinsic =10^(-6),
                       X=NULL, X_new=NULL, create_pdf_vario=FALSE, pdf_parameters=NULL,
  # ker.width.vario = 1.5, aggregation_kriging, method.analysis = 'Local mean',
                       metric_manifold, model_ts, vario_model, distance)

{
  if(aggregation_mean != "Weighted" && aggregation_mean != "Equal") stop("aggregation_mean expected to be: Equal or Weighted")

  if(ker.width.intrinsic==0 & aggregation_mean=="Weighted") warning("Aggregation for the mean will use equal weights since ker.width.intrinsic=0")

  # This function implements RDD-OOK procedure, by using the functions RDD_OOK_boot_man and RDD_OOK_aggr_man
  if(K==1)
  {
    resBootstrap=RDD_OOK_boot_man_mixed(data_coords=data_coords, data_val=data_val, K=K, grid=grid, nk_min=nk_min, B=1,
                                  # spdist=spdist,
                                  suppressMes=suppressMes,
                                  ker.width.intrinsic = 0,
                                  # mesh=mesh,
                                  graph.distance.complete= graph.distance.complete,
                                  # assign.matrix=assign.matrix, no.assg.grid=no.assg.grid,
                                  data.grid.distance = data.grid.distance,
                                  # is.observed=is.observed, border.length= border.length,
                                  # p=p,num.signif.entries=num.signif.entries,
                                  metric_ts=metric_ts,
                                  metric_manifold = metric_manifold)
    # tol=tol, max_it = max_it, n_h=n_h, tolerance_intrinsic=tolerance_intrinsic, X=X, X_new=X_new, X_tot=X_tot, plot=plot,
    # ker.width.vario = 0, method.analysis = method.analysis, metric_manifold = metric_manifold, metric_ts = metric_ts, model_ts = model_ts,
    # vario_model = vario_model, distance = distance


    # resAggregated=simplify2array(resBootstrap$fpred)[,,1]
    resAggregated = resBootstrap$fmean[[1]]
  }
  if(K>1)
  {
    resBootstrap=RDD_OOK_boot_man_mixed(data_coords=data_coords, data_val=data_val, K=K, grid=grid, nk_min=nk_min, B=B,
                                  # spdist=spdist,
                                  suppressMes=suppressMes,
                                  ker.width.intrinsic = ker.width.intrinsic,
                                  # mesh=mesh,
                                  graph.distance.complete= graph.distance.complete,
                                  # assign.matrix=assign.matrix, no.assg.grid=no.assg.grid,
                                  data.grid.distance = data.grid.distance,
                                  # is.observed=is.observed, border.length= border.length,
                                  # p=p,num.signif.entries=num.signif.entries,
                                  metric_ts=metric_ts,
                                  metric_manifold = metric_manifold)
    # tol=tol, max_it = max_it, n_h=n_h, tolerance_intrinsic =tolerance_intrinsic, X=X, X_new=X_new, X_tot=X_tot, plot=plot,
    # ker.width.vario = ker.width.vario, method.analysis = method.analysis,
    # metric_manifold = metric_manifold, metric_ts = metric_ts, model_ts = model_ts,
    # vario_model = vario_model, distance = distance
    if (aggregation_mean== "Equal") ker.width.intrinsic=0
    resAggregated=RDD_OOK_aggr_man(fOKBV = resBootstrap$fmean, weights_intrinsic = resBootstrap$kervalues_mean,
                                   ker.width.intrinsic=  ker.width.intrinsic) # p=p,  num.signif.entries = num.signif.entries

  }

  fmean = resAggregated

  model_pred = model_kriging_mixed (data_manifold = data_val, coords = data_coords, X = NULL, Sigma_data = fmean[1:N_samples], metric_manifold = metric_manifold,
                              model_ts = model_ts, vario_model = vario_model, # metric_ts = "Frobenius",
                              n_h=15, distance = distance, data_dist_mat=graph.distance.complete, data_grid_dist_mat=data.grid.distance, max_it = 100, tolerance = 1e-6, # weight_vario = NULL,
                              # weight_intrinsic = NULL, tolerance_intrinsic = 1e-6,
                              max_sill = NULL, max_a = NULL,
                              new_coords=prediction_grid, Sigma_new = fmean, X_new = NULL, create_pdf_vario = create_pdf_vario, pdf_parameters=pdf_parameters, suppressMes = FALSE)

  return(list(resBootstrap = resBootstrap,
              resAggregated = resAggregated,
              model_pred = model_pred))
}
