#' Create a GLS model
#'
#' @param data_manifold list or array [\code{p,p,N}] of \code{N} symmetric positive definite matrices of dimension \code{p*p}
#' @param coords \code{N*2} or \code{N*3} matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to
#' be provided in signed decimal degrees
#' @param X matrix (N rows and unrestricted number of columns) of additional covariates for the tangent space model, possibly NULL
#' @param Sigma \code{p*p} matrix representing the tangent point. If NULL the tangent point is computed as the intrinsic mean
#' of \code{data_manifold}
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot", "Correlation"
#' @param metric_ts metric used on the tangent space. It must be chosen among "Frobenius", "FrobeniusScaled", "Correlation"
#' @param model_ts type of model fitted on the tangent space. It must be chosen among "Intercept", "Coord1", "Coord2", "Additive"
#' @param vario_model type of variogram fitted. It must be chosen among "Gaussian", "Spherical", "Exponential"
#' @param n_h number of bins in the emprical variogram
#' @param distance type of distance between coordinates. It must be either "Eucldist" or "Geodist"
#' @param data_dist_mat Matrix of dimension \code{N*N} of distances between data points. If not provided it is computed using \code{distance}
#' @param max_it max number of iterations for the main loop
#' @param tolerance tolerance for the main loop
#' @param weight_intrinsic vector of length \code{N} to weight the locations in the computation of the intrinsic mean. If NULL
#' a vector of ones is used. Not needed if Sigma is provided
#' @param tolerance_intrinsic tolerance for the computation of the intrinsic mean. Not needed if Sigma is provided
#' @param max_sill maximum value allowed for \code{sill} in the fitted variogram. If NULL it is defined as \code{1.15*max(emp_vario_values)}
#' @param max_a maximum value for \code{a} in the fitted variogram. If NULL it is defined as \code{1.15*h_max}
#' @param param_weighted_vario List of 7 elements to be provided to consider Kernel weights for the variogram:
#' \code{weight_vario} (vector of length \code{N_tot} to weight the locations in the computation of the empirical variogram),
#' \code{distance_matrix_tot} (\code{N_tot*N_tot} matrix of distances between the locations),
#' \code{data_manifold_tot}
#' (list or array [\code{p,p,N_tot}] of \code{N_tot} symmetric positive definite matrices of dimension \code{p*p}),
#' \code{coords_tot} (\code{N_tot*2} or \code{N_tot*3} matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to
#' be provided in signed decimal degrees),
#' \code{X_tot} (matrix with N_tot rows and unrestricted number of columns, of additional covariates for the tangent space model. Possibly NULL),
#' \code{h_max} (maximum value of distance for which the variogram is computed)
#' \code{indexes_model} (indexes corresponding to \code{coords} in \code{coords_tot}). Required only in the case
#' \code{metric_manifold = "Correlation"}
#' @param plot boolean. If \code{TRUE} the empirical and fitted variograms are plotted
#' @param suppressMes boolean. If \code{TRUE} warning messagges are not printed
#' @param weight_extrinsic vector of length \code{N} to weight the locations in the computation of the extrinsic mean. If NULL
#' weight_intrinsic are used. Needed only if Sigma is not provided and \code{metric_manifold== "Correlation"}
#' @param tolerance_map_cor tolerance to use in the maps. Required only if
#' \code{metric_manifold== "Correlation"}
#' @return A list with the following fields:
#' \item{\code{beta}}{ vector of the beta matrices of the fitted model}
#' \item{\code{gamma_matrix}}{ \code{N*N} covariogram matrix}
#' \item{\code{residuals}}{ vector of the \code{N} residual matrices}
#' \item{\code{emp_vario_values}}{ vector of empircal variogram values in correspondence of \code{h_vec}}
#' \item{\code{h_vec}}{ vector of positions at which the empirical variogram is computed}
#' \item{\code{fitted_par_vario}}{ estimates of \emph{nugget}, \emph{sill-nugget} and \emph{practical range}}
#' \item{\code{iterations}}{ number of iterations of the main loop}
#' \item{\code{Sigma}}{ tangent point}
#' @description Given the coordinates and corresponding manifold values, this function creates a GLS model on the tangent space.
#' @details The manifold values are mapped on the tangent space and then a GLS model is fitted to them. A first estimate of the \code{beta} coefficients
#' is obtained assuming spatially uncorrelated errors. Then, in the main the loop, new estimates of the \code{beta} are obtained as a result of a
#' weighted least square problem where the weight matrix is the inverse of \code{gamma_matrix}. The residuals \code{(residuals = data_ts - fitted)}
#' are updated accordingly. The parameters of the variogram fitted to the residuals (and used in the evaluation of the \code{gamma_matrix}) are
#' computed using Gauss-Newton with backtrack method to solve the associated non-linear least square problem. The stopping criteria is based on the
#' absolute value of the variogram residuals' norm if \code{ker.width.vario=0}, while it is based on its increment otherwise.
#' @references D. Pigoli, A. Menafoglio & P. Secchi (2016):
#' Kriging prediction for manifold-valued random fields.
#' Journal of Multivariate Analysis, 145, 117-131.
#' @examples
#' data_manifold_model <- Manifoldgstat::rCov
#' coords_model <- Manifoldgstat::rGrid
#' Sigma <- matrix(c(2,1,1,1), 2,2)
#' model = model_GLS(data_manifold = data_manifold_model, coords = coords_model, Sigma = Sigma,
#'                  metric_manifold = "Frobenius", metric_ts = "Frobenius", model_ts = "Coord1",
#'                  vario_model = "Spherical", n_h = 15, distance = "Eucldist", max_it = 100,
#'                  tolerance = 1e-7, plot = TRUE)
#' @useDynLib Manifoldgstat
#' @export

model_GLS = function(data_manifold, coords, X = NULL, Sigma = NULL, metric_manifold = "Frobenius",
                     metric_ts = "Frobenius", model_ts = "Additive", vario_model = "Gaussian",
                     n_h=15, distance = "Geodist", max_it = 100, tolerance = 1e-6,
                     weight_intrinsic = NULL, tolerance_intrinsic = 1e-6, max_sill=NULL, max_a=NULL,
                     param_weighted_vario = NULL, plot = FALSE, suppressMes = FALSE, weight_extrinsic=NULL, tolerance_map_cor=1e-6){
  if ((metric_manifold=="Correlation" && metric_ts !="Correlation")
      || (metric_manifold!="Correlation" && metric_ts =="Correlation"))
    stop("Either metric_manifold and metric_ts are both Correlation, or none of them")

    coords = as.matrix(coords)
    N = dim(coords)[1]
    if(is.null(distance)) {
      if (is.null(data_dist_mat))
        stop("If distance is NULL data_dist_mat must be provided")
      else {
        # Controllo dimensioni matrici
        if(dim(data_dist_mat)[1]!=N || dim(data_dist_mat)[2]!=N) stop("data_dist_mat must be an N*N matrix")
      }
    }
    else {
      if (!is.null(data_dist_mat))
        warning("Since distance is not NULL parameter data_dist_mat will be discarded")
      if ( distance == "Geodist" & dim(coords)[2] != 2){
          stop("Geodist requires two coordinates")
      }
    }

  if(is.array(data_manifold)) {
    data_manifold = alply(data_manifold,3)
  }
  if(length(data_manifold) != dim(coords)[1]){
    stop("Dimension of data_manifold and coords must agree")
  }

  if (metric_manifold=="Correlation" && (diag(data_manifold[[1]]) !=rep(1,dim(data_manifold[[1]])[1])))
    stop ("Manifold data must be correlation matrices")


  if(!is.null(X)) {
    X = as.matrix(X)
    check = (dim(X)[1] == dim(coords)[1])
    if(!check) stop("X and coords must have the same number of rows")
  }

  if(is.null(Sigma)){
    if(is.null(weight_intrinsic)) weight_intrinsic = rep(1, length(data_manifold))
  }
  else{
    if(metric_manifold == "Correlation" && (diag(Sigma) != rep(1, dim(Sigma)[1])))
      stop("Sigma must be a correlation matrix")
  }


  # if(metric_manifold=="Correlation" && is.null(weight_extrinsic)) {weight_extrinsic = weight_intrinsic}
  if(is.null(weight_extrinsic)) {weight_extrinsic = weight_intrinsic}

  if(!is.null(param_weighted_vario)){
    if(is.array(param_weighted_vario$data_manifold_tot)) {
      param_weighted_vario$data_manifold_tot = alply(param_weighted_vario$data_manifold_tot,3)
    }
    param_weighted_vario$coords_tot = as.matrix(param_weighted_vario$coords_tot)
    N_tot = length(param_weighted_vario$weight_vario)

    if ( (dim(param_weighted_vario$coords_tot)[1] != N_tot) ||
         length(param_weighted_vario$data_manifold_tot) != N_tot ||
         dim(param_weighted_vario$distance_matrix_tot)[1] != N_tot ||
         dim(param_weighted_vario$distance_matrix_tot)[2] != N_tot){
      stop("Dimensions of weight_vario, coords_tot, data_manifold_tot and distance_matrix_tot must agree")
    }

    if(!is.null(param_weighted_vario$X_tot)) {
      param_weighted_vario$X_tot = as.matrix(param_weighted_vario$X_tot)
      check = (dim(param_weighted_vario$X_tot)[1] == N_tot && dim(param_weighted_vario$X_tot)[2]==dim(X)[2])
      if(!check) stop("X_tot must have the same number of rows of coords_tot and the same number of columns of X")
    }

    if(length(param_weighted_vario) != 7) stop("Param_weight_vario must be a list with length 7")

    result =.Call("get_model",data_manifold, coords,X, Sigma, distance, data_dist_mat, metric_manifold, metric_ts, model_ts, vario_model,
                  n_h, max_it, tolerance, max_sill, max_a, param_weighted_vario$weight_vario, param_weighted_vario$distance_matrix_tot,
                  param_weighted_vario$data_manifold_tot, param_weighted_vario$coords_tot, param_weighted_vario$X_tot,
                  param_weighted_vario$h_max, param_weighted_vario$indexes_model, weight_intrinsic, tolerance_intrinsic, weight_extrinsic, suppressMes, tolerance_map_cor)
  }

  else {
    result =.Call("get_model",data_manifold, coords,X, Sigma, distance, data_dist_mat, metric_manifold, metric_ts, model_ts, vario_model,
                  n_h, max_it, tolerance, max_sill, max_a, weight_vario = NULL, distance_matrix_tot = NULL, data_manifold_tot = NULL,
                  coords_tot = NULL, X_tot = NULL, h_max = NULL, indexes_model=NULL, weight_intrinsic, tolerance_intrinsic, weight_extrinsic, suppressMes, tolerance_map_cor)

  }

  empirical_variogram = list(emp_vario_values = result$emp_vario_values, h = result$h_vec)
  fitted_variogram = list(fit_vario_values = result$fit_vario_values, hh = result$hh)

  if(plot){
    plot_variogram(empirical_variogram = empirical_variogram, fitted_variogram = fitted_variogram, model = vario_model,
                   distance = distance)
  }
  result_list = result[-c(2,3)]
  return (result_list)
}
