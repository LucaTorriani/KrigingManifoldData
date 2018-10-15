#' Create a GLS model and directly perform kriging
#' 
#' @param data_manifold list or array [\code{p,p,N}] of \code{N} symmetric positive definite matrices of dimension \code{p*p}
#' @param coords \code{N*2} or \code{N*3} matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to
#' be provided in signed decimal degrees
#' @param X matrix (N rows and unrestricted number of columns) of additional covariates for the tangent space model, possibly NULL
#' @param Sigma \code{p*p} matrix representing the tangent point. If NULL the tangent point is computed as the intrinsic mean
#' of \code{data_manifold}
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot"
#' @param metric_ts metric used on the tangent space. It must be either "Frobenius" or "FrobeniusScaled"
#' @param model_ts type of model fitted on the tangent space. It must be chosen among "Intercept", "Coord1", "Coord2", "Additive"
#' @param vario_model type of variogram fitted. It must be chosen among "Gaussian", "Spherical", "Exponential"
#' @param n_h number of bins in the emprical variogram
#' @param distance type of distance between coordinates. It must be either "Eucldist" or "Geodist"
#' @param max_it max number of iterations for the main loop
#' @param tolerance tolerance for the main loop
#' @param weight_intrinsic vector of length \code{N} to weight the locations in the computation of the intrinsic mean. If NULL
#' a vector of ones is used. Not needed if Sigma is provided
#' @param tolerance_intrinsic tolerance for the computation of the intrinsic mean. Not needed if Sigma is provided
#' @param param_weighted_vario List of 7 elements to be provided to consider Kernel weights for the variogram: 
#' \code{weight_vario} (vector of length \code{N_tot} to weight the locations in the computation of the empirical variogram), 
#' \code{distance_matrix_tot} (\code{N_tot*N_tot} matrix of distances between the locations), 
#' \code{data_tspace_tot} (\code{N_tot*((p*(p+1))/2)} matrix where the i-th row represents projection on the tangent space of the i-th manifold data. It can be computed using .Call("map2_tangent_space")), 
#' \code{coords_tot} (\code{N_tot*2} or \code{N_tot*3} matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to
#' be provided in signed decimal degrees), 
#' \code{X_tot} (matrix with N_tot rows and unrestricted number of columns, of additional covariates for the tangent space model. Possibly NULL), 
#' \code{h_max} (maximum value of distance for which the variogram is computed)
#' \code{indexes_model} (indexes corresponding to \code{coords} in \code{coords_tot})
#' @param new_coords matrix of coordinates for the new locations where to perform kriging
#' @param X_new matrix (with the same number of rows of \code{new_coords}) of additional covariates for the new locations, possibly NULL
#' @param plot boolean. If \code{TRUE} the empirical and fitted variograms are plotted
#' @return list with the following fields:
#' \item{\code{beta}}{ vector of the beta matrices of the fitted model}
#' \item{\code{gamma_matrix}}{ \code{N*N} covariogram matrix}
#' \item{\code{residuals}}{ vector of the \code{N} residual matrices}
#' \item{\code{emp_vario_values}}{ vector of empircal variogram values in correspondence of \code{h_vec}}
#' \item{\code{h_vec}}{ vector of positions at which the empirical variogram is computed}
#' \item{\code{fitted_par_vario}}{ estimates of \emph{nugget}, \emph{sill-nugget} and \emph{practical range}}
#' \item{\code{iterations}}{ number of iterations of the main loop}
#' \item{\code{Sigma}}{ tangent point}
#' \item{\code{prediction}}{ vector of matrices predicted at the new locations}
#' @description Given the coordinates and corresponding manifold values, this function firstly creates a GLS model on the tangent space, and then
#' it performs kriging on the new locations.
#' @details The manifold values are mapped on the tangent space and then a GLS model is fitted to them. A first estimate of the beta coefficients
#' is obtained assuming spatially uncorrelated errors. Then, in the main the loop, new estimates of the beta are obtained as a result of a
#' weighted least square problem where the weight matrix is the inverse of \code{gamma_matrix}. The residuals \code{(residuals = data_ts - fitted)}
#' are updated accordingly. The parameters of the variogram fitted to the residuals (and used in the evaluation of the \code{gamma_matrix}) are
#' computed using Gauss-Newton with backtrack method to solve the associated non-linear least square problem.
#' Once the model is computed, simple kriging on the tangent space is performed in correspondence of the new locations and eventually
#' the estimates are mapped to the manifold.
#' @references D. Pigoli, A. Menafoglio & P. Secchi (2016):
#' Kriging prediction for manifold-valued random fields.
#' Journal of Multivariate Analysis, 145, 117-131.
#' @examples
#' data_manifold_tot <- Manifoldgstat::fieldCov
#' data_manifold_model <- Manifoldgstat::rCov
#' coords_model <- Manifoldgstat::rGrid
#' coords_tot <- Manifoldgstat::gridCov
#' Sigma <- matrix(c(2,1,1,1), 2,2)
#'
#' result = model_kriging (data_manifold = data_manifold_model, coords = coords_model, Sigma = Sigma, metric_manifold = "Frobenius",
#'                         metric_ts = "Frobenius", model_ts = "Coord1", vario_model = "Spherical", n_h = 15, distance = "Eucldist",
#'                         max_it = 100, tolerance = 10e-7,new_coords = coords_model)
#' result_tot = model_kriging (data_manifold = data_manifold_model, coords = coords_model, Sigma = Sigma, metric_manifold = "Frobenius",
#'                             metric_ts = "Frobenius",, model_ts = "Coord1", vario_model = "Spherical", n_h = 15, distance = "Eucldist",
#'                             max_it = 100, tolerance = 10e-7, new_coords = coords_tot, plot = FALSE)
#' x.min=min(coords_tot[,1])
#' x.max=max(coords_tot[,1])
#' y.min=min(coords_tot[,2])
#' y.max=max(coords_tot[,2])
#' dimgrid=dim(coords_tot)[1]
#' radius = 0.02
#'
#' par(cex=1.25)
#' plot(0,0, asp=1, col=fields::tim.colors(100), ylim=c(y.min,y.max), xlim=c(x.min, x.max), pch='', xlab='', ylab='', main = "Real Values")
#' for(i in 1:dimgrid)
#' { if(i %% 3 == 0) { car::ellipse(c(coords_tot[i,1],coords_tot[i,2]) , data_manifold_tot[,,i],radius=radius, center.cex=.5, col='navyblue')}}
#' rect(x.min, y.min, x.max, y.max)
#'
#' for(i in 1:250)
#' { car::ellipse(c(coords_model[i,1],coords_model[i,2]) , data_manifold_model[,,i],radius=radius, center.cex=.5, col='green')}
#' rect(x.min, y.min, x.max, y.max)
#'
#' par(cex=1.25)
#' plot(0,0, asp=1, col=fields::tim.colors(100), ylim=c(y.min,y.max),xlim=c(x.min, x.max), pch='', xlab='', ylab='',main = "Predicted values")
#' for(i in 1:dimgrid)
#' { if(i %% 3 == 0) { car::ellipse(c(coords_tot[i,1],coords_tot[i,2]) , (result_tot$prediction[[i]]),radius=radius, center.cex=.5, col='navyblue' )}}
#' rect(x.min, y.min, x.max, y.max)
#'
#' for(i in 1:250)
#' { car::ellipse(c(rGrid[i,1],rGrid[i,2]) , (result$prediction[[i]]),radius=radius, center.cex=.5, col='red')}
#' rect(x.min, y.min, x.max, y.max)
#' @useDynLib Manifoldgstat
#' @export
#' 
model_kriging = function(data_manifold, coords,  X = NULL, Sigma, metric_manifold = "Frobenius",
                             metric_ts = "Frobenius", model_ts = "Additive", vario_model = "Gaussian",
                             n_h=15, distance = "Geodist", max_it = 100, tolerance = 1e-6, weight_intrinsic = NULL, 
                             tolerance_intrinsic = 1e-6, param_weighted_vario = NULL, new_coords, X_new = NULL, plot = TRUE){
  
  
  if ( distance == "Geodist" & dim(coords)[2] != 2){
    stop("Geodist requires two coordinates")
  }
  coords = as.matrix(coords)
  new_coords = as.matrix(new_coords)
  
  if(!is.null(X)) {
    X = as.matrix(X)
    check = (dim(X)[1] == dim(coords)[1])
    if(!check) stop("X and coords must have the same number of rows")
    if(is.null(X_new)) stop("X and X_new must have the same number of columns")
    else {
      X_new = as.matrix(X_new)
      check = (dim(X_new)[1] == dim(new_coords)[1])
      if(!check) stop("X_new and new_coords must have the same number of rows")
      if (dim(X)[2]!=dim(X_new)[2]) stop("X and X_new must have the same number of columns")
    }
  }
  else {
    if (!is.null(X_new)) stop("X and X_new must have the same number of columns")
  }
  
  if( is.array(data_manifold)){
    data_manifold = alply(data_manifold,3)
  }
  
  if(length(data_manifold) != dim(coords)[1]){
    stop("Dimension of data_manifold and coords must agree")
  }
  
  if(is.null(Sigma)){
    if(is.null(weight_intrinsic)) weight_intrinsic = rep(1, length(data_manifold))
  }
  
  # controllare che else faccia riferimento a if precedente
  
  if(!is.null(param_weighted_vario)){
    param_weighted_vario$coords_tot = as.matrix(param_weighted_vario$coords_tot)
    N_tot = length(param_weighted_vario$weight_vario)
    
    if ( (dim(param_weighted_vario$coords_tot)[1] != N_tot) || 
         dim(param_weighted_vario$data_tspace_tot)[1] != N_tot ||
         dim(param_weighted_vario$distance_matrix_tot)[1] != N_tot ||
         dim(param_weighted_vario$distance_matrix_tot)[2] != N_tot){
      stop("Dimensions of weight_vario, coords_tot, data_tspace_tot and distance_matrix_tot must agree")
    } 
    
    if(!is.null(param_weighted_vario$X_tot)) {
      param_weighted_vario$X_tot = as.matrix(param_weighted_vario$X_tot)
      check = (dim(param_weighted_vario$X_tot)[1] == N_tot && dim(param_weighted_vario$X_tot)[2]==dim(X)[2])
      if(!check) stop("X_tot must have the same number of rows of coords_tot and the same number of columns of X")
    }
    
    if(length(param_weighted_vario) != 7) stop("Param_weighter_vario must be a list with length 7")
    
    result =.Call("get_model_and_kriging",data_manifold, coords,X, Sigma, distance, metric_manifold, metric_ts, model_ts, vario_model,
                  n_h, max_it, tolerance, param_weighted_vario$weight_vario, param_weighted_vario$distance_matrix_tot, 
                  param_weighted_vario$data_tspace_tot, param_weighted_vario$coords_tot, param_weighted_vario$X_tot, 
                  param_weighted_vario$h_max, param_weighted_vario$indexes_model, weight_intrinsic, tolerance_intrinsic, new_coords, X_new )
  }
  
  else {
    result =.Call("get_model_and_kriging",data_manifold, coords,X, Sigma, distance, metric_manifold, metric_ts, model_ts, vario_model,
                  n_h, max_it, tolerance, weight_vario = NULL, distance_matrix_tot = NULL, data_tspace_tot = NULL, 
                  coords_tot = NULL, X_tot = NULL, h_max = NULL, indexes_model = NULL, weight_intrinsic, tolerance_intrinsic, new_coords, X_new)
    
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