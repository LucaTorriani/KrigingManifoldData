#' Perform main routine for mixed_RDD
#' @param data_manifold array [\code{p,p,N}] of \code{N} symmetric positive definite matrices of dimension \code{p*p}
#' @param coords \code{N*2} or \code{N*3} matrix of [lat,long], [x,y] or [x,y,z] coordinates.
#' @param X matrix (N rows and unrestricted number of columns) of additional covariates for the tangent space model, possibly NULL
#' @param Sigma_data List of the \code{N} fictional tangent points in correspondence with the data used to build the model
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot"
#' @param model_ts type of model fitted on the tangent space. It must be chosen among "Intercept", "Coord1", "Coord2", "Additive"
#' @param vario_model type of variogram fitted. It must be chosen among "Gaussian", "Spherical", "Exponential"
#' @param n_h number of bins in the empirical variogram
#' @param distance type of distance between coordinates. It must be either "Eucldist" or "Geodist"
#' @param data_dist_mat \code{N*N} distance matrix (the [i,j] element is the length of the shortest path between points i and j)
#' @param data_grid_dist_mat \code{N*M} distance matrix between locations where the datum has been observed and locations where
#' the datum has to be predicted
#' @param max_it maximum number of iterations for the main loop of \code{model_kriging}
#' @param tolerance tolerance for the main loop of \code{model_kriging}
#' @param max_sill max value allowed for \code{sill} in the fitted variogram. If NULL it is defined as \code{1.15*max(emp_vario_values)}
#' @param max_a maximum value for \code{a} in the fitted variogram. If NULL it is defined as \code{1.15*h_max}
#' @param new_coords prediction grid, i.e. \code{M*2} or \code{M*3} matrix of coordinates where to predict
#' @param Sigma_new List of the \code{M} fictional tangent points in correspondence with the locations where we want to predict
#' @param X_new matrix (with the same number of rows of \code{new_coords}) of additional covariates for the new locations, possibly NULL
#' @param create_pdf_vario boolean. If \code{TRUE} the empirical and fitted variograms are plotted in a pdf file
#' @param pdf_parameters list with the fields \code{test_nr} and \code{sample_draw}. Additional parameters to name the pdf
#' @param suppressMes \{\code{TRUE}, \code{FALSE}\} controls the level of interaction and warnings given
#' @return it returns a list with the following fields
#' \item{\code{beta}}{ vector of the beta matrices of the fitted model}
#' \item{\code{gamma_matrix}}{ \code{N*N} covariogram matrix}
#' \item{\code{residuals}}{ vector of the \code{N} residual matrices}
#' \item{\code{emp_vario_values}}{ vector of empircal variogram values in correspondence of \code{h_vec}}
#' \item{\code{h_vec}}{ vector of positions at which the empirical variogram is computed}
#' \item{\code{fitted_par_vario}}{ estimates of \emph{nugget}, \emph{sill-nugget} and \emph{practical range}}
#' \item{\code{iterations}}{ number of iterations of the main loop}
#' \item{\code{prediction}}{ vector of matrices predicted at the new locations}
#' @details ...
#' @description ...
#' @useDynLib Manifoldgstat
model_kriging_mixed = function(data_manifold, coords, X = NULL, Sigma_data, metric_manifold = "Frobenius",
                         model_ts = "Additive", vario_model = "Gaussian", # metric_ts = "Frobenius",
                         n_h=15, distance = NULL, data_dist_mat=NULL,
                         data_grid_dist_mat=NULL, max_it = 100, tolerance = 1e-6, # weight_vario = NULL,
                         # weight_intrinsic = NULL, tolerance_intrinsic = 1e-6,
                         max_sill=NULL, max_a=NULL,
                         new_coords, Sigma_new, X_new = NULL, create_pdf_vario = TRUE, pdf_parameters=NULL, suppressMes = FALSE){
  coords = as.matrix(coords)
  new_coords = as.matrix(new_coords)
  N = dim(coords)[1]
  M = dim(new_coords)[1]
  if(is.null(distance)) {
    if ((is.null(data_grid_dist_mat)+is.null(data_dist_mat))!=0)
      stop("If distance is NULL data_dist_mat and data_grid_dist_mat must be provided")
    else {
      # Controllo dimensioni matrici
      if(dim(data_dist_mat)[1]!=N || dim(data_dist_mat)[2]!=N) stop("data_dist_mat must be an N*N matrix")
      if(dim(data_grid_dist_mat)[1]!=N || dim(data_grid_dist_mat)[2]!=M) stop("data_dist_mat must be an N*M matrix")
    }
  }
  else {
    if ((is.null(data_grid_dist_mat)+is.null(data_dist_mat))!=2)
      warning("Since distance is not NULL parameters data_dist_mat and data_grid_dist_mat will be discarded")
    if ( distance == "Geodist" & dim(coords)[2] != 2){
        stop("Geodist requires two coordinates")
    }
  }


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

  if(length(data_manifold) != N){
    stop("Dimension of data_manifold and coords must agree")
  }

  result =.Call("get_model_and_kriging_mixed",data_manifold, coords,X, Sigma_data, distance, data_dist_mat, data_grid_dist_mat,  metric_manifold, model_ts, vario_model, # metric_ts
                n_h, max_it, tolerance, max_sill, max_a, new_coords, Sigma_new, X_new, suppressMes ) # weight_vario, weight_intrinsic, tolerance_intrinsic,

  empirical_variogram = list(emp_vario_values = result$emp_vario_values, h = result$h_vec)
  fitted_variogram = list(fit_vario_values = result$fit_vario_values, hh = result$hh)

  if(create_pdf_vario){
    if (is.null(pdf_parameters)) pdf("Variogram-Method-SingleCell.pdf", width=14, height=7)
    else pdf(paste0("Variogram-Method-Mixed-Test_nr-", pdf_parameters$test_nr, "-K-", K, "-Sample_draw-", pdf_parameters$sample_draw,".pdf"), width=14, height=7)
    plot_variogram(empirical_variogram = empirical_variogram, fitted_variogram = fitted_variogram, model = vario_model,
                   distance = distance)
    dev.off()
  }

  result_list = result[-c(2,3)]
  return (result_list)
}
