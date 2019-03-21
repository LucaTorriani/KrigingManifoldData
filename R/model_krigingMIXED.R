#' @param data_manifold list or array [\code{p,p,N}] of \code{N} symmetric positive definite matrices of dimension \code{nxn}
#' @param coords \code{N*2} or \code{N*3} matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to
#' be provided in signed decimal degrees
#' @param X matrix (N rows and unrestricted number of columns) of additional covariates for the tangent space model, possibly NULL
#' @param Sigma_data List of \code{N} matrices of dimension \code{p*p} representing the tangent points in correspondence of the \code{coords}
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot"
#' @param model_ts type of model fitted on the tangent space. It must be chosen among "Intercept", "Coord1", "Coord2", "Additive"
#' @param vario_model type of variogram fitted. It must be chosen among "Gaussian", "Spherical", "Exponential"
#' @param n_h number of bins in the emprical variogram
#' @param distance type of distance between coordinates. It must be either "Eucldist" or "Geodist"
#' @param data_dist_mat Matrix of dimension \code{N*N} of distances between data points. If not provided it is computed using \code{distance}
#' @param data_grid_dist_mat Matrix of dimension \code{N*M} of distances between data points and grid points. If not provided it is computed using \code{distance}
#' @param max_it max number of iterations for the main loop
#' @param tolerance tolerance for the main loop
#' @param max_sill maximum value allowed for \code{sill} in the fitted variogram. If NULL it is defined as \code{1.15*max(emp_vario_values)}
#' @param max_a maximum value for \code{a} in the fitted variogram. If NULL it is defined as \code{1.15*h_max}
#' @param new_coords matrix of coordinates for the \code{M} new locations where to perform kriging
#' @param Sigma_new List of \code{M} matrices of dimension \code{p*p} representing the tangent points in correspondence of the \code{new_coords}
#' @param X_new matrix (with the same number of rows of \code{new_coords}) of additional covariates for the new locations, possibly NULL
#' @param create_pdf_vario boolean. If \code{TRUE} the empirical and fitted variograms are plotted in a pdf file
#' @param pdf_parameters list with the fields \code{test_nr}, \code{K} and \code{sample_draw}. Additional parameters to name the pdf
#' @param suppressMes boolean. If \code{TRUE} warning messagges are not printed
#' @return list with the following fields:
#' \item{\code{beta}}{ vector of the beta matrices of the fitted model}
#' \item{\code{gamma_matrix}}{ \code{N*N} covariogram matrix}
#' \item{\code{residuals}}{ vector of the \code{N} residual matrices}
#' \item{\code{emp_vario_values}}{ vector of empircal variogram values in correspondence of \code{h_vec}}
#' \item{\code{h_vec}}{ vector of positions at which the empirical variogram is computed}
#' \item{\code{fitted_par_vario}}{ estimates of \emph{nugget}, \emph{sill-nugget} and \emph{practical range}}
#' \item{\code{iterations}}{ number of iterations of the main loop}
#' \item{\code{prediction}}{ vector of matrices predicted at the new locations}
#' @description Given the coordinates and corresponding manifold values, this function firstly creates a GLS model on the tangent space, and then
#' performs kriging on the new locations.
#' @details For all \code{i in 1:N}, \code{data_manifold[,,i]} is mapped to the space tangent in \code{Sigma_data[,,i]}. Those values are, on their turn,
#' parallely transported to the common Hilbert space tangent in the identity matrix, where a GLS model is fitted to them. A first estimate of the beta coefficients
#' is obtained assuming spatially uncorrelated errors. Then, in the main the loop, new estimates of the beta are obtained as a result of a
#' weighted least square problem where the weight matrix is the inverse of \code{gamma_matrix}. The residuals \code{(residuals = data_ts - fitted)}
#' are updated accordingly. The parameters of the variogram fitted to the residuals (and used in the evaluation of the \code{gamma_matrix}) are
#' computed using Gauss-Newton with backtrack method to solve the associated non-linear least square problem.
#' Once the model is computed, simple kriging on the tangent space is performed in correspondence of all the new locations \code{new_coords[j,]}, with \code{j in 1:M}.
#' Then each estimate is parallely transported to the space tangent in \code{Sigma_new[,,j]} and finally mapped to the manifold.
#' @references D. Pigoli, A. Menafoglio & P. Secchi (2016):
#' Kriging prediction for manifold-valued random fields.
#' Journal of Multivariate Analysis, 145, 117-131.
#'
#' O Yair, M Ben-Chen, R Talmon (2018)
#' Parallel transport on the cone manifold of spd matrices for domain adaptation.
#' ArXiv Preprint, 1807.10479
#' @examples
#' data_manifold_tot <- Manifoldgstat::fieldCov
#' data_manifold_model <- Manifoldgstat::rCov
#' coords_model <- Manifoldgstat::rGrid
#' coords_tot <- Manifoldgstat::gridCov
#' Sigma <- matrix(c(2,1,1,1), 2,2)
#'
#' result = model_kriging_mixed (data_manifold = data_manifold_model, coords = coords_model, Sigma = Sigma, metric_manifold = "Frobenius",
#'                         metric_ts = "Frobenius", model_ts = "Coord1", vario_model = "Spherical", n_h = 15, distance = "Eucldist",
#'                         max_it = 100, tolerance = 10e-7,new_coords = coords_model)
#' result_tot = model_kriging_mixed (data_manifold = data_manifold_model, coords = coords_model, Sigma = Sigma, metric_manifold = "Frobenius",
#'                             metric_ts = "Frobenius",, model_ts = "Coord1", vario_model = "Spherical", n_h = 15, distance = "Eucldist",
#'                             max_it = 100, tolerance = 10e-7, new_coords = coords_tot, create_pdf_vario = FALSE)
#'
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
