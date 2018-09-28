#' Kriging prediction given the model
#'
#' @param GLS_model the object returned by \code{model_GLS}, or a list containing the fields:
#' \code{Sigma} (tangent point), \code{beta} (vector of the beta matrices of the fitted model),
#' \code{gamma_matrix} (\code{N*N} covariogram matrix), \code{residuals} (vector of the \code{N} residual matrices)
#' \code{fitted_par_vario} (estimates of \emph{nugget}, \emph{sill-nugget} and \emph{practical range})
#' @param coords \code{N*2} or \code{N*3} matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to
#' be provided in signed decimal degrees
#' @param new_coords matrix of coordinates for the new locations where to perform kriging
#' @param model_ts type of model fitted on the tangent space. It must be chosen among "Intercept", "Coord1", "Coord2", "Additive"
#' @param vario_model type of variogram fitted. It must be chosen among "Gaussian", "Spherical", "Exponential"
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot"
#' @param X_new matrix (with the same number of rows of \code{new_coords}) of additional covariates for the new locations, possibly NULL
#' @param distance type of distance between coordinates. It must be either "Eucldist" or "Geodist"
#' @return A list with a single field:
#' \item{\code{prediction}}{vector of matrices predicted at the new locations}
#' @description Given the GLS model kriging prediction on new location is performed.
#' @details The model provided is used to perform simple kriging on the tangent space in correspondence of the new locations.
#' The estimates are then mapped to the manifold to produce the actual prediction.
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
#' model = model_GLS(data_manifold = data_manifold_model, coords = coords_model, Sigma = Sigma,
#'                    metric_manifold = "Frobenius", metric_ts = "Frobenius", model_ts = "Coord1",
#'                    vario_model = "Spherical", n_h = 15, distance = "Eucldist", max_it = 100,
#'                    tolerance = 1e-7, plot = TRUE)
#' result = kriging (GLS_model = model, coords = coords_model, new_coords = coords_model, model_ts="Coord1",
#'                   vario_model= "Spherical", metric_manifold = "Frobenius", distance="Eucldist")
#' result_tot = kriging (GLS_model = model, coords = coords_model, new_coords = coords_tot, model_ts="Coord1",
#'                       vario_model= "Spherical", metric_manifold = "Frobenius", distance="Eucldist")
#'
#' x.min=min(coords_tot[,1])
#' x.max=max(coords_tot[,1])
#' y.min=min(coords_tot[,2])
#' y.max=max(coords_tot[,2])
#' dimgrid=dim(coords_tot)[1]
#' radius = 0.02
#' library(fields)
#'
#' par(cex=1.25)
#' plot(0,0, asp=1, col=tim.colors(100), ylim=c(y.min,y.max), xlim=c(x.min, x.max), pch='', xlab='', ylab='', main = "Real Values")
#' for(i in 1:dimgrid)
#' { if(i %% 3 == 0) { car::ellipse(c(coords_tot[i,1],coords_tot[i,2]), data_manifold_tot[,,i], radius=radius, center.cex=.5, col='navyblue')}}
#' rect(x.min, y.min, x.max, y.max)
#'
#' for(i in 1:250)
#' { car::ellipse(c(coords_model[i,1],coords_model[i,2]), data_manifold_model[,,i], radius=radius, center.cex=.5, col='green')}
#' rect(x.min, y.min, x.max, y.max)
#'
#' par(cex=1.25)
#' plot(0,0, asp=1, col=tim.colors(100), ylim=c(y.min,y.max),xlim=c(x.min, x.max), pch='', xlab='', ylab='',main = "Predicted values")
#' for(i in 1:dimgrid)
#' { if(i %% 3 == 0) { car::ellipse(c(coords_tot[i,1],coords_tot[i,2]), (result_tot$prediction[[i]]), radius=radius, center.cex=.5, col='navyblue' )}}
#' rect(x.min, y.min, x.max, y.max)
#'
#' for(i in 1:250)
#' { car::ellipse(c(coords_model[i,1],coords_model[i,2]), (result$prediction[[i]]), radius=radius, center.cex=.5, col='red')}
#' rect(x.min, y.min, x.max, y.max)
#' @useDynLib Manifoldgstat
#' @export
#'

kriging = function(GLS_model, coords, new_coords, model_ts= "Additive", vario_model="Gaussian",
                   metric_manifold="Frobenius",X_new = NULL, distance = "Geodist") {
  coords = as.matrix(coords)
  new_coords = as.matrix(new_coords)
  if(!is.null(X_new)) {X_new = as.matrix(X_new)}

  result = .Call("get_kriging", coords, new_coords, GLS_model$Sigma, distance, metric_manifold, model_ts, vario_model,
                 GLS_model$beta, GLS_model$gamma_matrix, GLS_model$fitted_par_vario, GLS_model$residuals, X_new)

  return (result)
}
