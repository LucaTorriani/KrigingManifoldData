#' Model given the data
#' @useDynLib Manifoldgstat
#' @export
#'
model_GLS = function(data_manifold, coords, X = NULL, Sigma = NULL, metric_manifold="Frobenius",
                                 metric_ts = "Frobenius", model_ts="additive", vario_model="Gaussian",
                                 n_h=15, distance="Geodist", max_it = 100, tolerance = 1e-6, weight_vario=NULL,
                                 weights_intrinsic = NULL, tolerance_intrinsic = 1e-6, plot = TRUE){

  if ( distance == "Geodist" & dim(coords)[2] != 2){
    stop("Geodist without two coordinates")
  }

  if( is.array(data_manifold_model)){
    data_manifold = alply(data_manifold,3)
  }

  if(!is.null(X)) {X = as.matrix(X)}
  if(is.null(Sigma)){
    if(is.null(weights_intrinsic)) weights_intrinsic = rep(1, length(data_manifold))
  }

  coords = as.matrix(coords)
  result =.Call("get_model",data_manifold, coords,X, Sigma, distance, metric_manifold, metric_ts, model_ts, vario_model,
                        n_h, max_it, tolerance, weight_vario, weights_intrinsic, tolerance_intrinsic )


  empirical_variogram = list(emp_vario_values = result$emp_vario_values, h = result$h_vec)
  fitted_variogram = list(fit_vario_values = result$fit_vario_values, hh = result$hh)

  if(plot){
  plot_variogram(empirical_variogram = empirical_variogram, fitted_variogram = fitted_variogram, model = vario_model,
                distance = distance)
  }
  #class(result) <- "modelGLS"

  return (result)
}
