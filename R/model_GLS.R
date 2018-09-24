#' Model given the data
#' @useDynLib Manifoldgstat
#' @export
#'
model_GLS = function(data_manifold, coords,X = NULL, Sigma, metric_manifold="Frobenius",
                                 metric_ts = "Frobenius", model_ts="additive", vario_model="Gaussian",
                                 n_h=15, distance="Geodist", max_it = 100, tolerance = 10^(-4), weight=NULL, plot = TRUE){

  if ( distance == "Geodist" & dim(coords)[2] != 2){
    stop("Geodist without two coordinates")
  }
  if(!is.null(X)) {X = as.matrix(X)}
  coords = as.matrix(coords)
  result =.Call("get_model",data_manifold, coords,X, Sigma, distance, metric_manifold, metric_ts, model_ts, vario_model,
                        n_h, max_it, tolerance,weight )


  empirical_variogram = list(emp_vario_values = result$emp_vario_values, h = result$h_vec)
  fitted_variogram = list(fit_vario_values = result$fit_vario_values, hh = result$hh)

  if(plot){
  plot_variogram(empirical_variogram = empirical_variogram, fitted_variogram = fitted_variogram, model = vario_model,
                distance = distance)
  }
  #class(result) <- "modelGLS"

  return (result)
}
