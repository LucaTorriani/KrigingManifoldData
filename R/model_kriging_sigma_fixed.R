model_kriging_sigma_fixed = function(data_manifold, coords,X = NULL, Sigma, metric_manifold="Frobenius",
                                 metric_ts = "Frobenius", model_ts="additive", vario_model="Gaussian",
                                 n_h=15, distance="Geodist", max_it = 100, tolerance = 10^(-4), weight=NULL,
                                 new_coords, X_new = NULL){

  if ( distance == "Geodist" & dim(coords)[2] != 2){
    stop("Geodist without two coordinates")
  }
  if(!is.null(X)) {X = as.matrix(X)}
  if(!is.null(X_new)) {X_new = as.matrix(X_new)}

  coords = as.matrix(coords)
  new_coords = as.matrix(new_coords)

  result =.Call("get_model_and_kriging",data_manifold, coords,X, Sigma, distance, metric_manifold, metric_ts, model_ts, vario_model,
                n_h, max_it, tolerance, weight, new_coords, X_new )

  empirical_variogram = list(emp_vario_values = result$emp_vario_values, h = result$h_vec)
  fitted_variogram = list(fit_vario_values = result$fit_vario_values, hh = result$hh)

  plot_variogram(empirical_variogram = empirical_variogram, fitted_variogram = fitted_variogram, model = vario_model,
                 distance = distance)

  return (list(beta_opt = result$beta, gamma_matrix = result$gamma_matrix, residuals = result$residuals,
               par = result$vario_parameters, iter = result$iterations, prediction = result$prediction))
}
