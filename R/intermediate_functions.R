dyn.load("/vagrant/KrigingManifoldData/src/interface_function.so")
library("Rcpp")
library("RcppEigen")

plot_variogram = function (empirical_variogram, fitted_variogram, model, distance) {
  hh = fitted_variogram$hh
  plot(hh[2:length(hh)],fitted_variogram$fit_vario_values[2:length(hh)],  col = 'blue', type = 'l', 
       ylim = c(0,1.15*max(empirical_variogram$emp_vario_values)), ylab = model, xlab = distance)
  points(empirical_variogram$h, empirical_variogram$emp_vario_values, pch = 4, col = 'blue')
  
}

model_GLS_sigma_fixed = function(data_manifold, coords,X = NULL, Sigma, metric_manifold="Frobenius", 
                                 metric_ts = "Frobenius", model_ts="additive", vario_model="Gaussian", 
                                 n_h=15, distance="Geodist", max_it = 100, tolerance = 10^(-4), weight=NULL){
  
  if ( distance == "Geodist" & dim(coords)[2] != 2){
    stop("Geodist without two coordinates")
  }
  if(!is.null(X)) {X = as.matrix(X)}
  coords = as.matrix(coords)
  result =.Call("get_model",data_manifold, coords,X, Sigma, distance, metric_manifold, metric_ts, model_ts, vario_model, 
                        n_h, max_it, tolerance,weight )

  
  empirical_variogram = list(emp_vario_values = result$emp_vario_values, h = result$h_vec)
  fitted_variogram = list(fit_vario_values = result$fit_vario_values, hh = result$hh)
  
  fitted_par_vario = result$vario_parameters
  beta = result$beta
  W = result$gamma_matrix
  residuals = result$residuals
  iter = result$iterations

  plot_variogram(empirical_variogram = empirical_variogram, fitted_variogram = fitted_variogram, model = vario_model,
                distance = distance)

  return (list(beta_opt = beta, gamma_matrix = W, residuals = residuals, par = fitted_par_vario, iter = iter))
}

kriging = function(GLS_model, Sigma, coords, new_coords, model_ts= "additive", vario_model="Gaussian",
                   metric_manifold="Frobenius",X_new = NULL, distance = "Geodist") {
  coords = as.matrix(coords)
  new_coords = as.matrix(new_coords)
  if(!is.null(X_new)) {X_new = as.matrix(X_new)}

  result = .Call("get_kriging", coords, new_coords, Sigma, distance, metric_manifold, model_ts, vario_model,
                 GLS_model$beta_opt, GLS_model$gamma_matrix, GLS_model$fitted_par_vario, GLS_model$residuals, X_new)
  # print(result)

  return (result)
}

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
