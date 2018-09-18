dyn.load("/vagrant/KrigingManifoldData/src/interface_function.so")
library("Rcpp")
library("RcppEigen")

plot_variogram = function(empirical_variogram, fitted_par_vario, model, distance){
  x11()
  xx<-seq(0,max(empirical_variogram$h), by = 0.01)
  if(model == 'Gaussian'){
    plot(xx[2:length(xx)],gauss_vario(fitted_par_vario,xx[2:length(xx)]),  col = 'blue', type = 'l', 
         ylim = c(0,1.15*max(empirical_variogram$emp_vario_values)), ylab = "GaussVariogram", xlab = distance)
    points(empirical_variogram$h, empirical_variogram$emp_vario_values, pch = 4, col = 'blue')
  }
  else if(model == 'Exponential'){
    plot(xx[2:length(xx)],exp_vario(fitted_par_vario,xx[2:length(xx)]),  col = 'blue', type = 'l', 
         ylim = c(0,1.15*max(empirical_variogram$emp_vario_values)), ylab = "ExpVariogram", xlab = distance)
    points(empirical_variogram$h, empirical_variogram$emp_vario_values, pch = 4, col = 'blue')
  }
  else if (model == 'Spherical'){
    plot(xx[2:length(xx)],sph_vario(fitted_par_vario,xx[2:length(xx)]),  col = 'blue', type = 'l', 
         ylim = c(0,1.15*max(empirical_variogram$emp_vario_values)), ylab = "SphVariogram", xlab = distance)
    points(empirical_variogram$h, empirical_variogram$emp_vario_values, pch = 4, col = 'blue')
  }
  else {
    stop('Model not available')
  }
}


model_GLS_sigma_fixed = function(data_manifold, coords,X = NULL, Sigma, metric_manifold="Frobenius", 
                                 metric_ts = "Frobenius", model_ts="additive", vario_model="Gaussian", 
                                 n_h=15, distance="Geodist", max_it = 100, tolerance = 10^(-4), weight=NULL){
  
  if ( distance == "Geodist" & dim(coords)[2] != 2){
    stop("Geodist without two coordinates")
  }
  
  if(!is.null(X)) X = as.matrix(X)
  coords = as.matrix(coords)
  result =.Call("get_model",data_manifold, coords,X, Sigma, distance, metric_manifold, metric_ts, model_ts, 
                        n_h, max_it, tolerance,weight )
  
  # empirical_variogram = list(result$emp_vario_values, result$h_vec)
  # fitted_par_vario = result$vario_parameters
  # beta = result$beta
  # W = result$gamma_matrix
  # residuals = result$residuals
  # 
  # 
  # plot_variogram(empirical_variogram = empirical_variogram, fitted_par_vario = fitted_par_vario, model = vario_model,
  #               distance = distance)
  # 
  # return (list(beta_opt = beta, gamma_matrix = W, residuals = residuals, par = fitted_par_vario))
  

}

