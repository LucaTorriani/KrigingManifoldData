source('map_functions.R')
source('distance_functions.R')
source('functions_structure.R')


### Arguments:***NEW*** Function modified on 15th march
# - data_manifold: data on the manifold
# - coords: coordinates 
# - X: dataframe of covariates (possibly NULL)
# - metric_manifold: metric used on the manifold ("Log_euclidean" or "Square_root" or "Frobenius")
# - metric_ts: metric used on the tangent space ("Frobenius" or "Scaled Frobenius")
# - model_ts: model on the tangent space ("coord1" or "coord2" or "additive")
# - vario_model: variogram model ("Gaussian" or "Exponential" or "Spherical")
# - n_h: number of points at which the variogram must be computed
# - distance: distane used ("Geodist" or "Eucldist")
# - max_it: maximum number of iterations to compute the betas 
# - tolerance: tolerance used as stopping criterion in the computation of the betas
# *NEW* argument weight
### Value: list (Sigma = tangent point, beta_opt = estimated beta, gamma_matrix = covariogram, residuals = residuals,
###              fitted_par_vario = list of the estimated parameters of the variogram)
model_GLS_sigma_fixed = function(data_manifold, coords,X = NULL, Sigma, metric_manifold="Frobenius", metric_ts = "Frobenius", 
                                 model_ts="additive", vario_model="Gaussian", n_h=15, distance="Geodist", max_it = 100, 
                                 tolerance = 10^(-4), weight=NULL){
  if (length(dim(data_manifold))==3) {
    data_manifold = matrixArray_to_matrix(data_manifold)
  }
  if (length(dim(data_manifold))!=2 & length(dim(data_manifold))!=3) {
    stop("Data manifold must be either a matrix or a array of matrices")
  }
  
  if(is.null(X)){
    X = matrix(0,dim(coords)[1],0)
  }
  if (distance == "Geodist") {
    vec_station_distance = compute_Geodist(coords)
  }
  else if (distance == "Eucldist") {
    vec_station_distance = compute_Eucldist(coords)
  }
  else{
    stop("Distance not available")
  }
  
  max_dist = max(vec_station_distance)
  W =diag(dim(data_manifold)[1])
  
  data_ts = t(apply(data_manifold, 1, logarithmic_map_vec, Sigma=Sigma, metric_manifold=metric_manifold)) # Dati sul piano tangente
  
  model = compute_beta(data_ts, coords, X,metric_manifold, model_ts, W)
  beta_new = model$coeff
  fit_values = model$fit_values
  iter = 0
  tol = tolerance + 1
  while (iter < max_it && tol >tolerance){
    Res = data_ts - fit_values
    empirical_variogram = emp_vario(Res, coords, vec_station_distance, n_h, metric_manifold, metric_ts, distance, Sigma, weight) # *NEW* argument weight
    fitted_par_vario = get_par_fitted_vario(vario_model, empirical_variogram, max_dist = max_dist)
    W = get_gamma_matrix(vario_model,fitted_par_vario,vec_station_distance, dim(data_ts)[1])
    beta_old = beta_new
    model = compute_beta(data_ts, coords, X,metric_manifold, model_ts, W)
    beta_new = model$coeff
    fit_values = model$fit_values
    iter = iter+1
    tol = 0
    for (i in dim(beta_new)[1]){
      if (metric_ts == "Frobenius") {
        tol = tol + tplane_distance (vec_to_matrix(beta_new[i,]), vec_to_matrix(beta_old[i,]), metric_ts)
      }
      else if (metric_ts =="Scaled_Frobenius") {
        tol = tol + tplane_distance (vec_to_matrix(beta_new[i,]), vec_to_matrix(beta_old[i,]), metric_ts, Sigma)
      }
      else {
        stop ("Tangent space metric not available")
      }
    }
    
  }
  #plot_variogram(empirical_variogram = empirical_variogram, fitted_par_vario = fitted_par_vario, model = vario_model, 
  #               distance = distance)
  return (list(Sigma = Sigma, beta_opt = beta_new, gamma_matrix = W, residuals = Res, par = fitted_par_vario))
}

