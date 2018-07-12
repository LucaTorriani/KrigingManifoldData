### Arguments:
# - GLS_model: list (Sigma_opt = tangent point, beta_opt = beta of the linear model on the tangent space, gamma_matrix = covariogram,
#                    residuals = residuals, par = parameters of the varigoram). All these parameters are returned by the function
#                    model_GLS_sigma_fixed.
# - coords: coordinates used to build the model
# - new_coords: coordinates at which the prediction must be computed
# - model_ts: model on the tangent space ("coord1" or "coord2" or "additive")
# - vario_model: variogram model ("Gaussian" or "Exponential" or "Spherical")
# - metric_manifold: metric used on the manifold ("Log_euclidean" or "Square_root" or "Frobenius")
# - X: dataframe of covariates (possibly NULL)
# - n_h: number of points at which the variogram must be computed
# - distance: distane used ("Geodist" or "Eucldist")
### Value: matrix in which row i represents the prediction corresponding to new_coords[i,]. To get the array of matrices use 
###        matrix_to_matrixArray         
kriging = function(GLS_model, coords, new_coords, model_ts= "additive", vario_model="Gaussian",
                   metric_manifold="Frobenius",X_new = NULL, distance = "Geodist"){
  if(is.null(X_new)){
    X_new = matrix(0,dim(new_coords)[1],0)
  }
  beta = GLS_model$beta_opt
  Sigma = GLS_model$Sigma_opt
  gamma_matrix = GLS_model$gamma_matrix
  residuals = GLS_model$residuals
  par = GLS_model$par
  
  new_coords = (as.matrix(new_coords))
  if (dim(new_coords)[2]!= 2){
    new_coords = t(new_coords)
  }
  design_matrix = create_design_matrix(new_coords, X_new, model_ts)
  m = dim(new_coords)[1] # Numero di nuove posizione in cui fare predizione
  prediction = matrix(0,nrow = m, ncol = dim(residuals)[2])
  for(i in 1:m){
    ci = compute_ci(coords, new_coords[i,], vario_model, par, distance)
    lambda_vector = as.vector(solve(gamma_matrix,ci))  
    first_two_terms = colSums(beta*as.vector(design_matrix[i,]))
    second_term = colSums(lambda_vector*residuals)
    prediction[i,] = matrix_to_vec(exponential_map(vec_to_matrix(first_two_terms + second_term), Sigma, metric_manifold))
  }
  return(prediction)
}