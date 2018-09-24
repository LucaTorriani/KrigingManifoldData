#' Kriging
#' @useDynLib KrigingManifoldData
#' @export
#'

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
