#' Plot empirical and fitted variogram
#' @param empirical_variogram A list containing the two following fields: 
#' - h_vec: vector of positions at which the empirical variogram is computed
#' - emp_vario_values: vector of empircal variogram values in correspondence of \code{h_vec}
#' @param fitted_variogram A list containing the two following fields: 
#' - hh: dense vector of positions at which \code{fit_vario_values} is computed
#' - fit_vario_values: Vector of fitted variogram values in correspondence of \code{hh}
#' @param model Type of variogram used for fitting (it will be reported on the y-axis). It can be "Gaussian", "Spherical" or "Exponential"
#' @param distance Type of distance used to compute \code{h_vec} (it will be reported on the x-axis). It must be either "Eucldist" or "Geodist" 
#' @useDynLib Manifoldgstat
plot_variogram = function (empirical_variogram, fitted_variogram, model, distance) {
  hh = fitted_variogram$hh
  plot(hh[2:length(hh)],fitted_variogram$fit_vario_values[2:length(hh)],  col = 'blue', type = 'l',
       ylim = c(0,1.15*max(empirical_variogram$emp_vario_values)), ylab = model, xlab = distance)
  points(empirical_variogram$h, empirical_variogram$emp_vario_values, pch = 4, col = 'blue')
}
