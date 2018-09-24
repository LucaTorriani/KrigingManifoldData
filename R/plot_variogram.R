#' Plot the variogram
#' @useDynLib KrigingManifoldData
#' @export
#'
plot_variogram = function (empirical_variogram, fitted_variogram, model, distance) {
  hh = fitted_variogram$hh
  plot(hh[2:length(hh)],fitted_variogram$fit_vario_values[2:length(hh)],  col = 'blue', type = 'l',
       ylim = c(0,1.15*max(empirical_variogram$emp_vario_values)), ylab = model, xlab = distance)
  points(empirical_variogram$h, empirical_variogram$emp_vario_values, pch = 4, col = 'blue')

}
