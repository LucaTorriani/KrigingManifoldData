#' Plot kernel  
#' @param data_coords coordinates of the data
#' @param id the index of the row of \code{data_coords} that will be used as center
#' @param xmax the maximum value for the x-coordinate (the minimum is 0)
#' @param ymax the maximum value for the y-coordinate (the minimum is 0)
#' @param m number of points on the grid in horizontal direction
#' @param n number of points on the grid in vertical direction
#' @param ker.width kernel width
#' @useDynLib Manifoldgstat
#' @export
plot_ker_rect = function(data_coords, id, xmax, ymax, m, n, ker.width) {
  dense_coeff = 2
  xmm_dense = seq(0,xmax,length=dense_coeff*m)
  ynn_dense = seq(0,ymax,length=dense_coeff*n)
  dense_grid = expand.grid(xmm_dense,ynn_dense)
  KER=kerfn(newdata=dense_grid,
            center=c(as.numeric(data_coords[id, 1:2]),id), ker.type = 'Gau',
            param = ker.width)
  
  ker.mat = matrix(KER, ncol=dense_coeff*n, nrow=dense_coeff*m, byrow = F) # Check 
  
  image.plot(xmm_dense, ynn_dense, ker.mat, col = tim.colors(M), asp=1)
  points(data_coords[id, 1:2], pch=20)
}