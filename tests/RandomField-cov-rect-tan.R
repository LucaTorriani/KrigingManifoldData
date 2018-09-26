load("/vagrant/KrigingManifoldData/data/RandomField-Cov-rect-tan.Rdata")  
load("/vagrant/KrigingManifoldData/data/RandomSample-Cov-rect-tan.Rdata")
coords_tot = gridCov
data_manifold_tot = fieldCov

dyn.load("/vagrant/KrigingManifoldData/src/interface_function.so")

# source ("/vagrant/KrigingManifoldData/R/intermediate_functions.R")

### Model parameters
metric_manifold = "Frobenius"      # "LogEuclidean" or "SquareRoot" or "Frobenius"
metric_ts = "Frobenius"            # "Frobenius" or "FrobeniusScaled"
vario_model="Spherical"            # "Gaussian" or "Spherical" or "Exponential"

distance = "Eucldist"              # "Eucldist" or "Geodist"
model_ts = "Coord1"                # "Intercept" or "Coord1" or "Coord2" or "Additive"

### Data manifold: matrix (Nx3) or array of matrices (Nxpxp or pxpxN)   (N= number of observations)
# data_manifold_model = data_manifold_tot  
library(plyr)

data_manifold_model = alply(rCov,3)

# data_manifold_matrix = matrixArray_to_matrix(data_manifold_model)
# data_manifold_array = matrix_to_matrixArray (data_manifold_matrix, p=2)

### Tangent point
# Sigma = matrix(c(15,11,11,27),2,2)
Sigma = matrix(c(2,1,1,1), 2,2) 
# Sigma2 = intrinsic_mean(data_manifold_tot, metric_manifold, metric_ts)  # DA IMPLEMENTARE IN C++
# Sigma3 = intrinsic_mean(data_manifold_model, metric_manifold, metric_ts) 
# Sigma = metodo Pigoli Secchi

### Coordinates where the data are observed
#coords_model = gridCov  # Total grid
coords_model = rGrid    # Sample of the grid

##############################################################################
################################## MODEL #####################################
##############################################################################

model_GLS = function(data_manifold, coords, X = NULL, Sigma = NULL, metric_manifold="Frobenius",
                     metric_ts = "Frobenius", model_ts="additive", vario_model="Gaussian",
                     n_h=15, distance="Geodist", max_it = 100, tolerance = 1e-6, weight_vario=NULL,
                     weights_intrinsic = NULL, tolerance_intrinsic = 1e-6, plot = TRUE){
  
  if ( distance == "Geodist" & dim(coords)[2] != 2){
    stop("Geodist without two coordinates")
  }
  
  if( is.array(data_manifold_model)){
    data_manifold = alply(data_manifold,3)
  }
  
  if(!is.null(X)) {X = as.matrix(X)}
  if(is.null(Sigma)){
    if(is.null(weights_intrinsic)) weights_intrinsic = rep(1, length(data_manifold))
  }
  
  coords = as.matrix(coords)
  result =.Call("get_model",data_manifold, coords,X, Sigma, distance, metric_manifold, metric_ts, model_ts, vario_model,
                n_h, max_it, tolerance, weight_vario, weights_intrinsic, tolerance_intrinsic )
  
  
  empirical_variogram = list(emp_vario_values = result$emp_vario_values, h = result$h_vec)
  fitted_variogram = list(fit_vario_values = result$fit_vario_values, hh = result$hh)
  
  if(plot){
    plot_variogram(empirical_variogram = empirical_variogram, fitted_variogram = fitted_variogram, model = vario_model,
                   distance = distance)
  }
  #class(result) <- "modelGLS"
  
  return (result)
}

plot_variogram = function (empirical_variogram, fitted_variogram, model, distance) {
  hh = fitted_variogram$hh
  plot(hh[2:length(hh)],fitted_variogram$fit_vario_values[2:length(hh)],  col = 'blue', type = 'l',
       ylim = c(0,1.15*max(empirical_variogram$emp_vario_values)), ylab = model, xlab = distance)
  points(empirical_variogram$h, empirical_variogram$emp_vario_values, pch = 4, col = 'blue')
  
}


beta_gamma_opt = model_GLS(data_manifold = data_manifold_model, coords = coords_model,  
                                       Sigma = Sigma, distance = distance, metric_manifold = metric_manifold, 
                                       metric_ts = metric_ts, model_ts = model_ts,
                                       vario_model = vario_model, n_h = 15,  
                                       max_it = 100, tolerance = 10^(-4))
# print(beta_gamma_opt$beta)
# print(beta_gamma_opt$vario_parameters)
# print(beta_gamma_opt$iterations)


# Rcpp::List result = Rcpp::List::create(Rcpp::Named("beta") = beta_vec_matrices,
#                                        Rcpp::Named("fit_vario_values") = fit_vario_values,
#                                        Rcpp::Named("hh") = hh,
#                                        Rcpp::Named("gamma_matrix") = gamma_matrix,
#                                        Rcpp::Named("residuals") = resVec,
#                                        Rcpp::Named("emp_vario_values") = emp_vario.get_emp_vario_values(),
#                                        Rcpp::Named("h_vec") = emp_vario.get_hvec(),
#                                        Rcpp::Named("vario_parameters") = parameters,
#                                        Rcpp::Named("iterations") = num_iter,
#                                        Rcpp::Named("Sigma")= Sigma);

if (FALSE) {
  
model = list(beta_opt = beta_gamma_opt$beta, gamma_matrix = beta_gamma_opt$gamma_matrix,
             residuals = beta_gamma_opt$residuals, fitted_par_vario = beta_gamma_opt$par)



# ##############################################################################
# ################################ KRIGING #####################################
# ##############################################################################
# 
# set.seed(2165)
# sample_idx = sample(1:1000, 300)
# coords_krig = coords_tot[sample_idx,]
# data_manifold_krig = matrixArray_to_matrix(data_manifold_tot)[sample_idx, ]
# 
# prediction = kriging (GLS_model = model, Sigma = Sigma, coords = coords_model, new_coords = coords_krig, model_ts=model_ts,
#                       vario_model= vario_model, metric_manifold = metric_manifold, distance=distance)

##############################################################################
################################### Test #####################################
##############################################################################
coords_krig = coords_model
prediction_model = kriging (GLS_model = model, Sigma = Sigma, coords = coords_model, new_coords = coords_krig, model_ts=model_ts,
                            vario_model= vario_model, metric_manifold = metric_manifold, distance=distance)
coords_krig = coords_tot
prediction_tot = kriging (GLS_model = model, Sigma = Sigma, coords = coords_model, new_coords = coords_krig, model_ts=model_ts,
                          vario_model= vario_model, metric_manifold = metric_manifold, distance=distance)

x.min=min(gridCov[,1])
x.max=max(gridCov[,1])
y.min=min(gridCov[,2])
y.max=max(gridCov[,2])
dimgrid=dim(gridCov)[1]
radius = 0.02  # 0.008
library(fields)

par(cex=1.25)
plot(0,0, asp=1, col=tim.colors(100), ylim=c(y.min,y.max),
     xlim=c(x.min, x.max), pch='', xlab='', ylab='', main = "Real Values")
for(i in 1:dimgrid)
{
  if(i %% 3 == 0)
  {
    car::ellipse(c(gridCov[i,1],gridCov[i,2]) , fieldCov[,,i],
                 radius=radius, center.cex=.5, col='navyblue')
  }
}
rect(x.min, y.min, x.max, y.max)

for(i in 1:250)
{

  car::ellipse(c(rGrid[i,1],rGrid[i,2]) , rCov[,,i],
               radius=radius, center.cex=.5, col='green')

}
rect(x.min, y.min, x.max, y.max)

#######################################################################################
## Tutto il campo + 250 fittati

par(cex=1.25)
plot(0,0, asp=1, col=tim.colors(100), ylim=c(y.min,y.max),
     xlim=c(x.min, x.max), pch='', xlab='', ylab='',main = "Predicted values")
for(i in 1:dimgrid)
{
  if(i %% 3 == 0)
  {
    car::ellipse(c(gridCov[i,1],gridCov[i,2]) , (prediction_tot$prediction[[i]]),
                 radius=radius, center.cex=.5, col='navyblue' )
  }
}
rect(x.min, y.min, x.max, y.max)

for(i in 1:250)
{

  car::ellipse(c(rGrid[i,1],rGrid[i,2]) , (prediction_model$prediction[[i]]),
               radius=radius, center.cex=.5, col='red')

}
rect(x.min, y.min, x.max, y.max)

# ###############################################################################
# ############################## Test single value ##############################
# ###############################################################################
# index = 152 # Un numero da 1 a 250
# predicted_value = vec_to_matrix(prediction_model[index,])
# predicted_value
# true_value = data_manifold_model[,,index]
# true_value


##############################################################################
############################# MODEL and KRIGING ##############################
##############################################################################
coords_krig = coords_model
result = model_kriging_sigma_fixed(data_manifold = data_manifold_model, coords = coords_model,  
                                       Sigma = Sigma, distance = distance, metric_manifold = metric_manifold, 
                                       metric_ts = metric_ts, model_ts = model_ts,
                                       vario_model = vario_model, n_h = 15,  
                                       max_it = 100, tolerance = 10^(-4), new_coords = coords_krig)
print (result$beta)
print (result$par)

coords_krig = coords_tot
result1 = model_kriging_sigma_fixed(data_manifold = data_manifold_model, coords = coords_model,  
                                   Sigma = Sigma, distance = distance, metric_manifold = metric_manifold, 
                                   metric_ts = metric_ts, model_ts = model_ts,
                                   vario_model = vario_model, n_h = 15,  
                                   max_it = 100, tolerance = 10^(-4), new_coords = coords_krig)

x.min=min(gridCov[,1])
x.max=max(gridCov[,1])
y.min=min(gridCov[,2])
y.max=max(gridCov[,2])
dimgrid=dim(gridCov)[1]
radius = 0.02  # 0.008
library(fields)

par(cex=1.25)
plot(0,0, asp=1, col=tim.colors(100), ylim=c(y.min,y.max),
     xlim=c(x.min, x.max), pch='', xlab='', ylab='', main = "Real Values")
for(i in 1:dimgrid)
{
  if(i %% 3 == 0)
  {
    car::ellipse(c(gridCov[i,1],gridCov[i,2]) , fieldCov[,,i],
                 radius=radius, center.cex=.5, col='navyblue')
  }
}
rect(x.min, y.min, x.max, y.max)

for(i in 1:250)
{
  
  car::ellipse(c(rGrid[i,1],rGrid[i,2]) , rCov[,,i],
               radius=radius, center.cex=.5, col='green')
  
}
rect(x.min, y.min, x.max, y.max)

#######################################################################################
## Tutto il campo + 250 fittati

par(cex=1.25)
plot(0,0, asp=1, col=tim.colors(100), ylim=c(y.min,y.max),
     xlim=c(x.min, x.max), pch='', xlab='', ylab='',main = "Predicted values")
for(i in 1:dimgrid)
{
  if(i %% 3 == 0)
  {
    car::ellipse(c(gridCov[i,1],gridCov[i,2]) , (result1$prediction[[i]]),
                 radius=radius, center.cex=.5, col='navyblue' )
  }
}
rect(x.min, y.min, x.max, y.max)

for(i in 1:250)
{
  
  car::ellipse(c(rGrid[i,1],rGrid[i,2]) , (result$prediction[[i]]),
               radius=radius, center.cex=.5, col='red')
  
}
rect(x.min, y.min, x.max, y.max)
}
