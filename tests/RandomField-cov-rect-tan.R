load("/vagrant/KrigingManifoldData/data/RandomField-Cov-rect-tan.Rdata")
load("/vagrant/KrigingManifoldData/data/RandomSample-Cov-rect-tan.Rdata")
coords_tot = gridCov
data_manifold_tot = fieldCov

dyn.load("/vagrant/KrigingManifoldData/src/interface_function.so")

source ("/vagrant/KrigingManifoldData/R/kriging.R")
source ("/vagrant/KrigingManifoldData/R/model_GLS.R")
source ("/vagrant/KrigingManifoldData/R/model_kriging.R")
source ("/vagrant/KrigingManifoldData/R/plot_variogram.R")


### Model parameters
metric_manifold = "Frobenius"      # "LogEuclidean" or "SquareRoot" or "Frobenius"
metric_ts = "Frobenius"            # "Frobenius" or "FrobeniusScaled"
vario_model="Spherical"            # "Gaussian" or "Spherical" or "Exponential"

distance = "Eucldist"              # "Eucldist" or "Geodist"
model_ts = "Coord1"                # "Intercept" or "Coord1" or "Coord2" or "Additive"

### Data manifold: matrix (Nx3) or array of matrices (Nxpxp or pxpxN)   (N= number of observations)
# data_manifold_model = data_manifold_tot
library(plyr)

#data_manifold_model = alply(rCov,3)
data_manifold_model = rCov
### Tangent point
# Sigma = matrix(c(15,11,11,27),2,2)
Sigma = matrix(c(2,1,1,1), 2,2)
# Sigma2 = intrinsic_mean(data_manifold_tot, metric_manifold, metric_ts)  # DA IMPLEMENTARE IN C++
# Sigma3 = intrinsic_mean(data_manifold_model, metric_manifold, metric_ts)
# Sigma = metodo Pigoli Secchi

### Coordinates where the data are observed
#coords_model = gridCov  # Total grid
coords_model = rGrid    # Sample of the grid

model = model_GLS(data_manifold = data_manifold_model, coords = coords_model, Sigma = Sigma,
                                       distance = distance, metric_manifold = metric_manifold,
                                       metric_ts = metric_ts, model_ts = model_ts,
                                       vario_model = vario_model, n_h = 15,
                                       max_it = 100, tolerance = 10^(-4))
print(model$beta)
print(model$fitted_par_vario)
print(model$iterations)
print(model$Sigma)


#model = list(beta = model$beta, gamma_matrix = model$gamma_matrix,
#             residuals = model$residuals, fitted_par_vario = model$fitted_par_vario, Sigma = model$Sigma)


if(FALSE){
# ##############################################################################
# ################################ KRIGING #####################################
# ##############################################################################
#
# set.seed(2165)
# sample_idx = sample(1:1000, 300)
# coords_krig = coords_tot[sample_idx,]
# data_manifold_krig = matrixArray_to_matrix(data_manifold_tot)[sample_idx, ]
#
# prediction = kriging (GLS_model = model, coords = coords_model, new_coords = coords_krig, model_ts=model_ts,
#                       vario_model= vario_model, metric_manifold = metric_manifold, distance=distance)

##############################################################################
################################### Test #####################################
##############################################################################
coords_krig = coords_model
prediction_model = kriging (GLS_model = model, coords = coords_model, new_coords = coords_krig, model_ts=model_ts,
                            vario_model= vario_model, metric_manifold = metric_manifold, distance=distance)
coords_krig = coords_tot
prediction_tot = kriging (GLS_model = model, coords = coords_model, new_coords = coords_krig, model_ts=model_ts,
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
result = model_kriging(data_manifold = data_manifold_model, coords = coords_model,
                                       Sigma = Sigma, distance = distance, metric_manifold = metric_manifold,
                                       metric_ts = metric_ts, model_ts = model_ts,
                                       vario_model = vario_model, n_h = 15,
                                       max_it = 100, tolerance = 10^(-4), new_coords = coords_krig)
print (result$beta)
print (result$par)

coords_krig = coords_tot
result1 = model_kriging(data_manifold = data_manifold_model, coords = coords_model,
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
