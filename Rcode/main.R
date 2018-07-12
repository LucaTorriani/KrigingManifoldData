library(ellipse)
library(maps)
library(sp)
library(gstat)
library(rgdal)
library(quadprog)
library(fields)
library(MASS)
library(psych)
library(nloptr)
library(sfsmisc)
library(purrr)
library(proxy)

setwd ('C:/Users/User/Google Drive/Progetto PACS/Codici_puliti/Ultimi dati simulati Menafoglio')
setwd ('/Users/ILARIASARTORI/Google Drive/Progetto PACS/Codici_puliti/Ultimi dati simulati Menafoglio')

### Cambiare directory
load("RandomField-Cov-rect-tan.Rdata")  
load("RandomSample-Cov-rect-tan.Rdata")
coords_tot = gridCov
data_manifold_tot = fieldCov

setwd ('C:/Users/User/Google Drive/Progetto PACS/CODICE_R_given sigma')
setwd ('/Users/ILARIASARTORI/Google Drive/Progetto PACS/CODICE_R_given sigma')
source ("model_GLS_sigma_fixed.R")
source("helper_model_kriging.R")
source ("kriging.R")

### Model parameters
metric_manifold = "Frobenius"      # "Log_euclidean" or "Square_root" or "Frobenius"
metric_ts = "Frobenius"            # "Frobenius" or "Scaled_Frobenius"
vario_model="Spherical"            # "Gaussian" or "Spherical" or "Exponential"
distance = "Eucldist"              # "Eucldist" or "Geodist"
model_ts = "coord1"                # "intercept" or "coord1" or "coord2" or "additive"

### Data manifold: matrix (Nx3) or array of matrices (Nxpxp or pxpxN)   (N= number of observations)
# data_manifold_model = data_manifold_tot  
data_manifold_model = rCov
# data_manifold_matrix = matrixArray_to_matrix(data_manifold_model)
# data_manifold_array = matrix_to_matrixArray (data_manifold_matrix, p=2)

### Tangent point
Sigma = matrix(c(15,11,11,27),2,2)
Sigma = matrix(c(2,1,1,1), 2,2) 
Sigma2 = intrinsic_mean(data_manifold_tot, metric_manifold, metric_ts) 
Sigma3 = intrinsic_mean(data_manifold_model, metric_manifold, metric_ts) 
# Sigma = metodo Pigoli Secchi

### Coordinates where the data are observed
#coords_model = gridCov  # Total grid
coords_model = rGrid    # Sample of the grid

##############################################################################
################################## MODEL #####################################
##############################################################################

beta_gamma_opt = model_GLS_sigma_fixed(data_manifold = data_manifold_model, coords = coords_model, X = NULL, 
                                       Sigma = Sigma, metric_manifold = metric_manifold, 
                                       metric_ts = metric_ts, model_ts = model_ts,
                                       vario_model = vario_model, n_h = 15, distance = distance, 
                                       max_it = 100, tolerance = 10^(-4))

model = list(Sigma_opt = beta_gamma_opt$Sigma, beta_opt = beta_gamma_opt$beta, gamma_matrix = beta_gamma_opt$gamma_matrix,
             residuals = beta_gamma_opt$residuals, par = beta_gamma_opt$par)

##############################################################################
################################ KRIGING #####################################
##############################################################################

set.seed(2165)
sample_idx = sample(1:1000, 300)
coords_krig = coords_tot[sample_idx,]
# data_manifold_krig = matrixArray_to_matrix(data_manifold_tot)[sample_idx, ]

prediction = kriging (GLS_model = model, coords = coords_model, new_coords = coords_krig, model_ts=model_ts,
                      vario_model= vario_model, metric_manifold = metric_manifold, distance=distance)



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

x11()
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

x11()
par(cex=1.25)
plot(0,0, asp=1, col=tim.colors(100), ylim=c(y.min,y.max),
     xlim=c(x.min, x.max), pch='', xlab='', ylab='',main = "Predicted values")
for(i in 1:dimgrid)
{
  if(i %% 3 == 0)
  {
    car::ellipse(c(gridCov[i,1],gridCov[i,2]) , vec_to_matrix(prediction_tot[i,]),
                 radius=radius, center.cex=.5, col='navyblue' )
  }
}
rect(x.min, y.min, x.max, y.max)

for(i in 1:250)
{
  
  car::ellipse(c(rGrid[i,1],rGrid[i,2]) , vec_to_matrix(prediction_model[i,]),
               radius=radius, center.cex=.5, col='red')
  
}
rect(x.min, y.min, x.max, y.max)


###############################################################################
############################## Test single value ##############################
###############################################################################
index = 152 # Un numero da 1 a 250
predicted_value = vec_to_matrix(prediction_model[index,])
predicted_value
true_value = data_manifold_model[,,index]
true_value


