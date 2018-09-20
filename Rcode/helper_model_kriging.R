### Arguments:
# - x: 2x2xn matrix array
# - tol
# - weight
### Value: the Frechet mean of the n matrices (Gradient descend method)
intrinsic_mean <- function(x, metric_manifold, metric_ts, tol = 1e-6, weight = NULL)
{
  if(length(dim(x))==2) {
    n=1
  }
  if(length(dim(x))==3) {
    n=dim(x)[3]
  }
  if(n==1){
    return(as.matrix(x))
  }
  M_prec = Xi_prec = x[,,1]
  tau=1
  if(is.null(weight)) {
    weight=rep(1, n)
  }
  if (metric_ts == "Scaled_Frobenius") {
    while (norm_mat_frob_scaled(A = Xi_prec, Sigma = M_prec)^2 > tol){
      Xi = (1/sum(weight)) * weight[1] * logarithmic_map(sdp = x[,,1], Sigma = M_prec, metric_manifold)
      for(i in 2:n) {
        Xi = Xi + (1/sum(weight)) * weight[i] * logarithmic_map(sdp = x[,,i], Sigma = M_prec, metric_manifold)
      }
      M <- exponential_map(symm = tau*Xi, Sigma = M_prec, metric_manifold)
      if (norm_mat_frob_scaled(A = Xi, Sigma = M) > norm_mat_frob_scaled(A = Xi_prec, Sigma = M_prec)) {
        tau = tau/2;
        Xi = Xi_prec
      }
      M_prec <- M
      Xi_prec <- Xi
    }
  }
  else if (metric_ts == "Frobenius") {
    while (norm_mat_frob(A = Xi_prec)^2 > tol){
      Xi = (1/sum(weight)) * weight[1] * logarithmic_map(sdp = x[,,1], Sigma = M_prec, metric_manifold)
      for(i in 2:n) {
        Xi = Xi + (1/sum(weight)) * weight[i] * logarithmic_map(sdp = x[,,i], Sigma = M_prec, metric_manifold)
      }
      M <- exponential_map(symm = tau*Xi, Sigma = M_prec, metric_manifold)
      if (norm_mat_frob(A = Xi) > norm_mat_frob(A = Xi_prec)) {
        tau = tau/2;
        Xi = Xi_prec
      }
      M_prec <- M
      Xi_prec <- Xi
    }
  }
  else {
    stop("Tangent space metric not available")
  }

  return(matrix_to_vec(M))
}

### Arguments: ***NEW*** Function modified on 15th march
# - data_ts: data on the tangent space
# - coords: coordinates
# - X: dataframe of covariates  (possibly NULL)
# - metric_manifold: metric used on the manifold ("Log_euclidean" or "Square_root" or "Frobenius")
# - model_ts: model on the tangent space ("coord1" or "coord2" or "additive")
# - W: matrix of weigth in the GLS
### Value: list containing coefficients beta and fitted values
compute_beta = function (data_ts, coords, X, metric_manifold, model_ts, W){
  # Map manifold values on the tangent space

  beta = NULL
  fit_values = NULL

  if(model_ts=="intercept") # *NEW*
  {
    dataset = as.data.frame(cbind(rep(1, dim(coords)[1]), X )) # *NEW* The intercept is included as constant regressor
    nomi = c("Intercept")
    if(dim(X)[2]!=0){
      for(i in 1:(dim(X)[2])){
        nomi = c(nomi, paste('V',i, sep=""))
      }
    }

  }
  else if(model_ts=="coord1")
  {
    dataset = as.data.frame(cbind(rep(1, dim(coords)[1]), coords[,1], X ))
    nomi = c("Intercept", "coord1")
    if(dim(X)[2]!=0){
      for(i in 1:(dim(X)[2])){
        nomi = c(nomi, paste('V',i, sep=""))
      }
    }

  }
  else if(model_ts=="coord2")
  {
    dataset = as.data.frame(cbind(rep(1, dim(coords)[1]), coords[,2], X )) # *NEW* The intercept is included as constant regressor
    nomi = c("Intercept","coord2")
    if(dim(X)[2]!=0){
      for(i in 1:(dim(X)[2])){
        nomi = c(nomi, paste('V',i, sep=""))
      }
    }

  }
  else if(model_ts=="additive")
  {
    dataset = as.data.frame(cbind(rep(1, dim(coords)[1]), coords, X ))
    nomi = c("Intercept","coord1", "coord2")
    if(dim(X)[2]!=0){
      for(i in 1:(dim(X)[2])){
        nomi = c(nomi, paste('V',i, sep=""))
      }
    }
  }

  else{
    stop("Model not available")
  }
  for (i in 1:dim(data_ts)[2])
  {
    fit = lm.gls(data_ts[,i] ~ -1 +., data = dataset, W = W, inverse=TRUE)
    beta = cbind(beta, (fit$coeff))
    fit_values = cbind(fit_values, fit$fitted.values)
  }
  rownames(beta) = nomi
  return (list(coeff = beta, fit_values = fit_values))

}


### Arguments:  ***NEW*** Function modified on 15th march
# - Res: residuals
# - coords: coordinates
# - vec_station_distance: distance between stations (vector form)
# - n_h: number of points where the empirical variogram is computed
# - metric_manifold: metric used on the manifold ("Log_euclidean" or "Square_root" or "Frobenius")
# - metric_ts: metric used on the tangent space
# - distance: type of distance used (Eucldist/Geodist)
# - Sigma: tangent space
### Value: a list containing h, the corresponding values of the variogram and the cardinality of each h
emp_vario = function(Res, coords, vec_station_distance , n_h, metric_manifold, metric_ts, distance, Sigma, weight = NULL){
  N = dim(coords)[1]
  h_max = compute_hmax(coords,distance)
  delta_h = h_max/n_h
  d = seq(0, h_max, by = delta_h)  # vettore degli h+-delta(h)

  h = rep(0,length(d)-1)
  var_values = rep(0,length(d)-1)
  card_Nh = rep(0,length(d)-1) # Cardinalita' di N(h) per cui e' stato possibile calcolare il variogramma

  if (metric_ts == "Frobenius") {
    Sigma =NULL
  }
  if (is.null(weight)) {
    weight = rep(1, N)
  }
  W = weight %*% t(weight)

  for (l in 2:length(d)){
    m = c()
    indKer=NULL #*NEW*
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        dist_loc_ij = vec_station_distance[N*(i-1) - i*(i-1)/2 + j-i]
        if (dist_loc_ij >=d[l-1] & dist_loc_ij<d[l]){ # Se la stazione j e' distante da i circa h
          m = c(m,tplane_distance(vec_to_matrix(Res[i,]),vec_to_matrix(Res[j,]),
                                  metric_ts, Sigma)) # vettore delle distanze tra due "valori" sul piano tangente
          card_Nh[l-1] = card_Nh[l-1]+1 # cardinalita' di N(h)
          indKer = rbind(indKer, c(i,j)) # memorizzo indici *NEW*
        }
      }
    }
    if(length(m)>0){
      var_values[l-1] = sum(W[indKer]*m^2)/(2*sum(W[indKer])) # *NEW* Variogramma empirico PESATO per h fissato
    }
    h[l-1] = (d[l]+d[l-1])/2 # Vettore degli h (in alalternativa si puo' fare come media delle distanze tra location
    # nell'intervallino in questione)
  }
  lab = which(card_Nh > 0) # indice degli h tali per cui hai un valore empirico del variogramma
  var_values = var_values[lab] # Valore del variogramma per gli h per cui e' stato possibile calcolarlo
  h = h[lab] # Valore degli h per cui ? stato possibile calcolare il variogramma
  card_Nh = card_Nh[lab] # Cardinalit? di N(h) per cui ? stato possibile calcolare il variogramma
  out = list(h=h,emp_vario_values=var_values,card_Nh=card_Nh)
  return (out)
}


### Arguments:
# - coords: coordinates
# - distance: type of distance used (Eucldist/Geodist)
### Value: the maximum value of h
compute_hmax = function(coords, distance){
  coords_spatial = SpatialPoints(coords)
  rectangle = bbox(coords_spatial)
  if(distance == "Geodist"){
    dist_max = Geodist(rectangle[,1], rectangle[,2])
  }
  else if (distance == "Eucldist") {
    dist_max = dist(rbind(rectangle[,1], rectangle[,2]))
  }
  else {
    stop("Distance not available")
  }

  return (1/3*(dist_max))
}


### Arguments: ***NEW*** Function modified on 15th march
# - vario_model: type of the model of the variogram (Gaussian or Spherical or Exponential)
# - empirical_variogram: list returned by emp_vario
# - max_dist : maximum distance among stations
### Value: list of the estimated parameters (nugget, sigma2, sill) of the vario_model
get_par_fitted_vario = function(vario_model, empirical_variogram, max_dist){
  par0 = get_init_par(empirical_variogram, vario_model)
  par0 = unlist(par0)
  if(par0[1] == 0) par0[1] = 1e-6 #*NEW*
  if(par0[2] == 0) par0[2] = 1e-6 #*NEW*

  max_sill = 1.15 * max(empirical_variogram$emp_vario_values)
  if (vario_model == "Gaussian"){
    max_a = 1/3*max_dist
    par_fitted_vario = constrOptim (par0, gauss_err, ui = rbind(c(1,0,0), c(-1,-1,0), c(0,1,0), c(0,0,1), c(0,0,-1) ),
                                    ci = c(0, -max_sill, 0, 0, -max_a),
                                    method = "Nelder-Mead",
                                    h = empirical_variogram$h,
                                    emp_vario_values = empiricfcal_variogram$emp_vario_values)$par
  }
  else if (vario_model == "Exponential"){
    max_a = 1/3*max_dist
    par_fitted_vario = constrOptim (par0, exp_err, ui = rbind(c(1,0,0), c(-1,-1,0), c(0,1,0), c(0,0,1), c(0,0,-1) ),
                                    ci = c(0, -max_sill, 0, 0, -max_a),
                                    method = "Nelder-Mead",
                                    h = empirical_variogram$h,
                                    emp_vario_values = empirical_variogram$emp_vario_values)$par
  }
  else if (vario_model == "Spherical"){
    max_a = max_dist
    par_fitted_vario = constrOptim (par0, sph_err, ui = rbind(c(1,0,0), c(-1,-1,0), c(0,1,0), c(0,0,1), c(0,0,-1) ),
                                    ci = c(0, -max_sill, 0, 0, -max_a),
                                    method = "Nelder-Mead",
                                    h = empirical_variogram$h,
                                    emp_vario_values = empirical_variogram$emp_vario_values)$par
  }
  else {
    stop("Model not available")
  }
  return (par_fitted_vario)
}


### Arguments: ***NEW*** Function modified on 15th march
# - empirical_variogram: list returned by emp_vario
# - model: type of the model of the variogram (Gaussian or Spherical or Exponential)
### Value: list of the intial parameters (starting points of the optimization of the get_par_fitted_vario)
get_init_par= function (empirical_variogram, model){
  first_two = head(empirical_variogram$emp_vario_values,2)
  last_four = tail(empirical_variogram$emp_vario_values, 4)
  N_h_first_two = head(empirical_variogram$card_Nh,2)
  N_h_last_four = tail(empirical_variogram$card_Nh,4)

  nugget = get_weighted_median(first_two, N_h_first_two)
  sill = get_weighted_median(last_four, N_h_last_four)

  if(model == "Gaussian" || model == "Exponential"){
    distance_to_sill = abs(empirical_variogram$emp_vario_values - 0.95*sill)
    tol = 0.0505*sill
    a = empirical_variogram$h[min(which(distance_to_sill <= tol))]*(1/3)
  }
  else if(model == "Spherical"){
    distance_to_sill = abs(empirical_variogram$emp_vario_values - sill)
    tol = 0.01*sill
    a = empirical_variogram$h[min(which(distance_to_sill <= tol))]
  }
  else {
    stop("Model not available")
  }
  return (list(nugget = nugget, sigma2 = max(sill-nugget, nugget*1e-3), a = a)) #*NEW*
}

### Arguments:
# - values: vector of values
# - card: vector of the same length of values containing the weigths of each element of values
### Value: wieghted median of the values
get_weighted_median = function(values, card){
  N = sum(card)
  cumCard = cumsum(card)
  if(N%%2){
    idx = min(which(cumCard >= (N+1)/2))
  }
  else {
    idx = min(which(cumCard >= N/2))
  }
  return(values[idx])

}

### Arguments:
# - vario_model: type of the model of the variogram (Gaussian or Spherical or Exponential)
# - fitted_par_vario: list of the fitted parameters of the variogram
# - distances: vector of distances between stations
# - N: number of observarions
### Value: covariogram
get_gamma_matrix = function (vario_model,fitted_par_vario, distances, N){
  if (vario_model == "Gaussian"){
    gamma_vec = gauss_cov(fitted_par_vario,distances)
    c0 = gauss_cov(fitted_par_vario, 0)

  }
  else if (vario_model == "Exponential"){
    gamma_vec = exp_cov(fitted_par_vario,distances)
    c0 = exp_cov(fitted_par_vario, 0)

  }
  else if (vario_model == "Spherical"){
    gamma_vec = sph_cov(fitted_par_vario,distances)
    c0 = sph_cov(fitted_par_vario, 0)

  }
  else {
    stop("Model not available")
  }
  gamma_matrix = matrix(0,N,N)
  gamma_matrix[lower.tri(gamma_matrix)] = gamma_vec
  gamma_matrix = t(gamma_matrix) + gamma_matrix
  diag(gamma_matrix) = c0*rep(1,N)

  return (gamma_matrix)
}

### Arguments:
# - h: single distance at which the gaussian variogram must be computed
# - par: parameters of the variogram (tau2, sigma2, a)
### Value: value of the gaussian variogram (with parameters par) in h
gauss_vario_univ = function(h, par){
  tau2 = par[1]
  sigma2 = par[2]
  a = par[3]
  if (h == 0) {
    gauss_variogram = 0
  }
  else {
    gauss_variogram = (tau2)+(sigma2)*(1-exp(-(h^2)/(a^2))) # Valori del variogramma gaussiano fittato usando i parametri di par
  }
  return (gauss_variogram)
}

### Arguments:
# - par: parameters of the variogram (tau2, sigma2, a)
# - supp: vector of points (distances) at which the gaussian variogram must be computed
### Value: vector containing the gaussian variogram (with parameters par) valued at h
gauss_vario = function(par, supp){
  supp = as.list(supp)
  gauss_variogram = map_dbl(supp, gauss_vario_univ, par = par)
  return(gauss_variogram)
}

### Arguments:
# - par: parameters of the variogram
# - h: vector of distances at which the gaussian covariogram must be computed
### Value: vector of values of of the gaussian variogram (with parameters par) in h
gauss_cov<-function(par,h){
  c0 = par[2] + par[1]
  gauss_covario = c0 - gauss_vario (par, h)
  return (gauss_covario)
}

### Arguments:
# - par: parameters of the variogram (tau2, sigma2, a)
# - h: points at which the error must be computed must be computed
# - emp_vario_values: emprical variogram values
### Value: sum of square of the differences between the gaussian variogram and the empirical variogram
gauss_err<-function(par,h,emp_vario_values){
  gauss_variogram = gauss_vario(par,h)
  err = sum((gauss_variogram - emp_vario_values)^2) # errore = somma dei quadrati degli errori
  return (err)
}


### Arguments:
# - h: single distance at which the exponetial variogram
# - par: parameters of the variogram (tau2, sigma2, a)
### Value: value of the exponetial variogram (with parameters par) in h
exp_vario_univ = function(h, par){
  tau2 = par[1]
  sigma2 = par[2]
  a = par[3]
  if (h == 0) {
    exp_variogram = 0
  }
  else {
    exp_variogram = (tau2)+(sigma2)*(1-exp(-(h)/a))
  }

  return(exp_variogram)
}

### Arguments:
# - par: parameters of the variogram (tau2, sigma2, a)
# - supp: vector of points (distances) at which the exponential variogram must be computed
### Value: vector containing the exponential variogram (with parameters par) valued at h
exp_vario = function(par, supp){
  supp = as.list(supp)
  exp_variogram = map_dbl(supp, exp_vario_univ, par = par)
  return(exp_variogram)
}

### Arguments:
# - par: parameters of the variogram
# - h: vector of distances at which the exponential covariogram must be computed
### Value: vector of values of of the gaussian exponential (with parameters par) in h
exp_cov<-function(par,h){
  c0 = par[2] + par[1]
  exp_covario = c0 - exp_vario (par, h)

  return (exp_covario)
}

### Arguments:
# - par: parameters of the variogram (tau2, sigma2, a)
# - h: points at which the error must be computed
# - emp_vario_values: emprical variogram values
### Value: sum of square of the differences between the exponential variogram and the empirical variogram
exp_err<-function(par,h,emp_vario_values){
  exp_variogram = exp_vario(par,h)
  err = sum((exp_variogram-emp_vario_values)^2)
  return (err)
}

### Arguments:
# - h: single distance at which the spherical variogram must be computed
# - par: parameters of the variogram (tau2, sigma2, a)
### Value: value of the spherical variogram (with parameters par) in h
sph_vario_univ = function(h, par) {
  tau2 = par[1]
  sigma2 = par[2]
  a = par[3]

  if(h == 0){
    sph_variogram = 0
  }
  else if(h >= a){
    sph_variogram = tau2 + sigma2
  }
  else{
    sph_variogram = tau2 + sigma2*(3/2*h/a - 1/2*(h/a)^3)
  }
  return (sph_variogram)
}

### Arguments:
# - par: parameters of the variogram (tau2, sigma2, a)
# - supp: vector of points (distances) at which the spherical variogram must be computed
### Value: vector containing the spherical variogram (with parameters par) valued at h
sph_vario = function(par, supp){
  supp = as.list(supp)
  sph_variogram = map_dbl(supp, sph_vario_univ, par = par)
  return(sph_variogram)
}

### Arguments:
# - par: parameters of the variogram
# - h: vector of distances at which the spherical covariogram must be computed
### Value: vector of values of of the spherical (with parameters par) in h
sph_cov = function(par,h){
  c0 = par[2] + par[1]
  sph_covario = c0 - sph_vario (par, h)
  return (sph_covario)
}

### Arguments:
# - par: parameters of the variogram (tau2, sigma2, a)
# - h: points at which the error must be computed
# - emp_vario_values: emprical variogram values
### Value: sum of square of the differences between the spherical variogram and the empirical variogram
sph_err = function(par,h,emp_vario_values){
  sph_variogram = sph_vario(par, h)
  err = sum((sph_variogram-emp_vario_values)^2)
  return (err)
}

### Arguments:
# - empirical_variogram: vector of values of the empirical variogram
# - fitted_par_vario: list of the fitted parameters of the selected variogram
# - model: model of the variogram
# - distance: type of distance used (Eucldist or Geodist)
### Value: plot of the empirical variogram with the corresponding fitted variogram
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


#############################################################################################
################################# kriging ###################################################
#############################################################################################

### Arguments: ***NEW*** Function modified on 15th march
# - new_coords: coordinates at which the prediction must be computed
# - X_new: dataframe/matrix of covariates in the new points (possibly NULL)
# - model_ts: model used in the tangent space (coord1 or coord2 or additive)
### Value: design matrix
create_design_matrix = function(new_coords, X_new, model_ts){
  n = dim(X_new)[1]
  if(model_ts=="intercept") # *NEW* Include stationary
  {
    design_matrix = as.data.frame(cbind(rep(1,n), X_new))
  }
  else if(model_ts=="coord1")
  {
    design_matrix = as.data.frame(cbind(rep(1,n),new_coords[,1], X_new))
  }
  else if(model_ts=="coord2")
  {
    design_matrix = as.data.frame(cbind(rep(1,n),new_coords[,2], X_new))
  }
  else if(model_ts=="additive")
  {
    design_matrix = as.data.frame(cbind(rep(1,n),new_coords, X_new))
  }
  else{
    stop("Model not available")
  }
  return(as.matrix(design_matrix))
}


### Arguments:
# - coords: coords used to build the model
# - new_coords: coordinates at which the prediction must be computed
# - vario_model: model selected for the variogram ("Gaussian" or "Exponential" or "Spherical")
# - par: list of the estimated parameters of the variogram
### Value: vector c used to computer the vector of weights lambda
compute_ci = function(coords, new_coords, vario_model, par, distance){
  if (distance == "Geodist") {
    distances_vec = apply (coords, 1, Geodist, c2 = new_coords)
  }
  else if (distance == "Eucldist") {
    distances_vec = apply (coords, 1, Eucldist, c2 = new_coords)
  }
  else {
    stop ("Distance not available")
  }

  if (vario_model == "Gaussian") {
    c = gauss_cov(par, distances_vec)
  }
  else if (vario_model == "Exponential") {
    c = exp_cov (par, distances_vec)
  }
  else if (vario_model =="Spherical") {
    c = sph_cov (par, distances_vec)
  }
  else {
    stop("Model not available")
  }
  return (c)
}
