expMat = function(A)
{
  eig = eigen(A)
  expA = eig$vectors%*%diag(exp(eig$values))%*%t(eig$vectors)
  return(expA)
}

logMat = function(A)
{
  eig = eigen(A)
  logA = eig$vectors%*%diag(log(eig$values))%*%t(eig$vectors)
  return(logA)
}


sqrtMat = function(A)
{
  eig = eigen(A)
  sqrtA = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
  return(sqrtA)
}


### Arguments: 
# - sdp
# - Sigma
# - metric_manifold
### Value: matrix_to_vec(logarithmic_map(sdp, Sigma, metric_manifold))
logarithmic_map_vec = function (sdp, Sigma, metric_manifold) {
  return (matrix_to_vec (logarithmic_map(sdp, Sigma, metric_manifold)))
}

### Arguments: 
# - sdp: pxp matrix or matrix_to_vec(pxp matrix) (in PD(p))
# - Sigma: Tangent point
# - metric_manifold: metric to be used on the manifold
### Value: The projection of sdp on the tangent space in Sigma
logarithmic_map = function(sdp, Sigma, metric_manifold){
  if(!is.matrix(sdp) || (is.matrix(sdp) & min(dim(sdp)) == 1)){
    sdp = vec_to_matrix(sdp) 
  }
  if(!is.matrix(Sigma) || (is.matrix(Sigma) & min(dim(Sigma)) == 1)){
    Sigma = vec_to_matrix(Sigma)
  }
  
  if(metric_manifold == "Frobenius")
  {
    sqrt_Sigma = sqrtMat(Sigma)
    inv_sqrt_Sigma = solve(sqrt_Sigma)
    
    return (sqrt_Sigma%*%logMat(inv_sqrt_Sigma%*%sdp%*%inv_sqrt_Sigma)%*%sqrt_Sigma)
  }
  else if (metric_manifold == "Log_euclidean")
  {
    return (logMat(sdp)-logMat(Sigma))
  }
  else if (metric_manifold == "Square_root")
  {
    return (sqrtMat(sdp)-sqrtMat(Sigma))
  }
  else {
    stop("Manifold metric not available")
  }
}

### Arguments: 
# - symm: pxp matrix or matrix_to_vec(pxp matrix) (in Sym(p))
# - Sigma: Tangent point
# - metric_manifold: metric to be used on the manifold
### Value: The projection of symm on the manifold
exponential_map = function(symm, Sigma, metric_manifold){
  if(!is.matrix(symm) || (is.matrix(symm) & min(dim(symm)) == 1)){
    symm = vec_to_matrix(symm) 
  }
  if(!is.matrix(Sigma) || (is.matrix(Sigma) & min(dim(Sigma)) == 1)){
    Sigma = vec_to_matrix(Sigma)
  }
  if(metric_manifold == "Frobenius")
  {
    sqrt_Sigma = sqrtMat(Sigma)
    inv_sqrt_Sigma = solve(sqrt_Sigma)
    return (sqrt_Sigma%*%expMat(inv_sqrt_Sigma%*%symm%*%inv_sqrt_Sigma)%*%sqrt_Sigma)
  }
  else if (metric_manifold == "Square_root")
  {
    tmp = sqrtMat(Sigma)+symm
    return (t(tmp)%*%tmp)
  }
  else if (metric_manifold == "Log_euclidean")
  {
    tmp = logMat(Sigma)+symm
    return (t(tmp)%*%tmp)
  }
  else {
    stop("Manifold metric not available")
  }
}
