### Arguments
# - c1: (Latitude, Longitude) of point 1
# - c2: (Latitude, Longitude) of point 2
### Value: Geodetic distance between point 1 and point 2
Geodist=function(c1,c2){
  R=6371
  r=(pi/2)/90
  tmp = sin(c1[1]*r)*sin(c2[1]*r)+cos(c1[1]*r)*cos(c2[1]*r)*cos(abs(-c2[2]*r+c1[2]*r))
  if( abs(tmp-1) < 2*.Machine$double.eps){
    return(0)
  }
  else if(abs(tmp+1) < 2*.Machine$double.eps){
    return(pi*R)
  }
  else{
    d=as.double(R*acos(tmp))
    return (d)   
  }
}

### Arguments
# - coords: nx2 matrix. The i-th row indicates (Latitute, Longitude) of point i
### Value: vector d of geodietic distances among the n points. The distance between point i and point j corresponds to d[n*(i-1) - i*(i-1)/2 + j-i]
compute_Geodist = function(coords){
  d = as.vector(proxy::dist(coords, Geodist))
  return(d)
}

### Arguments
# - c1: (x, y) of point 1
# - c2: (x, y) of point 2
### Value: Euclidean distance between point 1 and point 2
Eucldist = function (c1, c2) {
  return (dist(rbind(c1,c2)))
}

### Arguments
# - coords: nx2 matrix. The i-th row indicates (x, y) of point i
### Value: vector d of euclidean distances among the n points. The distance between point i and point j corresponds to d[n*(i-1) - i*(i-1)/2 + j-i]
compute_Eucldist = function(coords){
  d = dist(coords)
  return (d)
}

### Arguments
# - A: matrix whose norm has to be evaluated
# - Sigma: Tangent point
### Value: Frobenius-scaled norm of A
norm_mat_frob_scaled = function(A, Sigma) {
  invsigma = solve(Sigma)
  return(sqrt(sum(diag(invsigma%*%A%*%invsigma%*%A))))
}

### Arguments
# - A: matrix whose norm has to be evaluated
# - Sigma: Tangent point
### Value: Frobenius norm of A
norm_mat_frob = function(A) {
  return(sqrt(tplane_product(A,A)))
}

### Arguments:
# - P1: matrix_1
# - P2: matrix_2
### Value: Frobenius product between P1 and P2
tplane_product = function(P1, P2){
    return (tr(t(P1)%*%P2))
}

### Arguments:
# - P1: matrix_1 (in Sym(p))
# - P2: matrix_2 (in Sym(p))
# - metric_ts: distance to consider (in the tangent space)
# - Sigma (default NULL): Sigma to use in case of metric_ts = "Scaled_Frobenius"
### Value: Distance between P1 and P2, according to the metric indicated by metric_ts
tplane_distance = function(P1, P2, metric_ts, Sigma = NULL){
  if (metric_ts == "Frobenius") {
    return (sqrt(tplane_product(P1-P2, P1-P2)))
  }
  else if (metric_ts == "Scaled_Frobenius") {
    if (is.null(Sigma)) {
      stop ("To use the scaled Frobenius metric Sigma must be provided")
    }
    return (norm_mat_frob_scaled(P1-P2, Sigma))
  }
  else {
    stop ("Tangent space metric not available")
  }
  return (0)
}

### Arguments:
# - P1: matrix_to_vec(matrix_1) (in PD(p))
# - P2: matrix_to_vec(matrix_2) (in PD(p))
# - metric_: distance to use on the manifold
### Value: Distance between P1 and P2, according to the metric indicated by metric_manifold
manifold_distance = function(P1,P2,metric_manifold = "Frobenius"){ 
  P1 = vec_to_matrix(P1)
  P2 = vec_to_matrix(P2)
  if(metric_manifold == "Frobenius")
  {
    eigenvalues = eigen(solve(P1)%*%P2)$values
    return (sqrt(sum(log(eigenvalues)^2)))
  }
  else if (metric_manifold == "Log_euclidean")
  {
    return (norm(logMat(P1)-logMat(P2),"F"))
  }
  else if(metric_manifold == "Square_root")
  {
    return (norm(sqrtMat(P1)-sqrtMat(P2),"F"))
  }
  else
  {
    stop("Distance chosen is not available") 
  }
}

