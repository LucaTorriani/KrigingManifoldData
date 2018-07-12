### Arguments: 
# - dataframe
### Value: covariance matrix of dataframe
create_covariance_matrix = function (dataframe){
  return (cov(dataframe))
}

### Arguments:
# - D: nxn matrix (in Sym(p))
### Value: n*(n+1)/2 vector with the elements of the upper triangual part of D (diagonal included), saved row-wise
matrix_to_vec = function(D){
  d = rep(0,dim(D)[1]*(dim(D)[1]+1)/2) 
  k = 1
  for (i in 1:dim(D)[1]){
    for (j in i:dim(D)[1]){
      d[k] = D[i,j]
      k = k+1
    }
  }
  return  (d)
}

### Arguments:
# - d: n*(n+1)/2 vector 
### Value: nxn matrix (in Sym(p)). The symmetric matrix obtained putting, row-wise, in the upper triangular part (diagonal included) the elements of d
vec_to_matrix = function(d){
  n = (sqrt(8*length(d)+1)-1)/2
  D = matrix(0,n,n)
  k = 1
  for (i in 1:n){
    for (j in i:n){
      D[i,j] = as.numeric(d[k])
      D[j,i] = as.numeric(d[k])
      k = k+1
    }
  }
  return (D)
}

### Arguments:
# - arr: array of n pxp matrices (nxpxp or pxpxn)
### Value: list of
# - dim: n
# - pos: position of n in dim(arr) (i.e. 1 if arr is nxpxp, 3 if arr is pxpxn)
find_dim_array = function(arr){
  if(dim(arr)[1]==dim(arr)[2]){
    return(list(dim = dim(arr)[3], pos = 3))
  }
  else{
    return(list(dim = dim(arr)[1], pos = 1))
  }
}

### Arguments:
# - arr: array of n pxp symmetric matrices (nxpxp or pxpxn)
### Value: nx (p*(p+1)/2) matrix. The i-th row corresponds to matrix_to_vec(i-th matrix of the array)
matrixArray_to_matrix = function(arr){
  A = NULL
  dim = find_dim_array(arr)$dim
  pos = find_dim_array(arr)$pos
  for(i in 1:dim){
    if(pos == 1){
      A = rbind(A,matrix_to_vec(arr[i,,]))
    }
    else{
      A = rbind(A,matrix_to_vec(arr[,,i]))
      
    }
  }
  return(A)
}

### Arguments:
# - matr: nx(p*(p+1)/2) matrix
# - p: dimension of the matrices 
### Value: array of n pxp symmetric matrices (pxpxn). The i-th matrix is vec_to_matrix(i-th row of matr) 
matrix_to_matrixArray = function(matr, p){
  arr = array(dim = c(p,p, dim(matr)[1]))
  for(i in 1:dim(matr)[1]){
    arr[,,i] = vec_to_matrix(matr[i,])
  }
  return(arr)
}

