
################################ RDD #####################################

create.rdd = function(K, method.rdd='Voronoi', data_coords,
                      # border.length, spdist,
                      graph.distance,
                      # mesh,
                      nk_min,grid,
                      # is.observed, graph.distance.complete, assign.matrix,
                      data.grid.distance, suppressMes=T)
{
  if(method.rdd!='Voronoi')
  {
    print('Error: only Voronoi RDD implemented so far')
    return(-1)
  }

  N_samples = dim(data_coords)[1]

  card_min=0
  first=T
  while(card_min < nk_min)  # If at least a neighborhood has too few points, repeat the RDD
  {
    if(first==F & suppressMes==F)
      print("Repeating sampling")
    ind = sort(sample(1:N_samples, size=K)) # Indices of the data matrix corresponding to the centers
    centers = cbind(data_coords$x[ind],data_coords$y[ind],ind) # Extract centers

    # Assign data to centers
    assign = rep(0, N_samples)
    for(i in 1:N_samples){
      # For each datum, we compute the distance between the datum and each center
      # d = datapoints2centers.distance(x=cbind(data$x[i],data$y[i],i),
      #                                 centers = centers, method = spdist,
      #                                 graph.distance = graph.distance)
      d = graph.distance[ind,i]
      assign[i] = assign2center(d)
    }
    first=F
    card_min=min(table(assign))
  }

  if(is.null(grid))
  {
    list.ret = list(assign=assign, centers=centers)
    return(list.ret)
  }
  ngrid=dim(grid)[1]

  ### ------------- 1.2 Assign grid points to centers

  # data.grid.distance = loccoords(coords = cbind(data$x,data$y), locations = grid)
  # # coords = nsub*2 matrix with the coordinates of the observed data
  # # locations = ngrid*2 matrix with the coordinates od the prediction locations
  # # data.grid.distance is a nsub*ngrid matrix and the (i,j) element is the distance beteween the i-th observed point and the j-th prediction location
  # # if data.grid.distance(i,j)==0 -----> the j-th prediction point is an already observed value
  #
  # already.observed = which(data.grid.distance == 0, arr.ind = TRUE)
  # is.observed = rep(0,ngrid)
  # # If the i-th grid point is already observed ---> is.observed[i] contains the index of the point in the data matrix
  # # If the i-th grid point is unobserved       ---> is.observed[i] = 0
  # is.observed[already.observed[,2]]=already.observed[,1]
  #
  # Besides the coordinates of the grid point we need two other informations:
  # The index of the grid point in the grid matrix and
  # The index of the grid point in the data matrix if it's an already observed point

  # grid.new = data.frame(cbind(grid, 1:ngrid, is.observed))
  # colnames(grid.new)=c('x','y','grid.position','is.observed')

  # assigng = factor(rep(0, ngrid), levels=c(0,1:K))
  assigng = rep(0, ngrid)  # *NEW*

  # graph.distance.grid.centers : K*ngrid matrix (the element (k,i) is the distance on the graph between the k-th center and the i-th grid point)
  # (This matrix will be used to compute the kernel value if ker.width.intrinsic > 0)

  # graph.distance.grid.centers = matrix(NA, K, ngrid) # kx411
  graph.distance.grid.centers = data.grid.distance[ind,]
  # triangles = mesh$T

  for(i in 1:ngrid){
    # d = gridpoints2centers.distance(x=grid.new[i,], centers = centers,
    #                                 method = spdist,
    #                                 graph.distance = graph.distance.complete,
    #                                 assign.matrix = assign.matrix,
    #                                 data.grid.distance = data.grid.distance,
    #                                 triangles = triangles)
    # graph.distance.grid.centers[,i]=d
    d = data.grid.distance[ind,i]
    assigng[i]=assign2center(d)
  }

  # nk.assigng = rep(0, K)
  # for(k in 1:K)
  #   nk.assigng[k] = length(which(assigng==k))

  gridk = list()
  for(k in 1:K)
    gridk[[k]]=grid[assigng==k,]

  list.ret = list(assign=assign, centers=centers,
                  assigng=assigng, gridk=gridk,
                  graph.distance.grid.centers=graph.distance.grid.centers)
  return(list.ret)

}

# distance.on.graph.constrained = function(data, mesh, distance)
# {
#   # INPUT
#   # data = data matrix (nsub*3 dim)
#   # mesh = mesh obtained with Delaunay triangulation (with triangulate)
#   # distance = distance used to compute the edges' weights. "Eucldist" or "Geodist"
#   # OUTPUT
#   # graph.distance = nsub*nsub symmetric matrix (element (i,j) is the shortest length path between i and j on the graph defined by the mesh)
#
#   # Edges list : edge.num * 2 matrix (each row contains the indeces of the nodes i and j defining the edge e = {ij})
#   edge.list = mesh$E
#   edge.num = dim(edge.list)[1]
#
#   # Construct the graph (undirected)
#   graph = graph_from_edgelist(edge.list, directed = FALSE)
#
#   # Given two linked nodes i and j the weigth of edge (i,j) is:
#   # w(i,j) = euclidian distance (i, j)
#   edge.weights = rep(0, edge.num)
#
#   for(l in 1:edge.num){
#     node.i = edge.list[l,1]
#     node.j = edge.list[l,2]
#     point.i = mesh$P[node.i,1:2]
#     point.j = mesh$P[node.j,1:2]
#     if(distance == "Eucldist") edge.weights[l] = dist(rbind(point.i,point.j), method = 'euclidean')
#     else if (distance == "Geodist") edge.weights[l] = Geodist(point.i,point.j)
#   }
#
#   graph.distance = distances(graph, v = V(graph), to = V(graph),
#                              mode = c("all"), weights = edge.weights,
#                              algorithm = c("automatic"))
#   return(graph.distance)
# }

# grid2triangle.constrained = function(grid.matrix, data.matrix, mesh)
# {
#   # INPUT:
#   # grid : ngrid*2 matrix (the i-th row contains the coordinates of the i-th grid point)
#   # data : nsub*2 matrix (the i-th row contains the coordinates of the i-th data point)
#   # mesh : mesh obtaines with Delauny Triangulation (using create.MESH.2D function)
#   # OUTPUT:
#   # assign.matrix : ngrid*num.triangles (the (i,j) element = TRUE if the point grid[i,] is inside the j-th triangle, FALSE otherwise)
#   # no.assg.complete :indices of grid points that the function does not assign to any triangle
#   triangles = mesh$T
#   num.triangles = dim(triangles)[1]
#
#   ngrid = dim(grid.matrix)[1]
#   assign.matrix = matrix(NA, ngrid, num.triangles)
#
#   for(j in 1:num.triangles){
#     bnd.triangles = triangles[j,]
#     assign.matrix[,j] = in.out(bnd = cbind(mesh$P[bnd.triangles,1],mesh$P[bnd.triangles,2]),
#                                x = cbind(as.numeric(grid.matrix[,1]),as.numeric(grid.matrix[,2]))) # togliere punti veri?
#   }
#
#   check.complete = rep(NA, ngrid)
#   for(i in 1:ngrid){
#     check.complete[i] = length(which(assign.matrix[i,]==TRUE))
#   }
#   no.assg.complete = which(check.complete == 0)
#
#   output = list()
#   output[[1]]=assign.matrix
#   output[[2]]=no.assg.complete
#
#   return(output)
# }

assign2center=function(distance.vector)
{
  # INPUT :
  # distance.vector : K-dimensional vector (distance.vector[k] = distance between a fixed point and the k-th center)
  # OUTPUT :
  # the index of the center to which the point is assigned
  assigned.center = which.min(distance.vector)
  # if(!is.empty(assigned.center)){
  if(length(assigned.center)!=0){
    return(assigned.center)}
  # If the distance.vector is a vector with only NA's, it means that it's a grid point which cannot be assigned to any center
  # if(is.empty(assigned.center)){
  if(length(assigned.center)==0){
    return(0)
  }
}

# datapoints2centers.distance = function(x,centers, method, graph.distance=NULL)
# {
#   # This function compute the distance between an observed point and each center
#   # INPUT :
#   # x : it's a 3-dimension vector (x[1] = x-coordinate,
#   #                                x[2] = y-coordinate,
#   #                                x[3] = it indicates the index of the row in the data matrix corresponding the observed point)
#   #
#   # centers : K*3 matrix (centers[k,1:2] = coordinates of the k-th center,
#   #                       centers[k, 3] = indicates the row index in the data matrix of the k-th center)
#   #
#   # method : 'euclidean' or 'graph.distance'
#   #
#   # graph.distance : nsub*nsub matrix (graph.distance(i,j) = shortest path between i and j on the graph defined by Delaunay's Traingulation)
#   x = as.numeric(x)
#   if(method=='euclidean'){
#     d = as.matrix(apply(centers[,1:2],1, Eucldist, x[1:2]))
#   }
#
#   if(method=='graph.distance'){
#     ncenters = dim(centers)[1]
#     d = rep(NA, ncenters)
#     for(i in 1:ncenters){
#       d[i] = graph.distance[x[3],centers[i,3]]}
#   }
#   return(d)
# }

# gridpoints2centers.distance = function(x, centers, method, graph.distance=NULL, assign.matrix = NULL, data.grid.distance = NULL, triangles = NULL)
# {
#   # This function compute the distance between a grid point and each center
#   # INPUT :
#   # x : it's a 4-dimension vector (x[1] = x-coordinate,
#   #                                x[2] = y-coordinate,
#   #                                x[3] = index of the row corresponding to the grid point in the grid matrix,
#   #                                x[4] = if the point is an UNOBSERVED locations is 0, otherwise it indicates the
#   #                                       index of the row of the data matrix corresponding to x[1:2] coodinates)
#   #
#   # centers : K*3 matrix (centers[k,1:2] = coordinates of the k-th center,
#   #                       centers[k, 3] = indicates the row index in the data matrix of the k-th center)
#   #
#   # method : 'euclidean' or 'graph.distance'
#   #
#   # graph.distance : nsub*nsub matrix (graph.distance(i,j) = shortest path between i and j on the graph defined by Delaunay's Traingulation)
#   #
#   # assign.matrix : ngrid*nsub matrix (assign.matrix[i,j] = 1 if the grid point i is assigned to the j triangle, 0 otherwise)
#   #
#   # data.grid.distance : nsub*ngrid matrix (data.grid.distance(i,j) is the euclidian distance between the observed i point and the grid point j)
#   #
#   # triangles : ntriangles*3 matrix (each row contains the 3 indices of the vertices defining the triangle)
#   # OUTPUT :
#   # d : K-dimensional vector (d[i] is the distance between the point x and the i-th center)
#
#   x = as.numeric(x)
#   ncenters = dim(centers)[1]
#
#   if(method=='euclidean'){
#     d = as.matrix(apply(centers[,1:2],1, Eucldist, x[1:2]))
#   }
#
#   if(method=='graph.distance'){
#     # If x[4] > 0 it means that the grid point where we want to do a prediction is an already observed location.
#     # Therefore it's a vertex of the graph, the distance between the point and each center is computed directly as
#     # shortest distance path on the constructed graph (is already computed and stored in graph.distance)
#     if(x[4] > 0){
#       d = rep(NA, ncenters)
#       for(i in 1:ncenters){
#         d[i] = graph.distance[x[4],centers[i,3]]}
#     }
#     # If x[4] == 0 it means that the grid points is an unobserved location
#     if(x[4] == 0){
#       # Index of the row in the triangles matrix corresponding to the triangle/s to which the point has been assigned
#       ind.triangle = which(assign.matrix[x[3],]==1)
#       # It's assigned to one and only one triangle if it's inside the triangle (not on the border)
#       # or it could be assigned to two adjacent triangles if it's on the common border of the triangles (it never happens)
#       if(length(ind.triangle) == 1){
#         # d.Pj = euclidian distance between the point x and the vertex j of the triangle to which x has been assigned
#         d.P1 = data.grid.distance[triangles[ind.triangle,][1],x[3]]
#         d.P2 = data.grid.distance[triangles[ind.triangle,][2],x[3]]
#         d.P3 = data.grid.distance[triangles[ind.triangle,][3],x[3]]
#         d = rep(NA, ncenters)
#         for(i in 1:ncenters){
#           d[i] = min(d.P1 + graph.distance[triangles[ind.triangle,][1],centers[i,3]],
#                      d.P2 + graph.distance[triangles[ind.triangle,][2],centers[i,3]],
#                      d.P3 + graph.distance[triangles[ind.triangle,][3],centers[i,3]])
#         }
#       }
#       if(length(ind.triangle) == 2){
#         # There are 4-vertices to evaluate d.Pj
#         ind.vertex = union(triangles[ind.triangle[1],],triangles[ind.triangle[2],])
#         d.P1 = data.grid.distance[ind.vertex[1],x[3]]
#         d.P2 = data.grid.distance[ind.vertex[2],x[3]]
#         d.P3 = data.grid.distance[ind.vertex[3],x[3]]
#         d.P4 = data.grid.distance[ind.vertex[4],x[3]]
#         d = rep(NA, ncenters)
#         for(i in 1:ncenters){
#           d[i] = min(d.P1 + graph.distance[ind.vertex[1],centers[i,3]],
#                      d.P2 + graph.distance[ind.vertex[2],centers[i,3]],
#                      d.P3 + graph.distance[ind.vertex[3],centers[i,3]],
#                      d.P4 + graph.distance[ind.vertex[4],centers[i,3]])
#         }
#       }
#       # If the grid point is unobserved it cannot happen that length(ind.triangle)>2, because otherwise it's a vertex of a triangle and
#       # therefore it's a data point
#
#       # There are some grid points (for example grid points outside the convex.hull) which are not assigned to any triangle.
#       # These grid points are not used to perform prediction and we don't prefict the response variable in these points
#       if(length(ind.triangle)==0)
#       {d=rep(NA, ncenters)}
#     }
#   }
#   return(d)
# }

kerfn= function(newdata,center,ker.type='Gau',param) # dist, distance.matrix = NULL
{
  # This function compute the value of the kernel for a point given a reference center
  # Input:
  # newdata = coordinates of the point. Used only if distance='euclidean'
  # center = coordinates of the reference center
  # dist = method to compute distances ('euclidean' or 'graph.distance')
  #
  # ker.type = type of kernel (only Gaussian kernel implemented so far)
  # param = parameters that define the kernel
  # distance.matrix = 1. if the point is an observed location (data point)
  #                      distance.matrix = graph.distance nsub*nsub
  #                                        (the (i,j) element is the length of the shortest path between the observed points i and j)
  #                   2. if the point is an unobserved location (grid point)
  #                      distance.matrix = graph.distance.grid.centers K*ngrid
  #                                        (the (k,l) element is the distance on the graph between the k-th center and the l-th grid point)
  #                  Used only if distance='graph.distance'
  # Output:
  # value of the kernel function
  center = as.numeric(center)
  d = as.matrix(apply(newdata,1, Eucldist, center[1:2]))

  # if(dist == 'euclidean'){
  #   d = as.matrix(apply(newdata,1, Eucldist, center[1:2]))
  # }
  # if(dist == 'graph.distance'){
  #   # d = k-th row of the distance matrix
  #   d = distance.matrix[as.integer(center[3]),]
  # }
  #
  if(ker.type!='Gau')
  {
    print('Error: only Gaussian kernel implemented so far (*Gau*)')
  }
  eps=param[1]
  return(exp(-1/(2*eps^2)*(d^2)))
}


########################### MATRIX MANIPULATION ##########################



## Arguments:
# - ll: list of n pxp symmetric matrices
## Value: nx (p*(p+1)/2) matrix. The i-th row corresponds to matrix_to_vec(i-th matrix of the array)
list_to_matrix = function(ll){
  # tmp = map(ll, matrix_to_vec)
  # return(map_df(ll,rbind))
  vec_list = map(ll, matrix_to_vec)
  return(vec_list %>% reduce(rbind))
}

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

######################## GEOGRAPHICAL DISTANCE ###########################

# Haversine formula
Geodist=function(c1,c2){
  R=6371
  r=(pi/2)/90
  tmp = sqrt( sin( (c2[1]-c1[1])*r/2 )^2 + cos(c1[1]*r) * cos(c2[1]*r) * sin( (c2[2] - c1[2])*r/2) ^2 )
  if( tmp > 1){
    tmp = 1
  }
  d=as.double(2*R*asin(tmp))
  return (d)

}

Eucldist = function (c1, c2) {
  return (dist(rbind(c1,c2)))
}
