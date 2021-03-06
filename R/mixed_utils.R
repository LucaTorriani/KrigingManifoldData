#' Return a given element of a list
#' @param lista A list
#' @param i The index of the element to extract
#' @return It returns the i-th element of \code{lista}
#' @useDynLib Manifoldgstat
return_ith_list_element = function(lista,i) {
  return(lista[[i]])
}

# RDD_OOK_aggr_man_mixed=function(fOKBV,weights_intrinsic, ker.width.intrinsic, N_samples=NULL) # p, num.signif.entries=1
# {
#   # This functions aggregates the results of bootstrap in OKBV
#   # Input:
#   # fOKBV = result of OKBV, as returned by the first element of output of OKBV()
#   # weight_intrinsic = weights to aggregate the results
#   # ker.width.intrinsic
#   # N_samples = sample size
#   # p = dimension of the matrices on the manifold
#   # num_signif_entries = (p*(p+1))/2
#   # Output:
#   # Average prediction (on the grid used to produce fOKBV)
#   #
#
#   # if(num.signif.entries==1)
#   # {
#   #   print('Please provide the grid length *n*')
#   #   return(-1)
#   # }
#
#   # fOKBV1=simplify2array(fOKBV)
#   # B = dim(fOKBV2)[2]
#
#   # if(dim(fOKBV1)[1]==1)
#   if(length(fOKBV[[1]])==1)
#   {
#     print('Please check source code for grid dimension == 1!')
#     return(-1)
#   }
#
#   weight = NULL
#   ngrid=length(fOKBV[[1]])
#   # fpred.ave = array(NA, dim=c(ngrid, num.signif.entries))
#   fpred.ave = list()
#
#   # id.na = which(is.na(fOKBV1[,1,]))
#
#   if(ker.width.intrinsic != 0)
#   {
#     print("Using weights for aggregation")
#     W.b = do.call(cbind,weights_intrinsic)
#     W.tot=apply(W.b,1,sum, na.rm=TRUE)
#     weight=W.b/W.tot # It divides each row of W.b by the corresponding element of the vector W.tot
#   }
#
#
#   # if(length(id.na)>0)
#   # {
#   #   for(i in (1:ngrid)[-id.na]) {
#   #     # fpred.ave[i,]= matrix_to_vec(intrinsic_mean(data=matrix_to_matrixArray(t(fOKBV1[i,,]),p = p),
#   #     #                                             metric_manifold = metric_manifold, metric_ts = metric_ts, weight_intrinsic = weight[i,]))
#   #     fpred.ave[[i]]= intrinsic_mean(data= map(fOKBV,return_ith_list_element, k=i),
#   #                                    metric_manifold = metric_manifold, metric_ts = metric_ts, weight_intrinsic = weight[i,])
#   #   }
#   #
#   #   return(fpred.ave)
#   #
#   # }
#   # if(length(id.na)==0)
#   # {
#   #   for(i in (1:ngrid)){
#   #     # fpred.ave[i,]= matrix_to_vec(intrinsic_mean(data=matrix_to_matrixArray(t(fOKBV1[i,,]),p = p),
#   #     #                                             metric_manifold = metric_manifold, metric_ts = metric_ts, weight_intrinsic = weight[i,]))
#   #     fpred.ave[[i]]= intrinsic_mean(data= map(fOKBV,return_ith_list_element, k=i),
#   #                                    metric_manifold = metric_manifold, metric_ts = metric_ts, weight_intrinsic = weight[i,])
#   #
#   #   }
#   #   return(fpred.ave)
#   #
#   # }
#   for(i in (1:ngrid)){
#     # fpred.ave[i,]= matrix_to_vec(intrinsic_mean(data=matrix_to_matrixArray(t(fOKBV1[i,,]),p = p),
#     #                                             metric_manifold = metric_manifold, metric_ts = metric_ts, weight_intrinsic = weight[i,]))
#     fpred.ave[[i]]= intrinsic_mean(data= map(fOKBV,return_ith_list_element, i=i),
#                                    metric_manifold = metric_manifold, metric_ts = metric_ts, weight_intrinsic = weight[i,])
#   }
#
#   return(fpred.ave)
# }


#' Main routine for mixed_RDD
#' @param data_coords \code{N*2} or \code{N*3} matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to
#' be provided in signed decimal degrees
#' @param data_val array [\code{p,p,N}] of \code{N} symmetric positive definite matrices of dimension \code{p*p}
#' @param K number of cells the domain is subdivided in
#' @param grid prediction grid, i.e. \code{M*2} or \code{M*3} matrix of coordinates where to predict
#' @param nk_min minimum number of observations within a cell
#' @param B number of \emph{divide} iterations to perform
#' @param suppressMes \{\code{TRUE}, \code{FALSE}\} controls the level of interaction and warnings given
#' @param ker.width.intrinsic parameter controlling the width of the Gaussian kernel for the computation of the local mean (if 0,
#' a "step kernel" is used, giving weight 1 to all the data within the cell and 0 to those outside of it)
#' @param graph.distance.complete \code{N*N} distance matrix (the [i,j] element is the length of the shortest path between points i and j)
#' @param data.grid.distance \code{N*M} distance matrix between locations where the datum has been observed and locations where
#' the datum has to be predicted
#' @param metric_ts metric used on the tangent space. It must be chosen among "Frobenius", "FrobeniusScaled", "Correlation"
#' @param vario_model type of variogram fitted. It must be chosen among "Gaussian", "Spherical", "Exponential"
#' @param tol tolerance for the main loop of \code{model_kriging}
#' @param max_it maximum number of iterations for the main loop of \code{model_kriging}
#' @param n_h number of bins in the empirical variogram
#' @param tolerance_intrinsic tolerance for the computation of the intrinsic mean
#' @param X matrix (N rows and unrestricted number of columns) of additional covariates for the tangent space model, possibly NULL
#' @param X_new matrix (with the same number of rows of \code{new_coords}) of additional covariates for the new locations, possibly NULL
#' @param metric_manifold metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot"
#' @param model_ts type of model fitted on the tangent space. It must be chosen among "Intercept", "Coord1", "Coord2", "Additive"
#' @param distance type of distance between coordinates. It must be either "Eucldist" or "Geodist"

#' @return it returns a list with the following fields
#' \itemize{
#'       \item \code{fmean} {list of length \code{B}. Each field contains the prediction (at iteration \code{b}) for each location, obtained
#'                             as the intrinsic mean of the data within the tile it belongs to}
#'       \item \code{kervalues_mean} {Weights used for aggregating \code{fmean}}
#'       }
#' @details ...
#' @description ...
#' @useDynLib Manifoldgstat
RDD_OOK_boot_man_mixed = function(data_coords, data_val, K, grid, nk_min, B,
                            # spdist,
                            suppressMes,
                            ker.width.intrinsic,
                            # mesh,
                            graph.distance.complete,
                            # assign.matrix, no.assg.grid,
                            data.grid.distance,
                            # is.observed, border.length,
                            # num.signif.entries,
                            metric_ts,
                            vario_model,  tol, max_it, # method.analysis,
                            n_h, tolerance_intrinsic, X, X_new, # X_tot, ker.width.vario,
                            metric_manifold, model_ts, distance)
{
  # This function implements the bootstrap step of Ordinary Kriging with Random Domain
  # Decompositions for object data
  # Input:
  # data    = data.frame containing the data (the first two columns contain the coordinates, from the third on it contains the data)
  # K       = number of neighborhood (i.e., centers) to sample at each iteration
  # grid    = the grid of prediction
  # nk_min  = minimum number of observations within a neighborhood (default: 1)
  # B       = number of bootstap iterations (default: 100)
  # vario_model   = parametric variogram model ("Gaussian" or "Spherical" or "Exponential")
  # spdist  = method to compute distances in space ('euclidean' or 'graph.distance')
  # suppressMes = {T, F} controls the level of interaction and warnings given
  # tol     = tolerance for the main loop of model_kriging
  # max_it  = maximum number of iterations for the main loop of model_kriging
  # n_h     = number of bins in the empirical variogram
  # tolerance_intrinsic = Tolerance
  # X       = Additional covariates for the locations used to create the model
  # X_new   = Additional covariates for the locations where to perform kriging
  # X_tot   = Additional covariates for all N_samples
  # ker.width.intrinsic = parameter controlling the width of the Gaussian kernel for the computation of the local mean (if 0, no kernel is used)
  # ker.width.vario = parameter controlling the width of the Gaussian kernel for the computation of the empirical variogram (if 0, no kernel is used)
  # mesh    = contains the parameters defining the computed Delaunay triagulation
  # graph.distance.complete = nsub*nsub matrix (the [i,j] element is the length of the shortest path between points i and j)
  # assign.matrix = ngrid*num.triangles (the (i,j) element = TRUE if the point grid[i,] is inside the j-th triangle, FALSE otherwise)
  # no.assg.grid = indices of grid points that the function does not assign to any triangle
  # data.grid.distance
  # is.observed = If the i-th grid point is already observed ---> is.observed[i] contains the index of the point in the data matrix
  #               If the i-th grid point is unobserved       ---> is.observed[i] = 0
  # border.length = number of border points
  # p       = dimension of the manifold matrices
  # num.signif.entries = (p*(p+1))/2
  # t_step
  # method.analysis = {Local mean, Kriging} which method should be used
  # Output:
  # if (method.analysis == 'Local mean')
  #     list(fpred, variofit, kervalues_mean, ), with
  #          fpred   = list of B components, each containing a grid realization
  #          kervalues_mean = list of B components, each containing to be used for the aggregation of the local means (if aggregation_mean=="Weighted")
  #          kervalues_krig = list of B components, each containing to be used for the aggregation of the kriging predictions (if aggregation_krig=="Weighted")
  #          variofit = list of B components, each containing the parameters of the variogram over the grid
  #                     (i.e., for each grid point, the parameter used to make the prediction at that points are given)
  # if (method.analysis == 'Kriging')
  #     list(fmean, fpred, variofit, kervalues_mean, ), with
  #          fmean = list of B components, each containing the tangent points
  #          fpred   = list of B components, each containing a grid realization
  #          kervalues_mean = list of B components, each containing to be used for the aggregation of the local means (if aggregation_mean=="Weighted")
  #          kervalues_krig = list of B components, each containing to be used for the aggregation of the kriging predictions (if aggregation_krig=="Weighted")
  #          variofit = list of B components, each containing the parameter of the variogram over the grid
  #                     (i.e., for each grid point, the parameter used to make the prediction at that points are given)

  if(suppressMes)  oldw <- getOption("warn");  options(warn = -1)
  # Set parameters and initialize lists
  #*lists and parameters
  kervalues_mean=fmean_grid=fmean_data = list() # kervalues_krig=fmean=vfit
  # nk=rep(0,K)
  #*data
  colnames(data_coords)=c('x','y')
  # colnames(data_val)= c(paste0('z',1:num.signif.entries))
  nsub = dim(data)[1]
  ngrid = dim(grid)[1]

  #*distances
  graph.distance = graph.distance.complete
  # idx.start_sample = border.length+1;
  # idx.end_sample = border.length+N_samples
  # indexes_samples = idx.start_sample:idx.end_sample

  # data_manifold_tot_array=matrix_to_matrixArray(data[,3:(2+num.signif.entries)], p=p)

  # #------------------ Bootstrap step
  for(b in 1:B) # Repeat over the B bootstrap replicates
  {
    pFstat.predictTrK=gridk=kergrid_mean=list() #=kergrid_krig=vfitk
    if(suppressMes==F)
      print(paste("Repetition B = ", b, sep=''))

    # Step 1: Extract centers
    ### ------------- 1. Create RDD
    rdd = create.rdd(K=K, method.rdd = 'Voronoi', data_coords = data_coords,
                     # border.length = border.length, spdist = spdist,
                     graph.distance = graph.distance,
                     # mesh = mesh,
                     nk_min = nk_min, grid = grid,
                     # is.observed = is.observed, graph.distance.complete = graph.distance.complete, assign.matrix = assign.matrix,
                     data.grid.distance = data.grid.distance, suppressMes = suppressMes)
    assign = rdd$assign
    centers = rdd$centers
    assigng = rdd$assigng
    nk = table(assign)
    gridk = rdd$gridk
    graph.distance.grid.centers = rdd$graph.distance.grid.centers

    # veclocmean = matrix(NA, nrow=K, ncol=num.signif.entries)
    # veclocmean = list()

    # fpred[[b]]=matrix(NA,ngrid,num.signif.entries); colnames(fpred[[b]])=c(paste("Z",1:num.signif.entries))
    fmean_grid[[b]]= list() 
    fmean_data[[b]]= list() # Lista di ngrid matrici p*p

    if(ker.width.intrinsic>0) {
      kervalues_mean[[b]]=matrix(NA,ngrid,1); colnames(kervalues_mean[[b]])="Ker.val.mean"
    }
    # if(ker.width.vario>0) {
    #   kervalues_krig[[b]]=matrix(NA,ngrid,1); colnames(kervalues_krig[[b]])="Ker.val.krig"
    # }

    for(k in 1:K)
    {
      ##### Centro
      center = centers[k,]

      ##### Pesi
      weight.intrinsic = NULL

      if(ker.width.intrinsic>0) #compute kergrid = weights for aggregation (tangent point and prediction)
      {
        # if(K==1) dist_mat = t(as.matrix(graph.distance.grid.centers[,which(assigng == k)]))
        # else dist_mat = graph.distance.grid.centers[,which(assigng == k)]

        kergrid_mean[[k]]=kerfn(newdata=coordinates(gridk[[k]])[,1:2],center=c(center[1],center[2],k), ker.type = 'Gau', param = ker.width.intrinsic) # dist = spdist, distance.matrix = dist_mat
        kervalues_mean[[b]][which(assigng==k),1] = kergrid_mean[[k]]

        weight.intrinsic = kerfn(newdata=data[assign==k,1:2], center=center, ker.type = 'Gau', param = ker.width.intrinsic)  # dist=spdist, distance.matrix = graph.distance[, which(assign==k)]
      }

      # if(ker.width.vario>0) #compute kergrid = weights for aggregation (tangent point and prediction)
      # {
      #   # if(K==1) dist_mat = t(as.matrix(graph.distance.grid.centers[,which(assigng == k)]))
      #   # else dist_mat = graph.distance.grid.centers[,which(assigng == k)]
      #
      #   kergrid_krig[[k]]=kerfn(newdata=coordinates(gridk[[k]])[,1:2],center=c(center[1],center[2],k), ker.type = 'Gau', param = ker.width.vario) # dist = spdist,  distance.matrix = dist_mat
      #   kervalues_krig[[b]][which(assigng==k),1] = kergrid_krig[[k]]
      #
      #   }

      ##### Dati
      if(nk[as.character(k)]>0)  #table(factor(assigng, levels=1:K))[as.character(k)]>0
      {
        # datak=data.frame(data[assign==k,]) # extract data in k-th neighb.
        datak=data_val[,,assign==k] # extract data in k-th neighb.
        # # coordinates.datak=datak[,1:2]      # extract coordinates
        # datamat = matrix_to_matrixArray(datak[,-c(1,2)],p = p)

        # --- Estimate local tangent point
        # Sigma = intrinsic_mean(data = datamat, metric_manifold = metric_manifold,
        #                        metric_ts = metric_ts, weight_intrinsic = weight.intrinsic)
        # veclocmean[k,] = matrix_to_vec(Sigma)
        Sigma = intrinsic_mean(data = datak, metric_manifold = metric_manifold,
                               metric_ts = metric_ts, weight_intrinsic = weight.intrinsic)
        # veclocmean[[k]] = Sigma
        # --- Assign the local tangent point to each grid point & store results
        # fpred[[b]][which(assigng==k),]=matrix(rep(veclocmean[k,],nk[as.character(k)]), nrow = nk[as.character(k)], ncol=num.signif.entries, byrow = T)

        fmean_grid[[b]][which(assigng==k)]= lapply(seq_len(nk[as.character(k)]), function(X) Sigma)
        fmean_data[[b]][which(assign==k)]= lapply(seq_len(nk[as.character(k)]), function(X) Sigma)

      }
    } # for su K

    # # Grid points excluded from the prediction (set by user)
    # if (length(no.assg.grid)!=0) {
    #   fpred[[b]][no.assg.grid,]=NA
    #   if(method.analysis == 'Kriging') {
    #     fmean[[b]][no.assg.grid,]=NA
    #     vfit[[b]][no.assg.grid,]=NA
    #   }
    #   if(ker.width.intrinsic>0) {
    #     kervalues_mean[[b]][no.assg.grid,]=NA
    #   }
    #   if(ker.width.vario>0) {
    #     kervalues_krig[[b]][no.assg.grid,]=NA
    #   }
    # }

  } # for su B
  # if (method.analysis == 'Local mean')

  list.ret = list(fmean_data=fmean_data, fmean_grid=fmean_grid,  kervalues_mean=kervalues_mean) # ,kervalues_krig=kervalues_krig, variofit=vfit

  if(suppressMes) options(warn = oldw)

  return(list.ret)
}
