#' Perform mixed_RDD
#' @useDynLib Manifoldgstat
#' @export
mixed_RDD = function(data_coords, data_val, K, grid, nk_min=1, B=100,
                       # spdist='euclidean',
                       suppressMes=F,
                       ker.width.intrinsic = 0,
                       # mesh,
                       graph.distance.complete,
                       # assign.matrix, no.assg.grid,
                       data.grid.distance,
                       # is.observed, border.length,
                       N_samples, p, num.signif.entries,
                       aggregation_mean, metric_ts)
  # tol=1e-12, max_it=100, n_h=15, tolerance_intrinsic =10^(-6), X=NULL, X_new=NULL, X_tot=NULL, plot=FALSE,
  # ker.width.vario = 1.5, aggregation_kriging, method.analysis = 'Local mean',
  # metric_manifold, metric_ts, model_ts, vario_model, distance

{
  if(aggregation_mean != "Weighted" && aggregation_mean != "Equal") stop("aggregation_mean expected to be: Equal or Weighted")

  if(ker.width.intrinsic==0 & aggregation_mean=="Weighted") warning("Aggregation for the mean will use equal weights since ker.width.intrinsic=0")

  # This function implements RDD-OOK procedure, by using the functions RDD_OOK_boot_man and RDD_OOK_aggr_man
  if(K==1)
  {
    resBootstrap=RDD_OOK_boot_man_mixed(data_coords=data_coords, data_val=data_val, K=K, grid=grid, nk_min=nk_min, B=1,
                                  # spdist=spdist,
                                  suppressMes=suppressMes,
                                  ker.width.intrinsic = 0,
                                  # mesh=mesh,
                                  graph.distance.complete= graph.distance.complete,
                                  # assign.matrix=assign.matrix, no.assg.grid=no.assg.grid,
                                  data.grid.distance = data.grid.distance,
                                  # is.observed=is.observed, border.length= border.length,
                                  p=p,num.signif.entries=num.signif.entries, metric_ts=metric_ts)
    # tol=tol, max_it = max_it, n_h=n_h, tolerance_intrinsic=tolerance_intrinsic, X=X, X_new=X_new, X_tot=X_tot, plot=plot,
    # ker.width.vario = 0, method.analysis = method.analysis, metric_manifold = metric_manifold, metric_ts = metric_ts, model_ts = model_ts,
    # vario_model = vario_model, distance = distance


    # resAggregated=simplify2array(resBootstrap$fpred)[,,1]
    resAggregated = resBootstrap$fmean[[1]]
  }
  if(K>1)
  {
    resBootstrap=RDD_OOK_boot_man_mixed(data_coords=data_coords, data_val=data_val, K=K, grid=grid, nk_min=nk_min, B=B,
                                  # spdist=spdist,
                                  suppressMes=suppressMes,
                                  ker.width.intrinsic = ker.width.intrinsic,
                                  # mesh=mesh,
                                  graph.distance.complete= graph.distance.complete,
                                  # assign.matrix=assign.matrix, no.assg.grid=no.assg.grid,
                                  data.grid.distance = data.grid.distance,
                                  # is.observed=is.observed, border.length= border.length,
                                  p=p,num.signif.entries=num.signif.entries, metric_ts=metric_ts)
    # tol=tol, max_it = max_it, n_h=n_h, tolerance_intrinsic =tolerance_intrinsic, X=X, X_new=X_new, X_tot=X_tot, plot=plot,
    # ker.width.vario = ker.width.vario, method.analysis = method.analysis,
    # metric_manifold = metric_manifold, metric_ts = metric_ts, model_ts = model_ts,
    # vario_model = vario_model, distance = distance
    if (aggregation_mean== "Equal") ker.width.intrinsic=0
    resAggregated=RDD_OOK_aggr_man_mixed(fOKBV = resBootstrap$fmean, weights_intrinsic = resBootstrap$kervalues_mean,
                                   ker.width.intrinsic=  ker.width.intrinsic, N_samples = N_samples, p=p, num.signif.entries = num.signif.entries )

  }

  fmean = resRDD_OOK$resAggregated

  model_pred = model_kriging_mixed (data_manifold = data_val, coords = data_coords, X = NULL, Sigma_data = fmean[1:N_samples], metric_manifold = metric_manifold,
                              model_ts = model_ts, vario_model = vario_model, # metric_ts = "Frobenius",
                              n_h=15, distance = distance, data_dist_mat=graph.distance.complete, data_grid_dist_mat=, max_it = 100, tolerance = 1e-6, # weight_vario = NULL,
                              # weight_intrinsic = NULL, tolerance_intrinsic = 1e-6,
                              max_sill = NULL, max_a = NULL,
                              new_coords=prediction_grid, Sigma_new = fmean, X_new = NULL, plot = TRUE, suppressMes = FALSE)

  return(list(resBootstrap = resBootstrap,
              resAggregated = resAggregated  ))
}
