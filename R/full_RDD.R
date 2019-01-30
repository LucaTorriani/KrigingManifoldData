#' Perform full_RDD
#' @useDynLib Manifoldgstat
#' @export
full_RDD = function(data, K, grid, nk_min=1, B=100,
                    # spdist='euclidean',
                    suppressMes=F,
                    tol=1e-12, max_it=100, n_h=15, tolerance_intrinsic =10^(-6), X=NULL, X_new=NULL, X_tot=NULL, plot=FALSE,
                    ker.width.intrinsic = 0, ker.width.vario = 1.5,
                    # mesh,
                    graph.distance.complete,
                    # assign.matrix, no.assg.grid,
                    data.grid.distance,
                    # is.observed, border.length,
                    N_samples, p, num.signif.entries,
                    aggregation_mean, aggregation_kriging, method.analysis = 'Local mean',
                    metric_manifold, metric_ts, model_ts,
                    vario_model, distance)

{
  if(aggregation_mean != "Weighted" && aggregation_mean != "Equal") stop("aggregation_mean expected to be: Equal or Weighted")
  if(aggregation_kriging != "Weighted" && aggregation_kriging != "Equal") stop("aggregation_mean expected to be: Equal or Weighted")

  if(ker.width.intrinsic==0 & aggregation_mean=="Weighted") warning("Aggregation for the mean will use equal weights since ker.width.intrinsic=0")
  if(ker.width.vario==0 & aggregation_kriging=="Weighted") warning("Aggregation for the kriging predictors will use equal weights since ker.width.vario=0")

  # This function implements RDD-OOK procedure, by using the functions RDD_OOK_boot_man and RDD_OOK_aggr_man
  if(K==1)
  {
    resBootstrap=RDD_OOK_boot_man_mixed(data=data, K=K, grid=grid, nk_min=nk_min, B=1,
                                  # spdist=spdist,
                                  suppressMes=suppressMes, tol=tol, max_it = max_it,
                                  n_h=n_h, tolerance_intrinsic=tolerance_intrinsic, X=X, X_new=X_new, X_tot=X_tot, plot=plot,
                                  ker.width.intrinsic = 0,
                                  ker.width.vario = 0,
                                  # mesh=mesh,
                                  graph.distance.complete= graph.distance.complete,
                                  # assign.matrix=assign.matrix, no.assg.grid=no.assg.grid,
                                  data.grid.distance = data.grid.distance,
                                  # is.observed=is.observed, border.length= border.length,
                                  p=p,num.signif.entries=num.signif.entries,
                                  method.analysis = method.analysis,
                                  metric_manifold = metric_manifold, metric_ts = metric_ts, model_ts = model_ts,
                                  vario_model = vario_model, distance = distance)


    resAggregated=simplify2array(resBootstrap$fpred)[,,1]
    if(method.analysis == 'Kriging'){
      resLocalMean=simplify2array(resBootstrap$fmean)[,,1]
      return(list(resBootstrap = resBootstrap,
                  resAggregated = resAggregated,
                  resLocalMean = resLocalMean))
    }
  }
  if(K>1)
  {
    resBootstrap=RDD_OOK_boot_man_mixed(data=data, K=K, grid=grid, nk_min=nk_min, B=B,
                                  # spdist=spdist,
                                  suppressMes=suppressMes, tol=tol, max_it = max_it,
                                  n_h=n_h, tolerance_intrinsic =tolerance_intrinsic, X=X, X_new=X_new, X_tot=X_tot, plot=plot,
                                  ker.width.intrinsic = ker.width.intrinsic,
                                  ker.width.vario = ker.width.vario,
                                  # mesh=mesh,
                                  graph.distance.complete= graph.distance.complete,
                                  # assign.matrix=assign.matrix, no.assg.grid=no.assg.grid,
                                  data.grid.distance = data.grid.distance,
                                  # is.observed=is.observed, border.length= border.length,
                                  p=p,num.signif.entries=num.signif.entries,
                                  method.analysis = method.analysis,
                                  metric_manifold = metric_manifold, metric_ts = metric_ts, model_ts = model_ts,
                                  vario_model = vario_model, distance = distance)

    if(method.analysis == 'Local mean') {
      if (aggregation_mean== "Equal") ker.width.intrinsic=0
      resAggregated=RDD_OOK_aggr_man_mixed(fOKBV = resBootstrap$fpred, weights_intrinsic = resBootstrap$kervalues_mean,
                                     ker.width.intrinsic=  ker.width.intrinsic, N_samples = N_samples, p=p, num.signif.entries = num.signif.entries )
    }


    if(method.analysis == 'Kriging')
    {
      if (aggregation_kriging== "Equal") ker.width.vario = 0
      resAggregated=RDD_OOK_aggr_man_mixed(fOKBV = resBootstrap$fpred, weights_intrinsic = resBootstrap$kervalues_krig,
                                     ker.width.intrinsic=  ker.width.vario, N_samples = N_samples, p=p, num.signif.entries = num.signif.entries )

      if (aggregation_mean== "Equal") ker.width.intrinsic = 0
      resLocalMean=RDD_OOK_aggr_man_mixed(fOKBV = resBootstrap$fmean, weights_intrinsic = resBootstrap$kervalues_mean,
                                    ker.width.intrinsic =  ker.width.intrinsic, N_samples = N_samples, p=p, num.signif.entries = num.signif.entries )


      return(list(resBootstrap = resBootstrap,
                  resAggregated = resAggregated,
                  resLocalMean = resLocalMean))
    }
  }
  return(list(resBootstrap = resBootstrap,
              resAggregated = resAggregated  ))
}
