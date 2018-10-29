#include "IntrinsicMean.hpp"
#include <iostream>
#include <Rcpp.h>

MatrixXd intrinsic_mean_C(const std::vector<MatrixXd>& data_manifold, std::string distance_Manifold_name, map_functions::logarithmicMap& logMap,
   map_functions::exponentialMap& expMap,  distances_tplane::DistanceTplane& distanceTplane, double tolerance, const Vec& weight_intrinsic, const Vec& weight_extrinsic){

    // Weights
    double sum_weight_intrinsic(weight_intrinsic.sum());

    if (distance_Manifold_name == "Correlation") {
      // Punto tangente
      Eigen::MatrixXd Sigma = matrix_manipulation::Chol_decomposition(extrinsic_mean(data_manifold, weight_extrinsic));
      unsigned int p = Sigma.rows();

      logMap.set_members(Sigma);
      expMap.set_members(Sigma);
      distanceTplane.set_members(Sigma);

      // CODE
      Eigen::MatrixXd Xk(p,p);
      Eigen::MatrixXd Xk_prec(p,p);

      Xk = logMap.map2tplane(data_manifold[0]);

      double tau(1.0);
      double tmp;
      double tolk;
      double tolk_prec(tolerance + 1);
      unsigned int N(data_manifold.size());

      size_t num_iter(0);

      while(tolk_prec > tolerance && num_iter < 100) {

        Xk_prec = Xk;
        tmp = distanceTplane.norm(Xk_prec);
        tolk_prec = tmp*tmp;

        Xk.setZero();
        for (size_t i=0; i<N; i++) {
          Xk = Xk + weight_intrinsic(i)* logMap.map2tplane(data_manifold[i]);
        }
        Xk = Xk/sum_weight_intrinsic;
        Sigma = expMap.map2manifold(tau*Xk);

        logMap.set_members(Sigma);
        expMap.set_members(Sigma);
        distanceTplane.set_members(Sigma);

        tmp = distanceTplane.norm(Xk);
        tolk = tmp*tmp;
        if (tolk > tolk_prec) {
          tau = tau/2;
          Xk = Xk_prec;
        }
        num_iter++;
      }
      if(num_iter == 100) Rcpp::warning("Reached max number of iterations in intrinsic_mean");

      return (Sigma);
    }
    else {
      // Punto tangente
      Eigen::MatrixXd Sigma((data_manifold)[0]);
      unsigned int p = Sigma.rows();

      logMap.set_members(Sigma);
      expMap.set_members(Sigma);
      distanceTplane.set_members(Sigma);

      // CODE
      Eigen::MatrixXd Xk(p,p);
      Eigen::MatrixXd Xk_prec(p,p);

      Xk = data_manifold[0];

      double tau(1.0);
      double tmp;
      double tolk;
      double tolk_prec(tolerance + 1);
      unsigned int N(data_manifold.size());

      size_t num_iter(0);

      while(tolk_prec > tolerance && num_iter < 100) {

        Xk_prec = Xk;
        tmp = distanceTplane.norm(Xk_prec);
        tolk_prec = tmp*tmp;

        Xk.setZero();
        for (size_t i=0; i<N; i++) {
          Xk = Xk + weight_intrinsic(i)* logMap.map2tplane(data_manifold[i]);
        }
        Xk = Xk/sum_weight_intrinsic;
        Sigma = expMap.map2manifold(tau*Xk);

        logMap.set_members(Sigma);
        expMap.set_members(Sigma);
        distanceTplane.set_members(Sigma);

        tmp = distanceTplane.norm(Xk);
        tolk = tmp*tmp;
        if (tolk > tolk_prec) {
          tau = tau/2;
          Xk = Xk_prec;
        }
        num_iter++;
      }
      if(num_iter == 100) Rcpp::warning("Reached max number of iterations in intrinsic_mean");

      return (Sigma);
    }

}

MatrixXd extrinsic_mean (const std::vector<MatrixXd>& data_manifold, const Vec& weight_extrinsic) {
  unsigned int N(data_manifold.size());
  unsigned int p(data_manifold[0].rows());
  MatrixXd result(p,p);
  for (size_t i =0; i<N; i++) {
    result += data_manifold[i]*weight_extrinsic(i);  // Assuming data_manifold are already Chol
  }
  result = result/N;
  double col_norm;
  for (size_t j =0; j<p; j++) {
    col_norm = (result.col(j)).norm();
    result.col(j) /= col_norm;
  }
  return(result.transpose() * result);
}
