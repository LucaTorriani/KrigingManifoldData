#include <iostream>
#include <vector>
#include <utility>
#include <memory>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "DistanceManifold.hpp"
#include "Helpers.hpp"
#include "HelpersFactory.hpp"


extern "C"{
  RcppExport SEXP distance_manifold (SEXP s_data1, SEXP s_data2, SEXP s_N1, SEXP s_N2, SEXP s_manifold_metric) {
      BEGIN_RCPP
      // Data1
      unsigned int N1(Rcpp::as<unsigned int> (s_N1));
      std::vector<Eigen::MatrixXd> data1(N1);
      for(size_t i=0; i<N1; i++){
        data1[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data1,i));
      }

      // Data2
      unsigned int N2(Rcpp::as<unsigned int> (s_N2));
      std::vector<Eigen::MatrixXd> data2(N1);
      for(size_t i=0; i<N2; i++){
        data2[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data2,i));
      }

      // Distance manifold
      std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, FrobeniusScaled)
      manifold_factory::ManifoldFactory& manifold_fac (manifold_factory::ManifoldFactory::Instance());
      std::unique_ptr<distances_manifold::DistanceManifold> theManifoldDist = manifold_fac.create(distance_Manifold_name);

      std::vector<double> dist_vec(N1);
      if (N2==1) {
        for(size_t i=0; i<N1; i++) {
          dist_vec[i] = theManifoldDist->compute_distance(data1[i],data2[0]);
        }
      }
      else {
        for (size_t i=0; i<N1; i++) {
          dist_vec[i] = theManifoldDist->compute_distance(data1[i],data2[i]);
        }
      }

      return Rcpp::wrap(dist_vec);
      END_RCPP
}

}
