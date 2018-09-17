#include <iostream>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "Coordinates.hpp"
#include "Helpers.hpp"
#include "Distance.hpp"
#include "Point.hpp"


extern "C"{

  RcppExport SEXP get_model (SEXP s_data_manifold, SEXP s_coordinates, SEXP s_X, SEXP s_Sigma,
    SEXP s_distance, SEXP s_manifold_metric, SEXP s_ts_metric, SEXP s_ts_model, SEXP s_vario_model, SEXP s_n_h
    SEXP s_max_it, SEXP s_tolerance, SEXP s_weight) {

      // Punto tangente
      Eigen::Map<Eigen::MatrixXd> Sigma(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_Sigma));
      unsigned int n = Sigma.rows();

      // Distance
      std::string distance_name = Rcpp::as<std::string> (s_distance) ; //(Geodist, Euclidean)
      distances::Distance distance(distance_name);

      // Distance tplane
      std::string distance_name = Rcpp::as<std::string> (s_ts_metric) ; //(Frobenius, FrobeniusScaled)
      distances_tplane::DistanceTplane distanceTplane (distance_Tplane_name, Sigma);

      // Distance manifold
      std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
      distances_manifold::DistanceManifold distanceManifold (distance_Manifold_name, Sigma);

      // Map functions
      map_functions::logarithmicMap logMap(distanceManifold);

      // Coordinates
      Eigen::Map<Eigen::MatrixXd> coords_mat(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_coordinates));
      Coordinates coords(coords_mat);

      // Distance Matrix
      SpMat distanceMatrix(N,N);
      distanceMatrix = distance.create_distance_matrix(coords);

      // Data manifold
      Rcpp::List list_data_manifold(data_manifold);
      size_t N = L.size();
      std::vector<Eigen::MatrixXd> data_manifold(N);
      for(size_t i=0; i<N; i++) data_manifold[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(L,i));

      // Data tangent space
      std::vector<Eigen::MatrixXd> data_tspace(N);
      for (size_t i=0; i<N; i++) {
        data_tspace[i] = logMap.map2tplane(data_manifold[i]);
      }
      Eigen::MatrixXd big_matrix_data_tspace(N, (n+1)*n/2);
      big_matrix_data_tspace = matrix_manipulation::VecMatrices2bigMatrix(data_tspace);

      // Emp vario
      unsigned int n_h (Rcpp::as<unsigned int>( s_n_h));
      variogram_evaluation::EmpiricalVariogram emp_vario(coords, distance, n_h, distanceTplane, distanceMatrix);

      // Fitted vario
      vario_factory::VariogramFactory & vf(vario_factory::VariogramFactory::Instance());
      std::string variogram_type = Rcpp::as<std::string> (s_vario_model); // (Gaussian, Exponential, Spherical) //IMPLEMENTAREEE
      std::unique_ptr<variogram_evaluation::FittedVariogram> the_variogram = vf.create(variogram_type);

      // Gamma matrix
      Eigen::MatrixXd gamma_matrix(Eigen::MatrixXd::Identity(N,N));

      // Design matrix
      design_matrix::registerDesignMatrices();
      design_matrix::DesignMatrixFactory& design_matrix_fac = design_matrix::DesignMatrixFactory::Instance();
      std::string variogram_type (Rcpp::as<std::string> (s_vario_model)); //(Intercept, Coord1, Coord2, Additive)
      std::unique_ptr<design_matrix::DesignMatrix> theDesign_matrix = design_matrix_fac.create(model_name);

      Eigen::MatrixXd design_matrix = theDesign_matrix->compute_design_matrix(coords);

      // Model
      unsigned int n_covariates(design_matrix.cols());
      model_fit::Model model(big_matrix_data_tspace, design_matrix, n);
      model.update_model(gamma_matrix);
      Eigen::MatrixXd residuals(N, ((n+1)*n)/2);
      Eigen::MatrixXd beta(N, ((n+1)*n)/2);
      Eigen::MatrixXd beta_old(N, ((n+1)*n)/2);
      beta = model.get_beta();
      std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
      beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, n);
      std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

      unsigned int num_iter(0);
      unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
      double tolerance(Rcpp::as<unsigned int> (s_tolerance));
      std::vector<MatrixXd> res(_N);

      double tol = tolerance+1;
      while (num_iter < max_iter && tol > tolerance) {
        resMatrix = model.get_residuals();
        residuals = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, _n);
        emp_vario.update_emp_vario(residuals);
        the_variogram->evaluate_par_fitted(emp_vario);
        gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix, N);
        beta_old_vec_matrices = beta_vec_matrices;
        model.update_model(gamma_matrix);
        beta = model.get_beta();
        beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, n);
        tol=0;
        for (size_t i=0; i<N; i++) {
          tol += distanceTplane.compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
        }
        num_iter++;
      }

    }

    Eigen::Vector3d parameters = the_variogram->get_parameters();


    List result = List::create(Named("beta") = beta_vec_matrices,
                           Named("gamma_matrix") = gamma_matrix,
                           Named("residuals") = residuals,
                           Named("vario_parameters") = parameters);

    return wrap(result);


}














  // DISTANCE

  // errore se GEODIST e più di 2 _coords

  // COORDINATES
  // latitudine prima coordinata
  // lat e long in gradi
