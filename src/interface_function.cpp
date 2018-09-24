#include <iostream>
#include <vector>
#include <utility>
#include <memory>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "Coordinates.hpp"
#include "DesignMatrix.hpp"
#include "DistanceManifold.hpp"
#include "DistanceTplane.hpp"
#include "Distance.hpp"
#include "EmpiricalVariogram.hpp"
#include "FittedVariogram.hpp"
#include "Helpers.hpp"
#include "HelpersFactory.hpp"
#include "MapFunctions.hpp"
#include "Model.hpp"


extern "C"{

// CREATE MODEL
  RcppExport SEXP get_model (SEXP s_data_manifold, SEXP s_coordinates, SEXP s_X, SEXP s_Sigma,
    SEXP s_distance, SEXP s_manifold_metric, SEXP s_ts_metric, SEXP s_ts_model, SEXP s_vario_model, SEXP s_n_h,
    SEXP s_max_it, SEXP s_tolerance, SEXP s_weight) {

      Rcpp::Nullable<Vec> weight(s_weight);
      Rcpp::Nullable<Eigen::MatrixXd> X(s_X);

      // Punto tangente
      Eigen::Map<Eigen::MatrixXd> Sigma(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_Sigma));
      unsigned int n = Sigma.rows();

      // Distance tplane
      std::string distance_Tplane_name = Rcpp::as<std::string> (s_ts_metric) ; //(Frobenius, FrobeniusScaled)
      distances_tplane::DistanceTplane distanceTplane (distance_Tplane_name, std::make_shared<const Eigen::MatrixXd> (Sigma));

      // Distance manifold
      std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
      distances_manifold::DistanceManifold distanceManifold (distance_Manifold_name, std::make_shared<const Eigen::MatrixXd> (Sigma));

      // Map functions
      map_functions::logarithmicMap logMap(distanceManifold);

      // Data manifold
      Rcpp::List list_data_manifold(s_data_manifold);
      size_t N = list_data_manifold.size();
      std::vector<Eigen::MatrixXd> data_manifold(N);
      for(size_t i=0; i<N; i++){
        data_manifold[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(list_data_manifold,i));
      }

      // Data tangent space
      std::vector<Eigen::MatrixXd> data_tspace(N);
      for (size_t i=0; i<N; i++) {
        data_tspace[i] = logMap.map2tplane(data_manifold[i]);
      }

      Eigen::MatrixXd big_matrix_data_tspace(N, (n+1)*n/2);
      big_matrix_data_tspace = matrix_manipulation::VecMatrices2bigMatrix(data_tspace);

      // Distance
      std::string distance_name = Rcpp::as<std::string> (s_distance) ; //(Geodist, Eucldist)
      distances::Distance distance(distance_name);

      // Coordinates
      const Eigen::Map<Eigen::MatrixXd> coords_mat(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_coordinates));
      Coordinates coords(std::make_shared<const Eigen::MatrixXd>(coords_mat));

      // Distance Matrix
      SpMat distanceMatrix(N,N);
      distanceMatrix = distance.create_distance_matrix(coords, N);

      // Emp vario
      unsigned int n_h (Rcpp::as<unsigned int>( s_n_h));
      variogram_evaluation::EmpiricalVariogram emp_vario(coords, distance, n_h, distanceTplane, std::make_shared<const SpMat> (distanceMatrix));

      if(weight.isNotNull()) {
        Eigen::Map<Vec> weight(Rcpp::as<Eigen::Map<Vec>> (s_weight));
        emp_vario.set_weight(weight);
      }

      // Fitted vario
      vario_factory::VariogramFactory & vf(vario_factory::VariogramFactory::Instance());
      std::string variogram_type (Rcpp::as<std::string> (s_vario_model));
      std::unique_ptr<variogram_evaluation::FittedVariogram> the_variogram = vf.create(variogram_type);

      // Gamma matrix
      Eigen::MatrixXd gamma_matrix(Eigen::MatrixXd::Identity(N,N));

      // Design matrix
      design_matrix::registerDesignMatrices();
      design_matrix::DesignMatrixFactory& design_matrix_fac = design_matrix::DesignMatrixFactory::Instance();
      std::string model_name (Rcpp::as<std::string> (s_ts_model)); // (Additive, Coord1, Coord2, Intercept)
      std::unique_ptr<design_matrix::DesignMatrix> theDesign_matrix = design_matrix_fac.create(model_name);

      Eigen::MatrixXd design_matrix;
      if(X.isNotNull()) {
        Eigen::Map<Eigen::MatrixXd> X(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_X));
        design_matrix = theDesign_matrix->compute_design_matrix(coords, X);
      }
      else design_matrix = theDesign_matrix->compute_design_matrix(coords);

      // Model
      unsigned int n_covariates(design_matrix.cols());
      model_fit::Model model(std::make_shared<const Eigen::MatrixXd> (big_matrix_data_tspace), std::make_shared<const Eigen::MatrixXd> (design_matrix), n);
      model.update_model(gamma_matrix);
      Eigen::MatrixXd resMatrix(N, ((n+1)*n)/2);
      Eigen::MatrixXd beta(n_covariates, ((n+1)*n)/2);
      Eigen::MatrixXd beta_old(n_covariates, ((n+1)*n)/2);
      beta = model.get_beta(); // OK fino a qui
      std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
      beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, n);

      std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

      unsigned int num_iter(0);
      unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
      double tolerance(Rcpp::as<double> (s_tolerance));
      std::vector<MatrixXd> resVec(N);

      double tol = tolerance+1;

      while (num_iter < max_iter && tol > tolerance) {
        resMatrix = model.get_residuals();
        resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, n);

        emp_vario.update_emp_vario(resVec);
        the_variogram -> evaluate_par_fitted(emp_vario);

        gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix, N);
        beta_old_vec_matrices = beta_vec_matrices;

        model.update_model(gamma_matrix);
        beta = model.get_beta();
        beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, n);

        tol=0.0;
        for (size_t i=0; i<n_covariates; i++) {
          tol += distanceTplane.compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
        }

        num_iter++;
      }

      unsigned int n_hh(1000);
      Vec hh(n_hh);
      std::vector<double> emp_vario_values(emp_vario.get_card_h());
      emp_vario_values = emp_vario.get_hvec();

      hh.setLinSpaced(n_hh, 0, *std::max_element(emp_vario_values.begin(), emp_vario_values.end()));

      Vec fit_vario_values = the_variogram->get_vario_vec(hh, n_hh);


      Eigen::Vector3d parameters = the_variogram->get_parameters();

      Rcpp::List result = Rcpp::List::create(Rcpp::Named("beta") = beta_vec_matrices,
                           Rcpp::Named("fit_vario_values") = fit_vario_values,
                           Rcpp::Named("hh") = hh,
                           Rcpp::Named("gamma_matrix") = gamma_matrix,
                           Rcpp::Named("residuals") = resVec,
                           Rcpp::Named("emp_vario_values") = emp_vario.get_emp_vario_values(),
                           Rcpp::Named("h_vec") = emp_vario.get_hvec(),
                           Rcpp::Named("vario_parameters") = parameters,
                           Rcpp::Named("iterations") = num_iter);

      return Rcpp::wrap(result);
  }


// KRIGING
  RcppExport SEXP get_kriging (SEXP s_coordinates, SEXP s_new_coordinates,  SEXP s_Sigma,
    SEXP s_distance, SEXP s_manifold_metric, SEXP s_ts_model, SEXP s_vario_model,
    SEXP s_beta, SEXP s_gamma_matrix, SEXP s_vario_parameters, SEXP s_residuals, SEXP s_X_new) {

    Rcpp::Nullable<Eigen::MatrixXd> X_new(s_X_new);


    // Punto tangente
    Eigen::Map<Eigen::MatrixXd> Sigma(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_Sigma));
    unsigned int n = Sigma.rows();

    // Distance
    std::string distance_name = Rcpp::as<std::string> (s_distance) ; //(Geodist, Euclidean)
    distances::Distance distance(distance_name);

    // Distance manifold
    std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
    distances_manifold::DistanceManifold distanceManifold (distance_Manifold_name, std::make_shared<const Eigen::MatrixXd> (Sigma));

    // Map functions
    map_functions::exponentialMap expMap(distanceManifold);

    // Old coordinates
    const Eigen::Map<Eigen::MatrixXd> coords_mat(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_coordinates));
    Coordinates coords(std::make_shared<const Eigen::MatrixXd>(coords_mat));
    unsigned int N = coords_mat.rows();

    // New coordinates
    const Eigen::Map<Eigen::MatrixXd> new_coords_mat(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_new_coordinates));
    Coordinates new_coords(std::make_shared<const Eigen::MatrixXd>(new_coords_mat));
    unsigned int M = new_coords_mat.rows();

    // New Design matrix
    design_matrix::registerDesignMatrices();
    design_matrix::DesignMatrixFactory& design_matrix_fac = design_matrix::DesignMatrixFactory::Instance();
    std::string model_name (Rcpp::as<std::string> (s_ts_model)); // (Additive, Coord1, Coord2, Intercept)
    std::unique_ptr<design_matrix::DesignMatrix> theDesign_matrix = design_matrix_fac.create(model_name);

    Eigen::MatrixXd new_design_matrix;
    if(X_new.isNotNull()) {
      Eigen::Map<Eigen::MatrixXd> X_new(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_X_new));
      new_design_matrix = theDesign_matrix->compute_design_matrix(new_coords, X_new);
    }
    else new_design_matrix = theDesign_matrix->compute_design_matrix(new_coords);

    // Fitted vario
    Eigen::Map<Eigen::VectorXd> parameters(Rcpp::as<Eigen::Map<Eigen::VectorXd>> (s_vario_parameters));
    // Andrebbe fatto Vector3d ma non funziona

    vario_factory::VariogramFactory & vf(vario_factory::VariogramFactory::Instance());
    std::string variogram_type (Rcpp::as<std::string> (s_vario_model)); // (Gaussian, Exponential, Spherical) //IMPLEMENTAREEE
    std::unique_ptr<variogram_evaluation::FittedVariogram> the_variogram = vf.create(variogram_type);
    the_variogram->set_parameters(parameters);

    // Gamma matrix
    Eigen::Map<Eigen::MatrixXd> gamma_matrix(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_gamma_matrix));

    // Beta (lista di matrici)
    Rcpp::List list_beta(s_beta);
    size_t num_cov = list_beta.size();
    std::vector<Eigen::MatrixXd> beta_vec(N);
    for(size_t i=0; i<num_cov; i++) beta_vec[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(list_beta,i));

    // Residui (lista di matrici)
    Rcpp::List list_residuals(s_residuals);
    std::vector<Eigen::MatrixXd> residuals_vec(N);
    for(size_t i=0; i<N; i++) residuals_vec[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(list_residuals,i));

    std::vector<double> distanceVector(N);
    Vec ci(N);
    Vec lambda_vec(N);

    Eigen::LDLT<Eigen::MatrixXd> solver(N);
    solver.compute(gamma_matrix);

    Eigen::MatrixXd tmp(n,n);
    auto weighted_sum_beta = [&beta_vec, &num_cov, &tmp, n] (const Vec& design_matrix_row) { tmp.setZero(n,n); for (size_t j=0; j<num_cov; j++) tmp = tmp + beta_vec[j]*design_matrix_row(j); return tmp;};
    auto weighted_sum_residuals = [&residuals_vec, &N, &tmp, n] (const Vec& lambda_vec) { tmp.setZero(n,n); for (size_t j=0; j<N; j++) tmp = tmp + residuals_vec[j]*lambda_vec(j); return tmp;};

    Eigen::MatrixXd tplane_prediction(n,n);
    std::vector<Eigen::MatrixXd> manifold_prediction(M);

    for (size_t i=0; i<M; i++) {
      distanceVector = distance.create_distance_vector(coords, new_coords_mat.row(i));
      ci = the_variogram->get_covario_vec(distanceVector, M);
      lambda_vec = solver.solve(ci);
      tplane_prediction = weighted_sum_beta(new_design_matrix.row(i)) + weighted_sum_residuals(lambda_vec);
      manifold_prediction[i] = expMap.map2manifold(tplane_prediction);
    }


    Rcpp::List result = Rcpp::List::create(Rcpp::Named("prediction") = manifold_prediction);  // Solo questo?

    return Rcpp::wrap(result);
  }


  // CREATE MODEL AND KRIGNG
    RcppExport SEXP get_model_and_kriging (SEXP s_data_manifold, SEXP s_coordinates, SEXP s_X, SEXP s_Sigma,
      SEXP s_distance, SEXP s_manifold_metric, SEXP s_ts_metric, SEXP s_ts_model, SEXP s_vario_model, SEXP s_n_h,
      SEXP s_max_it, SEXP s_tolerance, SEXP s_weight, SEXP s_new_coordinates, SEXP s_X_new) {

        Rcpp::Nullable<Vec> weight(s_weight);
        Rcpp::Nullable<Eigen::MatrixXd> X(s_X);
        Rcpp::Nullable<Eigen::MatrixXd> X_new(s_X_new);

        // Punto tangente
        Eigen::Map<Eigen::MatrixXd> Sigma(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_Sigma));
        unsigned int n = Sigma.rows();

        // Distance tplane
        std::string distance_Tplane_name = Rcpp::as<std::string> (s_ts_metric) ; //(Frobenius, FrobeniusScaled)
        distances_tplane::DistanceTplane distanceTplane (distance_Tplane_name, std::make_shared<const Eigen::MatrixXd> (Sigma));

        // Distance manifold
        std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
        distances_manifold::DistanceManifold distanceManifold (distance_Manifold_name, std::make_shared<const Eigen::MatrixXd> (Sigma));

        // Map functions
        map_functions::logarithmicMap logMap(distanceManifold);

        // Data manifold
        Rcpp::List list_data_manifold(s_data_manifold);
        size_t N = list_data_manifold.size();
        std::vector<Eigen::MatrixXd> data_manifold(N);
        for(size_t i=0; i<N; i++){
          data_manifold[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(list_data_manifold,i));
        }

        // Data tangent space
        std::vector<Eigen::MatrixXd> data_tspace(N);
        for (size_t i=0; i<N; i++) {
          data_tspace[i] = logMap.map2tplane(data_manifold[i]);
        }

        Eigen::MatrixXd big_matrix_data_tspace(N, (n+1)*n/2);
        big_matrix_data_tspace = matrix_manipulation::VecMatrices2bigMatrix(data_tspace);

        // Distance
        std::string distance_name = Rcpp::as<std::string> (s_distance) ; //(Geodist, Eucldist)
        distances::Distance distance(distance_name);

        // Coordinates
        const Eigen::Map<Eigen::MatrixXd> coords_mat(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_coordinates));
        Coordinates coords(std::make_shared<const Eigen::MatrixXd>(coords_mat));

        // Distance Matrix
        SpMat distanceMatrix(N,N);
        distanceMatrix = distance.create_distance_matrix(coords, N);

        // Emp vario
        unsigned int n_h (Rcpp::as<unsigned int>( s_n_h));
        variogram_evaluation::EmpiricalVariogram emp_vario(coords, distance, n_h, distanceTplane, std::make_shared<const SpMat> (distanceMatrix));

        if(weight.isNotNull()) {
          Eigen::Map<Vec> weight(Rcpp::as<Eigen::Map<Vec>> (s_weight));
          emp_vario.set_weight(weight);
        }

        // Fitted vario
        vario_factory::VariogramFactory & vf(vario_factory::VariogramFactory::Instance());
        std::string variogram_type (Rcpp::as<std::string> (s_vario_model));
        std::unique_ptr<variogram_evaluation::FittedVariogram> the_variogram = vf.create(variogram_type);

        // Gamma matrix
        Eigen::MatrixXd gamma_matrix(Eigen::MatrixXd::Identity(N,N));

        // Design matrix
        design_matrix::registerDesignMatrices();
        design_matrix::DesignMatrixFactory& design_matrix_fac = design_matrix::DesignMatrixFactory::Instance();
        std::string model_name (Rcpp::as<std::string> (s_ts_model)); // (Additive, Coord1, Coord2, Intercept)
        std::unique_ptr<design_matrix::DesignMatrix> theDesign_matrix = design_matrix_fac.create(model_name);

        Eigen::MatrixXd design_matrix;
        if(X.isNotNull()) {
          Eigen::Map<Eigen::MatrixXd> X(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_X));
          design_matrix = theDesign_matrix->compute_design_matrix(coords, X);
        }
        else design_matrix = theDesign_matrix->compute_design_matrix(coords);

        // Model
        unsigned int n_covariates(design_matrix.cols());
        model_fit::Model model(std::make_shared<const Eigen::MatrixXd> (big_matrix_data_tspace), std::make_shared<const Eigen::MatrixXd> (design_matrix), n);
        model.update_model(gamma_matrix);
        Eigen::MatrixXd resMatrix(N, ((n+1)*n)/2);
        Eigen::MatrixXd beta(n_covariates, ((n+1)*n)/2);
        Eigen::MatrixXd beta_old(n_covariates, ((n+1)*n)/2);
        beta = model.get_beta(); // OK fino a qui
        std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
        beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, n);

        std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

        unsigned int num_iter(0);
        unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
        double tolerance(Rcpp::as<double> (s_tolerance));
        std::vector<MatrixXd> resVec(N);

        double tol = tolerance+1;

        while (num_iter < max_iter && tol > tolerance) {
          resMatrix = model.get_residuals();
          resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, n);

          emp_vario.update_emp_vario(resVec);
          the_variogram -> evaluate_par_fitted(emp_vario);

          gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix, N);
          beta_old_vec_matrices = beta_vec_matrices;

          model.update_model(gamma_matrix);
          beta = model.get_beta();
          beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, n);

          tol=0.0;
          for (size_t i=0; i<n_covariates; i++) {
            tol += distanceTplane.compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
          }

          num_iter++;
        }

        unsigned int n_hh(1000);
        Vec hh(n_hh);
        std::vector<double> emp_vario_values(emp_vario.get_card_h());
        emp_vario_values = emp_vario.get_hvec();

        hh.setLinSpaced(n_hh, 0, *std::max_element(emp_vario_values.begin(), emp_vario_values.end()));

        Vec fit_vario_values = the_variogram->get_vario_vec(hh, n_hh);

        // KRIGING
        // Map functions
        map_functions::exponentialMap expMap(distanceManifold);

        // New coordinates
        const Eigen::Map<Eigen::MatrixXd> new_coords_mat(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_new_coordinates));
        Coordinates new_coords(std::make_shared<const Eigen::MatrixXd>(new_coords_mat));
        unsigned int M = new_coords_mat.rows();

        // New Design matrix
        Eigen::MatrixXd new_design_matrix;
        if(X_new.isNotNull()) {
          Eigen::Map<Eigen::MatrixXd> X_new(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_X_new));
          new_design_matrix = theDesign_matrix->compute_design_matrix(new_coords, X_new);
        }
        else new_design_matrix = theDesign_matrix->compute_design_matrix(new_coords);

        std::vector<double> distanceVector(N);
        Vec ci(N);
        Vec lambda_vec(N);

        Eigen::LDLT<Eigen::MatrixXd> solver(N);
        solver.compute(gamma_matrix);

        unsigned int num_cov(beta_vec_matrices.size());
        Eigen::MatrixXd tmp(n,n);
        auto weighted_sum_beta = [&beta_vec_matrices, &num_cov, &tmp, n] (const Vec& design_matrix_row) { tmp.setZero(n,n); for (size_t j=0; j<num_cov; j++) tmp = tmp + beta_vec_matrices[j]*design_matrix_row(j); return tmp;};
        auto weighted_sum_residuals = [&resVec, &N, &tmp, n] (const Vec& lambda_vec) { tmp.setZero(n,n); for (size_t j=0; j<N; j++) tmp = tmp + resVec[j]*lambda_vec(j); return tmp;};

        Eigen::MatrixXd tplane_prediction(n,n);
        std::vector<Eigen::MatrixXd> manifold_prediction(M);

        for (size_t i=0; i<M; i++) {
          distanceVector = distance.create_distance_vector(coords, new_coords_mat.row(i));
          ci = the_variogram->get_covario_vec(distanceVector, M);
          lambda_vec = solver.solve(ci);
          tplane_prediction = weighted_sum_beta(new_design_matrix.row(i)) + weighted_sum_residuals(lambda_vec);
          manifold_prediction[i] = expMap.map2manifold(tplane_prediction);
        }

        Rcpp::List result = Rcpp::List::create(Rcpp::Named("beta") = beta_vec_matrices,
                             Rcpp::Named("fit_vario_values") = fit_vario_values,
                             Rcpp::Named("hh") = hh,
                             Rcpp::Named("gamma_matrix") = gamma_matrix,
                             Rcpp::Named("residuals") = resVec,
                             Rcpp::Named("emp_vario_values") = emp_vario.get_emp_vario_values(),
                             Rcpp::Named("h_vec") = emp_vario.get_hvec(),
                             Rcpp::Named("vario_parameters") = the_variogram->get_parameters(),
                             Rcpp::Named("iterations") = num_iter,
                             Rcpp::Named("prediction") = manifold_prediction);

        return Rcpp::wrap(result);
    }

}













  // DISTANCE

  // errore se GEODIST e più di 2 _coords

  // COORDINATES
  // latitudine prima coordinata
  // lat e long in gradi
