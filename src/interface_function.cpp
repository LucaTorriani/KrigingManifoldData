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
#include "IntrinsicMean.hpp"

extern "C"{

// CREATE MODEL
  RcppExport SEXP get_model (SEXP s_data_manifold, SEXP s_coordinates, SEXP s_X, SEXP s_Sigma,
    SEXP s_distance, SEXP s_manifold_metric, SEXP s_ts_metric, SEXP s_ts_model, SEXP s_vario_model, SEXP s_n_h,
    SEXP s_max_it, SEXP s_tolerance,
    SEXP s_weight_vario, SEXP s_distance_matrix_tot, SEXP s_matrix_data_tspace_tot, SEXP s_coordinates_tot, SEXP s_X_tot, SEXP s_hmax, // RDD
    SEXP s_weight_intrinsic, SEXP s_tolerance_intrinsic) {  // Tangent point

      BEGIN_RCPP
      Rcpp::Nullable<Eigen::MatrixXd> X(s_X);
      Rcpp::Nullable<Eigen::MatrixXd> X_tot(s_X_tot);
      Rcpp::Nullable<Eigen::MatrixXd> Sigma_n(s_Sigma);
      Rcpp::Nullable<Vec> weight_vario(s_weight_vario);

      // Coordinates model
      std::shared_ptr<const Eigen::MatrixXd> coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_coordinates));
      Coordinates coords(coords_ptr);
      unsigned int N = coords.get_N_station();

      // Data manifold model
      std::vector<Eigen::MatrixXd> data_manifold(N);
      for(size_t i=0; i<N; i++){
        data_manifold[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold,i));
      }
      unsigned int p = data_manifold[0].rows();

      // Distance tplane
      std::string distance_Tplane_name = Rcpp::as<std::string> (s_ts_metric) ; //(Frobenius, FrobeniusScaled)
      tplane_factory::TplaneFactory& tplane_fac (tplane_factory::TplaneFactory::Instance());
      std::unique_ptr<distances_tplane::DistanceTplane> theTplaneDist = tplane_fac.create(distance_Tplane_name);

      // Map functions
      std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
      map_factory::LogMapFactory& logmap_fac (map_factory::LogMapFactory::Instance());
      std::unique_ptr<map_functions::logarithmicMap> theLogMap = logmap_fac.create(distance_Manifold_name);

      // Punto tangente
      Eigen::MatrixXd Sigma(p,p);
      if(Sigma_n.isNotNull()) {
        Sigma = Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_Sigma);
        theTplaneDist->set_members(Sigma);
        theLogMap->set_members(Sigma);
      }
      else {
        double tolerance_intrinsic(Rcpp::as<double> (s_tolerance_intrinsic));
        Eigen::Map<Vec> weights_intrinsic(Rcpp::as<Eigen::Map<Vec>> (s_weight_intrinsic));
        map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
        std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);
        Sigma = intrinsic_mean_C(data_manifold, *theLogMap, *theExpMap, *theTplaneDist, tolerance_intrinsic, weights_intrinsic);
      }

      // Distance
      distance_factory::DistanceFactory& distance_fac (distance_factory::DistanceFactory::Instance());
      std::string distance_name( Rcpp::as<std::string> (s_distance)) ; //(Geodist, Eucldist)
      std::unique_ptr<distances::Distance> theDistance = distance_fac.create(distance_name);

      // Distance Matrix
      std::shared_ptr<const MatrixXd> distanceMatrix_ptr = theDistance->create_distance_matrix(coords, N);

      // Fitted vario
      vario_factory::VariogramFactory & vf(vario_factory::VariogramFactory::Instance());
      std::string variogram_type (Rcpp::as<std::string> (s_vario_model));
      std::unique_ptr<variogram_evaluation::FittedVariogram> the_variogram = vf.create(variogram_type);

      // Gamma matrix
      Eigen::MatrixXd gamma_matrix(Eigen::MatrixXd::Identity(N,N));

      // Design matrix
      design_factory::DesignFactory& design_matrix_fac (design_factory::DesignFactory::Instance());
      std::string model_name (Rcpp::as<std::string> (s_ts_model)); // (Additive, Coord1, Coord2, Intercept)
      std::unique_ptr<design_matrix::DesignMatrix> theDesign_matrix = design_matrix_fac.create(model_name);

      std::shared_ptr<Eigen::MatrixXd> design_matrix_ptr;
      if(X.isNotNull()) {
        Eigen::Map<Eigen::MatrixXd> X(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_X));
        design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords, X));
      }
      else design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords));

      unsigned int n_covariates(design_matrix_ptr->cols());

      // KERNEL
      if(weight_vario.isNotNull()) {
        // Rcpp::Rcout<<"Using kernel weights for variogram"<<"\n";
        // Weight vario
        Eigen::Map<Vec> weight_vario(Rcpp::as<Eigen::Map<Vec>> (s_weight_vario));

        // Distance Matrix tot
        std::shared_ptr<const Eigen::MatrixXd> distanceMatrix_tot_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_distance_matrix_tot));

        // Coordinates tot
        std::shared_ptr<const Eigen::MatrixXd> coords_tot_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_coordinates_tot));
        Coordinates coords_tot(coords_tot_ptr);
        unsigned int N_tot = coords_tot.get_N_station();

        // Data tspace tot
        std::shared_ptr<const Eigen::MatrixXd> big_matrix_data_tspace_tot_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_matrix_data_tspace_tot));

        // Emp vario
        unsigned int n_h (Rcpp::as<unsigned int>(s_n_h));
        double hmax (Rcpp::as<double>(s_hmax));
        variogram_evaluation::EmpiricalVariogram emp_vario(distanceMatrix_tot_ptr, n_h, N_tot, weight_vario, hmax);

        // Design matrix tot
        std::shared_ptr<Eigen::MatrixXd> design_matrix_tot_ptr;
        if(X_tot.isNotNull()) {
          Eigen::Map<Eigen::MatrixXd> X_tot(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_X_tot));
          design_matrix_tot_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords_tot, X_tot));
        }
        else design_matrix_tot_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords_tot));

        // Model
        model_fit::Model model(big_matrix_data_tspace_tot_ptr, design_matrix_ptr, design_matrix_tot_ptr, p);

        Eigen::MatrixXd resMatrix(N_tot, ((p+1)*p)/2);
        std::vector<MatrixXd> resVec(N_tot);

        // DA QUI è TUTTO UGUALE
        model.update_model(gamma_matrix);

        Eigen::MatrixXd beta(n_covariates, ((p+1)*p)/2);
        Eigen::MatrixXd beta_old(n_covariates, ((p+1)*p)/2);

        beta = model.get_beta();
        std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
        beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, p);
        std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

        unsigned int num_iter(0);
        unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
        double tolerance(Rcpp::as<double> (s_tolerance));

        double tol = tolerance+1;

        while (num_iter < max_iter && tol > tolerance) {
          resMatrix = model.get_residuals();
          resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, p);

          emp_vario.update_emp_vario(resVec, *(theTplaneDist));
          the_variogram -> evaluate_par_fitted(emp_vario);

          gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix_ptr, N);
          beta_old_vec_matrices = beta_vec_matrices;

          model.update_model(gamma_matrix);
          beta = model.get_beta();
          beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, p);

          tol=0.0;
          for (size_t i=0; i<n_covariates; i++) {
            tol += theTplaneDist->compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
          }
          num_iter++;
        }
        if(num_iter == max_iter) Rcpp::warning("Reached max number of iterations");

        unsigned int n_hh(1000);
        Vec hh(n_hh);
        std::vector<double> h_vario_values(emp_vario.get_card_h());
        h_vario_values = emp_vario.get_hvec();

        hh.setLinSpaced(n_hh, 0, *std::max_element(h_vario_values.begin(), h_vario_values.end()));

        Vec fit_vario_values = the_variogram->get_vario_vec(hh, n_hh);

        Rcpp::List result = Rcpp::List::create(Rcpp::Named("beta") = beta_vec_matrices,
                             Rcpp::Named("fit_vario_values") = fit_vario_values,
                             Rcpp::Named("hh") = hh,
                             Rcpp::Named("gamma_matrix") = gamma_matrix,
                             Rcpp::Named("residuals") = resVec,
                             Rcpp::Named("emp_vario_values") = emp_vario.get_emp_vario_values(),
                             Rcpp::Named("h_vec") = h_vario_values,
                             Rcpp::Named("fitted_par_vario") = the_variogram->get_parameters(),
                             Rcpp::Named("iterations") = num_iter,
                             Rcpp::Named("Sigma")= Sigma);

        return Rcpp::wrap(result);
        // return Rcpp::wrap(1);
      }
      else {  // EQUAL WEIGHTS
        // Rcpp::Rcout<<"Using equal weights for variogram"<<"\n";

        // Data tangent space
        std::vector<Eigen::MatrixXd> data_tspace(N);
        for (size_t i=0; i<N; i++) {
          data_tspace[i] = theLogMap->map2tplane(data_manifold[i]);
        }

        std::shared_ptr<const Eigen::MatrixXd> big_matrix_data_tspace_ptr = std::make_shared<const Eigen::MatrixXd>(matrix_manipulation::VecMatrices2bigMatrix(data_tspace));
        data_manifold.clear();
        data_tspace.clear();

        // Emp vario
        unsigned int n_h (Rcpp::as<unsigned int>( s_n_h));
        variogram_evaluation::EmpiricalVariogram emp_vario(distanceMatrix_ptr, n_h, coords, *(theDistance));

        // Model
        model_fit::Model model(big_matrix_data_tspace_ptr, design_matrix_ptr, p);

        Eigen::MatrixXd resMatrix(N, ((p+1)*p)/2);
        std::vector<MatrixXd> resVec(N);

        // DA QUI è TUTTO UGUALE
        model.update_model(gamma_matrix);

        Eigen::MatrixXd beta(n_covariates, ((p+1)*p)/2);
        Eigen::MatrixXd beta_old(n_covariates, ((p+1)*p)/2);

        beta = model.get_beta();
        std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
        beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, p);
        std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

        unsigned int num_iter(0);
        unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
        double tolerance(Rcpp::as<double> (s_tolerance));

        double tol = tolerance+1;

        while (num_iter < max_iter && tol > tolerance) {
          resMatrix = model.get_residuals();
          resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, p);

          emp_vario.update_emp_vario(resVec, *(theTplaneDist));
          the_variogram -> evaluate_par_fitted(emp_vario);

          gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix_ptr, N);
          beta_old_vec_matrices = beta_vec_matrices;

          model.update_model(gamma_matrix);
          beta = model.get_beta();
          beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, p);

          tol=0.0;
          for (size_t i=0; i<n_covariates; i++) {
            tol += theTplaneDist->compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
          }
          num_iter++;
        }
        if(num_iter == max_iter) Rcpp::warning("Reached max number of iterations");

        unsigned int n_hh(1000);
        Vec hh(n_hh);
        std::vector<double> h_vario_values(emp_vario.get_card_h());
        h_vario_values = emp_vario.get_hvec();

        hh.setLinSpaced(n_hh, 0, *std::max_element(h_vario_values.begin(), h_vario_values.end()));

        Vec fit_vario_values = the_variogram->get_vario_vec(hh, n_hh);

        Rcpp::List result = Rcpp::List::create(Rcpp::Named("beta") = beta_vec_matrices,
                             Rcpp::Named("fit_vario_values") = fit_vario_values,
                             Rcpp::Named("hh") = hh,
                             Rcpp::Named("gamma_matrix") = gamma_matrix,
                             Rcpp::Named("residuals") = resVec,
                             Rcpp::Named("emp_vario_values") = emp_vario.get_emp_vario_values(),
                             Rcpp::Named("h_vec") = h_vario_values,
                             Rcpp::Named("fitted_par_vario") = the_variogram->get_parameters(),
                             Rcpp::Named("iterations") = num_iter,
                             Rcpp::Named("Sigma")= Sigma);

        return Rcpp::wrap(result);
      }
      END_RCPP
  }


// KRIGING
  RcppExport SEXP get_kriging (SEXP s_coordinates, SEXP s_new_coordinates,  SEXP s_Sigma,
    SEXP s_distance, SEXP s_manifold_metric, SEXP s_ts_model, SEXP s_vario_model,
    SEXP s_beta, SEXP s_gamma_matrix, SEXP s_vario_parameters, SEXP s_residuals, SEXP s_X_new) {

    BEGIN_RCPP

    Rcpp::Nullable<Eigen::MatrixXd> X_new(s_X_new);

    // Punto tangente
    Eigen::Map<Eigen::MatrixXd> Sigma(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_Sigma));
    unsigned int p = Sigma.rows();

    // Distance
    distance_factory::DistanceFactory& distance_fac (distance_factory::DistanceFactory::Instance());
    std::string distance_name( Rcpp::as<std::string> (s_distance)) ; //(Geodist, Eucldist)
    std::unique_ptr<distances::Distance> theDistance = distance_fac.create(distance_name);

    // Distance manifold SERVE???
    std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)

    // Map functions
    map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
    std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);
    theExpMap->set_members(Sigma);

    // Old coordinates
    std::shared_ptr<const Eigen::MatrixXd> coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_coordinates));
    unsigned int N = coords_ptr->rows();
    Coordinates coords(coords_ptr);

    // New coordinates
    std::shared_ptr<const Eigen::MatrixXd> new_coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_new_coordinates));
    unsigned int M = new_coords_ptr->rows();
    Coordinates new_coords(new_coords_ptr);

    // New Design matrix
    design_factory::DesignFactory& design_matrix_fac (design_factory::DesignFactory::Instance());
    std::string model_name (Rcpp::as<std::string> (s_ts_model)); // (Additive, Coord1, Coord2, Intercept)
    std::unique_ptr<design_matrix::DesignMatrix> theDesign_matrix = design_matrix_fac.create(model_name);

    std::shared_ptr<Eigen::MatrixXd> new_design_matrix_ptr;
    if(X_new.isNotNull()) {
      Eigen::Map<Eigen::MatrixXd> X_new(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_X_new));
      new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords, X_new));
    }
    else new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords));

    // Fitted vario
    Eigen::Map<Eigen::VectorXd> parameters(Rcpp::as<Eigen::Map<Eigen::VectorXd>> (s_vario_parameters));

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

    Eigen::MatrixXd tmp(p,p);
    auto weighted_sum_beta = [&beta_vec, &num_cov, &tmp, p] (const Vec& design_matrix_row) { tmp.setZero(p,p); for (size_t j=0; j<num_cov; j++) tmp = tmp + beta_vec[j]*design_matrix_row(j); return tmp;};
    auto weighted_sum_residuals = [&residuals_vec, &N, &tmp, p] (const Vec& lambda_vec) { tmp.setZero(p,p); for (size_t j=0; j<N; j++) tmp = tmp + residuals_vec[j]*lambda_vec(j); return tmp;};

    Eigen::MatrixXd tplane_prediction(p,p);
    std::vector<Eigen::MatrixXd> manifold_prediction(M);

    for (size_t i=0; i<M; i++) {
      distanceVector = theDistance->create_distance_vector(coords, new_coords_ptr->row(i));
      ci = the_variogram->get_covario_vec(distanceVector, N);
      lambda_vec = solver.solve(ci);
      tplane_prediction = weighted_sum_beta(new_design_matrix_ptr->row(i)) + weighted_sum_residuals(lambda_vec);
      manifold_prediction[i] = theExpMap->map2manifold(tplane_prediction);
    }

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("prediction") = manifold_prediction);  // Solo questo?

    return Rcpp::wrap(result);
    END_RCPP
  }


  // CREATE MODEL AND KRIGNG
  RcppExport SEXP get_model_and_kriging (SEXP s_data_manifold, SEXP s_coordinates, SEXP s_X, SEXP s_Sigma,
    SEXP s_distance, SEXP s_manifold_metric, SEXP s_ts_metric, SEXP s_ts_model, SEXP s_vario_model, SEXP s_n_h,
    SEXP s_max_it, SEXP s_tolerance,
    SEXP s_weight_vario, SEXP s_distance_matrix_tot, SEXP s_data_manifold_tot, SEXP s_coordinates_tot, SEXP s_X_tot, SEXP s_hmax, SEXP s_indexes_model, // RDD
    SEXP s_weight_intrinsic, SEXP s_tolerance_intrinsic,
    SEXP s_new_coordinates, SEXP s_X_new) {  // KRIGING

      BEGIN_RCPP
      Rcpp::Nullable<Eigen::MatrixXd> X(s_X);
      Rcpp::Nullable<Eigen::MatrixXd> X_tot(s_X_tot);
      Rcpp::Nullable<Eigen::MatrixXd> Sigma_n(s_Sigma);
      Rcpp::Nullable<Vec> weight_vario(s_weight_vario);
      Rcpp::Nullable<Eigen::MatrixXd> X_new(s_X_new);

      // Coordinates model
      std::shared_ptr<const Eigen::MatrixXd> coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_coordinates));
      Coordinates coords(coords_ptr);
      unsigned int N = coords.get_N_station();

      // Data manifold model
      std::vector<Eigen::MatrixXd> data_manifold(N);
      for(size_t i=0; i<N; i++){
        data_manifold[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold,i));
      }
      unsigned int p = data_manifold[0].rows();

      // Distance tplane
      std::string distance_Tplane_name = Rcpp::as<std::string> (s_ts_metric) ; //(Frobenius, FrobeniusScaled)
      tplane_factory::TplaneFactory& tplane_fac (tplane_factory::TplaneFactory::Instance());
      std::unique_ptr<distances_tplane::DistanceTplane> theTplaneDist = tplane_fac.create(distance_Tplane_name);

      // Map functions
      std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
      map_factory::LogMapFactory& logmap_fac (map_factory::LogMapFactory::Instance());
      std::unique_ptr<map_functions::logarithmicMap> theLogMap = logmap_fac.create(distance_Manifold_name);

      // Punto tangente
      Eigen::MatrixXd Sigma(p,p);
      if(Sigma_n.isNotNull()) {
        Sigma = Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_Sigma);
        theTplaneDist->set_members(Sigma);
        theLogMap->set_members(Sigma);
      }
      else {
        double tolerance_intrinsic(Rcpp::as<double> (s_tolerance_intrinsic));
        Eigen::Map<Vec> weights_intrinsic(Rcpp::as<Eigen::Map<Vec>> (s_weight_intrinsic));
        map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
        std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);
        Sigma = intrinsic_mean_C(data_manifold, *theLogMap, *theExpMap, *theTplaneDist, tolerance_intrinsic, weights_intrinsic);
      }

      // Distance
      distance_factory::DistanceFactory& distance_fac (distance_factory::DistanceFactory::Instance());
      std::string distance_name( Rcpp::as<std::string> (s_distance)) ; //(Geodist, Eucldist)
      std::unique_ptr<distances::Distance> theDistance = distance_fac.create(distance_name);

      // Distance Matrix
      std::shared_ptr<const MatrixXd> distanceMatrix_ptr = theDistance->create_distance_matrix(coords, N);

      // Fitted vario
      vario_factory::VariogramFactory & vf(vario_factory::VariogramFactory::Instance());
      std::string variogram_type (Rcpp::as<std::string> (s_vario_model));
      std::unique_ptr<variogram_evaluation::FittedVariogram> the_variogram = vf.create(variogram_type);

      // Gamma matrix
      Eigen::MatrixXd gamma_matrix(Eigen::MatrixXd::Identity(N,N));

      // Design matrix
      design_factory::DesignFactory& design_matrix_fac (design_factory::DesignFactory::Instance());
      std::string model_name (Rcpp::as<std::string> (s_ts_model)); // (Additive, Coord1, Coord2, Intercept)
      std::unique_ptr<design_matrix::DesignMatrix> theDesign_matrix = design_matrix_fac.create(model_name);

      std::shared_ptr<Eigen::MatrixXd> design_matrix_ptr;
      if(X.isNotNull()) {
        Eigen::Map<Eigen::MatrixXd> X(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_X));
        design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords, X));
      }
      else design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords));

      unsigned int n_covariates(design_matrix_ptr->cols());

      // KERNEL
      if(weight_vario.isNotNull()) {
        data_manifold.clear();
        // Weight vario
        Eigen::Map<Vec> weight_vario(Rcpp::as<Eigen::Map<Vec>> (s_weight_vario));

        // Distance Matrix tot
        std::shared_ptr<const Eigen::MatrixXd> distanceMatrix_tot_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_distance_matrix_tot));

        // Coordinates tot
        std::shared_ptr<const Eigen::MatrixXd> coords_tot_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_coordinates_tot));
        Coordinates coords_tot(coords_tot_ptr);
        unsigned int N_tot = coords_tot.get_N_station();

        // Data manifold tot
        std::vector<Eigen::MatrixXd> data_manifold_tot(N_tot);
        for(size_t i=0; i<N_tot; i++){
          data_manifold_tot[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold_tot,i));
        }

        // Data tangent space tot
        std::vector<Eigen::MatrixXd> data_tspace_tot(N_tot);
        for (size_t i=0; i<N_tot; i++) {
          data_tspace_tot[i] = theLogMap->map2tplane(data_manifold_tot[i]);
        }

        std::shared_ptr<const Eigen::MatrixXd> big_matrix_data_tspace_tot_ptr = std::make_shared<const Eigen::MatrixXd>(matrix_manipulation::VecMatrices2bigMatrix(data_tspace_tot));
        data_manifold_tot.clear();
        data_tspace_tot.clear();


        // Emp vario
        unsigned int n_h (Rcpp::as<unsigned int>(s_n_h));
        double hmax (Rcpp::as<double>(s_hmax));
        variogram_evaluation::EmpiricalVariogram emp_vario(distanceMatrix_tot_ptr, n_h, N_tot, weight_vario, hmax);

        // Design matrix tot
        std::shared_ptr<Eigen::MatrixXd> design_matrix_tot_ptr;
        if(X_tot.isNotNull()) {
          Eigen::Map<Eigen::MatrixXd> X_tot(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_X_tot));
          design_matrix_tot_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords_tot, X_tot));
        }
        else design_matrix_tot_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords_tot));

        // Model
        model_fit::Model model(big_matrix_data_tspace_tot_ptr, design_matrix_ptr, design_matrix_tot_ptr, p);

        Eigen::MatrixXd resMatrix(N_tot, ((p+1)*p)/2);
        std::vector<MatrixXd> resVec(N_tot);

        // DA QUI è TUTTO UGUALE tranne la parte ("Select residuals in the k-th cell")
        model.update_model(gamma_matrix);

        Eigen::MatrixXd beta(n_covariates, ((p+1)*p)/2);
        Eigen::MatrixXd beta_old(n_covariates, ((p+1)*p)/2);

        beta = model.get_beta();
        std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
        beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, p);
        std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

        unsigned int num_iter(0);
        unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
        double tolerance(Rcpp::as<double> (s_tolerance));

        double tol = tolerance+1;

        while (num_iter < max_iter && tol > tolerance) {
          resMatrix = model.get_residuals();
          resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, p);

          emp_vario.update_emp_vario(resVec, *(theTplaneDist));
          the_variogram -> evaluate_par_fitted(emp_vario);

          gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix_ptr, N);
          beta_old_vec_matrices = beta_vec_matrices;

          model.update_model(gamma_matrix);
          beta = model.get_beta();
          beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, p);

          tol=0.0;
          for (size_t i=0; i<n_covariates; i++) {
            tol += theTplaneDist->compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
          }
          num_iter++;
        }
        if(num_iter == max_iter) Rcpp::warning("Reached max number of iterations");

        unsigned int n_hh(1000);
        Vec hh(n_hh);
        std::vector<double> h_vario_values(emp_vario.get_card_h());
        h_vario_values = emp_vario.get_hvec();

        hh.setLinSpaced(n_hh, 0, *std::max_element(h_vario_values.begin(), h_vario_values.end()));

        Vec fit_vario_values = the_variogram->get_vario_vec(hh, n_hh);

        // KRIGING

        // Map functions
        map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
        std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);
        theExpMap->set_members(Sigma);

        // New coordinates
        std::shared_ptr<const Eigen::MatrixXd> new_coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_new_coordinates));
        unsigned int M = new_coords_ptr->rows();
        Coordinates new_coords(new_coords_ptr);

        // New Design matrix
        std::shared_ptr<Eigen::MatrixXd> new_design_matrix_ptr;
        if(X_new.isNotNull()) {
          Eigen::Map<Eigen::MatrixXd> X_new(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_X_new));
          new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords, X_new));
        }
        else new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords));

        std::vector<double> distanceVector(N);

        Vec ci(N);
        Vec lambda_vec(N);

        Eigen::LDLT<Eigen::MatrixXd> solver(N);
        solver.compute(gamma_matrix);

        unsigned int num_cov(beta_vec_matrices.size());
        Eigen::MatrixXd tmp(p,p);

        // Select residuals in the k-th cell
        std::vector<unsigned int> indexes_model(Rcpp::as<std::vector<unsigned int>> (s_indexes_model));

        std::vector<Eigen::MatrixXd> resVec_k(N);
        for (size_t ii=0; ii<N; ii++ ) {
            resVec_k[ii]=resVec[indexes_model[ii]-1];
        }
        //

        auto weighted_sum_beta = [&beta_vec_matrices, &num_cov, &tmp, p] (const Vec& design_matrix_row) { tmp.setZero(p,p); for (size_t j=0; j<num_cov; j++) tmp = tmp + beta_vec_matrices[j]*design_matrix_row(j); return tmp;};
        auto weighted_sum_residuals = [&resVec_k, &N, &tmp, p] (const Vec& lambda_vec) { tmp.setZero(p,p); for (size_t j=0; j<N; j++) tmp = tmp + resVec_k[j]*lambda_vec(j); return tmp;};

        Eigen::MatrixXd tplane_prediction(p,p);
        std::vector<Eigen::MatrixXd> manifold_prediction(M);

        for (size_t i=0; i<M; i++) {
          distanceVector = theDistance->create_distance_vector(coords, new_coords_ptr->row(i));
          ci = the_variogram->get_covario_vec(distanceVector, N);
          lambda_vec = solver.solve(ci);
          tplane_prediction = weighted_sum_beta(new_design_matrix_ptr->row(i)) + weighted_sum_residuals(lambda_vec);
          manifold_prediction[i] = theExpMap->map2manifold(tplane_prediction);
        }


        Rcpp::List result = Rcpp::List::create(Rcpp::Named("beta") = beta_vec_matrices,
                                 Rcpp::Named("fit_vario_values") = fit_vario_values,
                                 Rcpp::Named("hh") = hh,
                                 Rcpp::Named("gamma_matrix") = gamma_matrix,
                                 Rcpp::Named("residuals") = resVec,
                                 Rcpp::Named("emp_vario_values") = emp_vario.get_emp_vario_values(),
                                 Rcpp::Named("h_vec") = emp_vario.get_hvec(),
                                 Rcpp::Named("fitted_par_vario") = the_variogram->get_parameters(),
                                 Rcpp::Named("iterations") = num_iter,
                                 Rcpp::Named("Sigma") = Sigma,
                                 Rcpp::Named("prediction") = manifold_prediction);

        return Rcpp::wrap(result);

      }
      else {  // EQUAL WEIGHTS
        // Data tangent space
        std::vector<Eigen::MatrixXd> data_tspace(N);
        for (size_t i=0; i<N; i++) {
          data_tspace[i] = theLogMap->map2tplane(data_manifold[i]);
        }

        std::shared_ptr<const Eigen::MatrixXd> big_matrix_data_tspace_ptr = std::make_shared<const Eigen::MatrixXd>(matrix_manipulation::VecMatrices2bigMatrix(data_tspace));
        data_manifold.clear();
        data_tspace.clear();

        // Emp vario
        unsigned int n_h (Rcpp::as<unsigned int>( s_n_h));
        variogram_evaluation::EmpiricalVariogram emp_vario(distanceMatrix_ptr, n_h, coords, *(theDistance));

        // Model
        model_fit::Model model(big_matrix_data_tspace_ptr, design_matrix_ptr, p);

        Eigen::MatrixXd resMatrix(N, ((p+1)*p)/2);
        std::vector<MatrixXd> resVec(N);

        // DA QUI è TUTTO UGUALE
        model.update_model(gamma_matrix);

        Eigen::MatrixXd beta(n_covariates, ((p+1)*p)/2);
        Eigen::MatrixXd beta_old(n_covariates, ((p+1)*p)/2);

        beta = model.get_beta();
        std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
        beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, p);
        std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

        unsigned int num_iter(0);
        unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
        double tolerance(Rcpp::as<double> (s_tolerance));

        double tol = tolerance+1;

        while (num_iter < max_iter && tol > tolerance) {
          resMatrix = model.get_residuals();
          resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, p);

          emp_vario.update_emp_vario(resVec, *(theTplaneDist));
          the_variogram -> evaluate_par_fitted(emp_vario);

          gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix_ptr, N);
          beta_old_vec_matrices = beta_vec_matrices;

          model.update_model(gamma_matrix);
          beta = model.get_beta();
          beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, p);

          tol=0.0;
          for (size_t i=0; i<n_covariates; i++) {
            tol += theTplaneDist->compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
          }
          num_iter++;
        }
        if(num_iter == max_iter) Rcpp::warning("Reached max number of iterations");

        unsigned int n_hh(1000);
        Vec hh(n_hh);
        std::vector<double> h_vario_values(emp_vario.get_card_h());
        h_vario_values = emp_vario.get_hvec();

        hh.setLinSpaced(n_hh, 0, *std::max_element(h_vario_values.begin(), h_vario_values.end()));

        Vec fit_vario_values = the_variogram->get_vario_vec(hh, n_hh);

        // KRIGING

        // Map functions
        map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
        std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);
        theExpMap->set_members(Sigma);

        // New coordinates
        std::shared_ptr<const Eigen::MatrixXd> new_coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_new_coordinates));
        unsigned int M = new_coords_ptr->rows();
        Coordinates new_coords(new_coords_ptr);

        // New Design matrix
        std::shared_ptr<Eigen::MatrixXd> new_design_matrix_ptr;
        if(X_new.isNotNull()) {
          Eigen::Map<Eigen::MatrixXd> X_new(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_X_new));
          new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords, X_new));
        }
        else new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords));

        std::vector<double> distanceVector(N);

        Vec ci(N);
        Vec lambda_vec(N);

        Eigen::LDLT<Eigen::MatrixXd> solver(N);
        solver.compute(gamma_matrix);

        unsigned int num_cov(beta_vec_matrices.size());
        Eigen::MatrixXd tmp(p,p);
        auto weighted_sum_beta = [&beta_vec_matrices, &num_cov, &tmp, p] (const Vec& design_matrix_row) { tmp.setZero(p,p); for (size_t j=0; j<num_cov; j++) tmp = tmp + beta_vec_matrices[j]*design_matrix_row(j); return tmp;};
        auto weighted_sum_residuals = [&resVec, &N, &tmp, p] (const Vec& lambda_vec) { tmp.setZero(p,p); for (size_t j=0; j<N; j++) tmp = tmp + resVec[j]*lambda_vec(j); return tmp;};

        Eigen::MatrixXd tplane_prediction(p,p);
        std::vector<Eigen::MatrixXd> manifold_prediction(M);

        for (size_t i=0; i<M; i++) {
          distanceVector = theDistance->create_distance_vector(coords, new_coords_ptr->row(i));
          ci = the_variogram->get_covario_vec(distanceVector, N);
          lambda_vec = solver.solve(ci);
          tplane_prediction = weighted_sum_beta(new_design_matrix_ptr->row(i)) + weighted_sum_residuals(lambda_vec);
          manifold_prediction[i] = theExpMap->map2manifold(tplane_prediction);
        }


        Rcpp::List result = Rcpp::List::create(Rcpp::Named("beta") = beta_vec_matrices,
                                 Rcpp::Named("fit_vario_values") = fit_vario_values,
                                 Rcpp::Named("hh") = hh,
                                 Rcpp::Named("gamma_matrix") = gamma_matrix,
                                 Rcpp::Named("residuals") = resVec,
                                 Rcpp::Named("emp_vario_values") = emp_vario.get_emp_vario_values(),
                                 Rcpp::Named("h_vec") = emp_vario.get_hvec(),
                                 Rcpp::Named("fitted_par_vario") = the_variogram->get_parameters(),
                                 Rcpp::Named("iterations") = num_iter,
                                 Rcpp::Named("Sigma") = Sigma,
                                 Rcpp::Named("prediction") = manifold_prediction);

        return Rcpp::wrap(result);
      }
    END_RCPP
}


// INTRINSIC MEAN
RcppExport SEXP intrinsic_mean (SEXP s_data, SEXP s_N, SEXP s_manifold_metric, SEXP s_ts_metric,
  SEXP s_tolerance, SEXP s_weight) {
    BEGIN_RCPP
    // Data
    unsigned int N(Rcpp::as<unsigned int> (s_N));
    std::vector<Eigen::MatrixXd> data(N);
    for(size_t i=0; i<N; i++){
      data[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data,i));
    }
    unsigned int p = data[0].rows();

    // Distance tplane
    std::string distance_Tplane_name = Rcpp::as<std::string> (s_ts_metric) ; //(Frobenius, FrobeniusScaled)
    tplane_factory::TplaneFactory& tplane_fac (tplane_factory::TplaneFactory::Instance());
    std::unique_ptr<distances_tplane::DistanceTplane> theTplaneDist = tplane_fac.create(distance_Tplane_name);

    // Map functions
    std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
    map_factory::LogMapFactory& logmap_fac (map_factory::LogMapFactory::Instance());
    std::unique_ptr<map_functions::logarithmicMap> theLogMap = logmap_fac.create(distance_Manifold_name);
    map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
    std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);

    // Tolerance
    double tolerance (Rcpp::as<double> (s_tolerance));

    // Weights
    Eigen::Map<Vec> weight(Rcpp::as<Eigen::Map<Vec>> (s_weight));
    double sum_weight(weight.sum());

    // CODE
    Eigen::MatrixXd Result((data)[0]);

    theLogMap->set_members(Result);
    theExpMap->set_members(Result);
    theTplaneDist->set_members(Result);

    Eigen::MatrixXd Xk(p,p);
    Eigen::MatrixXd Xk_prec(p,p);

    Xk = data[0];

    double tau(1.0);
    double tmp;
    double tolk;
    double tolk_prec(tolerance + 1);

    size_t num_iter(0);

    while(tolk_prec > tolerance && num_iter < 100) {

      Xk_prec = Xk;
      tmp = theTplaneDist->norm(Xk_prec);
      tolk_prec = tmp*tmp;

      Xk.setZero();
      for (size_t i=0; i<N; i++) {
        Xk = Xk + weight(i)* theLogMap->map2tplane(data[i]);
      }
      Xk = Xk/sum_weight;
      Result = theExpMap->map2manifold(tau*Xk);

      theLogMap->set_members(Result);
      theExpMap->set_members(Result);
      theTplaneDist->set_members(Result);

      tmp = theTplaneDist->norm(Xk);
      tolk = tmp*tmp;
      if (tolk > tolk_prec) {
        tau = tau/2;
        Xk = Xk_prec;
      }
      num_iter++;
    }
    if(num_iter == 100) Rcpp::warning("Reached max number of iterations in intrinsic_mean");

    return Rcpp::wrap(Result);
    END_RCPP
}


// MAP TO TANGENT SPACE
// RcppExport SEXP map2_tangent_space (SEXP s_data_manifold, SEXP s_N, SEXP s_manifold_metric, SEXP s_Sigma) {
//      BEGIN_RCPP
//
//      // Data manifold model
//      unsigned int N(Rcpp::as<unsigned int>(s_N));
//
//      std::vector<Eigen::MatrixXd> data_manifold(N);
//      for(size_t i=0; i<N; i++){
//        data_manifold[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold,i));
//      }
//
//      // Map functions
//      std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
//      map_factory::LogMapFactory& logmap_fac (map_factory::LogMapFactory::Instance());
//      std::unique_ptr<map_functions::logarithmicMap> theLogMap = logmap_fac.create(distance_Manifold_name);
//      Eigen::Map<Eigen::MatrixXd> Sigma(Rcpp::as<Eigen::Map<Eigen::MatrixXd>> (s_Sigma));
//      theLogMap->set_members(Sigma);
//
//      // Data tangent space
//      std::vector<Eigen::MatrixXd> data_tspace(N);
//      for (size_t i=0; i<N; i++) {
//        data_tspace[i] = theLogMap->map2tplane(data_manifold[i]);
//      }
//      return Rcpp::wrap(matrix_manipulation::VecMatrices2bigMatrix(data_tspace));
//
//      END_RCPP
//   }

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
