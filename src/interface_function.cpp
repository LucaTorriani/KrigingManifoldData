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
#include "ParallelTransport.hpp"

extern "C"{

  // CREATE MODEL AND KRIGNG
  RcppExport SEXP get_model_and_kriging (SEXP s_data_manifold, SEXP s_coordinates, SEXP s_X, SEXP s_Sigma_data,
     SEXP s_distance, SEXP s_manifold_metric,  SEXP s_ts_model, SEXP s_vario_model, SEXP s_n_h, // SEXP s_ts_metric,
     SEXP s_max_it, SEXP s_tolerance, SEXP s_max_sill, SEXP s_max_a, // SEXP s_weight_vario, SEXP s_weight_intrinsic, SEXP s_tolerance_intrinsic,
     SEXP s_new_coordinates, SEXP s_Sigma_new, SEXP s_X_new, SEXP s_suppressMes) {

    BEGIN_RCPP
    Rcpp::Rcout << "Inizio " << "\n";

    // Rcpp::Nullable<Vec> weight_vario(s_weight_vario);
    Rcpp::Nullable<Eigen::MatrixXd> X(s_X);
    Rcpp::Nullable<Eigen::MatrixXd> X_new(s_X_new);
    Rcpp::Nullable<double> max_sill_n (s_max_sill);
    Rcpp::Nullable<double> max_a_n (s_max_a);

    // Coordinates
    std::shared_ptr<const Eigen::MatrixXd> coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_coordinates));
    Coordinates coords(coords_ptr);
    unsigned int N = coords.get_N_station();

    // Data manifold
    std::vector<Eigen::MatrixXd> data_manifold(N);
    for(size_t i=0; i<N; i++){
      data_manifold[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold,i));
    }

    unsigned int p = data_manifold[0].rows();
    Rcpp::Rcout << "QUI1 " << "\n";

    // Distance tplane
    // std::string distance_Tplane_name = Rcpp::as<std::string> (s_ts_metric) ; //(Frobenius, FrobeniusScaled)
    std::string distance_Tplane_name("Frobenius");
    tplane_factory::TplaneFactory& tplane_fac (tplane_factory::TplaneFactory::Instance());
    std::unique_ptr<distances_tplane::DistanceTplane> theTplaneDist = tplane_fac.create(distance_Tplane_name);

    // Map functions
    std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
    map_factory::LogMapFactory& logmap_fac (map_factory::LogMapFactory::Instance());
    std::unique_ptr<map_functions::logarithmicMap> theLogMap = logmap_fac.create(distance_Manifold_name);
    Rcpp::Rcout << "QUI2 " << "\n";

    // Tangent points
    Rcpp::List Sigma_data(s_Sigma_data);
    std::vector<Eigen::MatrixXd> Sigma_data_vec(N);
    for(size_t i=0; i<N; i++) Sigma_data_vec[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(Sigma_data,i));

    // Data tangent space
    std::vector<Eigen::MatrixXd> data_tspace(N);
    for (size_t i=0; i<N; i++) {
      theLogMap->set_members(Sigma_data_vec[i]);
      data_tspace[i] = parallel_transport::transport_to_TI(Sigma_data_vec[i], theLogMap->map2tplane(data_manifold[i]));
    }

    std::shared_ptr<const Eigen::MatrixXd> big_matrix_data_tspace_ptr = std::make_shared<const Eigen::MatrixXd>(matrix_manipulation::VecMatrices2bigMatrix(data_tspace));
    data_manifold.clear();
    data_tspace.clear();

    // Distance
    distance_factory::DistanceFactory& distance_fac (distance_factory::DistanceFactory::Instance());
    std::string distance_name( Rcpp::as<std::string> (s_distance)) ; //(Geodist, Eucldist)
    std::unique_ptr<distances::Distance> theDistance = distance_fac.create(distance_name);

    // Distance Matrix
    std::shared_ptr<const SpMat> distanceMatrix_ptr = theDistance->create_distance_matrix(coords, N);
    Rcpp::Rcout << "QUI3 " << "\n";

    // Emp vario
    unsigned int n_h (Rcpp::as<unsigned int>( s_n_h));
    variogram_evaluation::EmpiricalVariogram emp_vario(distanceMatrix_ptr, n_h, coords, *(theDistance));

    // if(weight_vario.isNotNull()) {
    //   Vec weight_vario(Rcpp::as<Vec> (s_weight_vario));
    //   emp_vario.set_weight(weight_vario);
    // }

    // Fitted vario parameters
    double max_a;
    double max_sill;
    if(max_a_n.isNotNull()) max_a = Rcpp::as<double> (s_max_a);
    if(max_sill_n.isNotNull()) max_sill= Rcpp::as<double> (s_max_sill);
    Rcpp::Rcout << "QUI4 " << "\n";

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
      Eigen::MatrixXd X(Rcpp::as<Eigen::MatrixXd> (s_X));
      design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords, X));
    }
    else design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords));
    Rcpp::Rcout << "QUI5 " << "\n";

    unsigned int n_covariates(design_matrix_ptr->cols());

    // Model
    model_fit::Model model(big_matrix_data_tspace_ptr, design_matrix_ptr, p);
    model.update_model(gamma_matrix);

    Eigen::MatrixXd resMatrix(N, ((p+1)*p)/2);
    Eigen::MatrixXd beta(n_covariates, ((p+1)*p)/2);
    Eigen::MatrixXd beta_old(n_covariates, ((p+1)*p)/2);
    Rcpp::Rcout << "QUI16 " << "\n";

    beta = model.get_beta();
    std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
    beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, p);
    std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

    unsigned int num_iter(0);
    unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
    double tolerance(Rcpp::as<double> (s_tolerance));
    std::vector<MatrixXd> resVec(N);

    double tol = tolerance+1;
    std::vector<double> emp_vario_values;
    Rcpp::Rcout << "QUI7 " << "\n";

    max_iter = 10;
    while (num_iter < max_iter && tol > tolerance) {
      resMatrix = model.get_residuals();
      resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, p);

      emp_vario.update_emp_vario(resVec, *(theTplaneDist));

      emp_vario_values = emp_vario.get_emp_vario_values();
      if(!max_sill_n.isNotNull()) max_sill = 1.15 * (*std::max_element(emp_vario_values.begin(), emp_vario_values.end()));
      if(!max_a_n.isNotNull()) max_a = 1.15 * emp_vario.get_hmax();
      the_variogram -> evaluate_par_fitted(emp_vario, max_sill, max_a);

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
    Rcpp::Rcout << "QUI8 " << "\n";

    Vec fit_parameters (the_variogram->get_parameters());

    bool suppressMes(Rcpp::as<bool> (s_suppressMes));
    if (!suppressMes) {
      if(fit_parameters(1)==max_sill-fit_parameters(0)) Rcpp::Rcout << "Parameter sill bounded from above" << "\n";
      if(fit_parameters(2)==max_a) Rcpp::Rcout << "Parameter a bounded from above" << "\n";
    }

    unsigned int n_hh(1000);
    Vec hh(n_hh);
    std::vector<double> h_vario_values(emp_vario.get_card_h());
    h_vario_values = emp_vario.get_hvec();

    hh.setLinSpaced(n_hh, 0, *std::max_element(h_vario_values.begin(), h_vario_values.end()));

    Vec fit_vario_values = the_variogram->get_vario_vec(hh, n_hh);

    // KRIGING
    Rcpp::Rcout << "Inizio kriging " << "\n";

    // New coordinates
    std::shared_ptr<const Eigen::MatrixXd> new_coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_new_coordinates));
    unsigned int M = new_coords_ptr->rows();
    Coordinates new_coords(new_coords_ptr);

    // New tangent points
    Rcpp::List Sigma_new(s_Sigma_new);
    std::vector<Eigen::MatrixXd> Sigma_new_vec(M);
    for(size_t i=0; i<M; i++) Sigma_new_vec[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(Sigma_new,i));

    // Map functions
    map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
    std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);

    // New Design matrix
    std::shared_ptr<Eigen::MatrixXd> new_design_matrix_ptr;
    if(X_new.isNotNull()) {
      Eigen::MatrixXd X_new(Rcpp::as<Eigen::MatrixXd> (s_X_new));
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
    Rcpp::Rcout << "Prima ciclo kriging " << "\n";

    for (size_t i=0; i<M; i++) {
      Rcpp::Rcout << "Ciclo kriging " << "\n";
      distanceVector = theDistance->create_distance_vector(coords, new_coords_ptr->row(i));
      ci = the_variogram->get_covario_vec(distanceVector, N);
      lambda_vec = solver.solve(ci);
      tplane_prediction = weighted_sum_beta(new_design_matrix_ptr->row(i)) + weighted_sum_residuals(lambda_vec);
      theExpMap->set_members(Sigma_new_vec[i]);
      manifold_prediction[i] = theExpMap->map2manifold(parallel_transport::transport_from_TI(Sigma_new_vec[i], tplane_prediction));
      if (i==13 || i==56 || i == 82 || i==235 || i==543) {
        Rcpp::Rcout << "tplane_prediction "<< "\n" << tplane_prediction << "\n";
        Rcpp::Rcout << "transported_from_TI "<< "\n" << parallel_transport::transport_from_TI(Sigma_new_vec[i], tplane_prediction) << "\n";
        Rcpp::Rcout << "manifold_prediction "<< "\n" << manifold_prediction[i] << "\n";
      }
    }

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("beta") = beta_vec_matrices,
                             Rcpp::Named("fit_vario_values") = fit_vario_values,
                             Rcpp::Named("hh") = hh,
                             Rcpp::Named("gamma_matrix") = gamma_matrix,
                             Rcpp::Named("residuals") = resVec,
                             Rcpp::Named("emp_vario_values") = emp_vario_values,
                             Rcpp::Named("h_vec") = h_vario_values,
                             Rcpp::Named("fitted_par_vario") = fit_parameters,
                             Rcpp::Named("iterations") = num_iter,
                             // Rcpp::Named("Sigma") = Sigma,
                             Rcpp::Named("prediction") = manifold_prediction);

    return Rcpp::wrap(result);
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
    Vec weight(Rcpp::as<Vec> (s_weight));
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

  // DISTANCE MANIFOLD
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
