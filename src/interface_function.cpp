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


/*! \file
 @brief Main functions to create the model and perform kriging, along with functions to compute the distance on the manifold and the intrinsic mean.
 */

extern "C"{
/*!
  @brief  Given the coordinates and corresponding manifold values, this function creates a GLS model on the tangent space.
  @details The manifold values are mapped on the tangent space and then a GLS model is fitted to them. A first estimate of the `beta` coefficients
  is obtained assuming spatially uncorrelated errors. Then, in the main the loop, new estimates of the `beta` are obtained as a result of a
  weighted least square problem where the weight matrix is the inverse of `gamma_matrix`. The residuals `(residuals = data_ts - fitted)`
  are updated accordingly. The parameters of the variogram fitted to the `residuals` (and used in the evaluation of the `gamma_matrix`) are
  computed using Gauss-Newton with backtrack method to solve the associated non-linear least square problem. The stopping criteria is based on the
  absolute value of the variogram residuals' norm if `ker.width.vario=0`, while it is based on its increment otherwise.
  @note
    Reference: "Kriging prediction for manifold-valued random fields." \n
    Authors: D. Pigoli, A. Menafoglio & P. Secchi (2016) \n
    Periodical: Journal of Multivariate Analysis, 145, 117-131.
  @param s_data_manifold list of \f$N\f$ symmetric positive definite matrices of dimension \f$\left(p*p\right)\f$
  @param s_coordinates \f$\left(N*2\right)\f$ or \f$\left(N*3\right)\f$ matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to be provided in signed decimal degrees
  @param s_X matrix Matrix with \f$N\f$ rows and unrestricted number of columns of additional covariates for the tangent space model, possibly `NULL`
  @param s_Sigma Matrix \f$\left(p*p\right)\f$ representing the tangent point. If `NULL` the tangent point is computed as the intrinsic mean of `s_data_manifold`
  @param s_distance Type of distance between coordinates. It must be either "Eucldist" or "Geodist"
  @param s_manifold_metric Metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot", "Correlation"
  @param s_ts_metric Metric used on the tangent space. It must be chosen among "Frobenius", "FrobeniusScaled", "Correlation"
  @param s_ts_model Type of model fitted on the tangent space. It must be chosen among "Intercept", "Coord1", "Coord2", "Additive"
  @param s_vario_model Type of variogram fitted. It must be chosen among "Gaussian", "Spherical", "Exponential"
  @param s_n_h Number of bins in the emprical variogram
  @param s_max_it Max number of iterations for the main loop
  @param s_tolerance Tolerance for the main loop
  @param s_max_sill Maximum value allowed for \e sill in the fitted variogram. If `NULL` it is defined as \f$1.15*\max(\mbox{emp\_vario\_values})\f$
  @param s_max_a Maximum value for \e a in the fitted variogram. If `NULL` it is defined as \f$1.15*\mbox{h\_max}\f$
  @param s_weight_vario Vector of length \f$N\_tot\f$ to weight the locations in the computation of the empirical variogram
  @param s_distance_matrix_tot Matrix \f$\left(N\_tot*N\_tot\right)\f$ of distances between the locations,
  @param s_data_manifold_tot List of \f$N\_tot\f$ symmetric positive definite matrices of dimension \f$\left(p*p\right)\f$
  @param s_coordinates_tot \f$\left(N\_tot*2\right)\f$ or \f$\left(N\_tot*3\right)\f$ matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to be provided in signed decimal degrees),
  @param s_X_tot Matrix with \f$N\_tot\f$ rows and unrestricted number of columns, of additional covariates for the tangent space model. Possibly `NULL`
  @param s_hmax Maximum value of distance for which the variogram is computed
  @param s_indexes_model Indexes corresponding to `coords` in `coords_tot`. Required only in the case `metric_manifold=="Correlation"`
  @param s_weight_intrinsic Vector of length \f$N\f$ to weight the locations in the computation of the intrinsic mean. If `NULL` a vector of ones is used. Not needed if `Sigma` is provided
  @param s_tolerance_intrinsic Tolerance for the computation of the intrinsic mean. Not needed if `Sigma` is provided
  @param s_weight_extrinsic Vector of length \f$N\f$ to weight the locations in the computation of the extrinsic mean. If `NULL` `weight_intrinsic` are used. Needed only if `Sigma` is not provided and `metric_manifold=="Correlation"`
  @param s_suppressMes Boolean. If `TRUE` warning messagges are not printed
  @param s_tolerance_map_cor Tolerance to use in the maps. Required only if `metric_manifold=="Correlation"`
  @return A list with the following fields:
   - `beta` Vector of the beta matrices of the fitted model
   - `fit_vario_values` Vector of fitted variogram values in correspondence of `hh`
   - `hh` Dense vector of positions at which `fit_vario_values` is computed
   - `gamma_matrix` Covariogram matrix \f$\left(N*N\right)\f$
   - `residuals` Vector of the \f$N\f$ residual matrices
   - `emp_vario_values` Vector of empircal variogram values in correspondence of  `h_vec`
   - `h_vec` Vector of positions at which the empirical variogram is computed
   - `fitted_par_vario` Estimates of \e nugget, \e sill-nugget and \e practical \e range
   - `iterations` Number of iterations of the main loop
   - `Sigma` Tangent point
*/
  RcppExport SEXP get_model (SEXP s_data_manifold, SEXP s_coordinates, SEXP s_X, SEXP s_Sigma,
    SEXP s_distance, SEXP s_data_dist_mat, SEXP s_manifold_metric, SEXP s_ts_metric, SEXP s_ts_model, SEXP s_vario_model, SEXP s_n_h,
    SEXP s_max_it, SEXP s_tolerance, SEXP s_max_sill, SEXP s_max_a,
    SEXP s_weight_vario, SEXP s_distance_matrix_tot, SEXP s_data_manifold_tot, SEXP s_coordinates_tot, SEXP s_X_tot,  SEXP s_indexes_model, // RDD
    SEXP s_weight_intrinsic, SEXP s_tolerance_intrinsic, SEXP s_weight_extrinsic,  // Tangent point
    SEXP s_suppressMes, SEXP s_tolerance_map_cor) {

      BEGIN_RCPP
      Rcpp::Nullable<Eigen::MatrixXd> X(s_X);
      Rcpp::Nullable<Eigen::MatrixXd> X_tot(s_X_tot);
      Rcpp::Nullable<Eigen::MatrixXd> Sigma_n(s_Sigma);
      Rcpp::Nullable<Vec> weight_vario(s_weight_vario);
      Rcpp::Nullable<double> max_sill_n (s_max_sill);
      Rcpp::Nullable<double> max_a_n (s_max_a);
      Rcpp::Nullable<std::string> distance_n (s_distance);

      // Coordinates model
      std::shared_ptr<const Eigen::MatrixXd> coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_coordinates));
      Coordinates coords(coords_ptr);
      unsigned int N = coords.get_N_station();

      // Map functions
      std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
      map_factory::LogMapFactory& logmap_fac (map_factory::LogMapFactory::Instance());
      std::unique_ptr<map_functions::logarithmicMap> theLogMap = logmap_fac.create(distance_Manifold_name);
      if (distance_Manifold_name == "Correlation") {
        double tolerance_map_cor(Rcpp::as<double> (s_tolerance_map_cor));
        theLogMap->set_tolerance(tolerance_map_cor);
      }

      // Data manifold model
      std::vector<Eigen::MatrixXd> data_manifold(N);
      if (distance_Manifold_name == "Correlation") {
        MatrixXd tmp;
        for(size_t i=0; i<N; i++){
          tmp = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold,i));
          data_manifold[i] = matrix_manipulation::Chol_decomposition(tmp);
        }
      }
      else {
        for(size_t i=0; i<N; i++){
          data_manifold[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold,i));
        }
      }
      unsigned int p = data_manifold[0].rows();

      // Distance tplane
      std::string distance_Tplane_name = Rcpp::as<std::string> (s_ts_metric) ; //(Frobenius, FrobeniusScaled)
      tplane_factory::TplaneFactory& tplane_fac (tplane_factory::TplaneFactory::Instance());
      std::unique_ptr<distances_tplane::DistanceTplane> theTplaneDist = tplane_fac.create(distance_Tplane_name);

      // Punto tangente
      Eigen::MatrixXd Sigma(p,p);
      if(Sigma_n.isNotNull()) {
        Sigma = Rcpp::as<Eigen::MatrixXd> (s_Sigma);
        if (distance_Manifold_name == "Correlation") { Sigma = matrix_manipulation::Chol_decomposition(Sigma); }
        theTplaneDist->set_members(Sigma);
        theLogMap->set_members(Sigma);
      }
      else {
        double tolerance_intrinsic(Rcpp::as<double> (s_tolerance_intrinsic));
        Vec weights_intrinsic(Rcpp::as<Vec> (s_weight_intrinsic));
        Vec weights_extrinsic(Rcpp::as<Vec> (s_weight_extrinsic));
        map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
        std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);
        Sigma = intrinsic_mean_C(data_manifold, distance_Manifold_name, *theLogMap, *theExpMap, *theTplaneDist, tolerance_intrinsic, weights_intrinsic, weights_extrinsic);
      }

      // Emp vario
      unsigned int n_h (Rcpp::as<unsigned int>(s_n_h));
      variogram_evaluation::EmpiricalVariogram emp_vario(n_h);

      // Distance
      std::unique_ptr<distances::Distance> theDistance; // Used only if distance is not NULL
      std::shared_ptr<const Eigen::MatrixXd> distanceMatrix_ptr;
      std::shared_ptr<const Eigen::MatrixXd> distanceDataGridMatrix_ptr; // Used only if distance is NULL

      if(distance_n.isNotNull()) {
        distance_factory::DistanceFactory& distance_fac (distance_factory::DistanceFactory::Instance());
        std::string distance_name( Rcpp::as<std::string> (s_distance)) ; //(Geodist, Eucldist)
        // std::unique_ptr<distances::Distance> theDistance = distance_fac.create(distance_name);
        theDistance = distance_fac.create(distance_name);

        // Distance Matrix
        // std::shared_ptr<const SpMat> distanceMatrix_ptr = theDistance->create_distance_matrix(coords, N);
        distanceMatrix_ptr = theDistance->create_distance_matrix(coords, N);

        // Emp vario
        // emp_vario.set_distance_and_h_max(distanceMatrix_ptr, coords, *(theDistance));
      }
      else {
        distanceMatrix_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_data_dist_mat));

        // Emp vario
        // emp_vario.set_distance_and_h_max(distanceMatrix_ptr, max_dist);
      }

      // Fitted vario parameters
      double max_a;
      double max_sill;
      if(max_a_n.isNotNull()) max_a = Rcpp::as<double> (s_max_a);
      if(max_sill_n.isNotNull()) max_sill= Rcpp::as<double> (s_max_sill);

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

      unsigned int n_covariates(design_matrix_ptr->cols());

      // KERNEL
      if(weight_vario.isNotNull()) {
        // Weight vario
        Vec weight_vario(Rcpp::as<Vec> (s_weight_vario));
        // Distance Matrix tot
        std::shared_ptr<const Eigen::MatrixXd> distanceMatrix_tot_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_distance_matrix_tot));

        // Coordinates tot
        std::shared_ptr<const Eigen::MatrixXd> coords_tot_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_coordinates_tot));
        Coordinates coords_tot(coords_tot_ptr);
        unsigned int N_tot = coords_tot.get_N_station();

        // Data manifold tot
        std::vector<unsigned int> indexes_model(Rcpp::as<std::vector<unsigned int>> (s_indexes_model));

        std::vector<Eigen::MatrixXd> data_manifold_tot(N_tot);
        size_t ii(0);
        if (distance_Manifold_name == "Correlation") {
          for(size_t i=0; i<N_tot; i++){
            if (ii < indexes_model.size() && i == (indexes_model[ii]-1)) {
              data_manifold_tot[i] = data_manifold[ii];
              ii++;
            }
            else {
              data_manifold_tot[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold_tot,i));
              data_manifold_tot[i] = matrix_manipulation::Chol_decomposition(data_manifold_tot[i]);
            }
          }
        }
        else {
          for(size_t i=0; i<N_tot; i++) { // Penso sia più veloce rileggerli tutti piuttosto che fare i vari if
             data_manifold_tot[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold_tot,i));
          }
        }


        // Data tangent space tot
        std::vector<Eigen::MatrixXd> data_tspace_tot(N_tot);
        for (size_t i=0; i<N_tot; i++) {
          data_tspace_tot[i] = theLogMap->map2tplane(data_manifold_tot[i]);
        }

        // Data tspace tot
        std::shared_ptr<const Eigen::MatrixXd> big_matrix_data_tspace_tot_ptr = std::make_shared<const Eigen::MatrixXd>(matrix_manipulation::VecMatrices2bigMatrix(data_tspace_tot));
        data_manifold.clear();
        data_manifold_tot.clear();
        data_tspace_tot.clear();

        // Emp vario
        emp_vario.set_weights(N_tot,weight_vario);
        emp_vario.set_distance_and_h_max(distanceMatrix_tot_ptr, distanceMatrix_tot_ptr->maxCoeff()); // Sia che distance sia NULL che altrimenti

        // Design matrix tot
        std::shared_ptr<Eigen::MatrixXd> design_matrix_tot_ptr;
        if(X_tot.isNotNull()) {
          Eigen::MatrixXd X_tot(Rcpp::as<Eigen::MatrixXd> (s_X_tot));
          design_matrix_tot_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords_tot, X_tot));
        }
        else design_matrix_tot_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords_tot));

        // Model
        model_fit::Model model(big_matrix_data_tspace_tot_ptr, design_matrix_ptr, design_matrix_tot_ptr, p, distance_Manifold_name);

        Eigen::MatrixXd resMatrix(N_tot, ((p+1)*p)/2);
        std::vector<MatrixXd> resVec(N_tot);

        // DA QUI è TUTTO UGUALE (tranne la selezione dei residui)
        model.update_model(gamma_matrix);

        Eigen::MatrixXd beta(n_covariates, ((p+1)*p)/2);
        Eigen::MatrixXd beta_old(n_covariates, ((p+1)*p)/2);

        beta = model.get_beta();
        std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
        beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, p, distance_Manifold_name);
        std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

        unsigned int num_iter(0);
        unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
        double tolerance(Rcpp::as<double> (s_tolerance));

        double tol = tolerance+1;
        std::vector<double> emp_vario_values;

        while (num_iter < max_iter && tol > tolerance) {
          resMatrix = model.get_residuals();
          resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, p, distance_Manifold_name);

          emp_vario.update_emp_vario(resVec, *(theTplaneDist));

          emp_vario_values = emp_vario.get_emp_vario_values();
          if(!max_sill_n.isNotNull()) max_sill = 1.15 * (*std::max_element(emp_vario_values.begin(), emp_vario_values.end()));
          if(!max_a_n.isNotNull()) max_a = 1.15 * emp_vario.get_hmax();
          the_variogram -> evaluate_par_fitted_W(emp_vario, max_sill, max_a);

          gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix_ptr, N);
          beta_old_vec_matrices = beta_vec_matrices;

          model.update_model(gamma_matrix);
          beta = model.get_beta();
          beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, p, distance_Manifold_name);

          tol=0.0;
          for (size_t i=0; i<n_covariates; i++) {
            tol += theTplaneDist->compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
          }
          num_iter++;
        }
        if(num_iter == max_iter) Rcpp::warning("Reached max number of iterations");
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

        if (distance_Manifold_name == "Correlation") { Sigma = Sigma.transpose() * Sigma; };

        // Select residuals in the k-th cell
        std::vector<Eigen::MatrixXd> resVec_k(N);
        for (size_t ii=0; ii<N; ii++ ) {
            resVec_k[ii]=resVec[indexes_model[ii]-1];
        }
        //

        Rcpp::List result = Rcpp::List::create( Rcpp::Named("beta") = beta_vec_matrices,
                             Rcpp::Named("fit_vario_values") = fit_vario_values,
                             Rcpp::Named("hh") = hh,
                             Rcpp::Named("gamma_matrix") = gamma_matrix,
                             Rcpp::Named("residuals") = resVec,
                             Rcpp::Named("emp_vario_values") = emp_vario_values,
                             Rcpp::Named("h_vec") = h_vario_values,
                             Rcpp::Named("fitted_par_vario") = fit_parameters,
                             Rcpp::Named("iterations") = num_iter,
                             Rcpp::Named("Sigma")= Sigma );

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
        Vec weight_vario(N);
        weight_vario.setOnes(N);
        emp_vario.set_weights(N, weight_vario);
        if(distance_n.isNotNull()) {
          emp_vario.set_distance_and_h_max(distanceMatrix_ptr, coords, *(theDistance));
        }
        else {
          emp_vario.set_distance_and_h_max(distanceMatrix_ptr, distanceMatrix_ptr->maxCoeff());
        }

        // Model
        model_fit::Model model(big_matrix_data_tspace_ptr, design_matrix_ptr, p, distance_Manifold_name);

        Eigen::MatrixXd resMatrix(N, ((p+1)*p)/2);
        std::vector<MatrixXd> resVec(N);

        // DA QUI è TUTTO UGUALE
        model.update_model(gamma_matrix);

        Eigen::MatrixXd beta(n_covariates, ((p+1)*p)/2);
        Eigen::MatrixXd beta_old(n_covariates, ((p+1)*p)/2);

        beta = model.get_beta();
        std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
        beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, p, distance_Manifold_name);
        std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

        unsigned int num_iter(0);
        unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
        double tolerance(Rcpp::as<double> (s_tolerance));

        double tol = tolerance+1;
        std::vector<double> emp_vario_values;

        while (num_iter < max_iter && tol > tolerance) {
          resMatrix = model.get_residuals();
          resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, p, distance_Manifold_name);

          emp_vario.update_emp_vario(resVec, *(theTplaneDist));

          emp_vario_values = emp_vario.get_emp_vario_values();
          if(!max_sill_n.isNotNull()) max_sill = 1.15 * (*std::max_element(emp_vario_values.begin(), emp_vario_values.end()));
          if(!max_a_n.isNotNull()) max_a = 1.15 * emp_vario.get_hmax();
          the_variogram -> evaluate_par_fitted_E(emp_vario, max_sill, max_a);

          gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix_ptr, N);
          beta_old_vec_matrices = beta_vec_matrices;

          model.update_model(gamma_matrix);
          beta = model.get_beta();
          beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, p, distance_Manifold_name);

          tol=0.0;
          for (size_t i=0; i<n_covariates; i++) {
            tol += theTplaneDist->compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
          }
          num_iter++;
        }
        if(num_iter == max_iter) Rcpp::warning("Reached max number of iterations");

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

        if (distance_Manifold_name == "Correlation") { Sigma = Sigma.transpose() * Sigma; };

        Rcpp::List result = Rcpp::List::create(Rcpp::Named("beta") = beta_vec_matrices,
                             Rcpp::Named("fit_vario_values") = fit_vario_values,
                             Rcpp::Named("hh") = hh,
                             Rcpp::Named("gamma_matrix") = gamma_matrix,
                             Rcpp::Named("residuals") = resVec,
                             Rcpp::Named("emp_vario_values") = emp_vario_values,
                             Rcpp::Named("h_vec") = h_vario_values,
                             Rcpp::Named("fitted_par_vario") = fit_parameters,
                             Rcpp::Named("iterations") = num_iter,
                             Rcpp::Named("Sigma")= Sigma);

        return Rcpp::wrap(result);
      }
      END_RCPP
  }

/*!
    @brief  Given the GLS model kriging prediction on new location is performed.
    @details The model provided is used to perform simple kriging on the tangent space in correspondence of the new locations. The estimates are then mapped to the manifold to produce the actual prediction.
    @note
      Reference: "Kriging prediction for manifold-valued random fields." \n
      Authors: D. Pigoli, A. Menafoglio & P. Secchi (2016) \n
      Periodical: Journal of Multivariate Analysis, 145, 117-131.
    @param s_coordinates \f$\left(N*2\right)\f$ or \f$\left(N*3\right)\f$ matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to be provided in signed decimal degrees
    @param s_new_coordinates \f$\left(N*2\right)\f$ or \f$\left(N*3\right)\f$ matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to be provided in signed decimal degrees
    @param s_Sigma Matrix \f$\left(p*p\right)\f$ representing the tangent point. If `NULL` the tangent point is computed as the intrinsic mean of `s_data_manifold`
    @param s_distance Type of distance between coordinates. It must be either "Eucldist" or "Geodist"
    @param s_manifold_metric Metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot", "Correlation"
    @param s_ts_model Type of model fitted on the tangent space. It must be chosen among "Intercept", "Coord1", "Coord2", "Additive"
    @param s_vario_model Type of variogram fitted. It must be chosen among "Gaussian", "Spherical", "Exponential"
    @param s_beta Vector of the beta matrices of the fitted model
    @param s_gamma_matrix Covariogram matrix \f$\left(N*N\right)\f$
    @param s_vario_parameters Estimates of \e nugget, \e sill-nugget and \e practical \e range
    @param s_residuals Vector of the \f$N\f$ residual matrices
    @param s_X_new Matrix (with the same number of rows of `s_new_coordinates`) of additional covariates for the new locations, possibly `NULL`
    @param s_tolerance_map_cor Tolerance to use in the maps. Required only if `metric_manifold=="Correlation"`
    @return A list with the following field:
     - `prediction` Vector of matrices predicted at the new locations
  */
RcppExport SEXP get_kriging (SEXP s_coordinates, SEXP s_new_coordinates,  SEXP s_Sigma,
    SEXP s_distance, SEXP s_data_grid_dist_mat, SEXP s_manifold_metric, SEXP s_ts_model, SEXP s_vario_model,
    SEXP s_beta, SEXP s_gamma_matrix, SEXP s_vario_parameters, SEXP s_residuals, SEXP s_X_new, SEXP s_tolerance_map_cor) {

    BEGIN_RCPP

    Rcpp::Nullable<Eigen::MatrixXd> X_new(s_X_new);
    Rcpp::Nullable<std::string> distance_n(s_distance);
    std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)

    // Punto tangente
    Eigen::MatrixXd Sigma(Rcpp::as<Eigen::MatrixXd> (s_Sigma));
    if (distance_Manifold_name == "Correlation") { Sigma = matrix_manipulation::Chol_decomposition(Sigma); }
    unsigned int p = Sigma.rows();

    // Distance
    std::unique_ptr<distances::Distance> theDistance; // Used only if distance is not NULL
    std::shared_ptr<const Eigen::MatrixXd> distanceDataGridMatrix_ptr; // Used only if distance is NULL

    if(distance_n.isNotNull()) {
      distance_factory::DistanceFactory& distance_fac (distance_factory::DistanceFactory::Instance());
      std::string distance_name( Rcpp::as<std::string> (s_distance)) ; //(Geodist, Eucldist)
      theDistance = distance_fac.create(distance_name);
    }
    else {
      distanceDataGridMatrix_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_data_grid_dist_mat));
    }

    // Map functions
    map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
    std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);
    theExpMap->set_members(Sigma);
    if (distance_Manifold_name == "Correlation") {
      double tolerance_map_cor(Rcpp::as<double> (s_tolerance_map_cor));
      theExpMap->set_tolerance(tolerance_map_cor);
    }

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
      Eigen::MatrixXd X_new(Rcpp::as<Eigen::MatrixXd> (s_X_new));
      new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords, X_new));
    }
    else new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords));

    // Fitted vario
    Eigen::VectorXd parameters(Rcpp::as<Eigen::VectorXd> (s_vario_parameters));

    vario_factory::VariogramFactory & vf(vario_factory::VariogramFactory::Instance());
    std::string variogram_type (Rcpp::as<std::string> (s_vario_model)); // (Gaussian, Exponential, Spherical)
    std::unique_ptr<variogram_evaluation::FittedVariogram> the_variogram = vf.create(variogram_type);
    the_variogram->set_parameters(parameters);

    // Gamma matrix
    Eigen::MatrixXd gamma_matrix(Rcpp::as<Eigen::MatrixXd> (s_gamma_matrix));

    // Beta (lista di matrici)
    Rcpp::List list_beta(s_beta);
    size_t num_cov = list_beta.size();
    std::vector<Eigen::MatrixXd> beta_vec(N);
    for(size_t i=0; i<num_cov; i++) beta_vec[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(list_beta,i));

    // Residui (lista di matrici)
    Rcpp::List list_residuals(s_residuals);
    std::vector<Eigen::MatrixXd> residuals_vec(N);
    for(size_t i=0; i<N; i++) residuals_vec[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(list_residuals,i));

    Vec distanceVector(N);
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
      if (distance_n.isNotNull()) distanceVector = theDistance->create_distance_vector(coords, new_coords_ptr->row(i));
      else distanceVector = distanceDataGridMatrix_ptr->col(i);
      ci = the_variogram->get_covario_vec(distanceVector, N);
      lambda_vec = solver.solve(ci);
      tplane_prediction = weighted_sum_beta(new_design_matrix_ptr->row(i)) + weighted_sum_residuals(lambda_vec);
      manifold_prediction[i] = theExpMap->map2manifold(tplane_prediction);
      if (distance_Manifold_name == "Correlation") { manifold_prediction[i] = manifold_prediction[i].transpose() * manifold_prediction[i]; }
    }

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("prediction") = manifold_prediction);  // Solo questo?

    return Rcpp::wrap(result);
    END_RCPP
  }

/*!
    @brief  Given the coordinates and corresponding manifold values, this function firstly creates a GLS model on the tangent space, and then it performs kriging on the new locations.
    @details The manifold values are mapped on the tangent space and then a GLS model is fitted to them. A first estimate of the `beta` coefficients
    is obtained assuming spatially uncorrelated errors. Then, in the main the loop, new estimates of the `beta` are obtained as a result of a
    weighted least square problem where the weight matrix is the inverse of `gamma_matrix`. The residuals `(residuals = data_ts - fitted)`
    are updated accordingly. The parameters of the variogram fitted to the `residuals` (and used in the evaluation of the `gamma_matrix`) are
    computed using Gauss-Newton with backtrack method to solve the associated non-linear least square problem. The stopping criteria is based on the
    absolute value of the variogram residuals' norm if `ker.width.vario=0`, while it is based on its increment otherwise.
    Once the model is computed, simple kriging on the tangent space is performed in correspondence of the new locations and eventually the estimates are mapped to the manifold.
    @note
      Reference: "Kriging prediction for manifold-valued random fields." \n
      Authors: D. Pigoli, A. Menafoglio & P. Secchi (2016)   \n
      Periodical: Journal of Multivariate Analysis, 145, 117-131.
    @param s_data_manifold list of \f$N\f$ symmetric positive definite matrices of dimension \f$\left(p*p\right)\f$
    @param s_coordinates \f$\left(N*2\right)\f$ or \f$\left(N*3\right)\f$ matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to be provided in signed decimal degrees
    @param s_X matrix Matrix with \f$N\f$ rows and unrestricted number of columns of additional covariates for the tangent space model, possibly `NULL`
    @param s_Sigma Matrix \f$\left(p*p\right)\f$ representing the tangent point. If `NULL` the tangent point is computed as the intrinsic mean of `s_data_manifold`
    @param s_distance Type of distance between coordinates. It must be either "Eucldist" or "Geodist"
    @param s_manifold_metric Metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot", "Correlation"
    @param s_ts_metric Metric used on the tangent space. It must be chosen among "Frobenius", "FrobeniusScaled", "Correlation"
    @param s_ts_model Type of model fitted on the tangent space. It must be chosen among "Intercept", "Coord1", "Coord2", "Additive"
    @param s_vario_model Type of variogram fitted. It must be chosen among "Gaussian", "Spherical", "Exponential"
    @param s_n_h Number of bins in the emprical variogram
    @param s_max_it Max number of iterations for the main loop
    @param s_tolerance Tolerance for the main loop
    @param s_max_sill Maximum value allowed for \e sill in the fitted variogram. If `NULL` it is defined as \f$1.15*\max(\text{emp\_vario\_values})\f$
    @param s_max_a Maximum value for \e a in the fitted variogram. If `NULL` it is defined as \f$1.15*\text{h\_max}\f$
    @param s_weight_vario Vector of length \f$N\_tot\f$ to weight the locations in the computation of the empirical variogram
    @param s_distance_matrix_tot Matrix \f$\left(N\_tot*N\_tot\right)\f$ of distances between the locations,
    @param s_data_manifold_tot List of \f$N\_tot\f$ symmetric positive definite matrices of dimension \f$\left(p*p\right)\f$
    @param s_coordinates_tot \f$\left(N\_tot*2\right)\f$ or \f$\left(N\_tot*3\right)\f$ matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to be provided in signed decimal degrees),
    @param s_X_tot Matrix with \f$N\_tot\f$ rows and unrestricted number of columns, of additional covariates for the tangent space model. Possibly `NULL`
    @param s_hmax Maximum value of distance for which the variogram is computed
    @param s_indexes_model Indexes corresponding to `coords` in `coords_tot`. Required only in the case `metric_manifold=="Correlation"`
    @param s_weight_intrinsic Vector of length \f$N\f$ to weight the locations in the computation of the intrinsic mean. If `NULL` a vector of ones is used. Not needed if `Sigma` is provided
    @param s_tolerance_intrinsic Tolerance for the computation of the intrinsic mean. Not needed if `Sigma` is provided
    @param s_weight_extrinsic Vector of length \f$N\f$ to weight the locations in the computation of the extrinsic mean. If `NULL` `weight_intrinsic` are used. Needed only if `Sigma` is not provided and `metric_manifold=="Correlation"`
    @param s_new_coordinates \f$\left(N*2\right)\f$ or \f$\left(N*3\right)\f$ matrix of [lat,long], [x,y] or [x,y,z] coordinates. [lat,long] are supposed to be provided in signed decimal degrees
    @param s_X_new Matrix (with the same number of rows of `s_new_coordinates`) of additional covariates for the new locations, possibly `NULL`
    @param s_suppressMes Boolean. If `TRUE` warning messagges are not printed
    @param s_tolerance_map_cor Tolerance to use in the maps. Required only if `metric_manifold=="Correlation"`
    @return A list with the following fields:
     - `beta` Vector of the beta matrices of the fitted model
     - `fit_vario_values` Vector of fitted variogram values in correspondence of `hh`
     - `hh` Dense vector of positions at which `fit_vario_values` is computed
     - `gamma_matrix` Covariogram matrix \f$\left(N*N\right)\f$
     - `residuals` Vector of the \f$N\f$ residual matrices
     - `emp_vario_values` Vector of empircal variogram values in correspondence of  `h_vec`
     - `h_vec` Vector of positions at which the empirical variogram is computed
     - `fitted_par_vario` Estimates of \e nugget, \e sill-nugget and \e practical \e range
     - `iterations` Number of iterations of the main loop
     - `Sigma` Tangent point
     - `prediction` Vector of matrices predicted at the new locations
*/
RcppExport SEXP get_model_and_kriging (SEXP s_data_manifold, SEXP s_coordinates, SEXP s_X, SEXP s_Sigma,
    SEXP s_distance, SEXP s_data_dist_mat, SEXP s_data_grid_dist_mat, SEXP s_manifold_metric, SEXP s_ts_metric, SEXP s_ts_model, SEXP s_vario_model, SEXP s_n_h,
    SEXP s_max_it, SEXP s_tolerance, SEXP s_max_sill, SEXP s_max_a,
    SEXP s_weight_vario, SEXP s_distance_matrix_tot, SEXP s_data_manifold_tot, SEXP s_coordinates_tot, SEXP s_X_tot, SEXP s_indexes_model, // RDD
    SEXP s_weight_intrinsic, SEXP s_tolerance_intrinsic, SEXP s_weight_extrinsic,
    SEXP s_new_coordinates, SEXP s_X_new,  // KRIGING
    SEXP s_suppressMes, SEXP s_tolerance_map_cor) {

      BEGIN_RCPP
      Rcpp::Nullable<Eigen::MatrixXd> X(s_X);
      Rcpp::Nullable<Eigen::MatrixXd> X_tot(s_X_tot);
      Rcpp::Nullable<Eigen::MatrixXd> Sigma_n(s_Sigma);
      Rcpp::Nullable<Vec> weight_vario(s_weight_vario);
      Rcpp::Nullable<Eigen::MatrixXd> X_new(s_X_new);
      Rcpp::Nullable<double> max_sill_n (s_max_sill);
      Rcpp::Nullable<double> max_a_n (s_max_a);
      Rcpp::Nullable<std::string> distance_n (s_distance);


      // Coordinates model
      std::shared_ptr<const Eigen::MatrixXd> coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_coordinates));
      Coordinates coords(coords_ptr);
      unsigned int N = coords.get_N_station();

      std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
      // Data manifold model
      std::vector<Eigen::MatrixXd> data_manifold(N);
      if (distance_Manifold_name == "Correlation") {
        MatrixXd tmp;
        for(size_t i=0; i<N; i++){
          tmp = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold,i));
          data_manifold[i] = matrix_manipulation::Chol_decomposition(tmp);
        }
      }
      else {
        for(size_t i=0; i<N; i++){
          data_manifold[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold,i));
        }
      }
      unsigned int p = data_manifold[0].rows();

      // Distance tplane
      std::string distance_Tplane_name = Rcpp::as<std::string> (s_ts_metric) ; //(Frobenius, FrobeniusScaled)
      tplane_factory::TplaneFactory& tplane_fac (tplane_factory::TplaneFactory::Instance());
      std::unique_ptr<distances_tplane::DistanceTplane> theTplaneDist = tplane_fac.create(distance_Tplane_name);

      // Map functions
      map_factory::LogMapFactory& logmap_fac (map_factory::LogMapFactory::Instance());
      std::unique_ptr<map_functions::logarithmicMap> theLogMap = logmap_fac.create(distance_Manifold_name);
      map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
      std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);
      if (distance_Manifold_name == "Correlation") {
        double tolerance_map_cor(Rcpp::as<double> (s_tolerance_map_cor));
        theLogMap->set_tolerance(tolerance_map_cor);
        theExpMap->set_tolerance(tolerance_map_cor);
      }

      // Punto tangente
      Eigen::MatrixXd Sigma(p,p);
      if(Sigma_n.isNotNull()) {
        Sigma = Rcpp::as<Eigen::MatrixXd> (s_Sigma);
        if (distance_Manifold_name == "Correlation") { Sigma = matrix_manipulation::Chol_decomposition(Sigma); }
        theTplaneDist->set_members(Sigma);
        theLogMap->set_members(Sigma);
        theExpMap->set_members(Sigma);
      }
      else {
        double tolerance_intrinsic(Rcpp::as<double> (s_tolerance_intrinsic));
        Vec weights_intrinsic(Rcpp::as<Vec> (s_weight_intrinsic));
        Vec weights_extrinsic(Rcpp::as<Vec> (s_weight_extrinsic));
        // map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
        // std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);
        Sigma = intrinsic_mean_C(data_manifold, distance_Manifold_name, *theLogMap, *theExpMap, *theTplaneDist, tolerance_intrinsic, weights_intrinsic, weights_extrinsic);
      }

      // Emp vario
      unsigned int n_h (Rcpp::as<unsigned int>(s_n_h));
      variogram_evaluation::EmpiricalVariogram emp_vario(n_h);

      // Distance
      std::unique_ptr<distances::Distance> theDistance; // Used only if distance is not NULL
      std::shared_ptr<const Eigen::MatrixXd> distanceMatrix_ptr;
      std::shared_ptr<const Eigen::MatrixXd> distanceDataGridMatrix_ptr; // Used only if distance is NULL

      if(distance_n.isNotNull()) {
        distance_factory::DistanceFactory& distance_fac (distance_factory::DistanceFactory::Instance());
        std::string distance_name( Rcpp::as<std::string> (s_distance)) ; //(Geodist, Eucldist)
        // std::unique_ptr<distances::Distance> theDistance = distance_fac.create(distance_name);
        theDistance = distance_fac.create(distance_name);

        // Distance Matrix
        // std::shared_ptr<const SpMat> distanceMatrix_ptr = theDistance->create_distance_matrix(coords, N);
        distanceMatrix_ptr = theDistance->create_distance_matrix(coords, N);

        // Emp vario
        // emp_vario.set_distance_and_h_max(distanceMatrix_ptr, coords, *(theDistance));
      }
      else {
        distanceMatrix_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_data_dist_mat));

        // Data grid dist mat
        distanceDataGridMatrix_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_data_grid_dist_mat));

        // Emp vario
        // emp_vario.set_distance_and_h_max(distanceMatrix_ptr, max_dist);
      }


      // Fitted vario parameters
      double max_a;
      double max_sill;
      if(max_a_n.isNotNull()) max_a = Rcpp::as<double> (s_max_a);
      if(max_sill_n.isNotNull()) max_sill= Rcpp::as<double> (s_max_sill);

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

      unsigned int n_covariates(design_matrix_ptr->cols());

      // KERNEL
      if(weight_vario.isNotNull()) {
        // Weight vario
        Vec weight_vario(Rcpp::as<Vec> (s_weight_vario));

        // Distance Matrix tot
        std::shared_ptr<const Eigen::MatrixXd> distanceMatrix_tot_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_distance_matrix_tot));

        // Coordinates tot
        std::shared_ptr<const Eigen::MatrixXd> coords_tot_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_coordinates_tot));
        Coordinates coords_tot(coords_tot_ptr);
        unsigned int N_tot = coords_tot.get_N_station();

        // Data manifold tot
        std::vector<unsigned int> indexes_model(Rcpp::as<std::vector<unsigned int>> (s_indexes_model));

        std::vector<Eigen::MatrixXd> data_manifold_tot(N_tot);
        size_t ii(0);
        if (distance_Manifold_name == "Correlation") {
          for(size_t i=0; i<N_tot; i++){
            if (ii < indexes_model.size() && i == (indexes_model[ii]-1)) {
              data_manifold_tot[i] = data_manifold[ii];
              ii++;
            }
            else {
              data_manifold_tot[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold_tot,i));
              data_manifold_tot[i] = matrix_manipulation::Chol_decomposition(data_manifold_tot[i]);
            }
          }
        }
        else {
          for(size_t i=0; i<N_tot; i++) { // Penso sia più veloce rileggerli tutti piuttosto che fare i vari if
             data_manifold_tot[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data_manifold_tot,i));
          }
        }

        // Data tangent space tot
        std::vector<Eigen::MatrixXd> data_tspace_tot(N_tot);
        for (size_t i=0; i<N_tot; i++) {
          data_tspace_tot[i] = theLogMap->map2tplane(data_manifold_tot[i]);
        }

        std::shared_ptr<const Eigen::MatrixXd> big_matrix_data_tspace_tot_ptr = std::make_shared<const Eigen::MatrixXd>(matrix_manipulation::VecMatrices2bigMatrix(data_tspace_tot));
        data_manifold.clear();
        data_manifold_tot.clear();
        data_tspace_tot.clear();

        emp_vario.set_weights(N_tot,weight_vario);
        emp_vario.set_distance_and_h_max(distanceMatrix_tot_ptr, distanceMatrix_tot_ptr->maxCoeff()); // Sia che distance sia NULL che altrimenti

        // Design matrix tot
        std::shared_ptr<Eigen::MatrixXd> design_matrix_tot_ptr;
        if(X_tot.isNotNull()) {
          Eigen::MatrixXd X_tot(Rcpp::as<Eigen::MatrixXd> (s_X_tot));
          design_matrix_tot_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords_tot, X_tot));
        }
        else design_matrix_tot_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(coords_tot));

        // Model
        model_fit::Model model(big_matrix_data_tspace_tot_ptr, design_matrix_ptr, design_matrix_tot_ptr, p, distance_Manifold_name);

        Eigen::MatrixXd resMatrix(N_tot, ((p+1)*p)/2);
        std::vector<MatrixXd> resVec(N_tot);

        // DA QUI è TUTTO UGUALE tranne la parte ("Select residuals in the k-th cell")
        model.update_model(gamma_matrix);

        Eigen::MatrixXd beta(n_covariates, ((p+1)*p)/2);
        Eigen::MatrixXd beta_old(n_covariates, ((p+1)*p)/2);

        beta = model.get_beta();

        std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
        beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, p, distance_Manifold_name);
        std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

        unsigned int num_iter(0);
        unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
        double tolerance(Rcpp::as<double> (s_tolerance));

        double tol = tolerance+1;
        std::vector<double> emp_vario_values;

        while (num_iter < max_iter && tol > tolerance) {
          resMatrix = model.get_residuals();
          resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, p, distance_Manifold_name);

          emp_vario.update_emp_vario(resVec, *(theTplaneDist));

          emp_vario_values = emp_vario.get_emp_vario_values();
          if(!max_sill_n.isNotNull()) max_sill = 1.15 * (*std::max_element(emp_vario_values.begin(), emp_vario_values.end()));
          if(!max_a_n.isNotNull()) max_a = 1.15 * emp_vario.get_hmax();
          the_variogram -> evaluate_par_fitted_W(emp_vario, max_sill, max_a);

          gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix_ptr, N);
          beta_old_vec_matrices = beta_vec_matrices;

          model.update_model(gamma_matrix);
          beta = model.get_beta();

          beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, p, distance_Manifold_name);

          tol=0.0;
          for (size_t i=0; i<n_covariates; i++) {
            tol += theTplaneDist->compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
          }
          num_iter++;
        }
        if(num_iter == max_iter) Rcpp::warning("Reached max number of iterations");

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

        // // Map functions
        // map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
        // std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);
        // theExpMap->set_members(Sigma);

        // New coordinates
        std::shared_ptr<const Eigen::MatrixXd> new_coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_new_coordinates));
        unsigned int M = new_coords_ptr->rows();
        Coordinates new_coords(new_coords_ptr);

        // New Design matrix
        std::shared_ptr<Eigen::MatrixXd> new_design_matrix_ptr;
        if(X_new.isNotNull()) {
          Eigen::MatrixXd X_new(Rcpp::as<Eigen::MatrixXd> (s_X_new));
          new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords, X_new));
        }
        else new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords));

        Vec distanceVector(N);

        Vec ci(N);
        Vec lambda_vec(N);

        Eigen::LDLT<Eigen::MatrixXd> solver(N);
        solver.compute(gamma_matrix);

        unsigned int num_cov(beta_vec_matrices.size());
        Eigen::MatrixXd tmp(p,p);

        // Select residuals in the k-th cell
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
          if (distance_n.isNotNull()) distanceVector = theDistance->create_distance_vector(coords, new_coords_ptr->row(i));
          else distanceVector = distanceDataGridMatrix_ptr->col(i);
          ci = the_variogram->get_covario_vec(distanceVector, N);
          lambda_vec = solver.solve(ci);
          tplane_prediction = weighted_sum_beta(new_design_matrix_ptr->row(i)) + weighted_sum_residuals(lambda_vec);
          manifold_prediction[i] = theExpMap->map2manifold(tplane_prediction);
          if (distance_Manifold_name == "Correlation") { manifold_prediction[i] = manifold_prediction[i].transpose() * manifold_prediction[i]; }
        }

        if (distance_Manifold_name == "Correlation") { Sigma = Sigma.transpose() * Sigma; };

        Rcpp::List result = Rcpp::List::create(Rcpp::Named("beta") = beta_vec_matrices,
                                 Rcpp::Named("fit_vario_values") = fit_vario_values,
                                 Rcpp::Named("hh") = hh,
                                 Rcpp::Named("gamma_matrix") = gamma_matrix,
                                 Rcpp::Named("residuals") = resVec_k,
                                 Rcpp::Named("emp_vario_values") = emp_vario_values,
                                 Rcpp::Named("h_vec") = h_vario_values,
                                 Rcpp::Named("fitted_par_vario") = fit_parameters,
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
        Vec weight_vario(N);
        weight_vario.setOnes(N);
        emp_vario.set_weights(N, weight_vario);
        if(distance_n.isNotNull()) {
          emp_vario.set_distance_and_h_max(distanceMatrix_ptr, coords, *(theDistance));
        }
        else {
          emp_vario.set_distance_and_h_max(distanceMatrix_ptr, distanceMatrix_ptr->maxCoeff());
        }

        // Model
        model_fit::Model model(big_matrix_data_tspace_ptr, design_matrix_ptr, p, distance_Manifold_name);

        Eigen::MatrixXd resMatrix(N, ((p+1)*p)/2);
        std::vector<MatrixXd> resVec(N);

        // DA QUI è TUTTO UGUALE
        model.update_model(gamma_matrix);

        Eigen::MatrixXd beta(n_covariates, ((p+1)*p)/2);
        Eigen::MatrixXd beta_old(n_covariates, ((p+1)*p)/2);

        beta = model.get_beta();
        std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
        beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, p, distance_Manifold_name);
        std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

        unsigned int num_iter(0);
        unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
        double tolerance(Rcpp::as<double> (s_tolerance));

        double tol = tolerance+1;
        std::vector<double> emp_vario_values;

        while (num_iter < max_iter && tol > tolerance) {
          resMatrix = model.get_residuals();
          resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, p, distance_Manifold_name);

          emp_vario.update_emp_vario(resVec, *(theTplaneDist));

          emp_vario_values = emp_vario.get_emp_vario_values();
          if(!max_sill_n.isNotNull()) max_sill = 1.15 * (*std::max_element(emp_vario_values.begin(), emp_vario_values.end()));
          if(!max_a_n.isNotNull()) max_a = 1.15 * emp_vario.get_hmax();
          the_variogram -> evaluate_par_fitted_E(emp_vario, max_sill, max_a);

          gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix_ptr, N);
          beta_old_vec_matrices = beta_vec_matrices;

          model.update_model(gamma_matrix);
          beta = model.get_beta();
          beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, p, distance_Manifold_name);

          tol=0.0;
          for (size_t i=0; i<n_covariates; i++) {
            tol += theTplaneDist->compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
          }
          num_iter++;
        }
        if(num_iter == max_iter) Rcpp::warning("Reached max number of iterations");

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

        // // Map functions
        // map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
        // std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);
        // theExpMap->set_members(Sigma);

        // New coordinates
        std::shared_ptr<const Eigen::MatrixXd> new_coords_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_new_coordinates));
        unsigned int M = new_coords_ptr->rows();
        Coordinates new_coords(new_coords_ptr);

        // New Design matrix
        std::shared_ptr<Eigen::MatrixXd> new_design_matrix_ptr;
        if(X_new.isNotNull()) {
          Eigen::MatrixXd X_new(Rcpp::as<Eigen::MatrixXd> (s_X_new));
          new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords, X_new));
        }
        else new_design_matrix_ptr = std::make_shared<Eigen::MatrixXd> (theDesign_matrix->compute_design_matrix(new_coords));

        Vec distanceVector(N);

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
          if (distance_n.isNotNull()) distanceVector = theDistance->create_distance_vector(coords, new_coords_ptr->row(i));
          else distanceVector = distanceDataGridMatrix_ptr->col(i);
          ci = the_variogram->get_covario_vec(distanceVector, N);
          lambda_vec = solver.solve(ci);
          tplane_prediction = weighted_sum_beta(new_design_matrix_ptr->row(i)) + weighted_sum_residuals(lambda_vec);
          manifold_prediction[i] = theExpMap->map2manifold(tplane_prediction);
          if (distance_Manifold_name == "Correlation") { manifold_prediction[i] = manifold_prediction[i].transpose() * manifold_prediction[i]; }
        }

        if (distance_Manifold_name == "Correlation") { Sigma = Sigma.transpose() * Sigma; };

        Rcpp::List result = Rcpp::List::create(Rcpp::Named("beta") = beta_vec_matrices,
                                 Rcpp::Named("fit_vario_values") = fit_vario_values,
                                 Rcpp::Named("hh") = hh,
                                 Rcpp::Named("gamma_matrix") = gamma_matrix,
                                 Rcpp::Named("residuals") = resVec,
                                 Rcpp::Named("emp_vario_values") = emp_vario_values,
                                 Rcpp::Named("h_vec") = h_vario_values,
                                 Rcpp::Named("fitted_par_vario") = fit_parameters,
                                 Rcpp::Named("iterations") = num_iter,
                                 Rcpp::Named("Sigma") = Sigma,
                                 Rcpp::Named("prediction") = manifold_prediction);

        return Rcpp::wrap(result);
      }
    END_RCPP
}

/*!
  @brief  Evaluate the intrinsic mean of a given set of symmetric positive definite matrices.
  @param s_data list of \f$N\f$ symmetric positive definite matrices of dimension \f$\left(p*p\right)\f$
  @param s_N Number of data. `N = s_data.size()`
  @param s_manifold_metric Metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot", "Correlation"
  @param s_ts_metric Metric used on the tangent space. It must be chosen among "Frobenius", "FrobeniusScaled", "Correlation"
  @param s_tolerance Tolerance for the computation of the intrinsic_mean
  @param s_weight_intrinsic Vector of length \f$N\f$ to weight the locations in the computation of the intrinsic mean. If `NULL` a vector of ones is used
  @param s_weight_extrinsic Vector of length \f$N\f$ to weight the locations in the computation of the extrinsic mean. If `NULL` `weight_intrinsic` are used
  @param s_tolerance_map_cor Tolerance to use in the maps. Required only if `metric_manifold=="Correlation"`
  @return A matrix representing the intrinsic mean of the `s_data`
*/
 RcppExport SEXP intrinsic_mean (SEXP s_data, SEXP s_N, SEXP s_manifold_metric, SEXP s_ts_metric,
  SEXP s_tolerance, SEXP s_weight_intrinsic, SEXP s_weight_extrinsic, SEXP s_tolerance_map_cor) {
    BEGIN_RCPP
    std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
    // Data
    unsigned int N(Rcpp::as<unsigned int> (s_N));
    std::vector<Eigen::MatrixXd> data(N);

    if (distance_Manifold_name == "Correlation") {
      for(size_t i=0; i<N; i++){
        data[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data,i));
        data[i] = matrix_manipulation::Chol_decomposition(data[i]);
      }
    }
    else {
      for(size_t i=0; i<N; i++){
        data[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data,i));
      }
    }
    unsigned int p = data[0].rows();

    // Distance tplane
    std::string distance_Tplane_name = Rcpp::as<std::string> (s_ts_metric) ; //(Frobenius, FrobeniusScaled)
    tplane_factory::TplaneFactory& tplane_fac (tplane_factory::TplaneFactory::Instance());
    std::unique_ptr<distances_tplane::DistanceTplane> theTplaneDist = tplane_fac.create(distance_Tplane_name);

    // Map functions
    map_factory::LogMapFactory& logmap_fac (map_factory::LogMapFactory::Instance());
    std::unique_ptr<map_functions::logarithmicMap> theLogMap = logmap_fac.create(distance_Manifold_name);
    map_factory::ExpMapFactory& expmap_fac (map_factory::ExpMapFactory::Instance());
    std::unique_ptr<map_functions::exponentialMap> theExpMap = expmap_fac.create(distance_Manifold_name);

    if (distance_Manifold_name == "Correlation") {
      double tolerance_map_cor(Rcpp::as<double> (s_tolerance_map_cor));
      theLogMap->set_tolerance(tolerance_map_cor);
      theExpMap->set_tolerance(tolerance_map_cor);
    }
    // Tolerance
    double tolerance (Rcpp::as<double> (s_tolerance));

    // Weights
    Vec weight_intrinsic(Rcpp::as<Vec> (s_weight_intrinsic));
    Vec weight_extrinsic(Rcpp::as<Vec> (s_weight_extrinsic));
    double sum_weight_intrinsic(weight_intrinsic.sum());

    if (distance_Manifold_name == "Correlation") {
      Eigen::MatrixXd Result = matrix_manipulation::Chol_decomposition(extrinsic_mean(data, weight_extrinsic));

      theLogMap->set_members(Result);
      theExpMap->set_members(Result);
      theTplaneDist->set_members(Result);

      Eigen::MatrixXd Xk(p,p);
      Eigen::MatrixXd Xk_prec(p,p);

      Xk = theLogMap->map2tplane(data[0]);

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
          Xk = Xk + weight_intrinsic(i)* theLogMap->map2tplane(data[i]);
        }
        Xk = Xk/sum_weight_intrinsic;
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

      return Rcpp::wrap(Result.transpose() * Result);
    }
    else {
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
          Xk = Xk + weight_intrinsic(i)* theLogMap->map2tplane(data[i]);
        }
        Xk = Xk/sum_weight_intrinsic;
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
    }

    END_RCPP
}

/*!
  @brief Compute the manifold distance between symmetric positive definite matrices.
  @details If \f$ N1==N2 \f$ then the result is a vector of length \f$ N1=N2 \f$ containing in position `i` the manifold distance beetween `data1[[i]]` and `data2[[i]]`.
  Instead if \f$ N2==1 \f$ and \f$ N1!=1 \f$ the result is a vector of length \f$ N1 \f$ containing in position `i` the manifold distance between `data1[[i]]`and `data2[[1]]`
  @param s_data1 list of \f$N1\f$ symmetric positive definite matrices of dimension \f$\left(p*p\right)\f$
  @param s_data2 list of \f$N2\f$ symmetric positive definite matrices of dimension \f$\left(p*p\right)\f$
  @param s_N1 Number of data1. `N1 = s_data1.size()`
  @param s_N2 Number of data2. `N2 = s_data2.size()`
  @param s_manifold_metric Metric used on the manifold. It must be chosen among "Frobenius", "LogEuclidean", "SquareRoot", "Correlation"
  @return A double or a vector of distances
*/
 RcppExport SEXP distance_manifold (SEXP s_data1, SEXP s_data2, SEXP s_N1, SEXP s_N2, SEXP s_manifold_metric) {
      BEGIN_RCPP

      std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, FrobeniusScaled)

      // Data1
      unsigned int N1(Rcpp::as<unsigned int> (s_N1));
      std::vector<Eigen::MatrixXd> data1(N1);
      if (distance_Manifold_name == "Correlation") {
        MatrixXd tmp;
        for(size_t i=0; i<N1; i++){
          tmp = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data1,i));
          data1[i] = matrix_manipulation::Chol_decomposition(tmp);
        }
      }
      else {
        for(size_t i=0; i<N1; i++){
          data1[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data1,i));
        }
      }

      // Data2
      unsigned int N2(Rcpp::as<unsigned int> (s_N2));
      std::vector<Eigen::MatrixXd> data2(N1);
      if (distance_Manifold_name == "Correlation") {
        MatrixXd tmp;
        for(size_t i=0; i<N2; i++){
          tmp = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data2,i));
          data2[i] = matrix_manipulation::Chol_decomposition(tmp);
        }

      }
      else {
        for(size_t i=0; i<N2; i++){
          data2[i] = Rcpp::as<Eigen::MatrixXd>(VECTOR_ELT(s_data2,i));
        }
      }

      // Distance manifold
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

// CREATE MODEL AND KRIGNG
  RcppExport SEXP get_model_and_kriging_mixed (SEXP s_data_manifold, SEXP s_coordinates, SEXP s_X, SEXP s_Sigma_data,
     SEXP s_distance, SEXP s_data_dist_mat, SEXP s_data_grid_dist_mat, SEXP s_manifold_metric,  SEXP s_ts_model, SEXP s_vario_model, SEXP s_n_h, // SEXP s_ts_metric,
     SEXP s_max_it, SEXP s_tolerance, SEXP s_max_sill, SEXP s_max_a, // SEXP s_weight_vario, SEXP s_weight_intrinsic, SEXP s_tolerance_intrinsic,
     SEXP s_new_coordinates, SEXP s_Sigma_new, SEXP s_X_new, SEXP s_suppressMes) {

    BEGIN_RCPP

    // Rcpp::Nullable<Vec> weight_vario(s_weight_vario);
    Rcpp::Nullable<Eigen::MatrixXd> X(s_X);
    Rcpp::Nullable<Eigen::MatrixXd> X_new(s_X_new);
    Rcpp::Nullable<double> max_sill_n (s_max_sill);
    Rcpp::Nullable<double> max_a_n (s_max_a);
    Rcpp::Nullable<std::string> distance_n(s_distance);

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

    // Distance tplane
    // std::string distance_Tplane_name = Rcpp::as<std::string> (s_ts_metric) ; //(Frobenius, FrobeniusScaled)
    std::string distance_Tplane_name("Frobenius"); // Always Frobenius, since we consider the plane tangent in the identity
    tplane_factory::TplaneFactory& tplane_fac (tplane_factory::TplaneFactory::Instance());
    std::unique_ptr<distances_tplane::DistanceTplane> theTplaneDist = tplane_fac.create(distance_Tplane_name);

    // Map functions
    std::string distance_Manifold_name = Rcpp::as<std::string> (s_manifold_metric) ; //(Frobenius, SquareRoot, LogEuclidean)
    map_factory::LogMapFactory& logmap_fac (map_factory::LogMapFactory::Instance());
    std::unique_ptr<map_functions::logarithmicMap> theLogMap = logmap_fac.create(distance_Manifold_name);

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

    // Emp vario
    unsigned int n_h (Rcpp::as<unsigned int>(s_n_h));
    variogram_evaluation::EmpiricalVariogram emp_vario(n_h);
    Vec weight_vario(N);
    weight_vario.setOnes(N);
    emp_vario.set_weights(N, weight_vario);

    // Distance
    std::unique_ptr<distances::Distance> theDistance; // Used only if distance is not NULL
    std::shared_ptr<const Eigen::MatrixXd> distanceMatrix_ptr;
    std::shared_ptr<const Eigen::MatrixXd> distanceDataGridMatrix_ptr; // Used only if distance is NULL

    if(distance_n.isNotNull()) {
      distance_factory::DistanceFactory& distance_fac (distance_factory::DistanceFactory::Instance());
      std::string distance_name( Rcpp::as<std::string> (s_distance)) ; //(Geodist, Eucldist)
      // std::unique_ptr<distances::Distance> theDistance = distance_fac.create(distance_name);
      theDistance = distance_fac.create(distance_name);

      // Distance Matrix
      // std::shared_ptr<const SpMat> distanceMatrix_ptr = theDistance->create_distance_matrix(coords, N);
      distanceMatrix_ptr = theDistance->create_distance_matrix(coords, N);

      // Emp vario
      emp_vario.set_distance_and_h_max(distanceMatrix_ptr, coords, *(theDistance));
    }
    else {
      // Data dist mat
      // Vec disance_vec(Rcpp::as<Vec> (s_data_dist_vec));
      // std::vector<TripType> tripletList;
      // tripletList.reserve((N*(N-1))/2);
      // for (size_t i=0; i<(N-1); i++ ) {
      //   for (size_t j=(i+1); j<N; j++ ) {
      //     tripletList.push_back( TripType(i,j,disance_vec(N*i-(i+1)*i/2 + j-i -1)) );
      //   }
      // }
      // SpMat distance_matrix(N, N);
      // distance_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
      // distanceMatrix_ptr = std::make_shared<const SpMat> (distance_matrix);
      std::shared_ptr<const Eigen::MatrixXd> distanceMatrix_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_data_dist_mat));


      // Data grid dist mat
      distanceDataGridMatrix_ptr = std::make_shared<const Eigen::MatrixXd> (Rcpp::as<Eigen::MatrixXd> (s_data_grid_dist_mat));

      // Emp vario
      emp_vario.set_distance_and_h_max(distanceMatrix_ptr, distanceMatrix_ptr->maxCoeff());
    }




    // if(weight_vario.isNotNull()) {
    //   Vec weight_vario(Rcpp::as<Vec> (s_weight_vario));
    //   emp_vario.set_weight(weight_vario);
    // }

    // Fitted vario parameters
    double max_a;
    double max_sill;
    if(max_a_n.isNotNull()) max_a = Rcpp::as<double> (s_max_a);
    if(max_sill_n.isNotNull()) max_sill= Rcpp::as<double> (s_max_sill);

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

    unsigned int n_covariates(design_matrix_ptr->cols());

    // Model
    model_fit::Model model(big_matrix_data_tspace_ptr, design_matrix_ptr, p, distance_Manifold_name);
    model.update_model(gamma_matrix);

    Eigen::MatrixXd resMatrix(N, ((p+1)*p)/2);
    Eigen::MatrixXd beta(n_covariates, ((p+1)*p)/2);
    Eigen::MatrixXd beta_old(n_covariates, ((p+1)*p)/2);

    beta = model.get_beta();
    std::vector<Eigen::MatrixXd> beta_vec_matrices(n_covariates);
    beta_vec_matrices= matrix_manipulation::bigMatrix2VecMatrices(beta, p, distance_Manifold_name);
    std::vector<Eigen::MatrixXd> beta_old_vec_matrices(n_covariates);

    unsigned int num_iter(0);
    unsigned int max_iter(Rcpp::as<unsigned int> (s_max_it));
    double tolerance(Rcpp::as<double> (s_tolerance));
    std::vector<MatrixXd> resVec(N);

    double tol = tolerance+1;
    std::vector<double> emp_vario_values;

    while (num_iter < max_iter && tol > tolerance) {
      resMatrix = model.get_residuals();
      resVec = matrix_manipulation::bigMatrix2VecMatrices(resMatrix, p, distance_Manifold_name);

      emp_vario.update_emp_vario(resVec, *(theTplaneDist));

      emp_vario_values = emp_vario.get_emp_vario_values();
      if(!max_sill_n.isNotNull()) max_sill = 1.15 * (*std::max_element(emp_vario_values.begin(), emp_vario_values.end()));
      if(!max_a_n.isNotNull()) max_a = 1.15 * emp_vario.get_hmax();
      the_variogram -> evaluate_par_fitted_E(emp_vario, max_sill, max_a);

      gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix_ptr, N);
      beta_old_vec_matrices = beta_vec_matrices;

      model.update_model(gamma_matrix);
      beta = model.get_beta();
      beta_vec_matrices = matrix_manipulation::bigMatrix2VecMatrices(beta, p, distance_Manifold_name);

      tol=0.0;
      for (size_t i=0; i<n_covariates; i++) {
        tol += theTplaneDist->compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
      }

      num_iter++;
    }
    if(num_iter == max_iter) Rcpp::warning("Reached max number of iterations");

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

    Vec distanceVector(N);

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
      if (distance_n.isNotNull()) distanceVector = theDistance->create_distance_vector(coords, new_coords_ptr->row(i));
      else distanceVector = distanceDataGridMatrix_ptr->col(i);
      ci = the_variogram->get_covario_vec(distanceVector, N);
      lambda_vec = solver.solve(ci);
      tplane_prediction = weighted_sum_beta(new_design_matrix_ptr->row(i)) + weighted_sum_residuals(lambda_vec);
      theExpMap->set_members(Sigma_new_vec[i]);
      manifold_prediction[i] = theExpMap->map2manifold(parallel_transport::transport_from_TI(Sigma_new_vec[i], tplane_prediction));
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

}
