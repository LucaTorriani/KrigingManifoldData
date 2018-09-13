#include<iostream>
#include<vector>
#include<utility>
#include<memory>
#include "Helpers.hpp"
#include "HelpersFactory.hpp"
#include "FittedVariogram.hpp"
#include "EmpiricalVariogram.hpp"

using namespace Eigen;

int main(){
  unsigned int N; // Number of stations
  unsigned int n; // Dimensione delle matrici sulla manifold

  MatrixXd tangent_point(n,n); // Matrice simmmetrica e definita positiva
  // tangent_point(0, 0) = 1;
  // tangent_point(0, 1) = 2;
  // tangent_point(1, 0) = 3;
  // tangent_point(1, 1) = 3;

  std::string distance_name("Euclidean"); //(Geodist, Euclidean)
  Distance distance(distance_name);

  std::string distance_Tplane_name("Frobenius"); //(Frobenius, FrobeniusScaled)
  DistanceTplane distanceTplane (distance_Tplane_name, tangent_point);

  std::string distance_Manifold_name("Frobenius"); //(Frobenius, SquareRoot, LogEuclidean)
  DistanceManifold distanceManifold (distance_Manifold_name, tangent_point);

  std::string log_map_name("Frobenius");  //(Frobenius, SquareRoot, LogEuclidean) //STESSO DI SOPRA?
  logarithmicMap logMap(distanceManifold);


  // FIXED DATA (Opzione1) passare e salvare come const reference nelle classi a cui servono
  //            (Opzione2) mettere tutte insieme in una classe e passare questa allo stesso modo
  std::vector<Point> point_coords(N);
  Coordinates coords(point_coords, distance);
  std::vector<MatrixXd> data_manifold(N);
  std::vector<MatrixXd> data_tspace(N);
  for (size_t i=0; i<N; i++) {
    data_tspace[i]=logMap.map2tplane(data_manifold[i]);
  }
  MatrixXd big_matrix_data_tspace(N, (n+1)*n/2);
  big_matrix_data_tspace = VecMatrices2bigMatrix(data_tspace, n);

  //

  // Emp vario
  unsigned int n_h;
  EmpiricalVariogram emp_vario(coords, distance, n_h, distanceTplane);

  // Fitted vario
  VariogramFactory & vf(VariogramFactory::Instance());
  std::string variogram_type;
  unique_ptr<FittedVariogram> the_variogram = vf.create(variogram_type);

  // Gamma matrix
  MatrixXd gamma_matrix(MatrixXd::Identity(N,N));

  // Model
  // Create design_matrix
  unsigned int n_covariates(design_matrix.cols());
  Model model(data_tspace, design_matrix);
  model.update_model(gamma_matrix);
  MatrixXd residuals(N, ((n+1)*n)/2);
  MatrixXd beta(N, ((n+1)*n)/2);
  MatrixXd beta_old(N, (n+1)*n)/2);
  beta = model.getBeta();
  beta_vec_matrices = bigMatrix2VecMatrices(beta);

  unsigned int num_iter(0);
  unsigned int max_iter;
  double tolerance;
  double tol = tolerance+1;
  while (num_iter < max_iter && tol > tolerance) {
    rediduals = model.get_residuals();
    emp_vario.update_emp_vario(rediduals);
    the_variogram->evaluate_par_fitted(emp_vario);
    gamma_matrix = the_variogram->compute_gamma_matrix();
    beta_old_vec_matrices = beta_vec_matrices;
    model.update_model(gamma_matrix);
    beta = model.getBeta();
    beta_vec_matrices = bigMatrix2VecMatrices(beta);
    tol=0;
    for (size_t i=0; i<N; i++) {
      tol += distanceTplane.compute_distance(beta_old_vec_matrices[i], beta_vec_matrices[i]);
    }
    num_iter++;
  }


 //  max_dist = max(vec_station_distance)
 // W =diag(dim(data_manifold)[1])
 //
 // data_ts = t(apply(data_manifold, 1, logarithmic_map_vec, Sigma=Sigma, metric_manifold=metric_manifold)) # Dati sul piano tangente
 //
 // model = compute_beta(data_ts, coords, X,metric_manifold, model_ts, W)
 // beta_new = model$coeff
 // fit_values = model$fit_values
 // iter = 0
 // tol = tolerance + 1
 // while (iter < max_it && tol >tolerance){
 //   Res = data_ts - fit_values
 //   empirical_variogram = emp_vario(Res, coords, vec_station_distance, n_h, metric_manifold, metric_ts, distance, Sigma, weight) # *NEW* argument weight
 //   fitted_par_vario = get_par_fitted_vario(vario_model, empirical_variogram, max_dist = max_dist)
 //   W = get_gamma_matrix(vario_model,fitted_par_vario,vec_station_distance, dim(data_ts)[1])
 //   beta_old = beta_new
 //   model = compute_beta(data_ts, coords, X,metric_manifold, model_ts, W)
 //   beta_new = model$coeff
 //   fit_values = model$fit_values
 //   iter = iter+1
 //   tol = 0
 //   for (i in dim(beta_new)[1]){
 //     if (metric_ts == "Frobenius") {
 //       tol = tol + tplane_distance (vec_to_matrix(beta_new[i,]), vec_to_matrix(beta_old[i,]), metric_ts)
 //     }
 //     else if (metric_ts =="Scaled_Frobenius") {
 //       tol = tol + tplane_distance (vec_to_matrix(beta_new[i,]), vec_to_matrix(beta_old[i,]), metric_ts, Sigma)
 //     }
 //     else {
 //       stop ("Tangent space metric not available")
 //     }
 //   }
 //
 // }






}
