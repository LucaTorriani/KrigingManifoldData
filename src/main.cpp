#include<iostream>
#include<vector>
#include<utility>
#include<memory>
#include "Coordinates.hpp"
#include "DesignMatrix.hpp"
#include "Distance_Manifold.hpp"
#include "Distance_Tplane.hpp"
#include "Distance.hpp"
#include "EmpiricalVariogram.hpp"
#include "FittedVariogram.hpp"
#include "Helpers.hpp"
#include "HelpersFactory.hpp"
#include "MapFunctions.hpp"
#include "Model.hpp"


using namespace Eigen;
using namespace distances;
using namespace distances_tplane;
using namespace distances_manifold;
using namespace map_functions;
using namespace variogram_evaluation;
using namespace vario_factory;
using namespace design_matrix;
using namespace model_fit;
using namespace matrix_manipulation;

int main(){
  unsigned int N(7); // Number of stations
  unsigned int n(2); // Dimension of the matrices on the manifold

  // Tangent point
  MatrixXd tangent_point(n,n); // Matrice simmmetrica e definita positiva
  tangent_point = MatrixXd::Random(n,n);
  tangent_point = tangent_point.transpose()*tangent_point;

  std::string distance_name("Euclidean"); //(Geodist, Euclidean)
  Distance distance(distance_name);

  std::string distance_Tplane_name("Frobenius"); //(Frobenius, FrobeniusScaled)
  DistanceTplane distanceTplane (distance_Tplane_name, tangent_point);

  std::string distance_Manifold_name("Frobenius"); //(Frobenius, SquareRoot, LogEuclidean)
  DistanceManifold distanceManifold (distance_Manifold_name, tangent_point);

  logarithmicMap logMap(distanceManifold);


  // FIXED DATA (Opzione1) passare e salvare come const reference nelle classi a cui servono
  //            (Opzione2) mettere tutte insieme in una classe e passare questa allo stesso modo
  // coords
  MatrixXd coords_mat(N, 2);
  coords_mat = MatrixXd::Random(N,2);
  Coordinates coords(coords_mat);

  // Distance Matrix
  SpMat distanceMatrix(N,N);
  distanceMatrix = distance.create_distance_matrix(coords);

  // Data manifold
  std::vector<MatrixXd> data_manifold(N);
  MatrixXd manifold_mat(n,n);    // Matrice simmmetrica e definita positiva
  for (size_t i=0; i<N;i++) {
    manifold_mat = MatrixXd::Random(n,n);
    manifold_mat = manifold_mat.transpose()*manifold_mat;
    data_manifold.push_back(manifold_mat);
  }

  std::vector<MatrixXd> data_tspace(N);
  for (size_t i=0; i<N; i++) {
    data_tspace[i]=logMap.map2tplane(data_manifold[i]);
  }
  MatrixXd big_matrix_data_tspace(N, (n+1)*n/2);
  big_matrix_data_tspace = VecMatrices2bigMatrix(data_tspace);

  // Emp vario
  unsigned int n_h(15);
  EmpiricalVariogram emp_vario(coords, distance, n_h, distanceTplane, distanceMatrix);

  // Fitted vario
  VariogramFactory & vf(VariogramFactory::Instance());
  std::string variogram_type("Gaussian"); // (Gaussian, Exponential, Spherical) //IMPLEMENTAREEE
  std::unique_ptr<FittedVariogram> the_variogram = vf.create(variogram_type);

  // Gamma matrix
  MatrixXd gamma_matrix(MatrixXd::Identity(N,N));

  // Design matrix
  registerDesignMatrices();
  DesignMatrixFactory& design_matrix_fac = DesignMatrixFactory::Instance();
  std::string model_name("Coord1"); //(Intercept, Coord1, Coord2, Additive)
  std::unique_ptr<DesignMatrix> theDesign_matrix = design_matrix_fac.create(model_name);

  MatrixXd design_matrix = theDesign_matrix->compute_design_matrix(coords);


  // Model
  unsigned int n_covariates(design_matrix.cols());
  Model model(big_matrix_data_tspace, design_matrix, n);
  model.update_model(gamma_matrix);
  MatrixXd residuals(N, ((n+1)*n)/2);
  MatrixXd beta(N, ((n+1)*n)/2);
  MatrixXd beta_old(N, ((n+1)*n)/2);
  beta = model.get_beta();
  std::vector<MatrixXd> beta_vec_matrices(n_covariates);
  beta_vec_matrices= bigMatrix2VecMatrices(beta, n);
  std::vector<MatrixXd> beta_old_vec_matrices(n_covariates);

  unsigned int num_iter(0);
  unsigned int max_iter(100);
  double tolerance(1e-4);
  double tol = tolerance+1;
  while (num_iter < max_iter && tol > tolerance) {
    residuals = model.get_residuals();
    emp_vario.update_emp_vario(residuals);
    the_variogram->evaluate_par_fitted(emp_vario);
    gamma_matrix = the_variogram->compute_gamma_matrix(distanceMatrix, N);
    beta_old_vec_matrices = beta_vec_matrices;
    model.update_model(gamma_matrix);
    beta = model.get_beta();
    beta_vec_matrices = bigMatrix2VecMatrices(beta, n);
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
