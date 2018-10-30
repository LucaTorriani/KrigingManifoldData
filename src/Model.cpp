#include "Model.hpp"
#include <Eigen/Eigenvalues>

using namespace model_fit;


void Model::update_model(const MatrixXd &gamma_matrix) {
  MatrixXd inv_gamma_matrix(_N, _N);
  MatrixXd Id(_N,_N);
  Id.setIdentity(_N,_N);
  LDLT<MatrixXd> ldlt(_N);
  ldlt.compute(gamma_matrix);
  inv_gamma_matrix= ldlt.solve(Id);

  MatrixXd A(_num_cov, _num_cov);
  A.triangularView<Lower>() = (_design_matrix_model->transpose()) *inv_gamma_matrix* (*(_design_matrix_model));
  A.triangularView<StrictlyUpper>() = A.transpose();

  LDLT<MatrixXd> solver(_num_cov);
  solver.compute(A);

  Vec b(_num_cov);
  Vec y(_N);

  MatrixXd tmp(_num_cov, _N);
  tmp = (_design_matrix_model->transpose())*inv_gamma_matrix;

  if (_distance_Manifold_name== "Correlation") {
    _beta_matrix.col(0).setZero(_num_cov);
    for(size_t l=1; l< _num_coeff; l++){
      y = _data_tspace->col(l);
      b = tmp*y;
      _beta_matrix.col(l) = solver.solve(b);
    }
  }
  else {
    for(size_t l=0; l< _num_coeff; l++){
      y = _data_tspace->col(l);
      b = tmp*y;
      _beta_matrix.col(l) = solver.solve(b);
    }
  }

  _fitted_values = (*(_design_matrix_tot))*_beta_matrix;
  _residuals = *(_data_tspace) - _fitted_values;
}

MatrixXd Model::get_residuals() const {
  return _residuals;
}

MatrixXd Model::get_fitted_values() const {
  return _fitted_values;
}

MatrixXd Model::get_beta() const{
  return _beta_matrix;
}
