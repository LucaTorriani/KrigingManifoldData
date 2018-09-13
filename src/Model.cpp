#include "Model.hpp"

using namespace model_fit;


void Model::update_model(const MatrixXd &gamma_matrix) {

  MatrixXd inv_gamma_matrix(_N, _N);
  MatrixXd Id(_N,_N);
  Id.setIdentity(_N,_N);
  inv_gamma_matrix= gamma_matrix.llt().solve(Id);

  MatrixXd A(_(_n*(_n+1))/2, (_n*(_n+1))/2);
  A = _design_matrix.transpose()*inv_gamma_matrix*_design_matrix;
  A.triangularView<StrictlyUpper>() = A.transpose();

  LeastSquaresConjugateGradient<MatrixXd> lscg;
  lscg.compute(A);

  Vec b((_n*(_n+1))/2);
  Vec y(_N);

  MatrixXd tmp((_n*(_n+1))/2), _N);
  tmp = _design_matrix.transpose()*inv_gamma_matrix;
  for(size_t l=0; l< (_n*(_n+1))/2); l++){
    y = _fitted_values.col(l);
    b = tmp*y;
    _beta_matrix.col(l) = lscg.solve(b);
  }
  _fitted_values = _design_matrix*_beta_matrix;
  _residuals = _data_tspace - _fitted_values;
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
