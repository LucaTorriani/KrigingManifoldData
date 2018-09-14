#ifndef _MODEL_HPP_
#define _MODEL_HPP_

#include <vector>
#include <utility>
#include <Eigen/IterativeLinearSolvers>
#include "Helpers.hpp"

namespace model_fit {
class Model {
private:
  const MatrixXd& _data_tspace;
  const MatrixXd& _design_matrix;
  const unsigned int _N; // Numero stazioni
  const unsigned int _n; // dimensione matrice manifold
  const unsigned int _num_cov; // numero covariate nel modello
  const unsigned int _num_coeff; // Number of coefficients in the upper triangular of a _nx_n
  MatrixXd _beta_matrix;

  MatrixXd _fitted_values;
  MatrixXd _residuals;
public:
  Model(const MatrixXd& data_tspace, const MatrixXd& design_matrix, unsigned int n): _data_tspace(data_tspace), _design_matrix(design_matrix),
      _N(_data_tspace.rows()), _n(n), _num_cov(_design_matrix.cols()), _num_coeff((_n*(_n+1))/2), _beta_matrix(_num_cov, _num_coeff){};
  void update_model(const MatrixXd&);  // Updates Beta, fitted e residuals
  MatrixXd get_beta() const;
  MatrixXd get_residuals() const;
  MatrixXd get_fitted_values() const;
};

};

#endif
