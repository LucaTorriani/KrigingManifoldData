#ifndef _MODEL_HPP_
#define _MODEL_HPP_

#include <vector>
#include <iostream>

#include <utility>
#include <Eigen/IterativeLinearSolvers>
#include <memory>
#include "Helpers.hpp"

namespace model_fit {
class Model {
private:
  const std::shared_ptr<const MatrixXd> _data_tspace;
  const std::shared_ptr<const MatrixXd> _design_matrix_model;
  const std::shared_ptr<const MatrixXd> _design_matrix_tot;

  const unsigned int _N; // Numero stazioni modello
  const unsigned int _p; // dimensione matrice manifold
  const unsigned int _num_cov; // numero covariate nel modello
  const unsigned int _num_coeff; // Number of coefficients in the upper triangular of a _nx_n
  MatrixXd _beta_matrix;

  MatrixXd _fitted_values;
  MatrixXd _residuals;
public:
  Model(const std::shared_ptr<const MatrixXd> data_tspace, const std::shared_ptr<const MatrixXd> design_matrix_model, unsigned int p):  // EQUAL WEIGHTS
        _data_tspace(data_tspace), _design_matrix_model(design_matrix_model), _design_matrix_tot(design_matrix_model),
       _N(design_matrix_model->rows()), _p(p), _num_cov(_design_matrix_model->cols()), _num_coeff((_p*(_p+1))/2), _beta_matrix(_num_cov, _num_coeff){};
  Model(const std::shared_ptr<const MatrixXd> data_tspace, const std::shared_ptr<const MatrixXd> design_matrix_model,  const std::shared_ptr<const MatrixXd> design_matrix_tot, unsigned int p): // KERNEL
        _data_tspace(data_tspace), _design_matrix_model(design_matrix_model), _design_matrix_tot(design_matrix_tot),
        _N(design_matrix_model->rows()), _p(p), _num_cov(_design_matrix_model->cols()), _num_coeff((_p*(_p+1))/2), _beta_matrix(_num_cov, _num_coeff){};

  void update_model(const MatrixXd&);  // Updates Beta, fitted e residuals
  MatrixXd get_beta() const;
  MatrixXd get_residuals() const;
  MatrixXd get_fitted_values() const;
};

}

#endif
