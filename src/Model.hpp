#ifndef _MODEL_HPP_
#define _MODEL_HPP_

#include <vector>
#include <utility>
#include "Helpers.hpp"

namespace model_fit {
class Model {
private:
  MatrixXd _beta_matrix;
  const MatrixXd& _data_tspace;
  const unsigned int _N; // Numero stazioni
  const unsigned int _n; // dimensione matrice manifold

  MatrixXd _fitted_values;
  MatrixXd _residuals;
  const MatrixXd& _design_matrix;
public:
  Model(const MatrixXd& data_tspace, const MatrixXd& design_matrix, unsigned int n): _data_tspace(data_tspace), _design_matrix(design_matrix),
      _N(_data_tspace.rows()), _n(n), _beta_matrix(_design_matrix.cols(), (_n*(_n+1))/2){};
  void update_model(MatrixXd);  // Updates Beta, fitted e residuals
  MatrixXd get_beta() const;
  MatrixXd get_residuals() const;
  MatrixXd get_fitted_values() const;
};

};

#endif
