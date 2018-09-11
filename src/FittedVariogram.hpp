#ifndef _FITTED_VARIOGRAM_HPP
#define _FITTED_VARIOGRAM_HPP

// #include "Point.hpp"
// #include "Helpers.hpp"
// #include "Distance_Tplane.hpp"
// #include "Distance.hpp"
// #include "Coordinates.hpp"
#include "EmpiricalVariogram.hpp"

namespace variogram_evaluation {
class FittedVariogram{
protected:
  double weighted_median (const std::vector<double> &, const std::vector<unsigned int> &);
  Vector3d _parameters; // tau2, sigma2, a
  virtual void get_init_par(const EmpiricalVariogram &) = 0;
  void backtrack(const Vec &,Vec &,Vec &,MatrixXd &,const std::vector<double> &, unsigned int, double, double, const Vec&);
  virtual MatrixXd compute_jacobian(const std::vector<double> &, unsigned int) const = 0;

public:
  double get_tau2() const;
  double get_sigma2() const;
  double get_a() const;
  void evaluate_par_fitted(const EmpiricalVariogram &);
  virtual double get_vario_univ(const double &) const = 0;
  virtual double get_covario_univ(const double &) const = 0;
  virtual Vec get_vario_vec(const std::vector<double> &, unsigned int) const = 0;
  // virtual SpMat get_gamma_matrix() const = 0;
};


class GaussVariogram : public FittedVariogram {
  void get_init_par(const EmpiricalVariogram &) override;
  MatrixXd compute_jacobian(const std::vector<double> &, unsigned int) const override;

public:
  double get_vario_univ(const double &) const override;
  double get_covario_univ(const double &) const override;
  Vec get_vario_vec(const std::vector<double> &, unsigned int) const override;
};
}

#endif
