#ifndef _FITTED_VARIOGRAM_HPP
#define _FITTED_VARIOGRAM_HPP

#include "EmpiricalVariogram.hpp"
#include <iostream>
#include <Rcpp.h>

namespace variogram_evaluation {
class FittedVariogram{
protected:
  Vector3d _parameters; // tau2, sigma2, a
  double weighted_median (const std::vector<double> &, const std::vector<unsigned int> &);
  virtual void get_init_par(const EmpiricalVariogram &) = 0;
  void backtrack(const Vector3d &,Vector3d &,Vec &,const std::vector<double> &, unsigned int, double, double, const Vec&, double, double);
  virtual MatrixXd compute_jacobian(const std::vector<double> &, unsigned int) const = 0;

public:
  double get_tau2() const;
  double get_sigma2() const;
  double get_a() const;
  void evaluate_par_fitted_E(const EmpiricalVariogram &, double, double);
  void evaluate_par_fitted_W(const EmpiricalVariogram &, double, double);
  virtual double get_vario_univ(const double &) const = 0;
  double get_covario_univ(const double &) const;
  Vec get_vario_vec(const std::vector<double> &, unsigned int) const;
  Vec get_vario_vec(const Vec &, unsigned int) const;
  MatrixXd compute_gamma_matrix(const std::shared_ptr<const MatrixXd>, unsigned int) const;
  Vector3d get_parameters() const;
  void set_parameters(const Vector3d&);
  Vec get_covario_vec(const std::vector<double> &, unsigned int) const;

};


class GaussVariogram : public FittedVariogram {
  void get_init_par(const EmpiricalVariogram &) override;
  MatrixXd compute_jacobian(const std::vector<double> &, unsigned int) const override;

public:
  double get_vario_univ(const double &) const override;
};


class ExpVariogram : public FittedVariogram {
  void get_init_par(const EmpiricalVariogram &) override;
  MatrixXd compute_jacobian(const std::vector<double> &, unsigned int) const override;

public:
  double get_vario_univ(const double &) const override;
};

class SphVariogram : public FittedVariogram {
  void get_init_par(const EmpiricalVariogram &) override;
  MatrixXd compute_jacobian(const std::vector<double> &, unsigned int) const override;

public:
  double get_vario_univ(const double &) const override;
};


}

#endif
