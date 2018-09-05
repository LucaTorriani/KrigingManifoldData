#ifndef _VARIOGRAM_HPP
#define _VARIOGRAM_HPP

#include "Helpers.cpp"

namespace VariogramEvaluation{

class EmpiricalVariogram {
  Vec _emp_vario_values;
  Vec _hvec;
  double compute_hmax() const;

public:
  Vec get_emp_vario_values () const;
  Vec get_hvec() const;
  EmpiricalVariogram(/*residui, n_h*, il resto sono static members*/);
}

class FittedVariogram{
protected:
  Vec parameters; // tau2, sigma2, a
  void get_init_par(const EmpiricalVariogram &)) = 0;
  double compute_jacobian(const Vec &) const = 0;
  void backtrack(const Vec &,Vec &,Vec &,Matrixd &,const Vec &, double, double, std::function<Vec()>);

public:
  double get_tau2() const;
  double get_sigma2() const;
  double get_a() const;
  void evaluate_par_fitted(const EmpiricalVariogram &);
  double get_vario_univ(const double &) const = 0;
  double get_covario_univ(const double &) const = 0;

  Vec get_vario_vec(const Vec &) const = 0;
  Vec get_covario_vec(const Vec &) const = 0;


}


class GaussVariogram : public FittedVariogram {

  void get_init_par() override;
  double compute_jacobian(const double &) const;

public:

  void update_par_fitted() override;
  double get_vario_univ(const double &) const override;
  double get_covario_univ(const double &) const override;

}

}
