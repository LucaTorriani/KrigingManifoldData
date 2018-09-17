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
  Vector3d _parameters; // tau2, sigma2, a
  // const SpMat& _distanceMatrix; // SERVE COSTRUTTORE CHE LA RICEVA COME PARAMETRO
  double weighted_median (const std::vector<double> &, const std::vector<unsigned int> &);
  virtual void get_init_par(const EmpiricalVariogram &) = 0;
  void backtrack(const Vector3d &,Vector3d &,Vec &,MatrixXd &,const std::vector<double> &, unsigned int, double, double, const Vec&);
  virtual MatrixXd compute_jacobian(const std::vector<double> &, unsigned int) const = 0;

public:
  // FittedVariogram(const SpMat& distanceMatrix, unsigned int N): _distanceMatrix(distanceMatrix), _N(N){};  // NON SI PUO' FARE CON I FACTORY TEMPLATE
  double get_tau2() const;
  double get_sigma2() const;
  double get_a() const;
  void evaluate_par_fitted(const EmpiricalVariogram &);
  virtual double get_vario_univ(const double &) const = 0;
  double get_covario_univ(const double &) const;
  Vec get_vario_vec(const std::vector<double> &, unsigned int) const;
  MatrixXd compute_gamma_matrix(const SpMat&, unsigned int) const;
  Vector3d get_parameters() const;
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
