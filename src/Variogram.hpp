#ifndef _VARIOGRAM_HPP
#define _VARIOGRAM_HPP

#include "Helpers.cpp"

namespace VariogramEvaluation{

class EmpiricalVariogram {
  //Variano
  std::vector<double> _emp_vario_values;
  std::vector<double> _hvec;
  std::vector<unsigned int> _N_hvec;
  // Costanti
  double _N;
  Vec _d; // Vettore h+-deltah
  double _hmax;
  const DistanceTplane _distanceTplane;
  const SpMat _distanceMatrix;
  Vec _weights;

  void compute_hmax(const Coordinates&, const Distance&);

public:
  Vec get_emp_vario_values () const;
  Vec get_N_hvec() const;
  Vec get_hvec() const;
  void update_emp_vario(/*residui*/);
}

class FittedVariogram{
protected:
  double weighted_median (const Vec &, const Vec &);
  Vector3d _parameters; // tau2, sigma2, a
  virtual void get_init_par(const EmpiricalVariogram &)) = 0;
  void backtrack(const Vec &,Vec &,Vec &,Matrixd &,const Vec &, double, double, std::function<Vec()>);
  virtual MatrixXd compute_jacobian(const Vec &) const = 0;

public:
  FittedVariogram();
  double get_tau2() const;
  double get_sigma2() const;
  double get_a() const;
  void evaluate_par_fitted(const EmpiricalVariogram &);
  virtual double get_vario_univ(const double &) const = 0;
  virtual double get_covario_univ(const double &) const = 0;
  virtual Vec get_vario_vec(const Vec &) const = 0;

}


class GaussVariogram : public FittedVariogram {
  void get_init_par(const EmpiricalVariogram &) override;
  MatrixXd compute_jacobian(const Vec &) const;

public:
  double get_vario_univ(const double &) const override;
  double get_covario_univ(const double &) const override;
  Vec get_vario_vec(const Vec &) const override;
}

}
