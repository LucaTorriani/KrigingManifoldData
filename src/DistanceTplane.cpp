#include "DistanceTplane.hpp"

#include <cmath>

using namespace distances_tplane;

double DistanceTplane::compute_distance(const MatrixXd& M1, const MatrixXd& M2) const{
  return (norm(M1-M2));
}

// FROBENIUS
double Frobenius::norm(const MatrixXd& M1) const{
  return (M1.norm());
}

void Frobenius::set_members(const MatrixXd& Sigma) {}

// FROBENIUS SCALED
double FrobeniusScaled::norm(const MatrixXd& M) const{
  MatrixXd tmp(_p, _p);
  tmp = _SigmaInv*M*_SigmaInv*M;
  return (sqrt(tmp.trace()));
}

void FrobeniusScaled::set_members(const MatrixXd& Sigma) {
  _p = Sigma.rows();
  Eigen::LDLT<MatrixXd> solver(_p);
  solver.compute(Sigma);
  MatrixXd Id(_p,_p);
  Id.setIdentity();
  _SigmaInv = solver.solve(Id);
}
