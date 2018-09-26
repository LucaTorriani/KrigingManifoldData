#include "DistanceTplane.hpp"

#include <cmath>

using namespace distances_tplane;

double DistanceTplane::compute_distance(const MatrixXd& M1, const MatrixXd& M2) const{
  return (DistanceTplane::norm(M1-M2));
}

// FROBENIUS
double Frobenius::norm(const MatrixXd& M1) const{
  return (M1.norm());
}

void Frobenius::initialize_members(const std::shared_ptr<const MatrixXd>) {}

// FROBENIUS SCALED
double FrobeniusScaled::norm(const MatrixXd& M) const{
  MatrixXd tmp(_n, _n);
  tmp = _SigmaInv*M*_SigmaInv*M;
  return (sqrt(tmp.trace()));
}

void FrobeniusScaled::initialize_members(const std::shared_ptr<const MatrixXd> Sigma) {
  _n = Sigma->rows();
  Eigen::LDLT<MatrixXd> solver(_n);
  solver.compute(*(Sigma));
  MatrixXd Id(_n,_n);
  Id.setIdentity();
  _SigmaInv = solver.solve(Id);
}
