#include "DistanceTplane.hpp"

#include <cmath>

using namespace distances_tplane;

// FROBENIUS
double Frobenius::norm(const MatrixXd& M1) const{
  return (M1.norm());
}

double Frobenius::operator()(const MatrixXd& M1, const MatrixXd& M2) const{
  return ((M1-M2).norm());
}

// FROBENIUS SCALED
double FrobeniusScaled::norm(const MatrixXd& M) const{

  unsigned int n(M.rows());
  Eigen::LDLT<MatrixXd> solver(n);
  solver.compute(*(_Sigma));
  MatrixXd Id(n,n);
  Id.setIdentity();
  MatrixXd SigmaInv(n,n);
  SigmaInv = solver.solve(Id);

  MatrixXd tmp(n, n);
  tmp = SigmaInv*M*SigmaInv*M;
  return (sqrt(tmp.trace()));
}

double FrobeniusScaled::operator()(const MatrixXd& M1, const MatrixXd& M2) const{
  return (FrobeniusScaled::norm((M1-M2)));
}


DistanceTplane::DistanceTplane(const std::string & distanceTplane, const std::shared_ptr<const MatrixXd> Sigma):_distanceTplane(distanceTplane){
  _distances.insert(std::pair<std::string, std::function<double(const MatrixXd&, const MatrixXd&)> >("Frobenius", Frobenius()));
  _distances.insert(std::pair<std::string, std::function<double(const MatrixXd&, const MatrixXd&)> >("FrobeniusScaled", FrobeniusScaled(Sigma)));

}

double DistanceTplane::compute_distance(const MatrixXd& M1, const MatrixXd& M2) const {
  double result = _distances.at(_distanceTplane)(M1, M2);
  return result;
}
