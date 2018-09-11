#include "Distance_Tplane.hpp"

#include <cmath>

using namespace distances_tplane;

// FROBENIUS

double Frobenius::norm(const SpMat& M1) const{
  SpMat tmp = M1.selfadjointView<Eigen::Lower>();
  return (tmp.norm());
}

double Frobenius::operator()(const SpMat& M1, const SpMat& M2) const{
  unsigned int n = M1.cols();
  SpMat tmp(n,n);
  tmp = (M1-M2).selfadjointView<Lower>();
  return (tmp.norm());
}

// FROBENIUS SCALED
double FrobeniusScaled::norm(const SpMat& M) const{

  SpMat MM = M.selfadjointView<Eigen::Lower>();
  Eigen::SimplicialLDLT<SpMat,Lower> solver(_Sigma);
  unsigned int n(_Sigma.rows());
  SpMat Id(n,n);
  Id.setIdentity();
  MatrixXd SigmaInv(n,n);
  SigmaInv = solver.solve(Id);
  MatrixXd tmp(n, MM.cols());
  tmp = SigmaInv*MM*SigmaInv*MM;
  return (sqrt(tmp.trace()));
}

double FrobeniusScaled::operator()(const SpMat& M1, const SpMat& M2) const{
  return (FrobeniusScaled::norm((M1-M2)));
}


DistanceTplane::DistanceTplane(const std::string & distanceTplane, const SpMat& Sigma):_distanceTplane(distanceTplane){
  _distances.insert(std::pair<std::string, std::function<double(const SpMat&, const SpMat&)> >("Frobenius", Frobenius()));
  _distances.insert(std::pair<std::string, std::function<double(const SpMat&, const SpMat&)> >("FrobeniusScaled", FrobeniusScaled(Sigma)));

}

double DistanceTplane::compute_distance(const SpMat& M1, const SpMat& M2) const {
  double result = _distances.at(_distanceTplane)(M1, M2);
  return result;
}
