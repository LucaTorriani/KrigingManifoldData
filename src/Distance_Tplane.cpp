#include "Distance_Tplane.hpp"

#include <cmath>

using namespace Eigen;
using namespace distances_tplane;

// FROBENIUS

double Frobenius::norm(const SpMat& M1){
  SpMat tmp = M1.selfadjointView<Eigen::Lower>();
  return (tmp.norm());
}

double Frobenius::tplane_dist(const SpMat& M1, const SpMat& M2) {
  return ((M1-M2).norm());
}

double Frobenius::operator()(const SpMat& M1, const SpMat& M2, const SpMat & Sigma = SpMat(1,1)) {
  return (Frobenius::tplane_dist(M1,M2));
}

// FROBENIUS SCALED
double FrobeniusScaled::norm(const SpMat& M, const SpMat& Sigma){

  auto MM = M.selfadjointView<Eigen::Lower>();
  Eigen::SimplicialLDLT<SpMat,Lower> solver(Sigma);
  unsigned int n(Sigma.rows());
  SpMat Id(n,n);
  Id.setIdentity();
  MatrixXd SigmaInv(n,n);
  SigmaInv = solver.solve(Id);
  MatrixXd tmp(n, MM.cols());
  tmp = MM*(MM*SigmaInv).transpose()*SigmaInv;
  return (sqrt(tmp.trace()));
}


double FrobeniusScaled::tplane_dist(const SpMat& M1, const SpMat& M2, const SpMat& Sigma){
  return (FrobeniusScaled::norm((M1-M2), Sigma));
}

double FrobeniusScaled::operator()(const SpMat& M1, const SpMat& M2, const SpMat& Sigma){
  return (FrobeniusScaled::tplane_dist(M1,M2, Sigma));
}


DistanceTplane::DistanceTplane(){
  distances.insert(std::pair<std::string, std::function<double(const SpMat&, const SpMat&,const SpMat&)>>("Frobenius", Frobenius()));
  distances.insert(std::pair<std::string, std::function<double(const SpMat&, const SpMat&, const SpMat&)>>("FrobeniusScaled", FrobeniusScaled()));

}

double DistanceTplane::compute_distance( const std::string & distance_type, const SpMat& M1, const SpMat& M2, const SpMat& Sigma = SpMat(1,1)) {
  double result = distances[distance_type](M1, M2, Sigma);
  return result;
}
