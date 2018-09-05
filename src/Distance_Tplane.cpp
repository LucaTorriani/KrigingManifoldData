#include "Distance_Tplane.hpp"

#include <cmath>


using namespace distances_tplane;

// FROBENIUS

static double Frobenius::norm(const SpMat& M1){
  auto tmp = M1.selfadjointView<Eigen::Lower>();
  return (tmp.norm());
}


static double Frobenius::operator()(const SpMat& M1, const SpMat& M2) const{
  return ((M1-M2).norm());
}

// FROBENIUS SCALED
static double norm(const SpMat& M, const SpMat& Sigma){

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

static double FrobeniusScaled::operator()(const SpMat& M1, const SpMat& M2, const SpMat& Sigma){
  return (FrobeniusScaled::norm((M1-M2), Sigma);
}



DistanceTplane::DistanceTplane(){
  distances.insert(std::pair<std::string, std::function<double(std::vector<double>, std::vector<double>)>>("Frobenius", Frobenius()));
  distances.insert(std::pair<std::string, std::function<double(std::vector<double>, std::vector<double>)>>("FrobeniusScaled", FrobeniusScaled()));

}
double DistanceManifold::compute_distanceconst SpMat& M1, const SpMat& M2, const std::string & distance_type){
  double result = dist[distance_type](M1, M2);
  return result;
}
