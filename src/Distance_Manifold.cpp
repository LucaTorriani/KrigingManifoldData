#include "Distance_Manifold.hpp"

#include <cmath>

#include <iostream>

using namespace distances_manifold;

// FROBENIUS
double Frobenius::manifold_distance(const SpMat& M1, const SpMat& M2 ){

    Eigen::SimplicialLDLT<SpMat,Lower> solver(M1);
    SpMat tmp = M2.selfadjointView<Lower>();
    Eigen::MatrixXd matrix_result(solver.solve(tmp));

    VectorXcd eigenvalues =  matrix_result.eigenvalues();

    double ssq = 0.0;
      for(auto i = 0;i < eigenvalues.size(); i++)  ssq += (std::log(std::abs(eigenvalues(i)))*std::log(std::abs(eigenvalues(i))));

    return (std::sqrt(ssq));
}

double Frobenius::operator()(const SpMat& M1, const SpMat& M2) {
  return (Frobenius::manifold_distance(M1,M2));
}

// LOGEUCLIDEAN
double LogEuclidean::manifold_distance(const SpMat& M1, const SpMat& M2){
  unsigned int n = M1.cols();
  SpMat tmp(n,n);
  tmp =  (matrix_manipulation::logMat(M1)-matrix_manipulation::logMat(M2)).selfadjointView<Lower>();
  return ( tmp.norm());
}

double LogEuclidean::operator()(const SpMat& M1, const SpMat& M2)  {
  return(LogEuclidean::manifold_distance(M1,M2));
}

// SQROOT
double SqRoot::manifold_distance(const SpMat& M1, const SpMat& M2){
  unsigned int n = M1.cols();
  SpMat tmp(n,n);
  tmp = (matrix_manipulation::sqrtMat(M1)-matrix_manipulation::sqrtMat(M2)).selfadjointView<Lower>();
  return ( tmp.norm());
}

double SqRoot::operator()(const SpMat& M1, const SpMat& M2) {
  return(SqRoot::manifold_distance(M1,M2));

}

// CLASS DISTANCE MANIFOLD

DistanceManifold::DistanceManifold(const std::string& distanceManifold, const SpMat& Sigma):_distanceManifold(distanceManifold),
  _Sigma(Sigma)

{
  distances.insert(std::pair<std::string, std::function<double(const SpMat&, const SpMat&)>> ("Frobenius", Frobenius()));
  distances.insert(std::pair<std::string, std::function<double(const SpMat&, const SpMat&)>>("LogEuclidean", LogEuclidean()));
  distances.insert(std::pair<std::string, std::function<double(const SpMat&, const SpMat&)>>("SquareRoot", SqRoot()));
}

double DistanceManifold::compute_distance (const SpMat& M1, const SpMat& M2) {
  double result = distances[_distanceManifold](M1, M2);
  return result;
}

SpMat DistanceManifold::get_Sigma() const {
  return _Sigma;
}

std::string DistanceManifold::get_distanceType() const{
  return _distanceManifold;
}
