#include "DistanceManifold.hpp"
#include <cmath>
#include <iostream>

using namespace distances_manifold;

// FROBENIUS
double Frobenius::operator()(const MatrixXd& M1, const MatrixXd& M2) const{
  unsigned int n(M1.cols());
  Eigen::LDLT<MatrixXd> solver(n); // Piu veloce specificando prima la dimensione
  solver.compute(M1);
  MatrixXd matrix_result(n,n);
  matrix_result = solver.solve(M2);

  VectorXcd eigenvalues =  matrix_result.eigenvalues();

  double ssq = 0.0;
  for(auto i = 0;i < eigenvalues.size(); i++)  ssq += (std::log(eigenvalues(i).real())*std::log(eigenvalues(i).real()));
  return (std::sqrt(ssq));
}

// LOGEUCLIDEAN
double LogEuclidean::operator()(const MatrixXd& M1, const MatrixXd& M2) const {
  unsigned int n = M1.cols();
  MatrixXd tmp(n,n);
  tmp =  (matrix_manipulation::logMat(M1)-matrix_manipulation::logMat(M2));
  return ( tmp.norm());}

// SQROOT
double SqRoot::operator()(const MatrixXd& M1, const MatrixXd& M2) const {
  unsigned int n = M1.cols();
  MatrixXd tmp(n,n);
  tmp = matrix_manipulation::sqrtMat(M1)-matrix_manipulation::sqrtMat(M2);
  return ( tmp.norm());
}

// CLASS DISTANCE MANIFOLD

DistanceManifold::DistanceManifold(const std::string& distanceManifold, const std::shared_ptr<const MatrixXd> Sigma):_distanceManifold(distanceManifold),
  _Sigma(Sigma)

{
  distances.insert(std::pair<std::string, std::function<double(const MatrixXd&, const MatrixXd&)>> ("Frobenius", Frobenius()));
  distances.insert(std::pair<std::string, std::function<double(const MatrixXd&, const MatrixXd&)>>("LogEuclidean", LogEuclidean()));
  distances.insert(std::pair<std::string, std::function<double(const MatrixXd&, const MatrixXd&)>>("SquareRoot", SqRoot()));
}

double DistanceManifold::compute_distance (const MatrixXd& M1, const MatrixXd& M2) const {
  double result = (distances.at(_distanceManifold))(M1, M2);
  return result;
}

const std::shared_ptr<const MatrixXd> DistanceManifold::get_Sigma() const {
  return _Sigma;
}

const std::string& DistanceManifold::get_distanceType() const{
  return _distanceManifold;
}
