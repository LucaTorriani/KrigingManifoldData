#include "DistanceManifold.hpp"
#include <cmath>
#include <iostream>

using namespace distances_manifold;

// FROBENIUS
double Frobenius::compute_distance(const MatrixXd& M1, const MatrixXd& M2) const{
  unsigned int n(M1.cols());
  Eigen::LDLT<MatrixXd> solver(n);
  solver.compute(M1);
  MatrixXd matrix_result(n,n);
  matrix_result = solver.solve(M2);

  VectorXcd eigenvalues =  matrix_result.eigenvalues();

  double ssq = 0.0;
  for(auto i = 0;i < eigenvalues.size(); i++)  ssq += (std::log(eigenvalues(i).real())*std::log(eigenvalues(i).real()));
  return (std::sqrt(ssq));
}

// LOGEUCLIDEAN
double LogEuclidean::compute_distance(const MatrixXd& M1, const MatrixXd& M2) const {
  unsigned int n = M1.cols();
  MatrixXd tmp(n,n);
  tmp =  (matrix_manipulation::logMat(M1)-matrix_manipulation::logMat(M2));
  return ( tmp.norm());}

// SQROOT
double SqRoot::compute_distance(const MatrixXd& M1, const MatrixXd& M2) const {
  unsigned int n = M1.cols();
  MatrixXd tmp(n,n);
  tmp = matrix_manipulation::sqrtMat(M1)-matrix_manipulation::sqrtMat(M2);
  return ( tmp.norm());
}
