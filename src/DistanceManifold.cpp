#include "DistanceManifold.hpp"
#include <cmath>
#include <iostream>

using namespace distances_manifold;

// FROBENIUS
double Frobenius::compute_distance(const MatrixXd& M1, const MatrixXd& M2) const{
  unsigned int p(M1.cols());
  Eigen::LDLT<MatrixXd> solver(p);
  solver.compute(M1);
  MatrixXd matrix_result(p,p);
  matrix_result = solver.solve(M2);

  VectorXcd eigenvalues =  matrix_result.eigenvalues();

  double ssq = 0.0;
  for(auto i = 0;i < eigenvalues.size(); i++)  ssq += (std::log(eigenvalues(i).real())*std::log(eigenvalues(i).real()));
  return (std::sqrt(ssq));
}

// LOGEUCLIDEAN
double LogEuclidean::compute_distance(const MatrixXd& M1, const MatrixXd& M2) const {
  unsigned int p = M1.cols();
  MatrixXd tmp(p,p);
  tmp =  (matrix_manipulation::logMat(M1)-matrix_manipulation::logMat(M2));
  return ( tmp.norm());
}

// SQROOT
double SqRoot::compute_distance(const MatrixXd& M1, const MatrixXd& M2) const {
  unsigned int p = M1.cols();
  MatrixXd tmp(p,p);
  tmp = matrix_manipulation::sqrtMat(M1)-matrix_manipulation::sqrtMat(M2);
  return ( tmp.norm());
}

// CORRELATION
double Chol::compute_distance(const MatrixXd& M1, const MatrixXd& M2) const {
  unsigned int p(M1.rows());
  // MatrixXd H1(p,p);    // Assuming data_manifold are already Chol
  // MatrixXd H2(p,p);
  // H1 = Chol_decomposition(M1);
  // H2 = Chol_decomposition(M2);
  double result(0);
  double tmp;
  for (size_t i=1; i<p; i++) {
    tmp = ((M1.col(i)).transpose() * (M2.col(i))).value();  // .value() casts a matrix(1,1) to a double
    tmp = acos(std::max( std::min(tmp, 1.0), -1.0));
    result += tmp*tmp;
  }
  return (sqrt(result));
}
