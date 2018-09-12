#include "Helpers.hpp"

// EXPMAT
MatrixXd matrix_manipulation::expMat(const MatrixXd& A) {
  unsigned int n(A.cols());
  EigenSolver<MatrixXd> eigensolver(n);
  eigensolver.compute(A);
  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  VectorXd expvalues(n);

  for (size_t i = 0; i < n; i++) expvalues(i) = exp(eigenvalues(i).real());

  MatrixXd result(n, n);
  result =  eigenvectors.real()*expvalues.asDiagonal()*eigenvectors.real().transpose();

  return result;
};

// LOMAT
MatrixXd matrix_manipulation::logMat(const MatrixXd& A) {

  unsigned int n(A.cols());
  EigenSolver<MatrixXd> eigensolver(n);
  eigensolver.compute(A);

  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  VectorXd logvalues(n);

  for (size_t i = 0; i < n; i++) logvalues(i) = log(eigenvalues(i).real());

  MatrixXd result(n, n);
  result =  eigenvectors.real()*logvalues.asDiagonal()*eigenvectors.real().transpose();

  return result;
};

// SQRTMAT
MatrixXd matrix_manipulation::sqrtMat(const MatrixXd& A) {

  unsigned int n(A.cols());
  EigenSolver<MatrixXd> eigensolver(n);
  eigensolver.compute(A);
  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  VectorXd sqrtvalues(n);
  for (size_t i = 0; i < n; i++) sqrtvalues(i) = sqrt(eigenvalues(i).real());

  MatrixXd result(n, n);
  result =  eigenvectors.real()*sqrtvalues.asDiagonal()*eigenvectors.real().transpose();

  return result;
};
