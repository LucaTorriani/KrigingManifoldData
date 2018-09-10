#include "Helpers.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
using namespace Eigen;

// EXPMAT
SpMat matrix_manipulation::expMat(const SpMat &A) {

  SelfAdjointEigenSolver<SpMat> eigensolver (A, ComputeEigenvectors);
  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  unsigned int n = eigenvectors.rows();

  VectorXd expvalues(n);

  for (size_t i = 0; i < n; i++) expvalues(i) = exp(eigenvalues(i).real());


  MatrixXd tmp(n, n);
  tmp =  eigenvectors.real()*expvalues.asDiagonal()*eigenvectors.real().transpose();

  std::vector<TripType> tripletList;
  tripletList.reserve(n*(n+1)/2);
  for (size_t j=0; j < n; j++ ) {
    for (size_t i=j; i < n; i++ ) {
      tripletList.push_back(TripType(i,j,tmp(i,j)));
    }
  }

  SpMat result(n, n);
  result.setFromTriplets(tripletList.begin(), tripletList.end());
  return result;
};

SpMat matrix_manipulation::expMat(const MatrixXd& A) {
  unsigned int n(A.cols());
  EigenSolver<MatrixXd> eigensolver(n);
  eigensolver.compute(A);
  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  VectorXd expvalues(n);

  for (size_t i = 0; i < n; i++) expvalues(i) = exp(eigenvalues(i).real());

  MatrixXd tmp(n, n);
  tmp =  eigenvectors.real()*expvalues.asDiagonal()*eigenvectors.real().transpose();

  std::vector<TripType> tripletList;
  tripletList.reserve(n*(n+1)/2);
  for (size_t j=0; j < n; j++ ) {
    for (size_t i=j; i < n; i++ ) {
      tripletList.push_back(TripType(i,j,tmp(i,j)));
    }
  }

  SpMat result(n, n);
  result.setFromTriplets(tripletList.begin(), tripletList.end());
  return result;
};


// LOMAT
SpMat matrix_manipulation::logMat(const SpMat &A) {

  SelfAdjointEigenSolver<SpMat> eigensolver (A, ComputeEigenvectors);
  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  unsigned int n = eigenvectors.rows();

  VectorXd logvalues(n);

  for (size_t i = 0; i < n; i++) logvalues(i) = log(eigenvalues(i).real());


  MatrixXd tmp(n, n);
  tmp =  eigenvectors.real()*logvalues.asDiagonal()*eigenvectors.real().transpose();

  std::vector<TripType> tripletList;
  tripletList.reserve(n*(n+1)/2);
  for (size_t j=0; j < n; j++ ) {
    for (size_t i=j; i < n; i++ ) {
      tripletList.push_back(TripType(i,j,tmp(i,j)));
    }
  }

  SpMat result(n, n);
  result.setFromTriplets(tripletList.begin(), tripletList.end());
  return result;
};

SpMat matrix_manipulation::logMat(const MatrixXd& A) {

  unsigned int n(A.cols());
  EigenSolver<MatrixXd> eigensolver(n);
  eigensolver.compute(A);

  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  VectorXd logvalues(n);

  for (size_t i = 0; i < n; i++) logvalues(i) = log(eigenvalues(i).real());


  MatrixXd tmp(n, n);
  tmp =  eigenvectors.real()*logvalues.asDiagonal()*eigenvectors.real().transpose();

  std::vector<TripType> tripletList;
  tripletList.reserve(n*(n+1)/2);
  for (size_t j=0; j < n; j++ ) {
    for (size_t i=j; i < n; i++ ) {
      tripletList.push_back(TripType(i,j,tmp(i,j)));
    }
  }

  SpMat result(n, n);
  result.setFromTriplets(tripletList.begin(), tripletList.end());
  return result;
};


// SQRTMAT
SpMat matrix_manipulation::sqrtMat(const SpMat &A) {

  SelfAdjointEigenSolver<SpMat> eigensolver (A, ComputeEigenvectors);
  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  unsigned int n = eigenvectors.rows();

  VectorXd sqrtvalues(n);
  for (size_t i = 0; i < n; i++) sqrtvalues(i) = sqrt(eigenvalues(i).real());


  MatrixXd tmp(n, n);
  tmp =  eigenvectors.real()*sqrtvalues.asDiagonal()*eigenvectors.real().transpose();

  std::vector<TripType> tripletList;
  tripletList.reserve(n*(n+1)/2);
  for (size_t j=0; j < n; j++ ) {
    for (size_t i=j; i < n; i++ ) {
      tripletList.push_back(TripType(i,j,tmp(i,j)));
    }
  }

  SpMat result(n, n);
  result.setFromTriplets(tripletList.begin(), tripletList.end());
  return result;
};


SpMat matrix_manipulation::sqrtMat(const MatrixXd &A) {

  unsigned int n(A.cols());
  EigenSolver<MatrixXd> eigensolver(n);
  eigensolver.compute(A);
  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  VectorXd sqrtvalues(n);
  for (size_t i = 0; i < n; i++) sqrtvalues(i) = sqrt(eigenvalues(i).real());

  MatrixXd tmp(n, n);
  tmp =  eigenvectors.real()*sqrtvalues.asDiagonal()*eigenvectors.real().transpose();

  std::vector<TripType> tripletList;
  tripletList.reserve(n*(n+1)/2);
  for (size_t j=0; j < n; j++ ) {
    for (size_t i=j; i < n; i++ ) {
      tripletList.push_back(TripType(i,j,tmp(i,j)));
    }
  }

  SpMat result(n, n);
  result.setFromTriplets(tripletList.begin(), tripletList.end());
  return result;
};
