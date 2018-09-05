#include "Helpers.hpp"

SpMat matrix_manipulation::expMat(const SpMat &A) {

  SelfAdjointEigenSolver<SpMat> eigensolver (A, ComputeEigenvectors);
  VectorXd eigenvalues =  eigensolver.eigenvalues();
  auto eigenvectors = eigensolver.eigenvectors();

  for (auto i = 0; i < eigenvalues.size(); i++) eigenvalues[i] = exp(eigenvalues[i]);

  unsigned int n = eigenvectors.rows();

  MatrixXd tmp(n, n);
  tmp =  eigenvectors*eigenvalues.asDiagonal()*eigenvectors.transpose();

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

SpMat matrix_manipulation::logMat(const SpMat &A) {

  SelfAdjointEigenSolver<SpMat> eigensolver (A, ComputeEigenvectors);
  VectorXd eigenvalues =  eigensolver.eigenvalues();
  auto eigenvectors = eigensolver.eigenvectors();

  for (auto i = 0; i < eigenvalues.size(); i++) eigenvalues[i] = log(eigenvalues[i]);

  unsigned int n = eigenvectors.rows();

  MatrixXd tmp(n, n);
  tmp =  eigenvectors*eigenvalues.asDiagonal()*eigenvectors.transpose();

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


SpMat matrix_manipulation::sqrtMat(const SpMat &A) {

  SelfAdjointEigenSolver<SpMat> eigensolver (A, ComputeEigenvectors);
  VectorXd eigenvalues =  eigensolver.eigenvalues();
  auto eigenvectors = eigensolver.eigenvectors();

  for (auto i = 0; i < eigenvalues.size(); i++) eigenvalues[i] = sqrt(eigenvalues[i]);

  unsigned int n = eigenvectors.rows();

  MatrixXd tmp(n, n);
  tmp =  eigenvectors*eigenvalues.asDiagonal()*eigenvectors.transpose();

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
