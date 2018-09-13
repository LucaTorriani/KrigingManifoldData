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
  result =  eigenvectors.real()*expvalues.asDiagonal()*eigenvectors.real().transpose(); // not tested

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

std::vector<MatrixXd> bigMatrix2VecMatrices(const MatrixXd& bigMatrix, unsigned int n){
  unsigned int N(bigMatrix.rows());
  std::vector<MatrixXd> result(N);
  unsigned int k;
  for(size_t l=0; l<N; l++){
    result[l].resize(n,n);
    k = 0;
    for(size_t i=0; i<n; i++){
      result[l](i,i) = bigMatrix(l,k);
      k++;
      for(size_t j=i; j<n; j++){
        result[l](i,j) = bigMatrix(l,k);
        result[l](j,i) = bigMatrix(l,k);
        k++;
      }
    }
  }
  return(result);
};

MatrixXd VecMatrices2bigMatrix(const std::vector<MatrixXd>& vecMatrices) {
  unsigned int N (vecMatrices.size());
  unsigned int n((vecMatrices[0]).rows());
  MatrixXd bigMatrix(N, ((n+1)*n)/2);
  unsigned int k;
  for (size_t l=0; l<N;l++) {
    k=0;
    for (size_t i=0; i<n; i++) {
      for(size_t j=i; j<n; j++) {
        bigMatrix(l,k) = (vecMatrices[l])(i,j);
        k++;
      }
    }
  }
  return bigMatrix;
};
