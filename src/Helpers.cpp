#include "Helpers.hpp"

using namespace matrix_manipulation;

// EXPMAT
MatrixXd matrix_manipulation::expMat(const MatrixXd& A) {
  unsigned int p(A.cols());
  EigenSolver<MatrixXd> eigensolver(p);
  eigensolver.compute(A);
  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  VectorXd expvalues(p);

  for (size_t i = 0; i < p; i++) expvalues(i) = std::exp(eigenvalues(i).real());

  MatrixXd result(p, p);
  result =  eigenvectors.real()*expvalues.asDiagonal()*eigenvectors.real().transpose();

  return result;
}

// LOGMAT
MatrixXd matrix_manipulation::logMat(const MatrixXd& A) {

  unsigned int p(A.cols());
  EigenSolver<MatrixXd> eigensolver(p);
  eigensolver.compute(A);

  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  VectorXd logvalues(p);

  for (size_t i = 0; i < p; i++) logvalues(i) = std::log(eigenvalues(i).real());

  MatrixXd result(p, p);
  result =  eigenvectors.real()*logvalues.asDiagonal()*eigenvectors.real().transpose();

  return result;
}

// SQRTMAT
MatrixXd matrix_manipulation::sqrtMat(const MatrixXd& A) {

  unsigned int p(A.cols());
  EigenSolver<MatrixXd> eigensolver(p);
  eigensolver.compute(A);
  VectorXcd eigenvalues =  eigensolver.eigenvalues();
  MatrixXcd eigenvectors = eigensolver.eigenvectors();

  VectorXd sqrtvalues(p);
  for (size_t i = 0; i < p; i++) sqrtvalues(i) = sqrt(eigenvalues(i).real());

  MatrixXd result(p, p);
  result =  eigenvectors.real()*sqrtvalues.asDiagonal()*eigenvectors.real().transpose();

  return result;
}

// BIG MATRIX -> VEC of MATRICES
std::vector<MatrixXd> matrix_manipulation::bigMatrix2VecMatrices(const MatrixXd& bigMatrix, unsigned int p){
  unsigned int N(bigMatrix.rows());
  std::vector<MatrixXd> result(N);
  unsigned int k;
  for(size_t l=0; l<N; l++){
    result[l].resize(p,p);
    k = 0;
    for(size_t i=0; i<p; i++){
      result[l](i,i) = bigMatrix(l,k);
      k++;
      for(size_t j=i+1; j<p; j++){
        result[l](i,j) = bigMatrix(l,k);
        result[l](j,i) = bigMatrix(l,k);
        k++;
      }
    }
  }
  return(result);
}

// VEC of MATRICES -> BIG MATRIX
MatrixXd matrix_manipulation::VecMatrices2bigMatrix(const std::vector<MatrixXd>& vecMatrices) {
  unsigned int N (vecMatrices.size());
  unsigned int p((vecMatrices[0]).rows());
  MatrixXd bigMatrix(N, ((p+1)*p)/2);
  unsigned int k;
  for (size_t l=0; l<N;l++) {
    k=0;
    for (size_t i=0; i<p; i++) {
      for(size_t j=i; j<p; j++) {
        bigMatrix(l,k) = (vecMatrices[l])(i,j);
        k++;
      }
    }
  }
  return bigMatrix;
}
