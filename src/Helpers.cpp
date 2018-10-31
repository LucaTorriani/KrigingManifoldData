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
std::vector<MatrixXd> matrix_manipulation::bigMatrix2VecMatrices(const MatrixXd& bigMatrix, unsigned int p, const std::string& distance_Manifold_name){
  unsigned int N(bigMatrix.rows());
  std::vector<MatrixXd> result(N);

  if (distance_Manifold_name == "Correlation") {
    unsigned int k;
    for(size_t l=0; l<N; l++){
      result[l].resize(p,p);
      result[l].setZero();
      k = 0;
      for(size_t i=0; i<p; i++){
        result[l](i,i) = bigMatrix(l,k);
        k++;
        for(size_t j=i+1; j<p; j++){
          result[l](i,j) = bigMatrix(l,k);
          k++;
        }
      }
    }
  }
  else {
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

// CHOLESKY FACTORIZATIOON for a positive semi-definite matrix
MatrixXd matrix_manipulation::Chol_semidef (const MatrixXd& M1) {
  unsigned int p(M1.rows());
  MatrixXd result(p,p);
  result.setZero(p,p);
  result(0,0)=1;
  double tmp;
  for (size_t i=1; i<p; i++) {
    for (size_t j=0; j<i; j++) {
      if (j==0) result(i,j) = M1(i,j);
      else result(i,j) = M1(i,j) - ((result.row(i) * result.row(j).transpose()).value())/result(j,j);
    }
    tmp=result.row(i).norm();
    result(i,i) = sqrt(1-tmp*tmp);
  }
  return(result.transpose());
}

// CHOLESKY DECOMPOSITION
MatrixXd matrix_manipulation::Chol_decomposition (const MatrixXd M1) {
  unsigned int p(M1.rows());
  MatrixXd H(p,p);
  // H.setZero(p,p);

  // LLT
  LLT<MatrixXd> llt(p);
  llt.compute(M1);
  if (llt.info() == Success) {
    // LLT
    H = llt.matrixU();
  }
  else {
    // PIGOLI
    H = Chol_semidef(M1);

    // LDLT
    // std::cout << "LDLT: " << std::endl;
    // LDLT<MatrixXd> ldlt(p);
    // ldlt.compute(H_cor);
    // // OPZIONE 1 Metto a zero (rischio di non ricostruire quella parte di matrice)
    // // Vec sqrt_diag_H = ldlt.vectorD();
    // // for (size_t i=0; i<p; i++) {
    // //   if(sqrt_diag_H(i) > 0) sqrt_diag_H(i) = sqrt(sqrt_diag_H(i));
    // //   else sqrt_diag_H(i) = 0;
    // // }
    //
    // // OPZIONE 2 Non metto a zero (rischio ottenere Nan)
    // Vec sqrt_diag_H = sqrt(ldlt.vectorD().array());
    //
    // MatrixXd U_H = ldlt.matrixU();
    // H = (U_H.array()).rowwise() * sqrt_diag_H.transpose().array();
  }
  return(H);
}
