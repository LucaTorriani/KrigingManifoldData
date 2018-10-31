#ifndef _HELPERS_HPP_
#define _HELPERS_HPP_

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <cmath>
#include <string>

typedef Eigen::VectorXd Vec;
typedef Eigen::Triplet<double> TripType;

using namespace Eigen;

namespace matrix_manipulation {

  MatrixXd expMat(const MatrixXd&);
  MatrixXd logMat(const MatrixXd&);
  MatrixXd sqrtMat(const MatrixXd&);

  std::vector<MatrixXd> bigMatrix2VecMatrices(const MatrixXd&, unsigned int, const std::string&);
  MatrixXd VecMatrices2bigMatrix(const std::vector<MatrixXd>&);
  MatrixXd Chol_semidef (const MatrixXd& M1);
  MatrixXd Chol_decomposition (const MatrixXd M1);

}


#endif
