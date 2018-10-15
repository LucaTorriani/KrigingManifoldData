#ifndef _HELPERS_HPP_
#define _HELPERS_HPP_

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <cmath>

typedef Eigen::VectorXd Vec;
typedef Eigen::Triplet<double> TripType;

using namespace Eigen;

namespace matrix_manipulation {

  MatrixXd expMat(const MatrixXd&);
  MatrixXd logMat(const MatrixXd&);
  MatrixXd sqrtMat(const MatrixXd&);

  std::vector<MatrixXd> bigMatrix2VecMatrices(const MatrixXd&, unsigned int);
  MatrixXd VecMatrices2bigMatrix(const std::vector<MatrixXd>&);

}


#endif
