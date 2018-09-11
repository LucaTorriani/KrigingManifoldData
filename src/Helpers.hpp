#ifndef _HELPERS_HPP_
#define _HELPERS_HPP_

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <cmath>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::VectorXd Vec;
typedef Eigen::Triplet<double> TripType;

using namespace Eigen;

namespace matrix_manipulation {

  SpMat expMat(const SpMat&);
  SpMat expMat(const Eigen::MatrixXd&);

  SpMat logMat(const SpMat&);
  SpMat logMat(const Eigen::MatrixXd&);

  SpMat sqrtMat(const SpMat&);
  SpMat sqrtMat(const Eigen::MatrixXd&);
};


#endif
