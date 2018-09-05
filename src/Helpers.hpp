#ifndef _HELPERS_HPP_
#define _HELPERS_HPP_

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::VectorXd Vec;
typedef Eigen::Triplet<double> TripType;

namespace matrix_manipulation {

  SpMat expMat(const SpMat&);
  SpMat logMat(const SpMat&);
  SpMat sqrtMat(const SpMat&)

}


#endif
