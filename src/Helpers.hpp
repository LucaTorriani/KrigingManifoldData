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

  MatrixXd expMat(const MatrixXd&);

  MatrixXd logMat(const MatrixXd&);

  MatrixXd sqrtMat(const Eigen::MatrixXd&);
};


#endif
