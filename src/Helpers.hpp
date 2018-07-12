#ifndef _HELPERS_HPP_
#define _HELPERS_HPP_

#include <Eigen/Sparse>
#include <Eigen/Dense>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::VectorXd Vec;
typedef Triplet<double> TripType;

enum checkDistance {
  EUCLIDEAN,
  GEODIST
};

#endif
