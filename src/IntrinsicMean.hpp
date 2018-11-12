#ifndef _INTRINSIC_MEAN_HPP_
#define _INTRINSIC_MEAN_HPP_

#include <Eigen/Dense>
#include <memory>
#include "DistanceTplane.hpp"
#include "Helpers.hpp"
#include "MapFunctions.hpp"


/*! \file
  @brief Functions to compute intrinsic and extrinsic mean for manifold data
*/

using namespace Eigen;
/*!
  @brief Compute the intrinsic mean of a set of manifold matrices
  @param data_manifold Vector of manifold matrices
  @param distance_Manifold_name Name of the metric on the manifold
  @param logMap Logarithmic map
  @param expMap Exponential map
  @param distanceTplane Distance on the tangent space
  @param tolerance Tolerance for the algorithm's loop
  @param weight_intrinsic Weights
  @param weight_extrinsic Weights for the computation of the extrinsic mean (used if distance_Manifold_name=="Correlation")
  @return Matrix identifying the intrinsic mean of data_manifold
*/
MatrixXd intrinsic_mean_C(const std::vector<MatrixXd>& data_manifold, std::string distance_Manifold_name,
                          map_functions::logarithmicMap& logMap, map_functions::exponentialMap& expMap,
                          distances_tplane::DistanceTplane& distanceTplane, double tolerance,
                          const Vec& weight_intrinsic, const Vec& weight_extrinsic);

/*!
  @brief Compute the extrinsic mean of a set of matrices in \f$ Chol\left(p\right) \f$
  @param data_manifold Vector of matrices in  \f$ Chol\left(p\right) \f$
  @param weight_extrinsic Weights
  @return Matrix identifying the extrinsic mean of data_manifold
*/
MatrixXd extrinsic_mean (const std::vector<MatrixXd>& data_manifold, const Vec& weight_extrinsic);


#endif // _INTRINSIC_MEAN_HPP_
