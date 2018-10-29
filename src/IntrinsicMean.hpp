#ifndef _INTRINSIC_MEAN_HPP_
#define _INTRINSIC_MEAN_HPP_

#include <Eigen/Dense>
#include <memory>
#include "DistanceTplane.hpp"
#include "Helpers.hpp"
#include "MapFunctions.hpp"

using namespace Eigen;

MatrixXd intrinsic_mean_C(const std::vector<MatrixXd>&, std::string, map_functions::logarithmicMap&,
                         map_functions::exponentialMap&,   distances_tplane::DistanceTplane&, double, const Vec&, const Vec&);

MatrixXd extrinsic_mean (const std::vector<MatrixXd>&, const Vec&);


#endif
