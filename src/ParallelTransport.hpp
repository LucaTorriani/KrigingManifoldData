#ifndef _PARALLEL_TRANSPORT_HPP_
#define _PARALLEL_TRANSPORT_HPP_

#include "Helpers.hpp"

using namespace Eigen;

namespace parallel_transport {
  MatrixXd trasport_to_TI(MatrixXd, MatrixXd);
  MatrixXd transport_from_TI(MatrixXd, MatrixXd);
}

#endif
