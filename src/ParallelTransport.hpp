#ifndef _PARALLEL_TRANSPORT_HPP_
#define _PARALLEL_TRANSPORT_HPP_

#include "Helpers.hpp"

using namespace Eigen;

/*! \file
  @brief Functions to implement parallel transport
*/

namespace parallel_transport {
  /*!
    @brief Transport a symmetric matrix `V` from the plane tangent to the manifold in `Sigma` to the plane tangent in the identity
    @note
      Reference: "Parallel transport on the cone manifold of spd matrices for domain adaptation." \n
      Authors: O. Yair, M. Ben-Chen, R. Talmon \n
      Periodical: IEEE Transactions on Signal Processing, 2019.
    @param Sigma Tangent point of the plane on which `V` lies
    @param V Symmetric matrix that must be transported
    @return The value obtained transporting `V` on the plane tangent in the identity
  */
  MatrixXd transport_to_TI(MatrixXd Sigma, MatrixXd V);
  /*!
    @brief Transport a symmetric matrix `V` from the plane tangent to the manifold in the identity to the plane tangent in `Sigma`
    @note
      Reference: "Parallel transport on the cone manifold of spd matrices for domain adaptation." \n
      Authors: O. Yair, M. Ben-Chen, R. Talmon \n
      Periodical: IEEE Transactions on Signal Processing, 2019.
    @param Sigma Tangent point of the plane where `V` must be transported
    @param V Symmetric matrix that must be transported
    @return The value obtained transporting `V` on the plane tangent in `Sigma`
  */
  MatrixXd transport_from_TI(MatrixXd Sigma, MatrixXd V);
}

#endif
