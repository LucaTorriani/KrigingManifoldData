#ifndef _DISTANCE_TPLANE_HPP_
#define _DISTANCE_TPLANE_HPP_

#include "Helpers.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>
#include <memory>

/*! \file
  @brief Classes to compute the distance on the tangent space according to its metric
*/

namespace distances_tplane{

  /*!
    @brief	Abstract class for the computation of the distance on the tangent space
  */
  class DistanceTplane{
  public:
    /*!
      @brief Compute the distance between two matrices on the tangent space
      @param M1 First matrix
      @param M2 Second matrix
      @return Distance
    */
    double compute_distance(const MatrixXd& M1, const MatrixXd& M2) const;

    /*!
      @brief Compute the norm of a matrix on the tangent space
      @param M1 Matrix
      @return Norm
    */
    virtual double norm(const MatrixXd& M1) const = 0;

    /*!
      @brief Set the members that will be used in the computation of the norms and distances, according to the metric on the tangent space
      @param Sigma Tangent point
    */
    virtual void set_members(const MatrixXd& Sigma) = 0;
    /*!
      @brief Destructor
    */
    virtual ~DistanceTplane() = default;
  };

  /*!
    @brief	Class for the computation of the distance on the tangent space when `ts_metric=="Frobenius"`
  */
 class Frobenius : public DistanceTplane{
 public:
   /*!
     @brief Destructor
   */
   ~Frobenius() = default;
   double norm(const MatrixXd & M1) const override;
   void set_members(const MatrixXd& Sigma) override;

 };

 /*!
   @brief	Class for the computation of the distance on the tangent space when `ts_metric=="FrobeniusScaled"`
 */
 class FrobeniusScaled : public DistanceTplane{
   /*! Inverse of the tangent point Sigma */
   MatrixXd _SigmaInv;
   /*! Dimension of the matrices on the tangent space */
   unsigned int _p;
 public:
   /*!
     @brief Destructor
   */
   ~FrobeniusScaled() = default;
   double norm(const MatrixXd & M1) const override;
   void set_members(const MatrixXd& Sigma) override;

};

/*!
  @brief	Class for the computation of the distance on the tangent space when `ts_metric=="Correlation"`
*/
class Chol : public DistanceTplane{
public:
  /*!
    @brief Destructor
  */
  ~Chol() = default;
  double norm(const MatrixXd & M1) const override;
  void set_members(const MatrixXd& Sigma) override;
};



}



#endif // _DISTANCE_TPLANE_HPP_
