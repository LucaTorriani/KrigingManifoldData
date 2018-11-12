#ifndef _DISTANCE_MANIFOLD_HPP_
#define _DISTANCE_MANIFOLD_HPP_

#include "Helpers.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>
#include <memory>

/*! \file
  @brief Classes to compute the distance on the manifold according to its metric
*/

using namespace Eigen;
namespace distances_manifold{

  /*!
    @brief	Abstract class for the computation of the distance on the manifold
  */
  class DistanceManifold{
  public:
    /*!
      @brief Compute the distance between two matrices on the manifold
      @param M1 First matrix
      @param M2 Second matrix
      @return Distance
    */
    virtual double compute_distance(const MatrixXd& M1, const MatrixXd& M2) const = 0;
    /*!
      @brief Destructor
    */
    virtual ~DistanceManifold() = default;
  };

  /*!
    @brief	Class for the computation of the distance on the manifold when \f$\texttt{manifold\_metric=="Frobenius"}\f$
  */
  class Frobenius : public DistanceManifold{
  public:
    double compute_distance(const MatrixXd& M1, const MatrixXd& M2) const override;
    /*!
      @brief Destructor
    */
    ~Frobenius() = default;
    };

    /*!
      @brief	Class for the computation of the distance on the manifold when \f$\texttt{manifold\_metric=="LogEuclidean"}\f$
    */
  class LogEuclidean : public DistanceManifold{
  public:
    double compute_distance(const MatrixXd& M1, const MatrixXd& M2) const override;
    /*!
      @brief Destructor
    */
    ~LogEuclidean() = default;
  };

  /*!
    @brief	Class for the computation of the distance on the manifold when \f$\texttt{manifold\_metric=="SqRoot"}\f$
  */
  class SqRoot : public DistanceManifold{
  public:
    double compute_distance(const MatrixXd& M1, const MatrixXd& M2) const override;
    /*!
      @brief Destructor
    */
    ~SqRoot() = default;
  };

  /*!
    @brief	Class for the computation of the distance on the manifold when \f$\texttt{manifold\_metric=="Chol"}\f$
    @note The data on the manifold must be correlation matrices
  */
  class Chol : public DistanceManifold{
  public:
    double compute_distance(const MatrixXd& M1, const MatrixXd& M2) const override;
    /*!
      @brief Destructor
    */
    ~Chol() = default;
  };

}



#endif
