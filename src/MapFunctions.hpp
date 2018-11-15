#ifndef _MAPFUNCTIONS_HPP_
#define _MAPFUNCTIONS_HPP_
#include <memory>

#include "DistanceManifold.hpp"

/*! \file
  @brief Classes to compute the exponential and logarithmic map according to the manifold metric
*/

namespace map_functions {

// LOGARITHMIC MAP
/*!
  @brief	Abstract class for the computation of the logarithmic map
*/
  class logarithmicMap{
  public:
    /*!
      @brief Destructor
    */
    virtual ~logarithmicMap() = default;
    /*!
      @brief Map a manifold matrix to the tangent space
      @param M Manifold matrix to map
      @return Tangent space matrix identifying the mapped data
    */
    virtual MatrixXd map2tplane(const MatrixXd& M) const = 0;
    /*!
      @brief Set class members
      @param Sigma Tangent point
    */
    virtual void set_members(const MatrixXd& Sigma) = 0;
    /*!
      @brief Set tolerance
      @param tolerance_map_cor Tolerance
    */
    virtual void set_tolerance(double tolerance_map_cor) = 0;
  };
  /*!
    @brief	Class for the computation of the logarithmic map when `manifold_metric=="Frobenius"`
  */
  class logMapFrob : public logarithmicMap{
    /*! Square root of the tangent point `Sigma` */
    MatrixXd  _sqrtSigma;
    /*! Inverse of the square root of the tangent point `Sigma` */
    MatrixXd  _sqrtSigmaInv;
  public:
    /*!
      @brief Destructor
    */
     ~logMapFrob() = default;
    MatrixXd map2tplane(const MatrixXd& M) const override;
    void set_members(const MatrixXd& Sigma) override;
    void set_tolerance(double tolerance_map_cor) override;
  };
  /*!
    @brief	Class for the computation of the logarithmic map when `manifold_metric=="LogEuclidean"`
  */
  class logMapLogEucl : public logarithmicMap{
    /*! Tangent point `Sigma` */
    MatrixXd _Sigma;
  public:
    /*!
      @brief Destructor
    */
    ~logMapLogEucl() = default;
    MatrixXd map2tplane(const MatrixXd& M) const override;
    void set_members(const MatrixXd& Sigma) override;
    void set_tolerance(double tolerance_map_cor) override;
  };
  /*!
    @brief	Class for the computation of the logarithmic map when `manifold_metric=="SqRoot"`
  */
  class logMapSqRoot : public logarithmicMap{
    /*! Tangent point `Sigma` */
    MatrixXd _Sigma;
  public:
    /*!
      @brief Destructor
    */
    ~logMapSqRoot() = default;
    MatrixXd map2tplane(const MatrixXd& M) const override;
    void set_members(const MatrixXd& Sigma) override;
    void set_tolerance(double tolerance_map_cor) override;
  };
  /*!
    @brief	Class for the computation of the logarithmic map when `manifold_metric=="Correlation"`
  */
  class logMapChol : public logarithmicMap{
    /*! Tangent point `Sigma` */
    MatrixXd _Sigma;
    /*! Tolerance on the norm of the columns to avoid `Nan` */
    double _tolerance_map_cor;
    /*!
      @brief Project a matrix in \f$ Chol(p) \f$ to the tangent space 
    */
    Vec proj2tspace(const Vec&, const Vec&) const;
  public:
    /*!
      @brief Destructor
    */
    ~logMapChol() = default;
    MatrixXd map2tplane(const MatrixXd& M) const override;
    void set_members(const MatrixXd& Sigma) override;
    void set_tolerance(double tolerance_map_cor) override;
  };


  // EXPONENTIAL MAP
  /*!
    @brief	Abstract class for the computation of the exponential map
  */
  class exponentialMap{
  public:
    /*!
      @brief Destructor
    */
    virtual ~exponentialMap() = default;
    /*!
      @brief Map a tangent space matrix to the manifold
      @param M Tangent space matrix to map
      @return Manifold matrix identifying the mapped data
    */
    virtual MatrixXd map2manifold(const MatrixXd& M) const = 0;
    /*!
      @brief Set class members
      @param Sigma Tangent point
    */
    virtual void set_members(const MatrixXd& Sigma) = 0;
    /*!
      @brief Set tolerance
      @param tolerance_map_cor Tolerance
    */
    virtual void set_tolerance(double tolerance_map_cor) = 0;
  };

  /*!
    @brief	Class for the computation of the exponential map when `manifold_metric=="Frobenius"`
  */
  class expMapFrob : public exponentialMap{
    /*! Square root of the tangent point `Sigma` */
    MatrixXd _sqrtSigma;
    /*! Inverse of the square root of the tangent point `Sigma` */
    MatrixXd _sqrtSigmaInv;
  public:
    /*!
      @brief Destructor
    */
    ~expMapFrob() = default;
    MatrixXd map2manifold(const MatrixXd& M) const override;
    void set_members(const MatrixXd& Sigma) override;
    void set_tolerance(double tolerance_map_cor) override;
  };
  /*!
    @brief	Class for the computation of the exponential map when `manifold_metric=="LogEuclidean"`
  */
  class expMapLogEucl : public exponentialMap{
    /*! Tangent point `Sigma` */
    MatrixXd _Sigma;
  public:
    /*!
      @brief Destructor
    */
    ~expMapLogEucl() = default;
    MatrixXd map2manifold(const MatrixXd& M ) const override;
    void set_members(const MatrixXd& Sigma) override;
    void set_tolerance(double tolerance_map_cor) override;
  };
  /*!
    @brief	Class for the computation of the exponential map when `manifold_metric=="SqRoot"`
  */
  class expMapSqRoot : public exponentialMap{
    /*! Tangent point `Sigma` */
    MatrixXd _Sigma;
  public:
    /*!
      @brief Destructor
    */
    ~expMapSqRoot() = default;
    MatrixXd map2manifold(const MatrixXd& M) const override;
    void set_members(const MatrixXd& Sigma) override;
    void set_tolerance(double tolerance_map_cor) override;
  };
  /*!
    @brief	Class for the computation of the exponential map when `manifold_metric=="Correlation"`
  */
  class expMapChol : public exponentialMap{
    /*! Tolerance on the norm of the columns to avoid `Nan` */
    double _tolerance_map_cor;
    /*! Tangent point `Sigma` */
    MatrixXd _Sigma;
  public:
    /*!
      @brief Destructor
    */
    ~expMapChol() = default;
    MatrixXd map2manifold(const MatrixXd& M) const override;
    void set_members(const MatrixXd& Sigma) override;
    void set_tolerance(double tolerance_map_cor) override;
  };


}

#endif  // _MAPFUNCTIONS_HPP_
