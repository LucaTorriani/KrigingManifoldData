#ifndef _MAPFUNCTIONS_HPP_
#define _MAPFUNCTIONS_HPP_
#include <memory>

#include "DistanceManifold.hpp"

namespace map_functions {

// LOGARITHMIC MAP
  class logarithmicMap{
  public:
    virtual ~logarithmicMap() = default;
    virtual MatrixXd map2tplane(const MatrixXd&) const = 0;
    virtual void set_members(const MatrixXd&) = 0;

  };

  class logMapFrob : public logarithmicMap{
    MatrixXd  _sqrtSigma;
    MatrixXd  _sqrtSigmaInv;
  public:
     ~logMapFrob() = default;
    MatrixXd map2tplane(const MatrixXd&) const override;
    void set_members(const MatrixXd&) override;
  };

  class logMapLogEucl : public logarithmicMap{
    MatrixXd _Sigma;
  public:
    ~logMapLogEucl() = default;
    MatrixXd map2tplane(const MatrixXd&) const override;
    void set_members(const MatrixXd&) override;

  };

  class logMapSqRoot : public logarithmicMap{
    MatrixXd _Sigma;
  public:
    ~logMapSqRoot() = default;
    MatrixXd map2tplane(const MatrixXd&) const override;
    void set_members(const MatrixXd&) override;
  };

  class logMapChol : public logarithmicMap{
    MatrixXd _Sigma;
    Vec proj2tspace(const Vec&, const Vec&) const;
  public:
    ~logMapChol() = default;
    MatrixXd map2tplane(const MatrixXd&) const override;
    void set_members(const MatrixXd&) override;
  };


  // EXPONENTIAL MAP

  class exponentialMap{
  public:
    virtual ~exponentialMap() = default;
    virtual MatrixXd map2manifold(const MatrixXd&) const = 0;
    virtual void set_members(const MatrixXd&) = 0;
  };


  class expMapFrob : public exponentialMap{
    MatrixXd _sqrtSigma;
    MatrixXd _sqrtSigmaInv;
  public:
    ~expMapFrob() = default;
    MatrixXd map2manifold(const MatrixXd&) const override;
    void set_members(const MatrixXd&) override;

  };

  class expMapLogEucl : public exponentialMap{
    MatrixXd _Sigma;
  public:
    ~expMapLogEucl() = default;
    MatrixXd map2manifold(const MatrixXd&) const override;
    void set_members(const MatrixXd&) override;
  };

  class expMapSqRoot : public exponentialMap{
    MatrixXd _Sigma;
  public:
    ~expMapSqRoot() = default;
    MatrixXd map2manifold(const MatrixXd&) const override;
    void set_members(const MatrixXd&) override;
  };

  class expMapChol : public exponentialMap{
    MatrixXd _Sigma;
  public:
    ~expMapChol() = default;
    MatrixXd map2manifold(const MatrixXd&) const override;
    void set_members(const MatrixXd&) override;
  };


}

#endif
