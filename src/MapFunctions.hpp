#ifndef _MAPFUNCTIONS_HPP_
#define _MAPFUNCTIONS_HPP_
#include <memory>

#include "DistanceManifold.hpp"

namespace map_functions {

// LOGARITHMIC MAP
  class logMapFrob{
    MatrixXd  _sqrtSigma;
    MatrixXd  _sqrtSigmaInv;
  public:
    logMapFrob(const MatrixXd&, const MatrixXd&);
    MatrixXd operator()(const MatrixXd&) const;
  };

  class logMapLogEucl{
    const std::shared_ptr<const MatrixXd> _Sigma;
  public:
    logMapLogEucl(const std::shared_ptr<const MatrixXd>);
    MatrixXd operator()(const MatrixXd&) const;
  };

  class logMapSqRoot{
    const std::shared_ptr<const MatrixXd> _Sigma;
  public:
    logMapSqRoot(const std::shared_ptr<const MatrixXd>);
    MatrixXd operator()(const MatrixXd&) const;
  };

  class logarithmicMap{
    const std::string& _distanceManifold;
    std::map<std::string, std::function<MatrixXd(const MatrixXd&)>> maps;
  public:
    logarithmicMap(const distances_manifold::DistanceManifold&);
    MatrixXd map2tplane(const MatrixXd&) const;
  };


  // EXPONENTIAL MAP MAP
  class expMapFrob{
    MatrixXd _sqrtSigma;
    MatrixXd _sqrtSigmaInv;
  public:
    expMapFrob(const MatrixXd&, const MatrixXd&);
    MatrixXd operator()(const MatrixXd&) const;
  };

  class expMapLogEucl{
    const std::shared_ptr<const MatrixXd> _Sigma;
  public:
    expMapLogEucl(const std::shared_ptr<const MatrixXd>);
    MatrixXd operator()(const MatrixXd&) const;
  };

  class expMapSqRoot{
    const std::shared_ptr<const MatrixXd> _Sigma;
  public:
    expMapSqRoot(const std::shared_ptr<const MatrixXd>);
    MatrixXd operator()(const MatrixXd&) const;
  };

  class exponentialMap{
    const std::string& _distanceManifold;
    std::map<std::string, std::function<MatrixXd(const MatrixXd&)>> maps;
  public:
    exponentialMap(const distances_manifold::DistanceManifold&);
    MatrixXd map2manifold(const MatrixXd&) const;
  };
};


#endif
