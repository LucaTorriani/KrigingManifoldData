#ifndef _MAPFUNCTIONS_HPP_
#define _MAPFUNCTIONS_HPP_

#include "Helpers.hpp"
#include "Distance_Manifold.hpp"

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
    MatrixXd _Sigma;
  public:
    logMapLogEucl(const MatrixXd&);
    MatrixXd operator()(const MatrixXd&) const;
  };

  class logMapSqRoot{
    MatrixXd _Sigma;
  public:
    logMapSqRoot(const MatrixXd&);
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
    MatrixXd _Sigma;
  public:
    expMapLogEucl(const MatrixXd&);
    MatrixXd operator()(const MatrixXd&) const;
  };

  class expMapSqRoot{
    MatrixXd _Sigma;
  public:
    expMapSqRoot(const MatrixXd&);
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
