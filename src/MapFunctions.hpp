#ifndef _MAPFUNCTIONS_HPP_
#define _MAPFUNCTIONS_HPP_

#include "Helpers.hpp"
#include "Distance_Manifold.hpp"

namespace map_functions {

// LOGARITHMIC MAP
  class logMapFrob{
    const MatrixXd _sqrtSigma;
    const MatrixXd _sqrtSigmaInv;
  public:
    logMapFrob(const MatrixXd&, const MatrixXd&);
    SpMat operator()(const SpMat&);
  };

  class logMapLogEucl{
    SpMat _Sigma;
  public:
    logMapLogEucl(const SpMat&);
    SpMat operator()(const SpMat&);
  };

  class logMapSqRoot{
    SpMat _Sigma;
  public:
    logMapSqRoot(const SpMat&);
    SpMat operator()(const SpMat&);
  };

  class logarithmicMap{
    const std::string _distanceManifold;
    std::map<std::string, std::function<SpMat(const SpMat&)>> maps;
  public:
    logarithmicMap(const distances_manifold::DistanceManifold&);
    SpMat map2tplane(const SpMat&);
  };


  // EXPONENTIAL MAP MAP
  class expMapFrob{
    const MatrixXd _sqrtSigma;
    const MatrixXd _sqrtSigmaInv;
  public:
    expMapFrob(const MatrixXd&, const MatrixXd&);
    SpMat operator()(const SpMat&);
  };

  class expMapLogEucl{
    SpMat _Sigma;
  public:
    expMapLogEucl(const SpMat&);
    SpMat operator()(const SpMat&);
  };

  class expMapSqRoot{
    SpMat _Sigma;
  public:
    expMapSqRoot(const SpMat&);
    SpMat operator()(const SpMat&);
  };

  class exponentialMap{
    const std::string _distanceManifold;
    std::map<std::string, std::function<SpMat(const SpMat&)>> maps;
  public:
    exponentialMap(const distances_manifold::DistanceManifold&);
    SpMat map2manifold(const SpMat&);
  };
};


#endif
