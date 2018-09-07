#ifndef _DISTANCE_MANIFOLD_HPP_
#define _DISTANCE_MANIFOLD_HPP_

#include "Helpers.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>

using namespace Eigen;
namespace distances_manifold{

class Frobenius{
public:
  static double manifold_distance(const SpMat&, const SpMat& );
  double operator()(const SpMat&, const SpMat& );
  };

class LogEuclidean{
public:
  static double manifold_distance(const SpMat&, const SpMat& );
  double operator()(const SpMat&, const SpMat& );
};

class SqRoot{
public:
  static double manifold_distance(const SpMat&, const SpMat& );
  double operator()(const SpMat&, const SpMat& );
};


class DistanceManifold{
  const std::string _distanceManifold;
  const SpMat _Sigma;
  std::map<std::string,std::function<double(const SpMat&, const SpMat&)>> distances;
public:
  DistanceManifold(const std::string&, const SpMat&);
  double compute_distance(const SpMat&, const SpMat&);
  SpMat get_Sigma() const;

};



}



#endif
