#ifndef _DISTANCE_MANIFOLD_HPP_
#define _DISTANCE_MANIFOLD_HPP_

#include "Point.hpp"
#include "Helpers.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>

namespace distances_manifold{

class Frobenius{
public:
  double operator()(const SpMat&, const SpMat& );
  };

class LogEuclidean{
public:
  double operator()(const SpMat&, const SpMat& );
};

class SqRoot{
public:
  double operator()(const SpMat&, const SpMat& );
};


class DistanceManifold{
  std::map<std::string,std::function<double(std::vector<double>, std::vector<double>)>> dist;
public:
  DistanceManifold();
  double compute_distance(const SpMat&, const SpMat&, const std::string &);

};



}



#endif
