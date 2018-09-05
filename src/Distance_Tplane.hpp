#ifndef _DISTANCE_MANIFOLD_HPP_
#define _DISTANCE_MANIFOLD_HPP_

#include "Point.hpp"
#include "Helpers.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>

namespace distances_tplane{

class Frobenius{
public:
  static double norm(const SpMat&) const; // Attenzione al nome
  static double operator()(const SpMat &, const SpMat & ) const;
};

class FrobeniusScaled{
public:
  static double norm (const SpMat &, const SpMat &) const; // Attenzione al nome
  static double operator()(const SpMat &, const SpMat& ) const;
};


class DistanceTplane{
  std::map<std::string,std::function<double(std::vector<double>, std::vector<double>)>> distances;
public:
  DistanceTplane();
  double compute_distance(const SpMat&, const SpMat&, const std::string &);

};



}



#endif
