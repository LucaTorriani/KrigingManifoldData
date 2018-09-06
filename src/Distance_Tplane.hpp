#ifndef _DISTANCE_TPLANE_HPP_
#define _DISTANCE_TPLANE_HPP_

#include "Helpers.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>

namespace distances_tplane{

class Frobenius{
public:
  static double norm(const SpMat&); // Attenzione al nome
  static double tplane_dist(const SpMat &, const SpMat &);
  double operator()(const SpMat &, const SpMat &, const SpMat & );
};

class FrobeniusScaled{
public:
  static double norm (const SpMat &, const SpMat &) ; // Attenzione al nome
  static double tplane_dist(const SpMat &, const SpMat &, const SpMat & );
  double operator()(const SpMat &, const SpMat&, const SpMat & ) ;
};


class DistanceTplane{
  std::map<std::string,std::function<double(const SpMat&, const SpMat&, const SpMat &)>> distances;
public:
  DistanceTplane();
  double compute_distance(const std::string &, const SpMat&, const SpMat&, const SpMat&);

};



}



#endif
