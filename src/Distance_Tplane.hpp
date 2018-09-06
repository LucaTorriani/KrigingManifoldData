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
  double norm(const SpMat&); // Attenzione al nome
  double operator()(const SpMat &, const SpMat &);
};

class FrobeniusScaled{
  const SpMat _Sigma;
public:
  FrobeniusScaled(const SpMat &Sigma):_Sigma(Sigma){};
  double norm (const SpMat &) ; // Attenzione al nome
  double operator()(const SpMat &, const SpMat&) ;
};


class DistanceTplane{
  std::map<std::string,std::function<double(const SpMat&, const SpMat&)>> distances;
  const std::string _distanceTplane;
public:
  DistanceTplane(const std::string &, const SpMat &);
  double compute_distance(const SpMat&, const SpMat&);
};



}



#endif
