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
  double norm(const MatrixXd&) const; // Attenzione al nome
  double operator()(const MatrixXd &, const MatrixXd &) const;
};

class FrobeniusScaled{
  const MatrixXd& _Sigma;
public:
  FrobeniusScaled(const MatrixXd &Sigma):_Sigma(Sigma){};
  double norm (const MatrixXd &) const; // Attenzione al nome
  double operator()(const MatrixXd &, const MatrixXd&) const;
};


class DistanceTplane{
  std::map<std::string, std::function<double(const MatrixXd&, const MatrixXd&)>> _distances;
  const std::string _distanceTplane;
public:
  DistanceTplane(const std::string &, const MatrixXd &);
  double compute_distance(const MatrixXd&, const MatrixXd&) const;
};



}



#endif
