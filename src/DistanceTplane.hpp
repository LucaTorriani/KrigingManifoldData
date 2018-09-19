#ifndef _DISTANCE_TPLANE_HPP_
#define _DISTANCE_TPLANE_HPP_

#include "Helpers.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>
#include <memory>

namespace distances_tplane{

class Frobenius{
public:
  double norm(const MatrixXd&) const; // Attenzione al nome
  double operator()(const MatrixXd &, const MatrixXd &) const;
};

class FrobeniusScaled{
  const std::shared_ptr<const MatrixXd> _Sigma;
public:
  FrobeniusScaled(const std::shared_ptr<const MatrixXd> Sigma):_Sigma(Sigma){};
  double norm (const MatrixXd &) const; // Attenzione al nome
  double operator()(const MatrixXd &, const MatrixXd&) const;
};


class DistanceTplane{
  std::map<std::string, std::function<double(const MatrixXd&, const MatrixXd&)>> _distances;
  const std::string _distanceTplane;
public:
  DistanceTplane(const std::string &, const std::shared_ptr<const MatrixXd>);
  double compute_distance(const MatrixXd&, const MatrixXd&) const;
};



}



#endif
