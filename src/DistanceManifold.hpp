#ifndef _DISTANCE_MANIFOLD_HPP_
#define _DISTANCE_MANIFOLD_HPP_

#include "Helpers.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>
#include <memory>


using namespace Eigen;
namespace distances_manifold{

class Frobenius{
public:
  double operator()(const MatrixXd&, const MatrixXd& ) const;
  };

class LogEuclidean{
public:
  double operator()(const MatrixXd&, const MatrixXd& ) const;
};

class SqRoot{
public:
  double operator()(const MatrixXd&, const MatrixXd& ) const;
};


class DistanceManifold{
  const std::string _distanceManifold;
  const std::shared_ptr<const MatrixXd> _Sigma;
  std::map<std::string,std::function<double(const MatrixXd&, const MatrixXd&)>> distances;
public:
  DistanceManifold(const std::string&, const std::shared_ptr<const MatrixXd>);
  double compute_distance(const MatrixXd&, const MatrixXd&) const; // Se vogliamo const togliamo static
  const std::shared_ptr<const MatrixXd> get_Sigma() const;
  const std::string& get_distanceType() const;
};



}



#endif
