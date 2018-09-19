#ifndef _DISTANCE_HPP_
#define _DISTANCE_HPP_

#include "Helpers.hpp"
#include "Coordinates.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>

namespace distances{

class EuclDist{
public:
      double operator()(const Vec&, const Vec& ) const;
  };

class GeoDist{
  static double constexpr Earth_R = 6371.0;
  static double constexpr eps_dbl = std::numeric_limits<double>::epsilon();

public:
  double operator()(const Vec&, const Vec&) const;
};

class Distance{
  std::string _distance_type;
  std::map<std::string,std::function<double(const Vec&, const Vec&)>> _dist;
public:
  Distance(const std::string &);
  double compute_distance(const Vec&, const Vec&) const;
  SpMat create_distance_matrix(const Coordinates &, unsigned int) const;
  std::vector<double> create_distance_vector(const Coordinates &, const Vec &) const;


};



}


// struct Distance{
//   static double eucl_dist(const Point&, const Point&);
//   static double geo_dist(const Point&, const Point&);
//   static SpMat create_distance_matrix(const std::vector<Point>&, checkDistance);
//   static const double Earth_R(6371.0);
//   static const double eps_dbl(std::numeric_limits<double>::epsilon());
// };
#endif
