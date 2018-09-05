#ifndef _DISTANCE_HPP_
#define _DISTANCE_HPP_

#include "Point.hpp"
#include "Helpers.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>

namespace distances{

class EuclDist{
public:
      static double operator()(const Point&, const Point& );
  };

class GeoDist{
  double constexpr Earth_R = 6371.0;
  double constexpr eps_dbl = std::numeric_limits<double>::epsilon;

public:
  static double operator()(const Point&, const Point&);
};

class Distance{
  std::map<std::string,std::function<double(std::vector<double>, std::vector<double>)>> dist;
public:
  Distance();
  double compute_distance(const Point&, const Point&, const std::string &);
  SpMat create_distance_matrix(const std::vector<Point> &, const std::string & );

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
