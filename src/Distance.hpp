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
      double operator()(const Point&, const Point& );
  };

class GeoDist{
  double constexpr Earth_R = 6371.0;
  double constexpr eps_dbl = std::numeric_limits<double>::epsilon;

public:
  double operator()(const Point&, const Point&);
};

class Distance{
  std::string _distance_type;
  std::map<std::string,std::function<double(const Point&, const Point&)>> dist;
public:
  Distance(const std::string &);
  double compute_distance(const Point&, const Point&);
  SpMat create_distance_matrix(const std::vector<Point> &);

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
