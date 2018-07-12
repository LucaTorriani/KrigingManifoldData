#ifndef _DISTANCE_HPP_
#define _DISTANCE_HPP_

#include "Point.hpp"
#include "Helpers.hpp"

struct Distance{
  static double eucl_dist(const Point&, const Point&);
  static double geo_dist(const Point&, const Point&);
  static SpMat create_distance_matrix(const std::vector<Point>&, checkDistance);
};

#endif
