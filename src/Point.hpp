#ifndef _POINT_HPP_
#define _POINT_HPP_

#include <iostream>
#include "Helpers.hpp"

class Point {
  Vec _coords;

public:
  friend Point operator-(const Point &, const Point &);
  friend bool check_dim(const Point &, const Point &);

  Point(const Vec & coords): _coords(coords){}; // Controllare come si passano elementi eigen
  unsigned get_dimension() const;
  Vec get_coords() const;
  double operator()(std::size_t i) const;
  double l2norm() const;
};

Point operator-(const Point &, const Point &);
bool check_dim(const Point &, const Point &);

#endif
