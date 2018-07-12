#ifndef _POINT_HPP_
#define _POINT_HPP_

#include <iostream>
#include "Helpers.hpp"

class Point {
  Vec _coords;

public:
  friend Point operator-(const Point &, const Point &);
  friend Point check_dim(const Point &, const Point &);

  Point(const Vec & coord): _coords(coords){}; // Controllare come si passano elementi eigen
  unsigned get_dimension() const;
  const Vec & get_coords() const;
  const double & operator[](std::size_t i);
  double l2norm() const;
};

friend Point operator-(const Point &, const Point &);
friend Point check_dim(const Point &, const Point &);

#endif
