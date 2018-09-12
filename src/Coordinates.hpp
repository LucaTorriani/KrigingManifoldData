#ifndef _COORDINATES_HPP_
#define _COORDINATES_HPP_

#include <iostream>
#include <string>

#include "Point.hpp"
#include "Helpers.hpp"
#include "Distance.hpp"

class Coordinates {
  std::vector<Point> _coords;
public:
  Coordinates(std::vector<Point>&, const distances::Distance &);
  unsigned int get_N_station() const;
  std::vector<Point> get_coords() const;
  unsigned int get_n_coords() const;
};

#endif

SpMat Sigma
