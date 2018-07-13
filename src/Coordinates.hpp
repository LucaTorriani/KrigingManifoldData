#ifndef _COORDINATES_HPP_
#define _COORDINATES_HPP_

#include <iostream>
#include <string>

#include "Point.hpp"
#include "Helpers.hpp"
#include "Distance.hpp"

class Coordinates {
  std::vector<Point> _coords;
  SpMat _distance_matrix;
  static checkDistance _distance_type;
public:
  Coordinates(const checkDistance, std::vector<Point>&);


};

#endif
