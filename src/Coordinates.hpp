#ifndef _COORDINATES_HPP_
#define _COORDINATES_HPP_

#include <iostream>
#include <string>

#include "Point.hpp"
#include "Helpers.hpp"
#include "Distance.hpp"

class Coordinates {
  std::vector<Point> _coords;
  SpMat _distanceMatrix;
public:
  Coordinates(std::vector<Point>&, const Distance &);
  std::vector<Point> get_coords const;
  SpMat get_distance_matrix const;

};

#endif
