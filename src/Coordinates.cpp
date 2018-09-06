#include "Coordinates.hpp"


Coordinates::Coordinates(std::vector<Point> &coords, const Distance& distance):
   _coords(coords){
    _distance_matrix = distance.compute_distance(_coords);
}

SpMat Coordinates::get_distance_matrix () {
  return _distanceMatrix;
}

std::vector<Point> Coordinates::get_coords () {
  return _coords;
}
