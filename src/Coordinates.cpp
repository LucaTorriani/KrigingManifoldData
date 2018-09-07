#include "Coordinates.hpp"


Coordinates::Coordinates(std::vector<Point> &coords, const Distance& distance):
   _coords(coords){
    _distance_matrix = distance.compute_distance(_coords);
}

SpMat Coordinates::get_distance_matrix () const {
  return _distanceMatrix;
}

std::vector<Point> Coordinates::get_coords () const {
  return _coords;
}

unsigned int Coordinates::get_n_coords() const {
  return _coords[0].get_dimension();
}
