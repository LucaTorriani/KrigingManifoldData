#include "Coordinates.hpp"


Coordinates::Coordinates(std::vector<Point> &coords ): _coords(coords){};


std::vector<Point> Coordinates::get_coords () const {
  return _coords;
}

unsigned int Coordinates::get_n_coords() const {
  return _coords[0].get_dimension();
}

unsigned int get_N_station() const{
  return _coords.size();
}
