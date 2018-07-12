#include "Coordinates.hpp"


Coordinates::Coordinates(const std::string distance_type, std::vector<Point> coords):
  _distance_type(distance_type), _coords(coords){
    _distance_matrix = Distance::create_distance_matrix(_coords, _distance_type );
}
