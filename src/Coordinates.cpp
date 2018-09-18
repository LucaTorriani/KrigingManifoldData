#include "Coordinates.hpp"


const MatrixXd& Coordinates::get_coords () const {
  std::cout << _coords << std::endl;
  return _coords;
}

unsigned int Coordinates::get_n_coords() const {
  return _coords.cols();
}

unsigned int Coordinates::get_N_station() const{
  return _coords.rows();
}
