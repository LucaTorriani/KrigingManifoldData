#include "Coordinates.hpp"


const std::shared_ptr<const MatrixXd> Coordinates::get_coords () const {
  return _coords;
}

unsigned int Coordinates::get_n_coords() const {
  return _coords->cols();
}

unsigned int Coordinates::get_N_station() const{
  return _coords->rows();
}
