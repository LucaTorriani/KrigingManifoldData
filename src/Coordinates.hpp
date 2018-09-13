#ifndef _COORDINATES_HPP_
#define _COORDINATES_HPP_

#include <iostream>
#include <string>

#include "Helpers.hpp"

class Coordinates {
  const MatrixXd& _coords;
public:
  Coordinates(const MatrixXd& coords): _coords(coords){};
  unsigned int get_N_station() const;
  const MatrixXd& get_coords() const;
  unsigned int get_n_coords() const;
};

#endif
