#ifndef _COORDINATES_HPP_
#define _COORDINATES_HPP_

#include <iostream>
#include <string>
#include <memory>


#include "Helpers.hpp"

class Coordinates {
  const std::shared_ptr<const MatrixXd> _coords;
public:
  Coordinates(const std::shared_ptr<const MatrixXd> coords): _coords(coords){};
  // Coordinates(const Coordinates&) = delete;
  // Coordinates& operator=(const Coordinates&) = delete;

  unsigned int get_N_station() const;
  const std::shared_ptr<const MatrixXd> get_coords() const;
  unsigned int get_n_coords() const;
};

#endif
