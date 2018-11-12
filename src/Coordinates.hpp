#ifndef _COORDINATES_HPP_
#define _COORDINATES_HPP_

#include <iostream>
#include <string>
#include <memory>


#include "Helpers.hpp"

/*! \file
  @brief Coordinates class
*/

/*!
  @brief Class to store the coordinates of the data points
*/

class Coordinates {
  /*! Matrix of the coordinates */
  const std::shared_ptr<const MatrixXd> _coords;
public:
  /*!
    @brief Constructor
    @param coords Matrix containg the coordinates of the data points
  */
  Coordinates(const std::shared_ptr<const MatrixXd> coords): _coords(coords){};

  /*!
    @brief Return the number of data points in the analysis
  */
  unsigned int get_N_station() const;

  /*!
    @brief Return the matrix of the coordinates
  */
  const std::shared_ptr<const MatrixXd> get_coords() const;

  /*!
    @brief Return the number of coordinates for each data point
  */
  unsigned int get_n_coords() const;
};

#endif  // _COORDINATES_HPP_
