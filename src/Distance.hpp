#ifndef _DISTANCE_HPP_
#define _DISTANCE_HPP_

#include "Helpers.hpp"
#include "Coordinates.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>

/*! \file
  @brief Classes to create the compute the distances between data locations
*/

namespace distances{
  /*!
    @brief	Abstract class for the computation of the distance between data locations
  */
  class Distance{
  public:

    /*!
      @brief Compute the distance between two locations
      @param P1 Vector of coordinates for the first location
      @param P2 Vector of coordinates for the second location
      @return Distance
    */
    virtual double compute_distance(const Vec& P1, const Vec& P2) const = 0;

    /*!
      @brief Compute the distance matrix among a set of locations
      @param coordinates Matrix of the coordinates of the locations
      @param N Number of locations
      @return Distance matrix
    */
    std::shared_ptr<const MatrixXd> create_distance_matrix(const Coordinates & coordinates, unsigned int N) const;

    /*!
      @brief Compute the vector of distances between a point and a set of locations
      @param coordinates Matrix of the coordinates of the locations
      @param new_coord Coordinates of the new location
      @return Distance vector
    */
    std::vector<double> create_distance_vector(const Coordinates & coordinates, const Vec & new_coord) const;
    /*!
      @brief Destructor
    */
    virtual ~Distance() = default;
  };
  /*!
    @brief	Class for the computation of the distance between data locations when \f$\texttt{distance=="Eucldist"}\f$
  */
  class EuclDist : public Distance{
  public:
    double compute_distance(const Vec& P1, const Vec& P2) const override;
    /*!
      @brief Destructor
    */
    ~EuclDist() = default;
  };
  /*!
    @brief	Class for the computation of the distance between data locations when \f$\texttt{distance=="Geodist"}\f$
  */
  class GeoDist : public Distance{
    /*! Earth radius */
    static double constexpr Earth_R = 6371.0;
    /*! Machine epsilon */
    static double constexpr eps_dbl = std::numeric_limits<double>::epsilon();
  public:
    double compute_distance(const Vec& P1, const Vec& P2) const override;
    /*!
      @brief Destructor
    */
    ~GeoDist() = default;
};


}


#endif // _DISTANCE_HPP_
