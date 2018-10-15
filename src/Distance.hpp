#ifndef _DISTANCE_HPP_
#define _DISTANCE_HPP_

#include "Helpers.hpp"
#include "Coordinates.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>

namespace distances{

  class Distance{

  public:
    virtual double compute_distance(const Vec&, const Vec&) const = 0;
    std::shared_ptr<const MatrixXd> create_distance_matrix(const Coordinates &, unsigned int) const;
    std::vector<double> create_distance_vector(const Coordinates &, const Vec &) const;

  };

  class EuclDist : public Distance{
  public:
    double compute_distance(const Vec&, const Vec&) const override;
  };

  class GeoDist : public Distance{
    static double constexpr Earth_R = 6371.0;
    static double constexpr eps_dbl = std::numeric_limits<double>::epsilon();

  public:
    double compute_distance(const Vec&, const Vec&) const override;
};


}


#endif
