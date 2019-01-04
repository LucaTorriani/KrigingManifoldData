#ifndef _DISTANCE_TPLANE_HPP_
#define _DISTANCE_TPLANE_HPP_

#include "Helpers.hpp"
#include <vector>
#include <utility>
#include <map>
#include <functional>
#include <memory>

namespace distances_tplane{

  class DistanceTplane{
  public:
    double compute_distance(const MatrixXd&, const MatrixXd&) const;
    virtual double norm(const MatrixXd&) const = 0;
    virtual void set_members(const MatrixXd&) = 0;
    virtual ~DistanceTplane() = default;
  };

 class Frobenius : public DistanceTplane{
 public:
   ~Frobenius() = default;
   double norm(const MatrixXd &) const override;
   void set_members(const MatrixXd&) override;

 };

//  class FrobeniusScaled : public DistanceTplane{
//    MatrixXd _SigmaInv;
//    unsigned int _p;
//  public:
//    ~FrobeniusScaled() = default;
//    double norm(const MatrixXd &) const override;
//    void set_members(const MatrixXd&) override;
//
// };



}



#endif
