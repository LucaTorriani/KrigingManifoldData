#ifndef _EMPIRICAL_VARIOGRAM_HPP
#define _EMPIRICAL_VARIOGRAM_HPP

#include "Helpers.hpp"
#include "DistanceTplane.hpp"
#include "Distance.hpp"
#include "Coordinates.hpp"

namespace variogram_evaluation{

class EmpiricalVariogram {
  const unsigned int _n_h;
  const unsigned int _N;
  const std::shared_ptr<const SpMat> _distanceMatrix;

  std::vector<double> _emp_vario_values;
  std::vector<double> _hvec;
  std::vector<unsigned int> _N_hvec;
  unsigned int _card_h;
  Vec _d; // Vettore h+-deltah
  double _hmax;
  Vec _weights;

  void compute_hmax(const Coordinates&, const distances::Distance&);

public:
  // EmpiricalVariogram()= default;
  EmpiricalVariogram (const std::shared_ptr<const SpMat>, unsigned int, const Coordinates&, const distances::Distance&);

  void update_emp_vario(const std::vector<MatrixXd>&, const distances_tplane::DistanceTplane &);
  // void set_weight(const Vec&);
  std::vector<double> get_emp_vario_values () const;
  std::vector<unsigned int> get_N_hvec() const;
  std::vector<double> get_hvec() const;
  unsigned int get_card_h() const;
  unsigned int get_N() const;
  double get_hmax() const;
};

}


#endif
