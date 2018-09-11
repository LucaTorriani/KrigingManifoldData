#ifndef _EMPIRICAL_VARIOGRAM_HPP
#define _EMPIRICAL_VARIOGRAM_HPP

#include "Helpers.hpp"
#include "Distance_Tplane.hpp"
#include "Distance.hpp"
#include "Coordinates.hpp"

namespace variogram_evaluation{

class EmpiricalVariogram {
  //Variano
  std::vector<double> _emp_vario_values;
  std::vector<double> _hvec;
  std::vector<unsigned int> _N_hvec;
  unsigned int _card_h;
  // Costanti
  unsigned int _N;
  Vec _d; // Vettore h+-deltah
  double _hmax;
  const distances_tplane::DistanceTplane _distanceTplane;
  const SpMat _distanceMatrix;
  Vec _weights;
  const unsigned int _n_h;

  void compute_hmax(const Coordinates&, const distances::Distance&);

public:
  EmpiricalVariogram()= default;
  EmpiricalVariogram (const Coordinates&, const distances::Distance&, unsigned int, const distances_tplane::DistanceTplane &, const Vec &);
  EmpiricalVariogram (const Coordinates&, const distances::Distance&, unsigned int, const distances_tplane::DistanceTplane &);
  std::vector<double> get_emp_vario_values () const;
  std::vector<unsigned int> get_N_hvec() const;
  std::vector<double> get_hvec() const;
  unsigned int get_card_h() const;
  void update_emp_vario(std::vector<SpMat>&);
  unsigned int get_N() const;
};

};


#endif
