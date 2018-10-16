#ifndef _EMPIRICAL_VARIOGRAM_HPP
#define _EMPIRICAL_VARIOGRAM_HPP

#include "Helpers.hpp"
#include "DistanceTplane.hpp"
#include "Distance.hpp"
#include "Coordinates.hpp"

namespace variogram_evaluation{

class EmpiricalVariogram {
  // Costanti
  const unsigned int _n_h;
  const unsigned int _N;
  const std::shared_ptr<const MatrixXd> _distanceMatrix;

  // Variano
  std::vector<double> _emp_vario_values;
  std::vector<double> _hvec;
  std::vector<unsigned int> _N_hvec;
  unsigned int _card_h;

  Vec _d; // Vettore h+-deltah
  double _hmax;
  // const distances_tplane::DistanceTplane & _distanceTplane;
  Vec _weights;

  void compute_hmax(const Coordinates&, const distances::Distance&);

public:
  // EmpiricalVariogram()= default;
  EmpiricalVariogram (const std::shared_ptr<const MatrixXd>, unsigned int,  // EQUAL WEIGHTS
                      const Coordinates&, const distances::Distance&);
  EmpiricalVariogram (const std::shared_ptr<const MatrixXd>, unsigned int, unsigned int,  // KERNEL
                      const Vec&, double);
  void update_emp_vario(const std::vector<MatrixXd>&, const distances_tplane::DistanceTplane &);
  std::vector<double> get_emp_vario_values () const;
  std::vector<unsigned int> get_N_hvec() const;
  std::vector<double> get_hvec() const;
  unsigned int get_card_h() const;
  unsigned int get_N() const;
};

}


#endif
