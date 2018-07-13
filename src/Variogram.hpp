#ifndef _VARIOGRAM_HPP
#define _VARIOGRAM_HPP

#include <vector>

namespace VariogramEvaluation{

class EmpiricalVariogram {
  std::vector<double> _vario_values;
  std::vector<double> _hvec;
  double compute_hmax() const;

public:
  EmpiricalVariogram(/*residui, n_h*, il resto sono static members*/);
}




}
