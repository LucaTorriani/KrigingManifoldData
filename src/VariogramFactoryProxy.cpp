#include<vector>
#include<utility>
#include "HelpersFactory.hpp"
#include "FittedVariogram.hpp"

namespace {
  using vario_factory::VariogramProxy;
  VariogramProxy<variogram_evaluation::GaussVariogram> GAU("Gaussian");
  //VariogramProxy<ExpVariogram> EXP("Exponential");
  //VariogramProxy<SphVariogram> SPH("Spherical");
}
