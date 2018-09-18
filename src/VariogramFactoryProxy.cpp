#include<vector>
#include<utility>
#include "HelpersFactory.hpp"
#include "FittedVariogram.hpp"

namespace {
  using vario_factory::VariogramProxy;
  VariogramProxy<variogram_evaluation::GaussVariogram> GAU("Gaussian");
  VariogramProxy<variogram_evaluation::ExpVariogram> EXP("Exponential");
  VariogramProxy<variogram_evaluation::SphVariogram> SPH("Spherical");
}
