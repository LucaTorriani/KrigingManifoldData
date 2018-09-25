#include<vector>
#include<utility>
#include "HelpersFactory.hpp"
// #include "FittedVariogram.hpp"
// #include "DesignMatrix.hpp"

namespace {
  using vario_factory::VariogramProxy;
  using design_factory::DesignProxy;

  VariogramProxy<variogram_evaluation::GaussVariogram> gau("Gaussian");
  VariogramProxy<variogram_evaluation::ExpVariogram> exp("Exponential");
  VariogramProxy<variogram_evaluation::SphVariogram> sph("Spherical");

  DesignProxy<design_matrix::InterceptDM> intercept("Intercept");
  DesignProxy<design_matrix::Coord1DM> coord1("Coord1");
  DesignProxy<design_matrix::Coord2DM> coord2("Coord2");
  DesignProxy<design_matrix::AdditiveDM> additive("Additive");

}
