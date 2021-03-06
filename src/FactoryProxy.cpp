#include<vector>
#include<utility>
#include "HelpersFactory.hpp"

namespace {
  using vario_factory::VariogramProxy;
  using design_factory::DesignProxy;
  using distance_factory::DistanceProxy;
  using map_factory::LogMapProxy;
  using map_factory::ExpMapProxy;
  using manifold_factory::ManifoldProxy;
  using tplane_factory::TplaneProxy;


  VariogramProxy<variogram_evaluation::GaussVariogram> gau("Gaussian");
  VariogramProxy<variogram_evaluation::ExpVariogram> exp("Exponential");
  VariogramProxy<variogram_evaluation::SphVariogram> sph("Spherical");

  DesignProxy<design_matrix::InterceptDM> intercept("Intercept");
  DesignProxy<design_matrix::Coord1DM> coord1("Coord1");
  DesignProxy<design_matrix::Coord2DM> coord2("Coord2");
  DesignProxy<design_matrix::AdditiveDM> additive("Additive");

  DistanceProxy<distances::EuclDist> eucldist("Eucldist");
  DistanceProxy<distances::GeoDist> geodist("Geodist");

  LogMapProxy<map_functions::logMapFrob> logfrob("Frobenius");
  LogMapProxy<map_functions::logMapLogEucl> loglogeucl("LogEuclidean");
  LogMapProxy<map_functions::logMapSqRoot> logsqroot("SquareRoot");
  LogMapProxy<map_functions::logMapChol> logchol("Correlation");


  ExpMapProxy<map_functions::expMapFrob> expfrob("Frobenius");
  ExpMapProxy<map_functions::expMapLogEucl> explogeucl("LogEuclidean");
  ExpMapProxy<map_functions::expMapSqRoot> expsqroot("SquareRoot");
  ExpMapProxy<map_functions::expMapChol> expchol("Correlation");

  ManifoldProxy<distances_manifold::Frobenius> manfrob("Frobenius");
  ManifoldProxy<distances_manifold::LogEuclidean> manlogeucl("LogEuclidean");
  ManifoldProxy<distances_manifold::SqRoot> mansqroot("SquareRoot");
  ManifoldProxy<distances_manifold::Chol> manchol("Correlation");

  TplaneProxy<distances_tplane::Frobenius> tplanfrob("Frobenius");
  TplaneProxy<distances_tplane::FrobeniusScaled> tplanfrobscal("FrobeniusScaled");
  TplaneProxy<distances_tplane::Chol> tplanchol("Correlation");
}
