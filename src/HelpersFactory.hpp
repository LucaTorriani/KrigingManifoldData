#ifndef _HELPERS_FACTORY_HPP_
#define _HELPERS_FACTORY_HPP_

#include "FittedVariogram.hpp"
#include "DesignMatrix.hpp"
#include "Distance.hpp"
#include "MapFunctions.hpp"

#include "Factory.hpp"
#include "Proxy.hpp"

namespace vario_factory {
  typedef generic_factory::Factory<variogram_evaluation::FittedVariogram, std::string> VariogramFactory;  // Use standard Builder

  template<typename ConcreteProduct>
  using VariogramProxy = generic_factory::Proxy<VariogramFactory,ConcreteProduct>;
}

namespace design_factory {
  typedef generic_factory::Factory<design_matrix::DesignMatrix, std::string> DesignFactory;  // Use standard Builder

  template<typename ConcreteProduct>
  using DesignProxy = generic_factory::Proxy<DesignFactory,ConcreteProduct>;
}

namespace distance_factory {
  typedef generic_factory::Factory<distances::Distance, std::string> DistanceFactory;  // Use standard Builder

  template<typename ConcreteProduct>
  using DistanceProxy = generic_factory::Proxy<DistanceFactory,ConcreteProduct>;
}

namespace map_factory {
  typedef generic_factory::Factory<map_functions::logarithmicMap, std::string> LogMapFactory;  // Use standard Builder
  typedef generic_factory::Factory<map_functions::exponentialMap, std::string> ExpMapFactory;  // Use standard Builder

  template<typename ConcreteProduct>
  using LogMapProxy = generic_factory::Proxy<LogMapFactory,ConcreteProduct>;

  template<typename ConcreteProduct>
  using ExpMapProxy = generic_factory::Proxy<ExpMapFactory,ConcreteProduct>;
}

#endif
