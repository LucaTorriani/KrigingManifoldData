#ifndef _HELPERS_FACTORY_HPP_
#define _HELPERS_FACTORY_HPP_

#include "FittedVariogram.hpp"
#include "DesignMatrix.hpp"
#include "Factory.hpp"
#include "Proxy.hpp"

namespace vario_factory {
  typedef generic_factory::Factory<variogram_evaluation::FittedVariogram, std::string> VariogramFactory;  // Use standard Builder

  template<typename ConcreteProduct>
  using VariogramProxy = generic_factory::Proxy<VariogramFactory,ConcreteProduct>;
};

namespace design_factory {
  typedef generic_factory::Factory<design_matrix::DesignMatrix, std::string> DesignFactory;  // Use standard Builder

  template<typename ConcreteProduct>
  using DesignProxy = generic_factory::Proxy<DesignFactory,ConcreteProduct>;
};

#endif
