#ifndef _HELPERS_FACTORY_HPP_
#define _HELPERS_FACTORY_HPP_

#include "FittedVariogram.hpp"
#include "Factory.hpp"
#include "Proxy.hpp"

namespace vario_factory {
  typedef generic_factory::Factory<variogram_evaluation::FittedVariogram, std::string> VariogramFactory;  // Use standard Builder

  template<typename ConcreteProduct>
  using VariogramProxy = generic_factory::Proxy<VariogramFactory,ConcreteProduct>;
};

#endif
