#ifndef _HELPERS_FACTORY_HPP_
#define _HELPERS_FACTORY_HPP_

#include "FittedVariogram.hpp"
#include "DesignMatrix.hpp"
#include "Distance.hpp"
#include "MapFunctions.hpp"
#include "DistanceManifold.hpp"
#include "DistanceTplane.hpp"

#include "Factory.hpp"
#include "Proxy.hpp"

/*! \file
  @brief Typedefs for the factory
*/

namespace vario_factory {
  /*! \var typedef generic_factory::Factory<variogram_evaluation::FittedVariogram, std::string> VariogramFactory
    @brief Factory for the fitted variogram
  */
  typedef generic_factory::Factory<variogram_evaluation::FittedVariogram, std::string> VariogramFactory;  // Use standard Builder

  //! Proxy for the fitted variogram
  template<typename ConcreteProduct>
  using VariogramProxy = generic_factory::Proxy<VariogramFactory,ConcreteProduct>;
}

namespace design_factory {
  /*! \var generic_factory::Factory<design_matrix::DesignMatrix, std::string> DesignFactory
    @brief Factory for the design matrix
  */
  typedef generic_factory::Factory<design_matrix::DesignMatrix, std::string> DesignFactory;  // Use standard Builder

  //! Proxy for the design matrix
  template<typename ConcreteProduct>
  using DesignProxy = generic_factory::Proxy<DesignFactory,ConcreteProduct>;
}

namespace distance_factory {
  /*! \var generic_factory::Factory<distances::Distance, std::string> DistanceFactory
    @brief Factory for the distance
  */
  typedef generic_factory::Factory<distances::Distance, std::string> DistanceFactory;  // Use standard Builder

  //! Proxy for the distance
  template<typename ConcreteProduct>
  using DistanceProxy = generic_factory::Proxy<DistanceFactory,ConcreteProduct>;
}

namespace map_factory {
  /*! \var generic_factory::Factory<map_functions::logarithmicMap, std::string> LogMapFactory
    @brief Factory for the logarithmic map
  */
  typedef generic_factory::Factory<map_functions::logarithmicMap, std::string> LogMapFactory;  // Use standard Builder
  /*! \var generic_factory::Factory<map_functions::exponentialMap, std::string> ExpMapFactory
    @brief Factory for the exponential map
  */
  typedef generic_factory::Factory<map_functions::exponentialMap, std::string> ExpMapFactory;  // Use standard Builder

  //! Proxy for the logarithmic map
  template<typename ConcreteProduct>
  using LogMapProxy = generic_factory::Proxy<LogMapFactory,ConcreteProduct>;

  //! Proxy for the exponential map
  template<typename ConcreteProduct>
  using ExpMapProxy = generic_factory::Proxy<ExpMapFactory,ConcreteProduct>;
}

namespace manifold_factory{
  /*! \var generic_factory::Factory<distances_manifold::DistanceManifold, std::string> ManifoldFactory
    @brief Factory for the distance on the manifold
  */
  typedef generic_factory::Factory<distances_manifold::DistanceManifold, std::string> ManifoldFactory;  // Use standard Builder

  //! Proxy for the distance on the manifold
  template<typename ConcreteProduct>
  using ManifoldProxy = generic_factory::Proxy<ManifoldFactory,ConcreteProduct>;
}

namespace tplane_factory{
  /*! \var generic_factory::Factory<distances_tplane::DistanceTplane, std::string> TplaneFactory
    @brief Factory for the distance on the tangent space
  */
  typedef generic_factory::Factory<distances_tplane::DistanceTplane, std::string> TplaneFactory;  // Use standard Builder

  //! Proxy for the distance on the tangent space
  template<typename ConcreteProduct>
  using TplaneProxy = generic_factory::Proxy<TplaneFactory,ConcreteProduct>;
}

#endif // _HELPERS_FACTORY_HPP_
