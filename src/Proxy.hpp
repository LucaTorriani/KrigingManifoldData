#ifndef GENERICPROXY_HPP_
#define GENERICPROXY_HPP_
#include <string>
#include <memory>
#include <iostream>
#include <type_traits>

/*! \file
  @brief Proxy class
*/

namespace generic_factory {

  /*! A simple proxy for registering into a factory.

    It provides the builder as static method
    and the automatic registration mechanism.

    \param Factory The type of the factory.
    \param ConcreteProduct Is the derived (concrete) type to be
    registered in the factory

    @note I have to use the default builder provided by the factory. No check is made to verify it
  */
  template <typename Factory, typename ConcreteProduct>
  class Proxy {
  public:
    /*! \var typedef typename  Factory::AbstractProduct_type AbstractProduct_type
      @brief The container for the rules.
    */
    typedef typename  Factory::AbstractProduct_type AbstractProduct_type;
    /*! \var typedef typename  Factory::Identifier_type Identifier_type
      @brief The identifier.
    */
    typedef typename  Factory::Identifier_type Identifier_type;
    /*! \var typedef typename  Factory::Builder_type Builder_type
      @brief The builder type.
    */
    typedef typename  Factory::Builder_type Builder_type;
    /*! \var typedef Factory Factory_type
      @brief The factory type.
    */
    typedef Factory Factory_type;

    //! The constructor does the registration.
    Proxy(Identifier_type const &);
    //! The builder. Must comply with the signature.
    static std::unique_ptr<AbstractProduct_type> Build(){
      return std::unique_ptr<AbstractProduct_type>(new ConcreteProduct());
    }

  private:
    //! Copy onstructor deleted since it is a Singleton
    Proxy(Proxy const &)=delete;
    //! Assignment operator deleted since it is a Singleton
    Proxy & operator=(Proxy const &)=delete;
  };


  template<typename F, typename C>
  Proxy<F,C>::Proxy(Identifier_type const & name) {
    // get the factory. First time creates it.
    Factory_type & factory(Factory_type::Instance());
    // Insert the builder.
    factory.add(name,&Proxy<F,C>::Build);
    // std::cout<<"Added "<< name << " to factory"<<std::endl;
  }
}

#endif // GENERICPROXY_HPP_
