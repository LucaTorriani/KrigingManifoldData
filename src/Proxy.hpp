#ifndef GENERICPROXY_HPP_
#define GENERICPROXY_HPP_
#include <string>
#include <memory>
#include <iostream>
#include <type_traits>
namespace generic_factory {


  template <typename Factory, typename ConcreteProduct>
  class Proxy {
  public:

    typedef typename  Factory::AbstractProduct_type AbstractProduct_type;
    typedef typename  Factory::Identifier_type Identifier_type;
    typedef typename  Factory::Builder_type Builder_type;
    typedef Factory Factory_type;

    Proxy(Identifier_type const &);
    static std::unique_ptr<AbstractProduct_type> Build(){
      return std::unique_ptr<AbstractProduct_type>(new ConcreteProduct());
    }

  private:
    Proxy(Proxy const &)=delete;
    Proxy & operator=(Proxy const &)=delete;
  };


  template<typename F, typename C>
  Proxy<F,C>::Proxy(Identifier_type const & name) {
    // get the factory. First time creates it.
    Factory_type & factory(Factory_type::Instance());
    // Insert the builder. The & is not needed.
    factory.add(name,&Proxy<F,C>::Build);
    std::cout<<"Added "<< name << " to factory"<<std::endl;
  }
}

#endif /* RULESHANDLER_HPP_ */
