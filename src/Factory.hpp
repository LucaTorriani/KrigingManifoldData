#ifndef H_FACTORY_HPP_
#define H_FACTORY_HPP_
#include <map>
#include <vector>
#include <memory>
#include <functional>
#include <stdexcept>
#include <type_traits>

namespace generic_factory{

  template <typename AbstractProduct, typename Identifier, typename Builder=std::function<std::unique_ptr<AbstractProduct> ()>>
  class Factory{

  public:
    using AbstractProduct_type = AbstractProduct;
    using  Identifier_type=Identifier;
    using  Builder_type=Builder;

    static Factory & Instance();
    std::unique_ptr<AbstractProduct> create(Identifier const & name) const;
    void add(Identifier const &, Builder_type const &);
    std::vector<Identifier> registered()const;
    void unregister(Identifier const & name){ _storage.erase(name);}
    ~Factory()=default;

  private:
    typedef std::map<Identifier,Builder_type> Container_type;
    Factory()=default;
    Factory(Factory const &)=delete;
    Factory & operator =(Factory const &)=delete;
    Container_type _storage;
  };



  template <typename AbstractProduct, typename Identifier, typename Builder>
  Factory<AbstractProduct,Identifier,Builder> &
  Factory<AbstractProduct,Identifier,Builder>::Instance() {
    static Factory theFactory;
    return theFactory;
  }


  template <typename AbstractProduct, typename Identifier, typename Builder>
  std::unique_ptr<AbstractProduct>
  Factory<AbstractProduct,Identifier,Builder>::create(Identifier const & name) const {

    auto f = _storage.find(name); //C++11
    if (f == _storage.end()) {
	     std::string out="Identifier " + name + " is not stored in the factory";
	      throw std::invalid_argument(out);
    }
    else {
	       return std::unique_ptr<AbstractProduct>(f->second());
    }
  }

  template <typename AbstractProduct, typename Identifier, typename Builder>
  void
  Factory<AbstractProduct,Identifier,Builder>::add(Identifier const & name, Builder_type const & func){

    auto f =  _storage.insert(std::make_pair(name, func));
    if (f.second == false)
      throw std::invalid_argument("Double registration in Factory");
  }


  template <typename AbstractProduct, typename Identifier, typename Builder>
  std::vector<Identifier>
  Factory<AbstractProduct,Identifier,Builder>::registered() const {
    std::vector<Identifier> tmp;
    tmp.reserve(_storage.size());
    for(auto i=_storage.begin(); i!=_storage.end();++i)
      tmp.push_back(i->first);
    return tmp;
  }

}


#endif
