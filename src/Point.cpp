#include "Point.hpp"

unsigned Point::get_dimension() const {
  return (_coords.size()); // Esiste??
}
const Vec & Point::get_coords() const {
  return (_coords);
}
const double & Point::operator()(std::size_t i) const {
  return (_coords(i));
}

double l2norm() const {
  return(_coords.norm());
}

Point operator-(const Point & P1, const Point & P2){
  if(!check_dim(P1, P2))
  {
    std::cerr << "Error operator-: points don't have the same dimension" << '\n';
  }
  Point P(P1._coords - P2._coords);
  return (P);

}

bool check_dim(const Point & P1, const Point & P2){
  return P1.get_dimensions() == P2.get_dimensions();
}
