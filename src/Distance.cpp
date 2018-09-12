#include "Distance.hpp"
#define _USE_MATH_DEFINES
#include <cmath>


using namespace distances;

double EuclDist::operator()(const Point& P1, const Point& P2) const{
  return ((P1-P2).l2norm());
}

// OLD FORMULA: great circle distance (commented 3/9/18)
// double Geodist::operator()(const Point& P1, const Point& P2){
//   // Great-circle distance
//   double coeff = M_PI_2/90;
//   double lat1 =  P1(1);
//   double long1 =  P1(2);
//   double lat2 =  P2(1);
//   double long2 =  P2(2);
//   double cos_angle = sin(lat1*coeff)*sin(lat2*coeff)+cos(lat1*coeff)*cos(lat2*coeff)*cos(cmath::abs(long2-long1)*coeff);
//
//   if (abs(cos_angle-1) < 2*eps_dbl)
//     return (0);
//   else if (abs(cos_angle+1) < 2*eps_dbl)
//     return (Earth_R * M_PI);
//   else
//     return (Earth_R*acos(cos_angle));
//
// }

// Haversine formula (To be tested)
double GeoDist::operator()(const Point& P1, const Point& P2) const{
  double coeff = M_PI_2/90;
  unsigned int n(P1.get_dimension());
  Vec P1_coords(n);
  Vec P2_coords(n);
  P1_coords = P1.get_coords();
  P2_coords = P2.get_coords();
  double lat1 =  P1_coords(1);
  double long1 =  P1_coords(2);
  double lat2 =  P2_coords(1);
  double long2 =  P2_coords(2);
  double sin_1(sin( (lat2-lat1)/2*coeff ));
  double sin_2((sin( (long2-long1)/2*coeff )));
  double sqrth(sqrt( sin_1*sin_1 + cos(lat1*coeff)*cos(lat2*coeff)* sin_2*sin_2 ));

  if (sqrth > 1)
    sqrth = 1;

  return (Earth_R*asin(sqrth));

}

Distance::Distance(const std::string& distance_type):_distance_type(distance_type){
  _dist.insert(std::pair<std::string, std::function<double(const Point&, const Point&)>>("Euclidean", EuclDist()));
  _dist.insert(std::pair<std::string, std::function<double(const Point&, const Point&)>>("Geodist", GeoDist()));

}
double Distance::compute_distance(const Point& P1, const Point& P2) const{
  return _dist.at(_distance_type)(P1, P2);
}

SpMat Distance::create_distance_matrix(const Coordinates & coordinates) const{
  size_t num_points = coordinates.get_N_station();
  std::vector<Point> coords(num_points);
  coords = coordinates.get_coords();
  
  std::vector<TripType> tripletList;
  tripletList.reserve(num_points*(num_points-1)/2);
  for (size_t i=0; i<(num_points-1); i++ ) {
    for (size_t j=(i+1); j<num_points; j++ ) {
      tripletList.push_back(TripType(i,j,compute_distance(coords[i], coords[j])));
    }
  }

  SpMat distance_matrix(num_points, num_points);
  distance_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

  return (distance_matrix);
}
