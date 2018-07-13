#include "Distance.hpp"
#define _USE_MATH_DEFINES
#include <cmath>


using namespace distances;

double EuclDist::operator()(const Point& P1, const Point& P2){
  return ((P1-P2).l2norm());
}

double Geodist::operator()(const Point& P1, const Point& P2){
  // Great-circle distance
  double coeff = M_PI_2/90;
  double lat1 =  P1(1);
  double long1 =  P1(2);
  double lat2 =  P2(1);
  double long2 =  P2(2);
  double cos_angle = sin(lat1*coeff)*sin(lat2*coeff)+cos(lat1*coeff)*cos(lat2*coeff)*cos(cmath::abs(long2-long1)*coeff);

  if (abs(cos_angle-1) < 2*eps_dbl)
    return (0);
  else if (abs(cos_angle+1) < 2*eps_dbl)
    return (Earth_R * M_PI);
  else
    return (Earth_R*acos(cos_angle));

}

Distance::Distance(){
  distances.insert(std::pair<std::string, std::function<double(std::vector<double>, std::vector<double>)>>("Euclidean", EuclDist()));
  distances.insert(std::pair<std::string, std::function<double(std::vector<double>, std::vector<double>)>>("Manhattan", Manhattan()));

}
double Distance::compute_distance(const Point& P1, const Point& P2, const std::string & distance_type){
  double result = dist[distance_type](P1, P2);
  return result;
}

SpMat Distance::create_distance_matrix(const std::vector<Point> & coords,const std::string & distance_type){
  size_t num_points = coords.size();

  std::vector<TripType> tripletList;
  tripletList.reserve(num_points*(num_points-1)/2);
  for (size_t i=0; i<(num_points-1); i++ ) {
    for (size_t j=(i+1); j<num_points; j++ ) {
      tripletList.push_back(TripType(i,j,compute_distance(coords[i], coords[j], distance_type)));
    }
  }

  SpMat distance_matrix(num_points, num_points);
  distance_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

  return (distance_matrix);
}
