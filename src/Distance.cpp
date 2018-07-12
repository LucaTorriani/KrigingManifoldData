#include "Distance.hpp"
#define _USE_MATH_DEFINES
#include <cmath>

static double Distance::eucl_dist(const Point& P1, const Point& P2){
  return ((P1-P2).l2norm());
}


static double Distance::geo_dist(const Point& P1, const Point& P2){
  // Great-circle distance
  double coeff = M_PI_2/90;
  double lat1 =  P1(1);
  double long1 =  P1(2);
  double lat2 =  P2(1);
  double long2 =  P2(2);
  double cos_angle = sin(lat1*coeff)*sin(lat2*coeff)+cos(lat1*coeff)*cos(lat2*coeff)*cos(cmath::abs(long2-long1)*coeff);

  if (abs(cos_angle-1) < 2*Distance::eps_dbl)
    return (0);
  else if (abs(cos_angle+1) < 2*Distance::eps_dbl)
    return (Distance::Earth_R * M_PI);
  else
    return (Distance::Earth_R*acos(cos_angle));

}



static SpMat Distance::create_distance_matrix(const &std::vector<Point> coords, checkDistance){
  size_t num_points = coords.size();

  std::vector<TripType> tripletList;
  tripletList.reserve(num_points*(num_points-1)/2);
  switch (checkDistance){
    case EUCLIDEAN: {
      for (size_t i=0; i<(num_points-1); i++ ) {
        for (size_t j=(i+1); j<num_points; j++ ) {
          tripletList.push_back(TripType(i,j,Distance::eucl_dist(coords[i], coords[j])));
        }
      }
    }
    break;
    case GEODIST: {
      for (size_t i=0; i<(num_points-1); i++ ) {
        for (size_t j=(i+1); j<num_points; j++ ) {
          tripletList.push_back(TripType(i,j,Distance::geo_dist(coords[i], coords[j])));
        }
      }
    }
    break;
  };

  SpMat distance_matrix(num_points, num_points);
  distance_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

  return (distance_matrix);
}
