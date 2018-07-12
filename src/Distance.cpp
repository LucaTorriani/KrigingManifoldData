#include "Distance.hpp"

static double Distance::eucl_dist(const Point& P1, const Point& P2){
  return ((P1-P2).l2norm());
}


static double Distance::geo_dist(const Point&, const Point&){


}



static SpMat Distance::create_distance_matrix(const &std::vector<Point> coords, checkDistance){


}
