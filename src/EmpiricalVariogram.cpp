#include "EmpiricalVariogram.hpp"
#include <cmath>

using namespace variogram_evaluation;

// EmpiricalVariogram
EmpiricalVariogram::EmpiricalVariogram (const Coordinates& coords, const distances::Distance& distance, unsigned int n_h, const distances_tplane::DistanceTplane & distanceTplane, const Vec & weights):
  _distanceTplane(distanceTplane), _distanceMatrix(coords.get_distance_matrix()), _weights(weights), _n_h(n_h) {
    compute_hmax(coords, distance);
    _d.resize(n_h +1);
    _d.setLinSpaced(n_h+1, 0, _hmax);
    _N = (coords.get_coords()).size();
}

EmpiricalVariogram::EmpiricalVariogram (const Coordinates& coords, const distances::Distance& distance, unsigned int n_h, const distances_tplane::DistanceTplane & distanceTplane):
  _distanceTplane(distanceTplane), _distanceMatrix(coords.get_distance_matrix()), _n_h(n_h) {
    compute_hmax(coords, distance);
    _d.resize(n_h +1);
    _d.setLinSpaced(n_h+1, 0, _hmax);
    _N = (coords.get_coords()).size();

    _weights.resize(_N);
    _weights.setOnes(_N);
}

unsigned int EmpiricalVariogram::get_N() const {
  return _N;
}

void EmpiricalVariogram::update_emp_vario(std::vector<SpMat>& res) {
  _emp_vario_values.clear();
  _hvec.clear();
  _N_hvec.clear();
  std::vector<double> w_ij;
  std::vector<double> tplanedist2_ij;
  unsigned int card_esitmate = (_N-1)*_N/_n_h;  // Numero diversi valori di distanze tra _N punti è (_N-1)*_N/2
  w_ij.reserve(card_esitmate);
  tplanedist2_ij.reserve(card_esitmate);
  _card_h = 0;

  for (size_t l=1; l<(_n_h+1); l++) {
    w_ij.clear();
    tplanedist2_ij.clear();
    for (size_t j =0; j<(_N-1); j++) {
      for (size_t i=(j+1); i<_N; i++) {
        if (_distanceMatrix.coeff(i,j) >= _d(l-1) && _distanceMatrix.coeff(i,j) <= _d(l)) {
          tplanedist2_ij.push_back( (_distanceTplane.compute_distance(res[i], res[j]))*(_distanceTplane.compute_distance(res[i], res[j])) );
          w_ij.push_back(_weights(i)*_weights(j));
        }
      }
    }
    unsigned int actual_size = tplanedist2_ij.size();
    if (actual_size>0) {
      _card_h += actual_size;
      _N_hvec.push_back(actual_size);
      double num_sum = 0;
      double denom_sum = 0;
      for (size_t i=0; i <actual_size; i++) {
        num_sum += tplanedist2_ij[i]*w_ij[i];
        denom_sum += w_ij[i];
      }
      _emp_vario_values.push_back(num_sum/(2*denom_sum));
      _hvec.push_back((_d(l)+_d(l-1))/2);
    }
  }
}

void EmpiricalVariogram::compute_hmax(const Coordinates& coords, const distances::Distance& distance) {
  unsigned int n_coords = coords.get_n_coords();
  Vec coords_min_point(n_coords);
  Vec coords_max_point(n_coords);
  std::vector<Point> vec_coords (coords.get_coords());
  coords_min_point = (vec_coords[0]).get_coords();
  coords_max_point = coords_min_point;

  for (size_t i= 1; i<vec_coords.size(); i++) {
    Vec point_coords(n_coords);
    point_coords = (vec_coords[i]).get_coords();
    for (size_t j=0; j< n_coords; j++) {

      if (point_coords(j) < coords_min_point(j))
        coords_min_point(j) = point_coords(j);

      if (point_coords(j) > coords_max_point(j))
        coords_max_point(j) = point_coords(j);
    }
  }

  Point min_point(coords_min_point);
  Point max_point(coords_max_point);

  _hmax = (1/3*distance.compute_distance(min_point, max_point));
}

unsigned int EmpiricalVariogram::get_card_h() const {
  return _card_h;
}

std::vector<double> EmpiricalVariogram::get_emp_vario_values () const {
  return _emp_vario_values;
};

std::vector<unsigned int> EmpiricalVariogram::get_N_hvec() const {
  return _N_hvec;
};

std::vector<double> EmpiricalVariogram::get_hvec() const {
  return _hvec;
};
