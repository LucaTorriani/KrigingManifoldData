#include "EmpiricalVariogram.hpp"
#include <cmath>

using namespace variogram_evaluation;

// EmpiricalVariogram
EmpiricalVariogram::EmpiricalVariogram (const Coordinates& coords, const distances::Distance& distance, unsigned int n_h, const distances_tplane::DistanceTplane & distanceTplane, const std::shared_ptr<const SpMat> distanceMatrix):
  _N(coords.get_N_station()), _n(coords.get_n_coords()), _distanceTplane(distanceTplane), _distanceMatrix(distanceMatrix), _n_h(n_h) {
    compute_hmax(coords, distance);
    _d.resize(n_h +1);
    _d.setLinSpaced(n_h+1, 0, _hmax);
    std::cout << "d " << _d << "\n" <<std::endl;
    std::cout << "h " << _hmax << "\n" <<std::endl;

    _weights.resize(_N);
    _weights.setOnes(_N);
}

unsigned int EmpiricalVariogram::get_N() const {
  return _N;
}


void EmpiricalVariogram::set_weight(const Vec& weights){
  _weights = weights;
}

void EmpiricalVariogram::update_emp_vario(const std::vector<MatrixXd>& res) {
  _emp_vario_values.clear();
  _hvec.clear();
  _N_hvec.clear();
  std::vector<double> w_ij;
  std::vector<double> tplanedist2_ij;
  unsigned int card_estimate = ((_N-1)*_N)/_n_h;  // Numero diversi valori di distanze tra _N punti è (_N-1)*_N/2
  w_ij.reserve(card_estimate);
  tplanedist2_ij.reserve(card_estimate);

  for (size_t l=1; l<(_n_h+1); l++) {
    w_ij.clear();
    tplanedist2_ij.clear();
    for (size_t i =0; i<(_N-1); i++) {

      for (size_t j=(i+1); j<_N; j++) {
        if (_distanceMatrix->coeff(i,j) >= _d(l-1) && _distanceMatrix->coeff(i,j) <= _d(l)) {
          double tmp(_distanceTplane.compute_distance(res[i], res[j]));
          tplanedist2_ij.push_back( tmp*tmp );
          w_ij.push_back(_weights(i)*_weights(j));
        }
      }
    }
    unsigned int actual_size = tplanedist2_ij.size();
    if (actual_size>0) {
      _N_hvec.push_back(actual_size);
      double num_sum = 0;
      double denom_sum = 0;
      for (size_t k=0; k <actual_size; k++) {
        num_sum += tplanedist2_ij[k]*w_ij[k];
        denom_sum += w_ij[k];
      }
      _emp_vario_values.push_back(num_sum/(2*denom_sum));
      _hvec.push_back((_d(l)+_d(l-1))/2);
    }
  }
  _card_h = _hvec.size();

}

void EmpiricalVariogram::compute_hmax(const Coordinates& coords, const distances::Distance& distance) {
  unsigned int n_coords = coords.get_n_coords();
  Vec min_point(n_coords);
  Vec max_point(n_coords);
  const std::shared_ptr<const MatrixXd> mat_coords (coords.get_coords());

  min_point = (mat_coords->colwise()).minCoeff();
  max_point = (mat_coords->colwise()).maxCoeff();
  _hmax = (1.0/3)*distance.compute_distance(min_point, max_point);
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
