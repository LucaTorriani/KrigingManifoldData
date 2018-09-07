#include "Variogram.hpp"
#include <cmath>

using namespace VariogramEvaluation;

// EmpiricalVariogram
EmpiricalVariogram::EmpiricalVariogram (const Coordinates& coords, const Distance& distance, unsigned int n_h, const DistanceTplane & distanceTplane, const Vec weights &):
  _distanceTplane(distanceTplane), _weights(weight) {
    std::vector<Point> local_coords(coords.get_coords());
    _distanceMatrix = coords.get_distance_matrix();
    compute_hmax(local_coords, distance);
    double delta_h = _hmax/n_h;
    _d.resize(n_h +1);
    _d.setLinSpaced(n_h+1, 0, _hmax);
    _N = local_coords.size();

    _emp_vario_values.resize(n_h);
    _hvec.resize(n_h);
    _N_hvec.resize(n_h);
}

EmpiricalVariogram::EmpiricalVariogram (const Coordinates& coords, const Distance& distance, unsigned int n_h, const DistanceTplane & distanceTplane):
  _distanceTplane(distanceTplane), _weights(weight) {
    std::vector<Point> local_coords(coords.get_coords());
    _distanceMatrix = coords.get_distance_matrix();
    compute_hmax(local_coords, distance);
    double delta_h = _hmax/n_h;
    _d.resize(n_h +1);
    _d.setLinSpaced(n_h+1, 0, _hmax);
    _N = local_coords.size();

    // _emp_vario_values.resize(n_h);
    // _hvec.resize(n_h);
    // _N_hvec.resize(n_h);

    _weights.resize(_N);
    _weights.setOnes(_N);
}

void EmpiricalVariogram::update_emp_vario(std::vector<SpMat> res) {
  _emp_vario_values.clear();
  _hvec.clear();
  _N_hvec.clear();

  for (auto l=1; l<(n_h+1); l++) {
    std::vector<double> w_ij;
    std::vector<double> tplanedist2_ij;
    for (auto j =0; i<(_N-1); j++) {
      for (auto i=(j+1); i<_N; i++) {
        if (_distanceMatrix(i,j) >= _d(l-1) && _distanceMatrix(i,j) <= _d(l)) {
          tplanedist2_ij.pushback((_distanceTplane.compute_distance(res(i), res(j)))^2);
          w_ij.pushback(weights(i)*weights(j));
        }
      }
    }
    unsigned int actual_size = tplanedist2_ij.size();
    if (actual_size>0) {
      _N_hvec.pushback(actual_size);
      double num_sum = 0;
      double denom_sum = 0;
      for (auto i=0; i <actual_size; i++) {
        num_sum += tplanedist2_ij[i]*w_ij[i];
        denom_sum += w_ij[i];
      }
      _emp_vario_values.pushbach(num_sum/(2*denom_sum));
      _hvec.pushback((_d(l)+_d(l-1))/2);
    }
  }
}

// FittedVariogram
double FittedVariogram::get_tau2() const {
  return _parameters(0)
}

double FittedVariogram::get_sigma2() const {
  return _parameters(1)
}

double FittedVariogram::get_a() const{
  return _parameters(2);
}

void FittedVariogram::evaluate_par_fitted(const EmpiricalVariogram & emp_vario){

  double c = 1e-4; //Valori presi da libro quarteroni
  double s = 0.25;
  double tol = 1e-4;
  bool converged = false;
  unsigned int iter = 0;
  unsigned int max_iter = 100;

  get_init_par(emp_vario);

  Vec emp_vario_values = emp_vario.get_emp_vario_values();
  Vec h_vec = emp_vario.get_hvec();
  auto compute_residual = [&h_vec, &emp_vario_values](){return (get_vario_vec(h_vec) - emp_vario_values)}

  Vec residuals = compute_residual();
  Vec new_residuals = residuals;

  MatrixXd J = compute_jacobian(h_vec);

  // GaussNewton
  Vec gk =  J.transpose()*residuals;

  while(!converged && iter < max_iter){

    iter++;
    Vec dir = J.fullPivHouseholderQr().solve(-residuals);

    residuals = new_residuals;

    backtrack(dir, gk, new_residuals, J, h_vec, c,s, compute_residual)

    converged = (abs(residuals.squaredNorm() - new_residuals.squaredNorm()) < tol);

  }


}

void FittedVariogram::backtrack(const Vec &dir,Vec &gk,Vec &res,Matrixd &J,const Vec & h_vec, double c, double s, std::function<Vec()> compute_residual){

  double alpha = 1;
  const double alphamin = 1.e-5;
  Vec residuals = res;
  double fk = residuals.squaredNorm();
  Vector3d parameters_k = _parameters;
  _parameters = parameters_k + alpha*dir;
  J = compute_jacobian(h_vec)
  Vec new_residuals =  compute_residual();
  while(new_residuals.squaredNorm() > fk + alpha*c*gk.transpose()*dir && alpha > alphamin){

    alpha = alpha*s;
    _parameters = parameters_k + alpha*dir;
    residuals = new_residuals;
    fk = residuals.squaredNorm();
    new_residuals = compute_residual();
    gk =  J.transpose()*residuals;
    J = compute_jacobian(h_vec);

  }

}

double FittedVariogram::weighted_median (const Vec & values, const VectorXui & card) {
  unsigned int N = card.sum();
  VectorXui cumsums = cumsum(card);
  size_t index = 0;
  if (N%2==1) {
    while (cumsums(index)< (N+1)/2) { index++;};
  }
  else {
    while (cumsums(index)< N/2) { index++;};
  }
  return values(index);
}

// GaussianVariogram
double GaussVariogram::get_vario_univ(const double & h) const {
  double vario_value;
  if (h==0) {
    vario_value = 0;
  }
  else {
    vario_value = _parameters(0) + _parameters(1)*(1-exp(-(h^2)/(_parameters(2)^2)));
  }
  return vario_value
}

double GaussVariogram::get_covario_univ(const double & h) const {
  return _parameters(0) + _parameters(1) - get_vario_univ(h)
}

Vec GaussVariogram::get_vario_vec(const Vec & h_vec) const {
  unsigned int card_h = h_vec.size();
  VectorXd vario_values (card_h);
  for (auto i=0; i<card_h; i++) {
    vario_values(i) = get_vario_univ(h_vec(i));
  }
  return vario_values
}

MatrixXd GaussVariogram::compute_jacobian(const Vec & h_vec) const {
  unsigned int card_h = h_vec.size();
  MatrixXd jacobian (card_h,3);

  for (auto i=0; i<card_h; i++) {
    double tmp = exp(-(h_vec(i))^2/(_parameters(2))^2)
    jacobian(i,1) = 1;
    jacobian(i,2) = 1-tmp;
    jacobian(i,3) = -2*(h_vec(i))^2*_parameters(1)*tmp/(_parameters(2))^3;
  }
  return jacobian;
}

void GaussVariogram::get_init_par(const EmpiricalVariogram & emp_vario) {
  Vec emp_vario_values = emp_vario.get_emp_vario_values();
  Vec N_hvec = emp_vario.get_N_hvec();
  Vec hvec = emp_vario.get_hvec();
  unsigned int card_h =  N_hvec.size();

  Vec first_two = emp_vario_values.head(2);
  Vec last_four = emp_vario_values.tail(4); // Oppure uso end-3:end?
  VectorXui N_h_first_two = N_hvec.head(2);
  VectorXui N_h_last_four = N_hvec.tail(4);

  double sill = weighted_median(last_four, N_h_last_four;
  _parameters(0) = weighted_median(first_two, N_h_first_two);
  _parameters(1) = max(sill-_parameters(0), _parameters(0)*1e-3);

  tol = 0.0505*sill;
  size_t i = 0;
  while (abs(emp_vario_values(i) -0.95*sill) > tol) {
    i++;
  }
  _parameters(2) = 1/3*hvec(i);


}
