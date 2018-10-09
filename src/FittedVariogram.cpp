#include "FittedVariogram.hpp"
#include <cmath>

using namespace variogram_evaluation;

// FittedVariogram
double FittedVariogram::get_tau2() const {
  return _parameters(0);
}

double FittedVariogram::get_sigma2() const {
  return _parameters(1);
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
  unsigned int card_h(emp_vario.get_card_h());
  std::vector<double> h_vec(card_h);
  h_vec = emp_vario.get_hvec();

  std::vector<double> emp_vario_val(emp_vario.get_emp_vario_values());
  Vec emp_vario_values = Eigen::Map<Vec, Eigen::Unaligned>(emp_vario_val.data(), emp_vario_val.size());
  // TRASFORMAZIONE STD::VECTOR<DOUBLE> --> VEC: https://stackoverflow.com/questions/17036818/initialise-eigenvector-with-stdvector

  Vec vario_residuals = get_vario_vec(h_vec, card_h)  - emp_vario_values;
  Vec new_vario_residuals = vario_residuals;

  MatrixXd J (card_h, 3);
  J = compute_jacobian(h_vec, card_h);

  // GaussNewton
  Vector3d gk(J.transpose()*vario_residuals);
  Matrix3d JJ;
  Vector3d bb;
  LDLT<Matrix3d> solver(3);
  Vector3d dir;

  while((!converged) && iter < max_iter){

    JJ = J.transpose()*J;
    solver.compute(JJ);
    bb = - J.transpose()*(vario_residuals);
    dir = solver.solve(bb);

    vario_residuals = new_vario_residuals;
    backtrack(dir, gk, new_vario_residuals, h_vec, card_h, c,s, emp_vario_values);
    J = compute_jacobian(h_vec, card_h);
    gk =  J.transpose()*new_vario_residuals;
    converged = (std::abs(vario_residuals.squaredNorm() - new_vario_residuals.squaredNorm()) < tol);
    iter++;
  }
}

void FittedVariogram::backtrack(const Vector3d &dir,Vector3d &gk, Vec &res,const std::vector<double> & h_vec, unsigned int card_h, double c, double s, const Vec& emp_vario_values){

  double alpha = 1;
  const double alphamin = 1.e-5;
  double fk = res.squaredNorm();

  Vector3d parameters_k = _parameters;
  _parameters = parameters_k + alpha*dir;
  if (_parameters(0) < 0) _parameters(0) = 1e-7;
  res = get_vario_vec(h_vec, card_h)  - emp_vario_values;

  while(res.squaredNorm() > fk + alpha*c*gk.transpose()*dir && alpha > alphamin){

    alpha = alpha*s;
    _parameters = parameters_k + alpha*dir;
    if (_parameters(0) < 0) _parameters(0) = 1e-7;

    res = get_vario_vec(h_vec, card_h) - emp_vario_values;

  }
}

double FittedVariogram::weighted_median (const std::vector<double> & values, const std::vector<unsigned int> & card) {
  unsigned int num_values = card.size();
  std::vector<int> cumsums (num_values);
  cumsums[0]= card[0];
  for (size_t i=1; i<num_values; i++) {
    cumsums[i]=card[i]+cumsums[i-1];
  }

  unsigned int N = cumsums[num_values-1];
  size_t index = 0;
  if (N%2==1) {
    while (cumsums[index]< (N+1)/2) { index++;};
  }
  else {
    while (cumsums[index]< N/2) { index++;};
  }
  return values[index];
}

Vec FittedVariogram::get_vario_vec(const std::vector<double> & h_vec, unsigned int card_h) const {
  Vec vario_values(card_h);
  for (size_t i=0; i<card_h; i++) {
    vario_values(i) = get_vario_univ(h_vec[i]);
  }
  return vario_values;
}

Vec FittedVariogram::get_vario_vec(const Vec & h_vec, unsigned int card_h) const {
  Vec vario_values(card_h);
  for (size_t i=0; i<card_h; i++) {
    vario_values(i) = get_vario_univ(h_vec(i));
  }
  return vario_values;
}

double FittedVariogram::get_covario_univ(const double & h) const {
  return _parameters(0) + _parameters(1) - get_vario_univ(h);
}

MatrixXd FittedVariogram::compute_gamma_matrix(const std::shared_ptr<const SpMat> distanceMatrix, unsigned int N) const{
  MatrixXd gamma_matrix(N,N);
  double c0 = get_covario_univ(0);

  for (size_t i=0; i<(N-1); i++ ) {
    gamma_matrix(i,i) = c0;
    for (size_t j=(i+1); j<N; j++ ) {
      gamma_matrix(i,j) = get_covario_univ(distanceMatrix->coeff(i,j));
      gamma_matrix(j,i) = gamma_matrix(i,j);
    }
  }
  gamma_matrix(N-1,N-1) = c0;

  return (gamma_matrix);
}

Vec FittedVariogram::get_covario_vec(const std::vector<double> & h_vec, unsigned int card_h) const {
  Vec covario_values(card_h);
  for (size_t i=0; i<card_h; i++) {
    covario_values(i) = get_covario_univ(h_vec[i]);
  }
  return covario_values;
}

void FittedVariogram::set_parameters(const Vector3d& parameters) {
  _parameters = parameters;
}

// GaussianVariogram
double GaussVariogram::get_vario_univ(const double & h) const {
  double vario_value;
  if (h==0) {
    vario_value = 0;
  }
  else {
    vario_value = _parameters(0) + _parameters(1)*(1-exp(-(h*h)/(_parameters(2)*_parameters(2))));
  };
  return vario_value;
}

MatrixXd GaussVariogram::compute_jacobian(const std::vector<double> & h_vec, unsigned int card_h) const {
  MatrixXd jacobian (card_h,3);

  for (size_t i=0; i<card_h; i++) {
    double tmp = exp(-(h_vec[i]*h_vec[i])/(_parameters(2)*_parameters(2)));
    jacobian(i,0) = 1.0;
    jacobian(i,1) = 1.0-tmp;
    jacobian(i,2) = -2*(h_vec[i]*h_vec[i])*_parameters(1)*tmp/(_parameters(2)*_parameters(2)*_parameters(2));
  }
  return jacobian;
}

void GaussVariogram::get_init_par(const EmpiricalVariogram & emp_vario) {
  std::vector<double> emp_vario_values(emp_vario.get_emp_vario_values());
  std::vector<unsigned int> N_hvec(emp_vario.get_N_hvec());
  std::vector<double> hvec(emp_vario.get_hvec());
  unsigned int card_h(emp_vario.get_card_h());

  if (card_h == 0) {
    throw std::invalid_argument("No empirical variogram values in the cell");
  }
  else if (card_h==1) {
    _parameters(0)= emp_vario_values[0];
    if (_parameters(0) == 0) _parameters(0) = 1e-6;
    _parameters(1)= _parameters(0);
    _parameters(2) = 1.0/3*hvec[0];  // ?
    Rcpp::warning("Only 1 empirical variogram value in the cell");
  }
  else if (card_h==2){
    double tmp = weighted_median(emp_vario_values, N_hvec);
    if(tmp==0) tmp = 1e-6;
    _parameters(0)= tmp;
    _parameters(1)= tmp;
    _parameters(2) = 1.0/3*hvec[0]; // ?
    Rcpp::warning("Only 2 empirical variogram values in the cell");
  }
  else if (card_h==3) {
    std::vector<double> first_two(2);
    first_two[0] = emp_vario_values[0];
    first_two[1] = emp_vario_values[1];
    std::vector<unsigned int> N_h_first_two(2);
    N_h_first_two[0] = N_hvec[0];
    N_h_first_two[1] = N_hvec[1];

    double sill = weighted_median(emp_vario_values, N_hvec);
    _parameters(0) = weighted_median(first_two, N_h_first_two);
    if (_parameters(0) == 0) _parameters(0) = 1e-6;
    _parameters(1) = std::max(sill-_parameters(0), _parameters(0)*1e-3);

    double tol = 0.0505*sill;
    size_t i = 0;
    while (std::abs(emp_vario_values[i]-0.95*sill) > tol) {
      i++;
    }
    _parameters(2) = 1.0/3*hvec[i];
    Rcpp::warning("Only 3 empirical variogram values in the cell");

  }
  else {
    std::vector<double> first_two(2);
    first_two[0] = emp_vario_values[0];
    first_two[1] = emp_vario_values[1];
    std::vector<unsigned int> N_h_first_two(2);
    N_h_first_two[0] = N_hvec[0];
    N_h_first_two[1] = N_hvec[1];

    std::vector<double> last_four(4);
    last_four[0] = emp_vario_values[card_h-4];
    last_four[1] = emp_vario_values[card_h-3];
    last_four[2] = emp_vario_values[card_h-2];
    last_four[3] = emp_vario_values[card_h-1];
    std::vector<unsigned int> N_h_last_four(4);
    N_h_last_four[0] = N_hvec[card_h-4];
    N_h_last_four[1] = N_hvec[card_h-3];
    N_h_last_four[2] = N_hvec[card_h-2];
    N_h_last_four[3] = N_hvec[card_h-1];

    double sill = weighted_median(last_four, N_h_last_four);
    _parameters(0) = weighted_median(first_two, N_h_first_two);
    if (_parameters(0) == 0) _parameters(0) = 1e-6;
    _parameters(1) = std::max(sill-_parameters(0), _parameters(0)*1e-3);

    double tol = 0.0505*sill;
    size_t i = 0;
    while (std::abs(emp_vario_values[i]-0.95*sill) > tol) {
      i++;
    }
    _parameters(2) = 1.0/3*hvec[i];
  }
}

// ExponentialVariogram
double ExpVariogram::get_vario_univ(const double & h) const {
  double vario_value;
  if (h==0) {
    vario_value = 0;
  }
  else {
    vario_value = _parameters(0) + _parameters(1)*(1-exp(-h/_parameters(2)));
  };
  return vario_value;
}

MatrixXd ExpVariogram::compute_jacobian(const std::vector<double> & h_vec, unsigned int card_h) const {
  MatrixXd jacobian (card_h,3);

  for (size_t i=0; i<card_h; i++) {
    double tmp = exp(-h_vec[i]/_parameters(2));
    jacobian(i,0) = 1.0;
    jacobian(i,1) = 1.0-tmp;
    jacobian(i,2) = -h_vec[i]* _parameters(1)*tmp/(_parameters(2)*_parameters(2));
  }
  return jacobian;
}

void ExpVariogram::get_init_par(const EmpiricalVariogram & emp_vario) {
  std::vector<double> emp_vario_values(emp_vario.get_emp_vario_values());
  std::vector<unsigned int> N_hvec(emp_vario.get_N_hvec());
  std::vector<double> hvec(emp_vario.get_hvec());
  unsigned int card_h(emp_vario.get_card_h());

  if (card_h == 0) {
    throw std::invalid_argument("No empirical variogram values in the cell");
  }
  else if (card_h==1) {
    _parameters(0)= emp_vario_values[0];
    if (_parameters(0) == 0) _parameters(0) = 1e-6;
    _parameters(1)= _parameters(0);
    _parameters(2) = 1.0/3*hvec[0];  // ?
    Rcpp::warning("Only 1 empirical variogram value in the cell");
  }
  else if (card_h==2){
    double tmp = weighted_median(emp_vario_values, N_hvec);
    if(tmp==0) tmp = 1e-6;
    _parameters(0)= tmp;
    _parameters(1)= tmp;
    _parameters(2) = 1.0/3*hvec[0]; // ?
    Rcpp::warning("Only 2 empirical variogram values in the cell");
  }
  else if (card_h==3) {
    std::vector<double> first_two(2);
    first_two[0] = emp_vario_values[0];
    first_two[1] = emp_vario_values[1];
    std::vector<unsigned int> N_h_first_two(2);
    N_h_first_two[0] = N_hvec[0];
    N_h_first_two[1] = N_hvec[1];

    double sill = weighted_median(emp_vario_values, N_hvec);
    _parameters(0) = weighted_median(first_two, N_h_first_two);
    if (_parameters(0) == 0) _parameters(0) = 1e-6;
    _parameters(1) = std::max(sill-_parameters(0), _parameters(0)*1e-3);

    double tol = 0.0505*sill;
    size_t i = 0;
    while (std::abs(emp_vario_values[i]-0.95*sill) > tol) {
      i++;
    }
    _parameters(2) = 1.0/3*hvec[i];
    Rcpp::warning("Only 3 empirical variogram values in the cell");
  }
  else {
    std::vector<double> first_two(2);
    first_two[0] = emp_vario_values[0];
    first_two[1] = emp_vario_values[1];
    std::vector<unsigned int> N_h_first_two(2);
    N_h_first_two[0] = N_hvec[0];
    N_h_first_two[1] = N_hvec[1];

    std::vector<double> last_four(4);
    last_four[0] = emp_vario_values[card_h-4];
    last_four[1] = emp_vario_values[card_h-3];
    last_four[2] = emp_vario_values[card_h-2];
    last_four[3] = emp_vario_values[card_h-1];
    std::vector<unsigned int> N_h_last_four(4);
    N_h_last_four[0] = N_hvec[card_h-4];
    N_h_last_four[1] = N_hvec[card_h-3];
    N_h_last_four[2] = N_hvec[card_h-2];
    N_h_last_four[3] = N_hvec[card_h-1];

    double sill = weighted_median(last_four, N_h_last_four);
    _parameters(0) = weighted_median(first_two, N_h_first_two);
    if (_parameters(0) == 0) _parameters(0) = 1e-6;
    _parameters(1) = std::max(sill-_parameters(0), _parameters(0)*1e-3);

    double tol = 0.0505*sill;
    size_t i = 0;
    while (std::abs(emp_vario_values[i]-0.95*sill) > tol) {
      i++;
    }
    _parameters(2) = 1.0/3*hvec[i];
  }
}


// SphericalVariogram
double SphVariogram::get_vario_univ(const double & h) const {
  double vario_value;
  if (h==0) {
    vario_value = 0;
  }
  else if (h>=_parameters(2)) {
    vario_value = _parameters(0) + _parameters(1);
  }
  else {
    double tmp(h/_parameters(2));
    vario_value = _parameters(0) + _parameters(1)*(3.0/2*tmp - 1.0/2*tmp*tmp*tmp);
  };
  return vario_value;
}

MatrixXd SphVariogram::compute_jacobian(const std::vector<double> & h_vec, unsigned int card_h) const {

  MatrixXd jacobian (card_h,3);
  for (size_t i=0; i<card_h; i++) {
    jacobian(i,0) = 1.0;
    if (h_vec[i] >= _parameters(2)) {
      jacobian(i,1) = 1.0;
      jacobian(i,2) = 0.0;
    }
    else {
      double tmp1(h_vec[i]/_parameters(2));
      double tmp2(tmp1/_parameters(2));
      jacobian(i,1) = (3.0/2*tmp1 - 1.0/2*tmp1*tmp1*tmp1);
      jacobian(i,2) = -3.0/2*_parameters(1)*(tmp2-tmp2*tmp2*h_vec[i]);
    }
  }

  return jacobian;
}

void SphVariogram::get_init_par(const EmpiricalVariogram & emp_vario) {
  std::vector<double> emp_vario_values(emp_vario.get_emp_vario_values());
  std::vector<unsigned int> N_hvec(emp_vario.get_N_hvec());
  std::vector<double> hvec(emp_vario.get_hvec());
  unsigned int card_h(emp_vario.get_card_h());

  if (card_h == 0) {
    throw std::invalid_argument("No empirical variogram values in the cell");
  }
  else if (card_h==1) {
    _parameters(0)= emp_vario_values[0];
    if (_parameters(0) == 0) _parameters(0) = 1e-6;
    _parameters(1)= _parameters(0);
    _parameters(2) = 1.0/3*hvec[0];  // ?
    Rcpp::warning("Only 1 empirical variogram value in the cell");
  }
  else if (card_h==2){
    double tmp = weighted_median(emp_vario_values, N_hvec);
    if(tmp==0) tmp = 1e-6;
    _parameters(0)= tmp;
    _parameters(1)= tmp;
    _parameters(2) = 1.0/3*hvec[0]; // ?
    Rcpp::warning("Only 2 empirical variogram values in the cell");
  }
  else if (card_h==3) {
    std::vector<double> first_two(2);
    first_two[0] = emp_vario_values[0];
    first_two[1] = emp_vario_values[1];
    std::vector<unsigned int> N_h_first_two(2);
    N_h_first_two[0] = N_hvec[0];
    N_h_first_two[1] = N_hvec[1];

    double sill = weighted_median(emp_vario_values, N_hvec);
    _parameters(0) = weighted_median(first_two, N_h_first_two);
    if (_parameters(0) == 0) _parameters(0) = 1e-6;
    _parameters(1) = std::max(sill-_parameters(0), _parameters(0)*1e-3);

    double tol = 0.0505*sill;
    size_t i = 0;
    while (std::abs(emp_vario_values[i]-0.95*sill) > tol) {
      i++;
    }
    _parameters(2) = 1.0/3*hvec[i];
    Rcpp::warning("Only 3 empirical variogram values in the cell");
  }
  else {
    std::vector<double> first_two(2);
    first_two[0] = emp_vario_values[0];
    first_two[1] = emp_vario_values[1];
    std::vector<unsigned int> N_h_first_two(2);
    N_h_first_two[0] = N_hvec[0];
    N_h_first_two[1] = N_hvec[1];

    std::vector<double> last_four(4);
    last_four[0] = emp_vario_values[card_h-4];
    last_four[1] = emp_vario_values[card_h-3];
    last_four[2] = emp_vario_values[card_h-2];
    last_four[3] = emp_vario_values[card_h-1];
    std::vector<unsigned int> N_h_last_four(4);
    N_h_last_four[0] = N_hvec[card_h-4];
    N_h_last_four[1] = N_hvec[card_h-3];
    N_h_last_four[2] = N_hvec[card_h-2];
    N_h_last_four[3] = N_hvec[card_h-1];

    double sill = weighted_median(last_four, N_h_last_four);
    _parameters(0) = weighted_median(first_two, N_h_first_two);
    if (_parameters(0) == 0) _parameters(0) = 1e-6;
    _parameters(1) = std::max(sill-_parameters(0), _parameters(0)*1e-3);

    double tol = 0.01*sill;
    size_t i = 0;
    while (std::abs(emp_vario_values[i]-sill) > tol) {
      i++;
    }
    _parameters(2) = hvec[i];
  }
}

Vector3d FittedVariogram::get_parameters() const{
  return _parameters;
}
