#include "MapFunctions.hpp"
#include<iostream>
using namespace map_functions;

// *** Logarithmic Map ***

//LOGMAPFROB
MatrixXd logMapFrob::map2tplane(const MatrixXd& M) const {
  unsigned int p(_sqrtSigmaInv.cols());

  MatrixXd prod(p,p);
  MatrixXd tmp(p,p);
  tmp  = _sqrtSigmaInv*M*_sqrtSigmaInv;
  prod = matrix_manipulation::logMat(tmp);

  MatrixXd result(p,p);
  result = _sqrtSigma*prod*_sqrtSigma;

  return result;
}

void logMapFrob::set_members(const MatrixXd& Sigma) {
  unsigned int p = Sigma.cols();
  _sqrtSigma =  matrix_manipulation::sqrtMat(Sigma);

  Eigen::LDLT<MatrixXd> solver(p);
  solver.compute(_sqrtSigma);
  MatrixXd Id(p,p);
  Id.setIdentity();
  _sqrtSigmaInv = solver.solve(Id);
}

//LOGMAPLOGEUCL
MatrixXd logMapLogEucl::map2tplane(const MatrixXd& M) const{
  return (matrix_manipulation::logMat(M) - matrix_manipulation::logMat(_Sigma));
}


void logMapLogEucl::set_members(const MatrixXd& Sigma) {
  _Sigma = Sigma;
}

//LOGMAPSQROOT
MatrixXd logMapSqRoot::map2tplane(const MatrixXd& M) const{
  return (matrix_manipulation::sqrtMat(M) - matrix_manipulation::sqrtMat(_Sigma));
}


void logMapSqRoot::set_members(const MatrixXd& Sigma) {
  _Sigma = Sigma;
}

// CORRELATION
MatrixXd logMapChol::map2tplane (const MatrixXd& H) const{
  unsigned int p(H.rows());
  MatrixXd result(p,p);
  result.setZero(p,p);
  result(0,0)=0;
  Vec H_vec;
  Vec Sigma_vec;
  Vec proj_diff;
  double proj_diff_norm;
  for (size_t i=1; i<p; i++) {
    H_vec = H.col(i);
    Sigma_vec = _Sigma.col(i);
    proj_diff = proj2tspace(H_vec-Sigma_vec, Sigma_vec);
    proj_diff_norm = proj_diff.norm();
    if (proj_diff_norm < 1e-10) result.col(i).setZero();
    else result.col(i) = acos((H_vec.transpose()*Sigma_vec).value()) * proj_diff/proj_diff_norm;
  }
  return result;
}

Vec logMapChol::proj2tspace(const Vec& x_vec, const Vec& sigma_vec) const {
  return (x_vec - ((x_vec.transpose()*sigma_vec).value())*sigma_vec);
}

void logMapChol::set_members(const MatrixXd& Sigma) {
  _Sigma = Sigma;
}

// *** Exponential Map ***

//EXPMAPFROB
MatrixXd expMapFrob::map2manifold(const MatrixXd& M) const{
  unsigned int p(_sqrtSigmaInv.cols());
  MatrixXd prod(p,p);
  MatrixXd tmp(p,p);
  tmp = _sqrtSigmaInv*M*_sqrtSigmaInv;
  prod = matrix_manipulation::expMat(tmp);

  MatrixXd result(p,p);
  result = _sqrtSigma*prod*_sqrtSigma;

  return result;
}

void expMapFrob::set_members(const MatrixXd& Sigma){

  unsigned int p = Sigma.cols();
  _sqrtSigma =  matrix_manipulation::sqrtMat(Sigma);

  Eigen::LDLT<MatrixXd> solver(p);
  solver.compute(_sqrtSigma);
  MatrixXd Id(p,p);
  Id.setIdentity();
  _sqrtSigmaInv = solver.solve(Id);
}

//EXPMAPLOEGEUCL
MatrixXd expMapLogEucl::map2manifold(const MatrixXd& M) const{
  unsigned int p(M.cols());

  MatrixXd tmp(p,p);
  tmp = matrix_manipulation::logMat(_Sigma) + M;
  MatrixXd result(p,p);
  result = tmp.transpose()*tmp;

  return (result);
}

void expMapLogEucl::set_members(const MatrixXd& Sigma){
  _Sigma = Sigma;
}

//LOGMAPSQROOT
MatrixXd expMapSqRoot::map2manifold(const MatrixXd& M) const{
  unsigned int p(M.cols());

  MatrixXd tmp(p,p);
  tmp = matrix_manipulation::sqrtMat(_Sigma) + M;
  MatrixXd result(p,p);
  result = tmp.transpose()*tmp;

  return (result);
}

void expMapSqRoot::set_members(const MatrixXd& Sigma){
  _Sigma = Sigma;
}

// CORRELATION
MatrixXd expMapChol::map2manifold (const MatrixXd& V) const{
  unsigned int p(V.rows());
  MatrixXd result(p,p);
  result.setZero(p,p);
  result(0,0)=1;
  for (size_t i=1; i<p; i++) {
    double col_norm (V.col(i).norm());
    if (col_norm < 1e-10) result.col(i) = _Sigma.col(i);
    else result.col(i) = cos(col_norm)* _Sigma.col(i) + sin(col_norm)*V.col(i)/col_norm;
  }
  return result;
}

void expMapChol::set_members(const MatrixXd& Sigma){
  _Sigma = Sigma;
}
