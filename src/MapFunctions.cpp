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
