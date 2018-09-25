#include "MapFunctions.hpp"
#include<iostream>
using namespace map_functions;

// *** Logarithmic Map ***

//LOGMAPFROB
MatrixXd logMapFrob::map2tplane(const MatrixXd& M) const {
  unsigned int n(_sqrtSigmaInv.cols());

  MatrixXd prod(n,n);
  MatrixXd tmp(n,n);
  tmp  = _sqrtSigmaInv*M*_sqrtSigmaInv;
  prod = matrix_manipulation::logMat(tmp);

  MatrixXd result(n,n);
  result = _sqrtSigma*prod*_sqrtSigma;

  return result;
}

void logMapFrob::initialize_members(const std::shared_ptr<const MatrixXd> Sigma) {
  unsigned int n = Sigma->cols();
  _sqrtSigma =  matrix_manipulation::sqrtMat(*(Sigma));

  Eigen::LDLT<MatrixXd> solver(n);
  solver.compute(_sqrtSigma);
  MatrixXd Id(n,n);
  Id.setIdentity();
  _sqrtSigmaInv = solver.solve(Id);
}

//LOGMAPLOGEUCL
MatrixXd logMapLogEucl::map2tplane(const MatrixXd& M) const{
  return (matrix_manipulation::logMat(M) - matrix_manipulation::logMat(*_Sigma));
}

void logMapLogEucl::initialize_members(const std::shared_ptr<const MatrixXd> Sigma) {
  _Sigma = Sigma;
}

//LOGMAPSQROOT
MatrixXd logMapSqRoot::map2tplane(const MatrixXd& M) const{
  return (matrix_manipulation::sqrtMat(M) - matrix_manipulation::sqrtMat(*_Sigma));
}

void logMapSqRoot::initialize_members(const std::shared_ptr<const MatrixXd> Sigma){
  _Sigma = Sigma;
}

// *** Exponential Map ***

//EXPMAPFROB
MatrixXd expMapFrob::map2manifold(const MatrixXd& M) const{
  unsigned int n(_sqrtSigmaInv.cols());
  MatrixXd prod(n,n);
  MatrixXd tmp(n,n);
  tmp = _sqrtSigmaInv*M*_sqrtSigmaInv;
  prod = matrix_manipulation::expMat(tmp);

  MatrixXd result(n,n);
  result = _sqrtSigma*prod*_sqrtSigma;

  return result;
}

void expMapFrob::initialize_members(const std::shared_ptr<const MatrixXd> Sigma){

  unsigned int n = Sigma->cols();
  _sqrtSigma =  matrix_manipulation::sqrtMat(*(Sigma));

  Eigen::LDLT<MatrixXd> solver(n);
  solver.compute(_sqrtSigma);
  MatrixXd Id(n,n);
  Id.setIdentity();
  _sqrtSigmaInv = solver.solve(Id);
}


//EXPMAPLOEGEUCL
MatrixXd expMapLogEucl::map2manifold(const MatrixXd& M) const{
  unsigned int n(M.cols());

  MatrixXd tmp(n,n);
  tmp = matrix_manipulation::logMat(*_Sigma) + M;
  MatrixXd result(n,n);
  result = tmp.transpose()*tmp;

  return (result);
}

void expMapLogEucl::initialize_members(const std::shared_ptr<const MatrixXd> Sigma){
  _Sigma = Sigma;
}

//LOGMAPSQROOT
MatrixXd expMapSqRoot::map2manifold(const MatrixXd& M) const{
  unsigned int n(M.cols());

  MatrixXd tmp(n,n);
  tmp = matrix_manipulation::sqrtMat(*_Sigma) + M;
  MatrixXd result(n,n);
  result = tmp.transpose()*tmp;

  return (result);
}

void expMapSqRoot::initialize_members(const std::shared_ptr<const MatrixXd> Sigma){
  _Sigma = Sigma;
}
