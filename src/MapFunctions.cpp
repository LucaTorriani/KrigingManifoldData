#include "MapFunctions.hpp"
#include<iostream>
using namespace map_functions;

// *** Logarithmic Map ***

//LOGMAPFROB
logMapFrob::logMapFrob(const MatrixXd& sqrtSigma, const MatrixXd& sqrtSigmaInv):_sqrtSigma(sqrtSigma),_sqrtSigmaInv(sqrtSigmaInv){};

MatrixXd logMapFrob::operator()(const MatrixXd& M) const {
  unsigned int n(_sqrtSigmaInv.cols());

  MatrixXd prod(n,n);
  prod = matrix_manipulation::logMat(_sqrtSigmaInv*M*_sqrtSigmaInv);

  MatrixXd result(n,n);
  result = _sqrtSigma*prod*_sqrtSigma;

  return result;
}

//LOGMAPLOEGEUCL
logMapLogEucl::logMapLogEucl(const MatrixXd& Sigma):_Sigma(Sigma){};

MatrixXd logMapLogEucl::operator()(const MatrixXd& M) const{
  return (matrix_manipulation::logMat(M) - matrix_manipulation::logMat(_Sigma));
}

//LOGMAPSQROOT

logMapSqRoot::logMapSqRoot(const MatrixXd& Sigma): _Sigma(Sigma){};

MatrixXd logMapSqRoot::operator()(const MatrixXd& M) const{
  return (matrix_manipulation::sqrtMat(M) - matrix_manipulation::sqrtMat(_Sigma));
}

//LOGARITHMICMAP

logarithmicMap::logarithmicMap(const distances_manifold::DistanceManifold& distanceManifoldObj): _distanceManifold(distanceManifoldObj.get_distanceType()) {

  MatrixXd Sigma(distanceManifoldObj.get_Sigma());

  unsigned int n = Sigma.cols();
  MatrixXd sqrtSigma(n,n);
  sqrtSigma =  matrix_manipulation::sqrtMat(Sigma);

  Eigen::LDLT<MatrixXd> solver(n); // Piu veloce specificando prima la dimensione
  solver.compute(sqrtSigma);
  MatrixXd Id(n,n);
  Id.setIdentity();
  MatrixXd sqrtSigmaInv(n,n);
  sqrtSigmaInv = solver.solve(Id);

  maps.insert(std::pair<std::string, std::function<MatrixXd(const MatrixXd&)>> ("Frobenius", logMapFrob(sqrtSigma, sqrtSigmaInv)));
  maps.insert(std::pair<std::string, std::function<MatrixXd(const MatrixXd&)>> ("SquareRoot", logMapSqRoot(Sigma)));
  maps.insert(std::pair<std::string, std::function<MatrixXd(const MatrixXd&)>> ("LogEuclidean", logMapLogEucl(Sigma)));
}

MatrixXd logarithmicMap::map2tplane(const MatrixXd& M) const{
  return maps.at(_distanceManifold)(M);
}


// *** Exponential Map ***

//EXPMAPFROB
expMapFrob::expMapFrob(const MatrixXd& sqrtSigma, const MatrixXd& sqrtSigmaInv):_sqrtSigma(sqrtSigma),_sqrtSigmaInv(sqrtSigmaInv){};

MatrixXd expMapFrob::operator()(const MatrixXd& M) const{
  unsigned int n(_sqrtSigmaInv.cols());
  MatrixXd prod(n,n);
  prod = matrix_manipulation::expMat(_sqrtSigmaInv*M*_sqrtSigmaInv);

  MatrixXd result(n,n);
  result = _sqrtSigma*prod*_sqrtSigma;

  return result;
}

//EXPMAPLOEGEUCL
expMapLogEucl::expMapLogEucl(const MatrixXd& Sigma):_Sigma(Sigma){};

MatrixXd expMapLogEucl::operator()(const MatrixXd& M) const{
  unsigned int n(M.cols());

  MatrixXd tmp(n,n);
  tmp = matrix_manipulation::logMat(_Sigma) + M;
  MatrixXd result(n,n);
  result = MatrixXd(tmp).transpose()*MatrixXd(tmp);

  return (result);
}

//LOGMAPSQROOT

expMapSqRoot::expMapSqRoot(const MatrixXd&  Sigma): _Sigma(Sigma){};

MatrixXd expMapSqRoot::operator()(const MatrixXd& M) const{
  unsigned int n(M.cols());

  MatrixXd tmp(n,n);
  tmp = matrix_manipulation::sqrtMat(_Sigma) + M;
  MatrixXd result(n,n);
  result = tmp.transpose()*tmp;

  return (result);
}

//LOGARITHMICMAP

exponentialMap::exponentialMap(const distances_manifold::DistanceManifold& distanceManifoldObj): _distanceManifold(distanceManifoldObj.get_distanceType()) {

  MatrixXd Sigma(distanceManifoldObj.get_Sigma());

  unsigned int n = Sigma.cols();
  MatrixXd sqrtSigma(n,n);
  sqrtSigma =  matrix_manipulation::sqrtMat(Sigma);

  Eigen::LDLT<MatrixXd> solver(n);
  solver.compute(sqrtSigma);
  MatrixXd Id(n,n);
  Id.setIdentity();
  MatrixXd sqrtSigmaInv(n,n);
  sqrtSigmaInv = solver.solve(Id);

  maps.insert(std::pair<std::string, std::function<MatrixXd(const MatrixXd&)>> ("Frobenius", expMapFrob(sqrtSigma, sqrtSigmaInv)));
  maps.insert(std::pair<std::string, std::function<MatrixXd(const MatrixXd&)>> ("SquareRoot", expMapSqRoot(Sigma)));
  maps.insert(std::pair<std::string, std::function<MatrixXd(const MatrixXd&)>> ("LogEuclidean", expMapLogEucl(Sigma)));
}

MatrixXd exponentialMap::map2manifold(const MatrixXd& M) const{
  return maps.at(_distanceManifold)(M);
}
