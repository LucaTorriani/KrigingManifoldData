#include "MapFunctions.hpp"
#include<iostream>
using namespace map_functions;

// *** Logarithmic Map ***

//LOGMAPFROB
logMapFrob::logMapFrob(const MatrixXd& sqrtSigma, const MatrixXd& sqrtSigmaInv):_sqrtSigma(sqrtSigma),_sqrtSigmaInv(sqrtSigmaInv){};

SpMat logMapFrob::operator()(const SpMat& M){
  unsigned int n(_sqrtSigmaInv.cols());
  SpMat MM (n,n);
  MM = M.selfadjointView<Lower>();

  SpMat prod(n,n);
  prod = matrix_manipulation::logMat(_sqrtSigmaInv*MM*_sqrtSigmaInv).selfadjointView<Lower>();

  MatrixXd tmp(n,n);
  tmp = _sqrtSigma*prod*_sqrtSigma;

  std::vector<TripType> tripletList;
  tripletList.reserve(n*(n+1)/2);
  for (size_t i=0; i < n; i++ ) {
    for (size_t j=0; j < i+1;j++ ) {
      tripletList.push_back(TripType(i,j,tmp(i,j)));
    }
  }

  SpMat result(n, n);
  result.setFromTriplets(tripletList.begin(), tripletList.end());

  return result;
}

//LOGMAPLOEGEUCL
logMapLogEucl::logMapLogEucl(const SpMat& Sigma):_Sigma(Sigma){};

SpMat logMapLogEucl::operator()(const SpMat& M){
  return (matrix_manipulation::logMat(M) - matrix_manipulation::logMat(_Sigma));
}

//LOGMAPSQROOT

logMapSqRoot::logMapSqRoot(const SpMat& Sigma): _Sigma(Sigma){};

SpMat logMapSqRoot::operator()(const SpMat& M){
  return (matrix_manipulation::sqrtMat(M) - matrix_manipulation::sqrtMat(_Sigma));
}

//LOGARITHMICMAP

logarithmicMap::logarithmicMap(const distances_manifold::DistanceManifold& distanceManifoldObj): _distanceManifold(distanceManifoldObj.get_distanceType()) {

  SpMat Sigma(distanceManifoldObj.get_Sigma());

  unsigned int n = Sigma.cols();
  SpMat sqrtSigmaSP(n,n);
  sqrtSigmaSP =  matrix_manipulation::sqrtMat(Sigma).selfadjointView<Lower>();

  MatrixXd sqrtSigma(n,n);
  sqrtSigma = MatrixXd(sqrtSigmaSP);

  // Computing inverse of sqrtSigma
  // Eigen::SimplicialLDLT<SpMat,Lower> solver(sqrtSigmaSP);
  // SpMat Id(n,n);
  // Id.setIdentity();
  // MatrixXd sqrtSigmaInv(n,n);
  // sqrtSigmaInv = solver.solve(Id);

  Eigen::LDLT<MatrixXd> solver(n); // Piu veloce specificando prima la dimensione
  solver.compute(sqrtSigma);
  MatrixXd Id(n,n);
  Id.setIdentity();
  MatrixXd sqrtSigmaInv(n,n);
  sqrtSigmaInv = solver.solve(Id);

  maps.insert(std::pair<std::string, std::function<SpMat(const SpMat&)>> ("Frobenius", logMapFrob(sqrtSigma, sqrtSigmaInv)));
  maps.insert(std::pair<std::string, std::function<SpMat(const SpMat&)>> ("SquareRoot", logMapSqRoot(Sigma)));
  maps.insert(std::pair<std::string, std::function<SpMat(const SpMat&)>> ("LogEuclidean", logMapLogEucl(Sigma)));
}

SpMat logarithmicMap::map2tplane(const SpMat& M){
  return maps[_distanceManifold](M);
}


// *** Exponential Map ***

//EXPMAPFROB
expMapFrob::expMapFrob(const MatrixXd& sqrtSigma, const MatrixXd& sqrtSigmaInv):_sqrtSigma(sqrtSigma),_sqrtSigmaInv(sqrtSigmaInv){};

SpMat expMapFrob::operator()(const SpMat& M){
  unsigned int n(_sqrtSigmaInv.cols());
  SpMat MM (n,n);
  MM = M.selfadjointView<Lower>();

  SpMat prod(n,n);
  prod = matrix_manipulation::expMat(_sqrtSigmaInv*MM*_sqrtSigmaInv).selfadjointView<Lower>();

  MatrixXd tmp(n,n);
  tmp = _sqrtSigma*prod*_sqrtSigma;

  std::vector<TripType> tripletList;
  tripletList.reserve(n*(n+1)/2);
  for (size_t i=0; i < n; i++ ) {
    for (size_t j=0; j < i+1;j++ ) {
      tripletList.push_back(TripType(i,j,tmp(i,j)));
    }
  }
  SpMat result(n, n);
  result.setFromTriplets(tripletList.begin(), tripletList.end());

  return result;
}

//EXPMAPLOEGEUCL
expMapLogEucl::expMapLogEucl(const SpMat& Sigma):_Sigma(Sigma){};

SpMat expMapLogEucl::operator()(const SpMat& M){
  unsigned int n(M.cols());
  SpMat tmpSP(n,n);
  tmpSP = (matrix_manipulation::logMat(_Sigma) + M).selfadjointView<Lower>();
  MatrixXd tmp(n,n);
  tmp = MatrixXd(tmpSP).transpose()*MatrixXd(tmpSP);

  std::vector<TripType> tripletList;
  tripletList.reserve(n*(n+1)/2);
  for (size_t i=0; i < n; i++ ) {
    for (size_t j=0; j < i+1;j++ ) {
      tripletList.push_back(TripType(i,j,tmp(i,j)));
    }
  }

  SpMat result(n, n);
  result.setFromTriplets(tripletList.begin(), tripletList.end());
  return (result);
}

//LOGMAPSQROOT

expMapSqRoot::expMapSqRoot(const SpMat&  Sigma): _Sigma(Sigma){};

SpMat expMapSqRoot::operator()(const SpMat& M){
  unsigned int n(M.cols());
  SpMat tmpSP(n,n);
  tmpSP = (matrix_manipulation::sqrtMat(_Sigma) + M).selfadjointView<Lower>();
  MatrixXd tmp(n,n);
  tmp = MatrixXd(tmpSP).transpose()*MatrixXd(tmpSP);
    std::vector<TripType> tripletList;
    tripletList.reserve(n*(n+1)/2);
    for (size_t i=0; i < n; i++ ) {
      for (size_t j=0; j < i+1;j++ ) {
        tripletList.push_back(TripType(i,j,tmp(i,j)));
      }
    }

    SpMat result(n, n);
    result.setFromTriplets(tripletList.begin(), tripletList.end());
    return (result);
}

//LOGARITHMICMAP

exponentialMap::exponentialMap(const distances_manifold::DistanceManifold& distanceManifoldObj): _distanceManifold(distanceManifoldObj.get_distanceType()) {

  SpMat Sigma(distanceManifoldObj.get_Sigma());

  unsigned int n = Sigma.cols();
  SpMat sqrtSigmaSP(n,n);
  sqrtSigmaSP =  matrix_manipulation::sqrtMat(Sigma).selfadjointView<Lower>();

  MatrixXd sqrtSigma(n,n);
  sqrtSigma = MatrixXd(sqrtSigmaSP);

  // Computing inverse of sqrtSigma
  // Eigen::SimplicialLDLT<SpMat,Lower> solver(sqrtSigmaSP);
  // SpMat Id(n,n);
  // Id.setIdentity();
  // MatrixXd sqrtSigmaInv(n,n);
  // sqrtSigmaInv = solver.solve(Id);

  Eigen::LDLT<MatrixXd> solver(n); // Piu veloce specificando prima la dimensione
  solver.compute(sqrtSigma);
  MatrixXd Id(n,n);
  Id.setIdentity();
  MatrixXd sqrtSigmaInv(n,n);
  sqrtSigmaInv = solver.solve(Id);

  maps.insert(std::pair<std::string, std::function<SpMat(const SpMat&)>> ("Frobenius", expMapFrob(sqrtSigma, sqrtSigmaInv)));
  maps.insert(std::pair<std::string, std::function<SpMat(const SpMat&)>> ("SquareRoot", expMapSqRoot(Sigma)));
  maps.insert(std::pair<std::string, std::function<SpMat(const SpMat&)>> ("LogEuclidean", expMapLogEucl(Sigma)));
}

SpMat exponentialMap::map2manifold(const SpMat& M){
  return maps[_distanceManifold](M);
}
