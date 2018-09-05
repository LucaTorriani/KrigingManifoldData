#include "Distance_Manifold.hpp"

#include <cmath>


using namespace distances_manifold;

// FROBENIUS
static double Frobenius::operator()(const SpMat& M1, const SpMat& M2) const{

  Eigen::SimplicialLDLT<SpMat,Lower> solver(M1);
  SpMat tmp = M2.selfadjointView<Lower>();
  Eigen::MatrixXd matrix_result(solver.solve(tmp));

  SelfAdjointEigenSolver<SpMat> eigensolver (A,EigenvaluesOnly);
  VectorXd eigenvalues =  eigensolver1.eigenvalues();

  double ssq = 0.0;
    for(auto i = 0;i < eigenvalues.size(); i++)  ssq += (std::log(std::abs(eigenvalues(i)))*std::log(std::abs(eigenvalues(i))));
  }
  return (std::sqrt(ssq));
}

// LOGEUCLIDEAN
static double LogEuclidean::operator()(const SpMat& M1, const SpMat& M2) const {
  return ( (matrix_manipulation::logMat(M1)-matrix_manipulation::logMat(M2)).norm())
}

// SQROOT
static double SqRoot::operator()(const SpMat& M1, const SpMat& M2) const{
  return ( (matrix_manipulation::sqrtMat(M1)-matrix_manipulation::sqrtMat(M2)).norm())

}

DistanceManifold::DistanceManifold(){
  distances.insert(std::pair<std::string, std::function<double(std::vector<double>, std::vector<double>)>>("Frobenius", Frobenius()));
  distances.insert(std::pair<std::string, std::function<double(std::vector<double>, std::vector<double>)>>("Log_euclidean", LogEuclidean()));
  distances.insert(std::pair<std::string, std::function<double(std::vector<double>, std::vector<double>)>>("Square_root", SqRoot()));

}
double DistanceManifold::compute_distance const (SpMat& M1, const SpMat& M2, const std::string & distance_type){
  double result = dist[distance_type](M1, M2);
  return result;
}