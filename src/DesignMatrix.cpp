#include "DesignMatrix.hpp"

using namespace design_matrix;

// **** InterceptDM ***
MatrixXd InterceptDM::compute_design_matrix(const Coordinates& coords) const {
  unsigned int N(coords.get_N_station());
  MatrixXd Z(N,1);
  Z.setOnes(N,1);
  return(Z);
}

MatrixXd InterceptDM::compute_design_matrix(const Coordinates& coords, const MatrixXd& X) const{
  const unsigned int N(coords.get_N_station());
  const unsigned int n(X.cols());

  MatrixXd Z(N,n+1);
  Z.block(0,0,N,1) = MatrixXd::Ones(N,1);
  Z.block(0,1, N,n) = X;

  return(Z);
}

// **** Coord1DM ***
MatrixXd Coord1DM::compute_design_matrix(const Coordinates& coords) const{
  const unsigned int N(coords.get_N_station());
  MatrixXd Z(N,2);
  Z.block(0,0,N,1) = MatrixXd::Ones(N,1);
  Z.block(0,1,N,1) = (coords.get_coords())->col(0);
  return(Z);
}

MatrixXd Coord1DM::compute_design_matrix(const Coordinates& coords, const MatrixXd& X) const {
  const unsigned int N(coords.get_N_station());
  const unsigned int n(X.cols());

  MatrixXd Z(N,n+2);
  Z.block(0,0,N,1) = MatrixXd::Ones(N,1);
  Z.block(0,1,N,1) = (coords.get_coords())->col(0);
  Z.block(0,2, N,n) = X;

  return(Z);
}

// **** Coord2DM ***
MatrixXd Coord2DM::compute_design_matrix(const Coordinates& coords) const{
  const unsigned int N(coords.get_N_station());
  MatrixXd Z(N,2);
  Z.block(0,0,N,1) = MatrixXd::Ones(N,1);
  Z.block(0,1,N,1) = (coords.get_coords())->col(1);
  return(Z);
}

MatrixXd Coord2DM::compute_design_matrix(const Coordinates& coords, const MatrixXd& X) const {
  const unsigned int N(coords.get_N_station());
  const unsigned int n(X.cols());

  MatrixXd Z(N,n+2);
  Z.block(0,0,N,1) = MatrixXd::Ones(N,1);
  Z.block(0,1,N,1) = (coords.get_coords())->col(1);
  Z.block(0,2, N,n) = X;

  return(Z);
}


// **** AdditiveDM ***
MatrixXd AdditiveDM::compute_design_matrix(const Coordinates& coords) const{
  const unsigned int N(coords.get_N_station());
  MatrixXd Z(N,3);
  Z.block(0,0,N,1) = MatrixXd::Ones(N,1);
  Z.block(0,1,N,2) = *(coords.get_coords());
  return(Z);
}

MatrixXd AdditiveDM::compute_design_matrix(const Coordinates& coords, const MatrixXd& X) const {
  const unsigned int N(coords.get_N_station());
  const unsigned int n(X.cols());

  MatrixXd Z(N,n+3);
  Z.block(0,0,N,1) = MatrixXd::Ones(N,1);
  Z.block(0,1,N,2) = *(coords.get_coords());
  Z.block(0,3, N,n) = X;

  return(Z);
}

// *** DesignMatrixFactory ***
DesignMatrixFactory& DesignMatrixFactory::Instance(){
  static DesignMatrixFactory theDesignMatrixFactory;
  return theDesignMatrixFactory;
}

std::unique_ptr<DesignMatrix> DesignMatrixFactory::create(const Identifier& model_name) const {
  auto f = _storage.find(model_name);
  return (f==_storage.end()) ? std::unique_ptr<DesignMatrix>():std::unique_ptr<DesignMatrix>(f->second());
}

void DesignMatrixFactory::add(const Identifier& identifier, const Builder& builder){
  _storage.insert(std::make_pair(identifier, builder));

}


// *** Registration ***
void design_matrix::registerDesignMatrices(){
  DesignMatrixFactory& design_matrices = DesignMatrixFactory::Instance();
  design_matrices.add("Intercept", InterceptDM());
  design_matrices.add("Coord1", Coord1DM());
  design_matrices.add("Coord2", Coord2DM());
  design_matrices.add("Additive", AdditiveDM());

};
