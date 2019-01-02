#include "ParallelTransport.hpp"

using namespace parallel_transport;

MatrixXd parallel_transport::transport_to_TI(MatrixXd Sigma, MatrixXd V) {
  unsigned int p = Sigma.rows();
  Eigen::LDLT<MatrixXd> solver(p);
  solver.compute(Sigma);
  MatrixXd Id(p,p);
  Id.setIdentity();
  MatrixXd SigmaInvSqrt(p,p);
  MatrixXd SigmaSqrt(p,p);
  MatrixXd tmp(p,p);

  SigmaInvSqrt = matrix_manipulation::sqrtMat(solver.solve(Id));
  SigmaSqrt = matrix_manipulation::sqrtMat(Sigma);
  tmp = SigmaSqrt* matrix_manipulation::expMat(0.5*SigmaInvSqrt*V*SigmaInvSqrt) * SigmaInvSqrt;

  return(tmp*Sigma*tmp.transpose());
}

MatrixXd parallel_transport::transport_from_TI(MatrixXd Lambda, MatrixXd V) {
  unsigned int p = Lambda.rows();
  MatrixXd tmp(p,p);
  tmp = matrix_manipulation::expMat(0.5*V);

  return(tmp*tmp.transpose());
}
