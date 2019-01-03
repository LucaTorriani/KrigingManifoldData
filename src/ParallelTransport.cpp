#include "ParallelTransport.hpp"

using namespace parallel_transport;

// MatrixXd parallel_transport::transport_to_TI(MatrixXd Sigma, MatrixXd V) {
//   unsigned int p = Sigma.rows();
//   Eigen::LDLT<MatrixXd> solver(p);
//   MatrixXd Id(p,p);
//   Id.setIdentity();
//   MatrixXd SigmaInvSqrt(p,p);
//   MatrixXd SigmaSqrt(p,p);
//   MatrixXd tmp(p,p);
//
//   SigmaSqrt = matrix_manipulation::sqrtMat(Sigma);
//   solver.compute(SigmaSqrt);
//   SigmaInvSqrt = solver.solve(Id);
//   tmp = SigmaSqrt* matrix_manipulation::expMat(0.5*SigmaInvSqrt*V*SigmaInvSqrt) * SigmaInvSqrt;
//
//   return(tmp*Sigma*tmp.transpose());
// }
//
// MatrixXd parallel_transport::transport_from_TI(MatrixXd Lambda, MatrixXd V) {
//   unsigned int p = Lambda.rows();
//   MatrixXd tmp(p,p);
//   tmp = matrix_manipulation::expMat(0.5*V);
//
//   return(tmp*tmp.transpose());
// }

MatrixXd parallel_transport::transport_to_TI(MatrixXd Sigma, MatrixXd V) {
  unsigned int p = Sigma.rows();
  Eigen::LDLT<MatrixXd> solver(p);
  MatrixXd Id(p,p);
  Id.setIdentity();
  MatrixXd SigmaInvSqrt(p,p);

  solver.compute(Sigma);
  SigmaInvSqrt = matrix_manipulation::sqrtMat(solver.solve(Id));

  return(SigmaInvSqrt*V*SigmaInvSqrt.transpose());
}

MatrixXd parallel_transport::transport_from_TI(MatrixXd Sigma, MatrixXd V) {
  unsigned int p = Sigma.rows();
  MatrixXd SigmaSqrt(p,p);
  SigmaSqrt = matrix_manipulation::sqrtMat(Sigma);

  return(SigmaSqrt*V*SigmaSqrt.transpose());
}
