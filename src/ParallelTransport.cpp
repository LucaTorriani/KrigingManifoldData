#include "ParallelTransport.hpp"

using namespace parallel_transport;

MatrixXd parallel_transport::trasport_to_TI(MatrixXd Sigma, MatrixXd V) {
  p = Sigma.rows();
  Eigen::LDLT<MatrixXd> solver(p);
  solver.compute(Sigma);
  MatrixXd Id(p,p);
  Id.setIdentity();
  MatrixXd SigmaInvSqrt(p,p);
  MatrixXd SigmaSqrt(p,p);
  MatrixXd tmp(p,p);

  SigmaInvSqrt = sqrtMat(solver.solve(Id));
  SigmaSqrt = sqrtMat(Sigma);
  tmp = SigmaSqrt* expMat(0.5*SigmaInvSqrt*V*SigmaInvSqrt) * SigmaInvSqrt;

  return(tmp*Sigma*tmp.transpose())
}

MatrixXd parallel_transport::trasport_to_TI(MatrixXd Lambda, MatrixXd V) {
  p = Lambda.rows();
  MatrixXd tmp(p,p);
  tmp = expMat(0.5*V);

  return(tmp*tmp.transpose())
}
