#include "Variogram.hpp"

using namespace VariogramEvaluation;

void FittedVariogram::evaluate_par_fitted(const EmpiricalVariogram & emp_vario){

  double c = 1e-4; //Valori presi da libro quarteroni
  double s = 0.25;
  double tol = 1e-4;
  bool converged = false;
  unsigned int iter = 0;
  unsigned int max_iter = 100;

  get_init_par(emp_vario);

  Vec emp_vario_values = emp_vario.get_emp_vario_values();
  Vec h_vec = emp_vario.get_hvec();
  auto compute_residual = [&h_vec, &emp_vario_values](){return (get_vario_vec(h_vec) - emp_vario_values)}

  Vec residuals = compute_residual();
  Vec new_residuals = residuals;

  MatrixXd J = compute_jacobian(h_vec);

  // GaussNewton
  Vec gk =  J.transpose()*residuals;

  while(!converged && iter < max_iter){

    iter++;
    Vec dir = J.fullPivHouseholderQr().solve(-residuals);

    residuals = new_residuals;

    backtrack(dir, gk, new_residuals, J, h_vec, c,s, compute_residual)

    converged = (abs(residuals.squaredNorm() - new_residuals.squaredNorm()) < tol);

  }


}

void FittedVariogram::backtrack(const Vec &dir,Vec &gk,Vec &res,Matrixd &J,const Vec & h_vec, double c, double s, std::function<Vec()> compute_residual){

  double alpha = 1;
  const double alphamin = 1.e-5;
  Vec residuals = res;
  double fk = residuals.squaredNorm();
  Vec parameters_k = parameters;
  parameters = parameters_k + alpha*dir;
  J = compute_jacobian(h_vec)
  Vec new_residuals =  compute_residual();
  while(new_residuals.squaredNorm() > fk + alpha*c*gk.transpose()*dir && alpha > alphamin){

    alpha = alpha*s;
    parameters = parameters_k + alpha*dir;
    residuals = new_residuals;
    fk = residuals.squaredNorm();
    new_residuals = compute_residual();
    gk =  J.transpose()*residuals;
    J = compute_jacobian(h_vec);

  }

}
