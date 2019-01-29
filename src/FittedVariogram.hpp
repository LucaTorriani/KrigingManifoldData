#ifndef _FITTED_VARIOGRAM_HPP
#define _FITTED_VARIOGRAM_HPP

#include "EmpiricalVariogram.hpp"
#include <iostream>
#include <Rcpp.h>

/*! \file
  @brief Classes to fit a model variogram (Gaussian, Exponential or Spherical) to an empirical one
*/

namespace variogram_evaluation {
  /*!
    @brief	Abstract class for computation and storage of the fitted variogram
  */
class FittedVariogram{
protected:
  /*! Vector storing the three parameters of the fitted varriogram \e (nugget, \e sill-nugget, \e practical \e range) */
  Vector3d _parameters; // tau2, sigma2, a
  /*!
    @brief Compute weighted median
    @param values Values whose median must be computed
    @param card Weights
    @return Median of `values`, weighted using `card`
  */
  double weighted_median (const std::vector<double> & values, const std::vector<unsigned int> & card);
  /*!
    @brief Initialize `_parameters`
    @details The `_parameters` are initialized as follows: \n
      \f$ \mbox{\_parameters}\left(0\right)=\text{weighted\_median} \left(\mbox{\_emp\_vario\_values.head} \left(2\right), \mbox{\_N\_hvec.head} \left(2\right)\right) \f$ \n
      \f$ \mbox{\_parameters}\left(1\right)=\text{weighted\_median} \left(\mbox{\_emp\_vario\_values.tail} \left(4\right), \mbox{\_N\_hvec.tail} \left(4\right)\right) - \mbox{\_parameters}\left(0\right) \f$  \n
      \f$ \mbox{\_parameters}\left(2\right)= \mbox{\_hvec}\left(i^\ast \right) \f$  \n where \f$i^\ast\f$ is the first \f$i \text{ s.t.} \left| \mbox{\_emp\_vario\_values}\left(i\right)- \left( \mbox{\_parameters}\left(0\right) + \mbox{\_parameters}\left(1\right) \right) \right| < \mbox{tol} \f$
    @param emp_vario Empirical variogram
  */
  virtual void get_init_par(const EmpiricalVariogram & emp_vario) = 0;
  /*!
    @brief Update `_parameters` moving along `dir`
  */
  void backtrack(const Vector3d &dir,Vector3d &gk, Vec &res,const std::vector<double> & h_vec, unsigned int card_h, double c, double s, const Vec& emp_vario_values, double max_sill, double max_a);
  /*!
    @brief Compute the jacobian of the variogram residuals (which coincides with the one of the model variogram), according to the variogram type
  */
  virtual MatrixXd compute_jacobian(const std::vector<double> & h_vec, unsigned int card_h) const = 0;
public:
  /*!
    @brief Return `_parameters(0)`, i.e. the \e nugget
  */
  double get_tau2() const;
  /*!
    @brief Return `_parameters(1)`, i.e. the \e sill-nugget
  */
  double get_sigma2() const;
  /*!
    @brief Return `_parameters(2)`, i.e. the \e practical \e range
  */
  double get_a() const;
  /*!
    @brief Compute the parameters of the fitted variogram
    @note Like ::variogram_evaluation::FittedVariogram::evaluate_par_fitted_W, but different stopping criteria. This function is used when equal weights are considered
    @details The parameters are computed using Gauss-Newton with backtrack method to solve the non-linear least square problem:
    \f{equation*}{
    \begin{aligned}
        & \underset{\mbox{\_parameters}}{\arg\min}
        & & \underset{h \in \mbox{\_hvec}} {\sum} \left( \gamma_m \left(\mbox{\_parameters}, h\right) - \widehat{\gamma} \left(h\right) \right)^2 \\
        & \text{subject to}
        & & 0 \leq \mbox{\_parameters}\left(0\right) \\
        &&& 0 \leq \mbox{\_parameters}\left(1\right) \leq \mbox{max\_sill} - \mbox{\_parameters}\left(0\right) \\
        &&& 0 \leq \mbox{\_parameters}\left(2\right) \leq \mbox{max\_a}
    \end{aligned}
    \f}
    The starting values for the `_parameters` are obtained through ::variogram_evaluation::FittedVariogram::get_init_par. \n
    The stopping criteria is based on the decrease of the error norm.
    @param emp_vario Empirical variogram
    @param max_sill Maximum value for the \e sill
    @param max_a Maximum value for  \e a
  */
  void evaluate_par_fitted_E(const EmpiricalVariogram & emp_vario, double max_sill, double max_a);
  /*!
    @brief Compute the parameters of the fitted variogram
    @note Like ::variogram_evaluation::FittedVariogram::evaluate_par_fitted_E, but different stopping criteria. This function is used when kernel weights are considered
    @details The parameters are computed using Gauss-Newton with backtrack method to solve the non-linear least square problem:
    \f{equation*}{
    \begin{aligned}
        & \underset{\text{\_parameters}}{\arg\min}
        & & \underset{h \in \mbox{\_hvec}} {\sum} \left( \gamma_m \left(\mbox{\_parameters}, h\right) - \widehat{\gamma} \left(h\right) \right)^2 \\
        & \text{subject to}
        & & 0 \leq \mbox{\_parameters}\left(0\right) \\
        &&& 0 \leq \mbox{\_parameters}\left(1\right) \leq \mbox{max\_sill} - \mbox{\_parameters}\left(0\right) \\
        &&& 0 \leq \mbox{\_parameters}\left(2\right) \leq \mbox{max\_a}
    \end{aligned}
    \f}
    The starting values for the `_parameters` are obtained through ::variogram_evaluation::FittedVariogram::get_init_par.  \n
    The stopping criteria is based on the difference in the decrease of the error norm between two consecutive iterations.
    @param emp_vario Empirical variogram
    @param max_sill Maximum value for the \e sill
    @param max_a Maximum value for  \e a
  */
  void evaluate_par_fitted_W(const EmpiricalVariogram & emp_vario, double max_sill, double max_a);
  /*!
    @brief Compute the value of the model variogram at a given distance, according to the variogram type
    @param h The distance where to evaluate the variogram
    @return Variogram value at distance h: \f$ \gamma_m \left(h\right) \f$
  */
  virtual double get_vario_univ(const double & h) const = 0;
  /*!
    @brief Compute the value of the model covariogram at a given distance, according to the variogram type
    @param h The distance where to evaluate the covariogram
    @return Covariogram value at distance \f$ h: C_m \left(h\right) = \left(\mbox{\_parameters}\left(0\right) + \mbox{\_parameters}\left(1\right) \right) - \gamma_m \left(h\right) \f$
  */
  double get_covario_univ(const double & h) const;
  /*!
    @brief Compute the values of the model variogram at a given vector of distances, according to the variogram type
    @param h_vec The distances where to evaluate the variogram
    @param card_h Number of distances where the variogram must be computed. \f$ card\_h = h\_vec.size() \f$
    @return Vector of variogram values at distances \f$ h \in  h\_vec \f$
  */
  Vec get_vario_vec(const std::vector<double> & h_vec, unsigned int card_h) const;
  /*!
    @brief Compute the values of the model variogram at a given vector of distances, according to the variogram type
    @param h_vec The distances where to evaluate the variogram
    @param card_h Number of distances where the variogram must be computed. \f$ card\_h = h\_vec.size() \f$
    @return Vector of variogram values at distances  \f$ h  \in  h\_vec \f$
  */
  Vec get_vario_vec(const Vec & h_vec, unsigned int card_h) const;
  /*!
    @brief Compute the covariogram matrix, according to the variogram type
    @param distanceMatrix Matrix of distances among the \f$N\f$ locations
    @param N Number of data locations
    @return Covariogram matrix
  */
  MatrixXd compute_gamma_matrix(const std::shared_ptr<const MatrixXd> distanceMatrix, unsigned int N) const;
  /*!
    @brief Return `_parameters`
  */
  Vector3d get_parameters() const;
  /*!
    @brief Set _parameters
    @param parameters Vector of parameters' values
  */
  void set_parameters(const Vector3d& parameters);
  /*!
    @brief Compute the values of the model covariogram at a given vector of distances, according to the variogram type
    @param h_vec The distances where to evaluate the covariogram
    @param card_h Number of distances where the covariogram must be computed. \f$ card\_h = h\_vec.size() \f$
    @return Vector of covariogram values at distances \f$ h \in  h\_vec \f$
  */
  Vec get_covario_vec(const Vec & h_vec, unsigned int card_h) const;
  /*!
    @brief Destructor
  */
  virtual ~FittedVariogram() = default;

};

/*!
  @brief	Class for computation and storage of the fitted variogram when `vario_model=="Gaussian"`
*/
class GaussVariogram : public FittedVariogram {
  void get_init_par(const EmpiricalVariogram &) override;
  MatrixXd compute_jacobian(const std::vector<double> &, unsigned int) const override;

public:
  double get_vario_univ(const double & h) const override;
  /*!
    @brief Destructor
  */
  ~GaussVariogram() = default;
};

/*!
  @brief	Class for computation and storage of the fitted variogram when `vario_model=="Exponential"`
*/
class ExpVariogram : public FittedVariogram {
  void get_init_par(const EmpiricalVariogram &) override;
  MatrixXd compute_jacobian(const std::vector<double> &, unsigned int) const override;

public:
  double get_vario_univ(const double & h) const override;
  /*!
    @brief Destructor
  */
  ~ExpVariogram() = default;
};

/*!
  @brief	Class for computation and storage of the fitted variogram when `vario_model=="Spherical"`
*/
class SphVariogram : public FittedVariogram {
  void get_init_par(const EmpiricalVariogram &) override;
  MatrixXd compute_jacobian(const std::vector<double> &, unsigned int) const override;

public:
  double get_vario_univ(const double & h) const override;
  /*!
    @brief Destructor
  */
  ~SphVariogram() = default;
};


}

#endif // _FITTED_VARIOGRAM_HPP
