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
  /*! Vector storing the three parameters of the fitted varriogram \f$ \left(\textit{nugget}, \textit{sill-nugget}, \textit{practical range} \right) \f$ */
  Vector3d _parameters; // tau2, sigma2, a
  /*!
    @brief Compute weighted median
    @param values Values whose median must be computed
    @param card Weights
    @return Median of \f$ \texttt{values} \f$, weighted using \f$ \texttt{card} \f$
  */
  double weighted_median (const std::vector<double> & values, const std::vector<unsigned int> & card);
  /*!
    @brief Initialize \f$\texttt{\_parameters}\f$
    @details The \f$\texttt{\_parameters}\f$ are initialized as follows: \n
      \f$ \text{\_parameters}\left(0\right)=\text{weighted\_median} \left(\text{\_emp\_vario\_values.head} \left(2\right), \text{\_N\_hvec.head} \left(2\right)\right) \f$ \n
      \f$ \text{\_parameters}\left(1\right)=\text{weighted\_median} \left(\text{\_emp\_vario\_values.tail} \left(4\right), \text{\_N\_hvec.tail} \left(4\right)\right) - \text{\_parameters}\left(0\right) \f$  \n
      \f$ \text{\_parameters}\left(2\right)= \text{\_hvec}\left(i^\ast \right) \f$  \n where \f$i^\ast\f$ is the first \f$i \text{ s.t.} \left| \text{\_emp\_vario\_values}\left(i\right)- \left( \text{\_parameters}\left(0\right) + \text{\_parameters}\left(1\right) \right) \right| < \text{tol} \f$
    @param emp_vario Empirical variogram
  */
  virtual void get_init_par(const EmpiricalVariogram & emp_vario) = 0;
  /*!
    @brief Update \f$\texttt{\_parameters}\f$ moving along \f$\texttt{dir}\f$
  */
  void backtrack(const Vector3d &dir,Vector3d &gk, Vec &res,const std::vector<double> & h_vec, unsigned int card_h, double c, double s, const Vec& emp_vario_values, double max_sill, double max_a);
  /*!
    @brief Compute the jacobian of the variogram residuals (which coincides with the one of the model variogram), according to the variogram type
  */
  virtual MatrixXd compute_jacobian(const std::vector<double> & h_vec, unsigned int card_h) const = 0;
public:
  /*!
    @brief Return \f$ \texttt{\_parameters}\left(0\right) \f$, i.e. the \f$ \textit{nugget}\f$
  */
  double get_tau2() const;
  /*!
    @brief Return \f$\texttt{\_parameters}\left(1\right)\f$, i.e. the \f$ \textit{sill-nugget}\f$
  */
  double get_sigma2() const;
  /*!
    @brief Return \f$\texttt{\_parameters}\left(2\right)\f$, i.e. the \f$ \textit{practical range}\f$
  */
  double get_a() const;
  /*!
    @brief Compute the parameters of the fitted variogram
    @note Like evaluate_par_fitted_W, but different stopping criteria. This function is used when equal weights are considered
    @details The parameters are computed using Gauss-Newton with backtrack method to solve the non-linear least square problem:
    \f{equation*}{
    \begin{aligned}
        & \underset{\text{\_parameters}}{\arg\min}
        & & \underset{h \in \_hvec} {\sum} \left( \gamma_m \left(\text{\_parameters}, h\right) - \widehat{\gamma} \left(h\right) \right)^2 \\
        & \text{subject to}
        & & 0 \leq \text{\_parameters}\left(0\right) \\
        &&& 0 \leq \text{\_parameters}\left(1\right) \leq \text{max\_sill} - \text{\_parameters}\left(0\right) \\
        &&& 0 \leq \text{\_parameters}\left(2\right) \leq \text{max\_a}
    \end{aligned}
    \f}
    The starting values for the \f$ \texttt{\_parameters} \f$ are obtained through ::variogram_evaluation::FittedVariogram::get_init_par. \n
    The stopping criteria is based on the decrease of the error norm.
    @param emp_vario Empirical variogram
    @param max_sill Maximum value for the \f$ \textit{sill} \f$
    @param max_a Maximum value for  \f$ \textit{a} \f$
  */
  void evaluate_par_fitted_E(const EmpiricalVariogram & emp_vario, double max_sill, double max_a);
  /*!
    @brief Compute the parameters of the fitted variogram
    @note Like evaluate_par_fitted_E, but different stopping criteria. This function is used when kernel weights are considered
    @details The parameters are computed using Gauss-Newton with backtrack method to solve the non-linear least square problem:
    \f{equation*}{
    \begin{aligned}
        & \underset{\text{\_parameters}}{\arg\min}
        & & \underset{h \in \_hvec} {\sum} \left( \gamma_m \left(\text{\_parameters}, h\right) - \widehat{\gamma} \left(h\right) \right)^2 \\
        & \text{subject to}
        & & 0 \leq \text{\_parameters}\left(0\right) \\
        &&& 0 \leq \text{\_parameters}\left(1\right) \leq \text{max\_sill} - \text{\_parameters}\left(0\right) \\
        &&& 0 \leq \text{\_parameters}\left(2\right) \leq \text{max\_a}
    \end{aligned}
    \f}
    The starting values for the \f$ \texttt{\_parameters} \f$ are obtained through ::variogram_evaluation::FittedVariogram::get_init_par.  \n
    The stopping criteria is based on the difference in the decrease of the error norm between two consecutive iterations.
    @param emp_vario Empirical variogram
    @param max_sill Maximum value for the \f$ \textit{sill} \f$
    @param max_a Maximum value for  \f$ \textit{a} \f$
  */
  void evaluate_par_fitted_W(const EmpiricalVariogram & emp_vario, double max_sill, double max_a);
  /*!
    @brief Compute the value of the model variogram at a given distance, according to the variogram type
    @param h Distance where to evaluate the variogram
    @return Variogram value at distance h: \f$ \gamma_m \left(h\right) \f$
  */
  virtual double get_vario_univ(const double & h) const = 0;
  /*!
    @brief Compute the value of the model covariogram at a given distance, according to the variogram type
    @param h Distance where to evaluate the covariogram
    @return Covariogram value at distance \f$ h: C_m \left(h\right) = \left(\text{\_parameters}\left(0\right) + \text{\_parameters}\left(1\right) \right) - \gamma_m \left(h\right) \f$
  */
  double get_covario_univ(const double & h) const;
  /*!
    @brief Compute the values of the model variogram at a given vector of distances, according to the variogram type
    @param h_vec Distances where to evaluate the variogram
    @param card_h Number of distances where the variogram must be computed. \f$ \text{card\_h} = \text{h\_vec.size()}  \f$
    @return Vector of variogram values at distances \f$h \in \text{\_h\_vec} \f$
  */
  Vec get_vario_vec(const std::vector<double> & h_vec, unsigned int card_h) const;
  /*!
    @brief Compute the values of the model variogram at a given vector of distances, according to the variogram type
    @param h_vec Distances where to evaluate the variogram
    @param card_h Number of distances where the variogram must be computed. \f$ \text{card\_h} = \text{h\_vec.size()}  \f$
    @return Vector of variogram values at distances \f$ h \in \text{h\_vec} \f$
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
    @brief Return \f$ \texttt{\_parameters} \f$
  */
  Vector3d get_parameters() const;
  /*!
    @brief Set _parameters
    @param parameters Vector of parameters' values
  */
  void set_parameters(const Vector3d& parameters);
  /*!
    @brief Compute the values of the model covariogram at a given vector of distances, according to the variogram type
    @param h_vec Distances where to evaluate the covariogram
    @param card_h Number of distances where the covariogram must be computed. \f$ \text{card\_h} = \text{h\_vec.size()}  \f$
    @return Vector of covariogram values at distances \f$ h \in \text{h\_vec} \f$
  */
  Vec get_covario_vec(const std::vector<double> & h_vec, unsigned int card_h) const;
  /*!
    @brief Destructor
  */
  virtual ~FittedVariogram() = default;

};

/*!
  @brief	Class for computation and storage of the fitted variogram when \f$\texttt{vario\_model=="Gaussian"}\f$
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
  @brief	Class for computation and storage of the fitted variogram when \f$\texttt{vario\_model=="Exponential"}\f$
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
  @brief	Class for computation and storage of the fitted variogram when \f$\texttt{vario\_model=="Spherical"}\f$
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
