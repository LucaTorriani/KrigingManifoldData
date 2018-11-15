#ifndef _MODEL_HPP_
#define _MODEL_HPP_

#include <vector>
#include <iostream>
#include <string>

#include <utility>
#include <Eigen/IterativeLinearSolvers>
#include <memory>
#include "Helpers.hpp"

/*! \file
  @brief Model class
*/

namespace model_fit {
  /*!
    @brief	Class to compute and store the linear model on the tangent space
  */
class Model {
private:
  /*! BigMatrix of the data on the tangent space */
  const std::shared_ptr<const MatrixXd> _data_tspace;
  /*! Design matrix for the data in the cell (used to compute the beta) */
  const std::shared_ptr<const MatrixXd> _design_matrix_model;
  /*! Design matrix for all the data in the domain (used ti compute the residuals) */
  const std::shared_ptr<const MatrixXd> _design_matrix_tot;
  /*! Name of the metric on the manifold */
  const std::string _distance_Manifold_name;

  /*! Number of stations in the cell */
  const unsigned int _N;
  /*! Dimension of the matrices on the manifold */
  const unsigned int _p;
  /*! Number of covariates in the model */
  const unsigned int _num_cov;
  /*! Number of significant entries in a symmetric \f$ \left(\text{\_p}*\text{\_p}\right) \f$ matrix.  \f$ \text{\_num\_coeff} = \frac{\text{\_p}*\left( \text{\_p}+1\right)}{2}  \f$ */
  const unsigned int _num_coeff;
  /*! \f$ \left(\text{\_num\_cov}*\text{\_num\_coeff}\right) \f$ matrix where the \f$i^{th}\f$ row contains the upper triangular part of \f$ \beta_{..i} \f$, the \f$i^{th}\f$ coefficient of the tangent space linear model  */
  MatrixXd _beta_matrix;
  /*! \f$ \left(\text{N\_tot}*\text{\_num\_coeff}\right) \f$ matrix where the \f$i^{th}\f$ row contains the upper triangular part of the matrix fitted in the \f$i^{th}\f$ location on the tangent space by the linear model. \f$ \mbox{\_fitted\_values} = \left(*\left(\mbox{\_design\_matrix\_tot}\right)\right)*\mbox{\_beta\_matrix} \f$ */
  MatrixXd _fitted_values;
  /*! \f$ \left(\text{N\_tot}*\text{\_num\_coeff}\right) \f$ matrix where the \f$i^{th}\f$ row contains the upper triangular part of the residual matrix in the \f$i^{th}\f$ location.  \f$  \mbox{\_residuals} = \left(*\left(\mbox{\_data\_tspace}\right)\right)-\mbox{\_fitted\_values} \f$  */
  MatrixXd _residuals;
public:
  /*!
    @brief Constructor
  */
  Model(const std::shared_ptr<const MatrixXd> data_tspace, const std::shared_ptr<const MatrixXd> design_matrix_model, unsigned int p, const std::string& distance_Manifold_name):  // EQUAL WEIGHTS
        _data_tspace(data_tspace), _design_matrix_model(design_matrix_model), _design_matrix_tot(design_matrix_model), _distance_Manifold_name(distance_Manifold_name),
       _N(design_matrix_model->rows()), _p(p), _num_cov(_design_matrix_model->cols()), _num_coeff((_p*(_p+1))/2), _beta_matrix(_num_cov, _num_coeff){};
  /*!
    @brief Constructor
  */
  Model(const std::shared_ptr<const MatrixXd> data_tspace, const std::shared_ptr<const MatrixXd> design_matrix_model,  const std::shared_ptr<const MatrixXd> design_matrix_tot, unsigned int p, const std::string& distance_Manifold_name): // KERNEL
        _data_tspace(data_tspace), _design_matrix_model(design_matrix_model), _design_matrix_tot(design_matrix_tot), _distance_Manifold_name(distance_Manifold_name),
        _N(design_matrix_model->rows()), _p(p), _num_cov(_design_matrix_model->cols()), _num_coeff((_p*(_p+1))/2), _beta_matrix(_num_cov, _num_coeff){};
  /*!
    @brief Update `_beta_matrix`, `_fitted_values` and `_residuals` according to the new covariogram matrix
    @param gamma_matrix Covariogram matrix
  */
  void update_model(const MatrixXd& gamma_matrix);
  /*!
    @brief Return `_beta_matrix`
  */
  MatrixXd get_beta() const;
  /*!
    @brief Return `_residuals`
  */
  MatrixXd get_residuals() const;
  /*!
    @brief Return `_fitted_values`
  */
  MatrixXd get_fitted_values() const;
};

}

#endif //_MODEL_HPP_
