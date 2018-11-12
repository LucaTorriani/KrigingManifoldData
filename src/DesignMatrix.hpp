#ifndef _DESIGN_MATRIX_HPP_
#define _DESIGN_MATRIX_HPP_

#include <memory>
#include<iostream>
#include<utility>
#include<map>
#include<string>
#include <Eigen/Dense>
#include "Coordinates.hpp"

using namespace Eigen;

/*! \file
  @brief Classes to create the design_matrix according to the model on the tangent space
*/


namespace design_matrix {

  /*!
    @brief	Abstract class for the computation of the design_matrix
  */
  class DesignMatrix {
  public:
    /*!
      @brief Compute the design matrix with no covariates besides the coordinates
      @param coords Matrix of the coordinates of the locations
      @return Design matrix
    */
    virtual MatrixXd compute_design_matrix(const Coordinates& coords) const = 0;

    /*!
      @brief Compute the design matrix with extra covariates besides the coordinates
      @param coords Matrix of the coordinates of the locations
      @param X Matrix of the additional covariates for the locations
      @return Design matrix
    */
    virtual MatrixXd compute_design_matrix(const Coordinates& coords, const MatrixXd& X) const =0;
    /*!
      @brief Destructor
    */
    virtual ~DesignMatrix() = default;
  };

  /*!
    @brief	Class for the computation of the design_matrix when \f$\texttt{model\_ts=="Intercept"}\f$
  */
  class InterceptDM : public DesignMatrix {
  public:
    MatrixXd compute_design_matrix(const Coordinates& coords) const override;
    MatrixXd compute_design_matrix(const Coordinates& coords, const MatrixXd& X) const override;
    // std::unique_ptr<DesignMatrix> operator()(){return (std::unique_ptr<DesignMatrix>(new InterceptDM));}
    /*!
      @brief Destructor
    */
    ~InterceptDM() = default;
  };

  /*!
    @brief	Class for the computation of the design_matrix when \f$\texttt{model\_ts=="Coord1"}\f$
  */
  class Coord1DM: public DesignMatrix {
  public:
    MatrixXd compute_design_matrix(const Coordinates& coords) const override;
    MatrixXd compute_design_matrix(const Coordinates& coords, const MatrixXd& X) const override;
    // std::unique_ptr<DesignMatrix> operator()(){return (std::unique_ptr<DesignMatrix>(new Coord1DM));}
    /*!
      @brief Destructor
    */
    ~Coord1DM() = default;
  };

  /*!
    @brief	Class for the computation of the design_matrix when \f$\texttt{model\_ts=="Coord2"}\f$
  */
  class Coord2DM: public DesignMatrix {
  public:
    MatrixXd compute_design_matrix(const Coordinates& coords) const override;
    MatrixXd compute_design_matrix(const Coordinates& coords, const MatrixXd& X) const override;
    // std::unique_ptr<DesignMatrix> operator()(){return (std::unique_ptr<DesignMatrix>(new Coord2DM));}
    /*!
      @brief Destructor
    */
    ~Coord2DM() = default;
  };

  /*!
    @brief	Class for the computation of the design_matrix when \f$\texttt{model\_ts=="Additive"}\f$
  */
  class AdditiveDM: public DesignMatrix {
  public:
    MatrixXd compute_design_matrix(const Coordinates& coords)const override;
    MatrixXd compute_design_matrix(const Coordinates& coords, const MatrixXd& X) const override;
    // std::unique_ptr<DesignMatrix> operator()(){return (std::unique_ptr<DesignMatrix>(new AdditiveDM));}
    /*!
      @brief Destructor
    */
    ~AdditiveDM() = default;
  };
}

#endif  // _DESIGN_MATRIX_HPP_
