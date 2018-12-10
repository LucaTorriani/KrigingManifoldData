#ifndef _HELPERS_HPP_
#define _HELPERS_HPP_

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <cmath>
#include <string>

/*! \file
  @brief Functions to manipulate matrices
*/

/*! \var typedef Eigen::VectorXd Vec
  @brief Eigen dynamic vector of doubles
*/
typedef Eigen::VectorXd Vec;
// typedef Eigen::Triplet<double> TripType;


using namespace Eigen;

namespace matrix_manipulation {
  /*!
    @brief Compute the exponential of a matrix
    @param A Matrix whose exponential must be computed
    @return Expontial of A
  */
  MatrixXd expMat(const MatrixXd& A);
  /*!
    @brief Compute the logarithm of a matrix
    @param A Matrix whose logarithm must be computed
    @return Logarithm of A
  */
  MatrixXd logMat(const MatrixXd& A);
  /*!
    @brief Compute the square root of a matrix
    @param A Matrix whose square root must be computed
    @return Square root of A
  */
  MatrixXd sqrtMat(const MatrixXd& A);

  /*!
    @brief Transform a \f$ \left(N, \frac{p*\left(p+1\right)}{2}\right)\f$ matrix (where each row represents a symmetric matrix) in a vector of length \f$ N \f$ of \f$ p*p\f$symmetric matrices
    @param bigMatrix Matrix whose rows store the values of the upper triangular parts of \f$ N p*p \f$ symmetric matrices
    @param p Dimension of the matrices
    @param distance_Manifold_name Name of the metric on the manifold to use. If `distance_Manifold_name=="Correlation"` then only the upper triangular parts of the matrices is filled
    @return Vector of length \f$ N \f$ of \f$ p*p\f$ symmetric matrices corresponding to the transformation of bigMatrix
  */
  std::vector<MatrixXd> bigMatrix2VecMatrices(const MatrixXd& bigMatrix, unsigned int p, const std::string& distance_Manifold_name);
  /*!
    @brief Transform a vector of length \f$ N \f$ of \f$ p*p\f$ symmetric matrices in a \f$ \left(N,\frac{p*\left(p+1\right)}{2}\right) \f$ matrix (where each row represents a symmetric matrix)
    @param vecMatrices Vector of symmetric matrices whose upper trinagular parts will be stored in the rows of the output matrix
    @return Matrix \f$ \left(\frac{p*\left(p+1\right)}{2}, N\right)\f$ correponsing to the transformation of vecMatrices
  */
  MatrixXd VecMatrices2bigMatrix(const std::vector<MatrixXd>& vecMatrices);
  /*!
    @brief Compute the Cholesky decomposition of a matrix not positive definite
    @param M1 Positive semidefinite matrix whose Cholesky decomposition must be computed
    @return Matrix in \f$ Chol\left(p\right) \f$, representing the upper trianglar matrix of A's Cholesky decomposition
  */
  MatrixXd Chol_semidef (const MatrixXd& M1);
  /*!
    @brief Compute the Cholesky decomposition of a matrix
    @param M1 Positive definite or semidefine matrix whose Cholesky decomposition must be computed
    @return Matrix in \f$ Chol\left(p\right) \f$, representing the upper trianglar matrix of A's Cholesky decomposition
  */
  MatrixXd Chol_decomposition (const MatrixXd M1);

}


#endif // _HELPERS_HPP_
