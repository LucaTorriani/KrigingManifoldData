#ifndef _EMPIRICAL_VARIOGRAM_HPP
#define _EMPIRICAL_VARIOGRAM_HPP

#include "Helpers.hpp"
#include "DistanceTplane.hpp"
#include "Distance.hpp"
#include "Coordinates.hpp"

/*! \file
  @brief Class to compute the empirical variogram
*/

namespace variogram_evaluation{

  /*!
    @brief	Class for computation and storage of the empirical variogram
  */
class EmpiricalVariogram {
  /*! Number of bins */
  const unsigned int _n_h;
  /*! Number of data points used to compute the empirical variogram */
  const unsigned int _N;
  /*! Matrix of distances among the N locations */
  const std::shared_ptr<const MatrixXd> _distanceMatrix;

  /*! Vector storing the variogram estimates. \f$mbox{\_emp_vario_values[i]} = \hat{\gamma}\left(\mbox{\_hvec[i]}\right)\f$  */
  std::vector<double> _emp_vario_values;
  /*! Vector storing the distances at which the variogram is estimated */
  std::vector<double> _hvec;
  /*! Vector storing the number of data used in the estimation of the corresponding empirical variogram value */
  std::vector<unsigned int> _N_hvec;
  /*! Number of distances for which the variogram has been evaluated. _card_h = _hvec.size() */
  unsigned int _card_h;
  /*! Vector of equispaced distances whose midpoints are the candidates to enter in _hvec */
  Vec _d; // Vettore h+-deltah
  /*! Maximum distance considered */
  double _hmax;
  /*! Vector of the weights for the N data */
  Vec _weights;
  /*!
    @brief Compute the maximum distance to be considered
  */
  void compute_hmax(const Coordinates&, const distances::Distance&);

public:
  /*!
    @brief Constructor
  */
  EmpiricalVariogram (const std::shared_ptr<const MatrixXd>, unsigned int,  // EQUAL WEIGHTS
                      const Coordinates&, const distances::Distance&);
  /*!
    @brief Constructor
  */
  EmpiricalVariogram (const std::shared_ptr<const MatrixXd>, unsigned int, unsigned int,  // KERNEL
                      const Vec&, double);

  /*!
    @brief Update \f$\texttt{\_emp\_vario\_values}\f$, \f$\texttt{\_hvec}\f$ and \f$\texttt{\_N\_hvec}\f$ according to the new residuals \f$ \Delta \left(\boldsymbol{s_{i}}\right) i=1,\cdot,\_N \f$.
    @details \f$\texttt{\_emp\_vario\_values}\f$ are estimated through the method of moments:
    \f[
      \widehat{\gamma}\left(h\right) = \frac{\underset{N\left(h\right)}{\sum} w\left(\boldsymbol{s_{i}}\right)w\left(\boldsymbol{s_{j}}\right)\left\Vert \Delta \left(\boldsymbol{s_{i}}\right)-\Delta \left(\boldsymbol{s_{j}}\right)\right\Vert ^{2} } {2\underset{N\left(h\right)}{\sum} w\left(\boldsymbol{s_{i}}\right)w\left(\boldsymbol{s_{j}}\right) }
    \f]
      where \f$ N\left(h\right) = \left\lbrace \left(\boldsymbol{s}_{i}, \boldsymbol{s}_{j} \in D\right)\, : \, h - \varDelta h \leq \left\Vert \boldsymbol{s}_{i} - \boldsymbol{s}_{j} \right\Vert \leq h - \varDelta h  \right\rbrace \f$
    @param res Vector of the N residual matrices
    @param distanceTplane Distance on the tangent space
  */
  void update_emp_vario(const std::vector<MatrixXd>& res, const distances_tplane::DistanceTplane & distanceTplane);
  /*!
    @brief Return \f$\texttt{\_emp\_vario\_values}\f$
  */
  std::vector<double> get_emp_vario_values () const;
  /*!
    @brief Return \f$\texttt{\_N\_hvec}\f$
  */
  std::vector<unsigned int> get_N_hvec() const;
  /*!
    @brief Return \f$\texttt{\_hvec}\f$
  */
  std::vector<double> get_hvec() const;
  /*!
    @brief Return \f$\texttt{\_card\_h}\f$
  */
  unsigned int get_card_h() const;
  /*!
    @brief Return \f$\texttt{\_N}\f$
  */
  unsigned int get_N() const;
  /*!
    @brief Return \f$\texttt{\_hmax}\f$
  */
  double get_hmax() const;
};

}


#endif // _EMPIRICAL_VARIOGRAM_HPP
