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
  unsigned int _N;
  /*! Matrix of distances among the \f$N\f$ locations */
  std::shared_ptr<const MatrixXd> _distanceMatrix;

  /*! Vector storing the variogram estimates. \f$\mbox{\_emp\_vario\_values[i]} = \hat{\gamma}\left(\mbox{\_hvec[i]}\right)\f$  */
  std::vector<double> _emp_vario_values;
  /*! Vector storing the distances at which the variogram is estimated */
  std::vector<double> _hvec;
  /*! Vector storing the number of data used in the estimation of the corresponding empirical variogram value */
  std::vector<unsigned int> _N_hvec;
  /*! Number of distances for which the variogram has been evaluated.  \f$\mbox{\_card\_h} = \mbox{\_hvec.size()}\f$  */
  unsigned int _card_h;
  /*! Vector of equispaced distances whose midpoints are the candidates to enter in \f$\_hvec\f$  */
  Vec _d; // Vettore h+-deltah
  /*! Maximum distance considered */
  double _hmax;
  /*! Vector of the weights for the \f$N\f$ data */
  Vec _weights;
  /*!
    @brief Compute the maximum distance to be considered
    @details `_hmax` is computed as one third of the diagonal of the bouding box for the data
    @param coords Data coordinates
    @param distance Geographical distance
  */
  void compute_hmax(const Coordinates& coords, const distances::Distance& distance);

public:
  /*!
    @brief Constructor
  */
  EmpiricalVariogram (unsigned int);

  /*!
    @brief Update `_emp_vario_values`, `_hvec` and `_N_hvec` according to the new residuals \f$ \Delta \left(\boldsymbol{s_{i}}\right) i=1,\ldots,\_N \f$.
    @details `_emp_vario_values` are estimated through the method of moments:
    \f[
      \widehat{\gamma}\left(h\right) = \frac{\underset{N\left(h\right)}{\sum} w\left(\boldsymbol{s_{i}}\right)w\left(\boldsymbol{s_{j}}\right)\left\Vert \Delta \left(\boldsymbol{s_{i}}\right)-\Delta \left(\boldsymbol{s_{j}}\right)\right\Vert ^{2} } {2\underset{N\left(h\right)}{\sum} w\left(\boldsymbol{s_{i}}\right)w\left(\boldsymbol{s_{j}}\right) }
    \f]
      where \f$ N\left(h\right) = \left\lbrace \left(\boldsymbol{s}_{i}, \boldsymbol{s}_{j} \in D\right)\, : \, h - \varDelta h \leq \left\Vert \boldsymbol{s}_{i} - \boldsymbol{s}_{j} \right\Vert \leq h + \varDelta h  \right\rbrace \f$
    @param res Vector of the \f$N\f$ residual matrices
    @param distanceTplane Distance on the tangent space
  */
  void update_emp_vario(const std::vector<MatrixXd>& res, const distances_tplane::DistanceTplane & distanceTplane);
  /*!
    @brief Return `_emp_vario_values`
  */
  std::vector<double> get_emp_vario_values () const;
  /*!
    @brief Return `_N_hvec`
  */
  std::vector<unsigned int> get_N_hvec() const;
  /*!
    @brief Return `_hvec`
  */
  std::vector<double> get_hvec() const;
  /*!
    @brief Return `_card_h`
  */
  unsigned int get_card_h() const;
  /*!
    @brief Return `_N`
  */
  unsigned int get_N() const;
  /*!
    @brief Return `_hmax`
  */
  double get_hmax() const;
  /*!
    @brief Set `_distanceMatrix`, `_d` and `_hmax`
    @details `_hmax` is computed using `compute_hmax`
    @param distanceMatrix Matrix of distances among the \f$N\f$ locations
    @param coords Data coordinates
    @param distance Geographical distance
  */
  void set_distance_and_h_max(const std::shared_ptr<const Eigen::MatrixXd> distanceMatrix, const Coordinates& coords, const distances::Distance& distance);
  /*!
    @brief Set `_distanceMatrix`, `_d` and `_hmax`
    @details `_hmax` is computed as one third of `max_dist`
    @param distanceMatrix Matrix of distances among the \f$N\f$ locations
    @param max_dist Maximum distance between data locations
  */
  void set_distance_and_h_max(const std::shared_ptr<const Eigen::MatrixXd> distanceMatrix, const double& max_dist);
  /*!
    @brief Set `_N` and `_weights`
    @param N Number \f$N\f$ of locations
    @param weights Weights for the \f$N\f$ locations used in the computation of the empirical variogram
  */
  void set_weights(unsigned int N, const Vec& weights);
};

}


#endif // _EMPIRICAL_VARIOGRAM_HPP
