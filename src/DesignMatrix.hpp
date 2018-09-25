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

namespace design_matrix {


  class DesignMatrix {

  public:
    virtual MatrixXd compute_design_matrix(const Coordinates&) const = 0;
    virtual MatrixXd compute_design_matrix(const Coordinates&, const MatrixXd&) const =0;

  };

  // class DesignMatrixFactory{
  //   typedef std::string Identifier;
  //   typedef std::function<std::unique_ptr<DesignMatrix>()> Builder;
  //   DesignMatrixFactory() = default;
  //   DesignMatrixFactory(const DesignMatrixFactory&) = delete;
  //   DesignMatrixFactory& operator=(const DesignMatrixFactory&) = delete;
  //   std::map<Identifier, Builder> _storage;
  //
  // public:
  //   static DesignMatrixFactory& Instance();
  //   std::unique_ptr<DesignMatrix> create(const Identifier&) const;
  //   void add(const Identifier&, const Builder&);
  //
  // };

  class InterceptDM : public DesignMatrix {
  public:
    MatrixXd compute_design_matrix(const Coordinates&) const override;
    MatrixXd compute_design_matrix(const Coordinates&, const MatrixXd&) const override;
    std::unique_ptr<DesignMatrix> operator()(){return (std::unique_ptr<DesignMatrix>(new InterceptDM));}
  };

  class Coord1DM: public DesignMatrix {
  public:
    MatrixXd compute_design_matrix(const Coordinates&) const override;
    MatrixXd compute_design_matrix(const Coordinates&, const MatrixXd&) const override;
    std::unique_ptr<DesignMatrix> operator()(){return (std::unique_ptr<DesignMatrix>(new Coord1DM));}
  };

  class Coord2DM: public DesignMatrix {
  public:
    MatrixXd compute_design_matrix(const Coordinates&) const override;
    MatrixXd compute_design_matrix(const Coordinates&, const MatrixXd&) const override;
    std::unique_ptr<DesignMatrix> operator()(){return (std::unique_ptr<DesignMatrix>(new Coord2DM));}
  };

  class AdditiveDM: public DesignMatrix {
  public:
    MatrixXd compute_design_matrix(const Coordinates&)const override;
    MatrixXd compute_design_matrix(const Coordinates&, const MatrixXd&) const override;
    std::unique_ptr<DesignMatrix> operator()(){return (std::unique_ptr<DesignMatrix>(new AdditiveDM));}
  };

  // void registerDesignMatrices();

};

#endif
