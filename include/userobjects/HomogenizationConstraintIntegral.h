//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementUserObject.h"

// Helpers common to the whole homogenization system
namespace HomogenizationConstants
{
  typedef std::vector<std::pair<unsigned int, unsigned int>> index_list;
  const std::map<bool, std::vector<index_list>> indices {
  {true, {
    {{0,0}},
    {{0,0},{1,1},{1,0},{0,1}},
    {{0,0},{1,0},{2,0},{0,1},{1,1},{2,1},{0,2},{1,2},{2,2}}
         }},
  {false, {
    {{0,0}},
    {{0,0},{1,1},{0,1}},
    {{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}}
          }}};
  const std::map<bool, std::vector<unsigned int>> required {
    {true, {1, 4, 9}},
    {false, {1, 3, 6}}};
  enum class ConstraintType { Stress, Strain };

  inline ConstraintType map_string(std::string input) {
    if ((input == "strain") || (input == "Strain"))
      return ConstraintType::Strain;
    else if ((input == "stress") || (input == "Stress"))
      return ConstraintType::Stress;
    else
      mooseError("Constraint type must be either stress or strain");
  }
}

class HomogenizationConstraintIntegral : public ElementUserObject
{
 public:
  static InputParameters validParams();

  HomogenizationConstraintIntegral(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

  virtual const RankTwoTensor & getResidual() const;
  virtual const RankFourTensor & getJacobian() const;

 protected:
  virtual RankTwoTensor computeResidual();
  virtual RankFourTensor computeJacobian();


  Real sdStressJacobian(unsigned int i, unsigned int j, 
                        unsigned int a, unsigned int b);
  Real sdStrainJacobian(unsigned int i, unsigned int j,
                        unsigned int a, unsigned int b);
  Real ldStressJacobian(unsigned int i, unsigned int j, 
                        unsigned int a, unsigned int b);
  Real ldStrainJacobian(unsigned int i, unsigned int j,
                        unsigned int a, unsigned int b);

  const bool _ld;
  unsigned int _ndisp;
  unsigned int _ncomps;

  const MaterialProperty<RankTwoTensor> &_stress;
  const MaterialProperty<RankFourTensor> &_material_jacobian;
  const MaterialProperty<RankTwoTensor> &_F;
  const MaterialProperty<RankTwoTensor> &_PK1;
  const MaterialProperty<Real> &_J;
  const MaterialProperty<RankTwoTensor> & _invF;
  const MaterialProperty<RankTwoTensor> & _df;

  std::vector<const Function*> _targets;
  std::vector<HomogenizationConstants::ConstraintType> _ctypes;

  unsigned int _qp;
  unsigned int _h;
  unsigned int _hh;

  HomogenizationConstants::index_list _indices;

  RankTwoTensor _residual;
  RankFourTensor _jacobian;
};
