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

namespace HomogenizationConstants
{
  typedef std::vector<std::pair<unsigned int, unsigned int>> index_list;
  const std::map<bool, std::vector<index_list>> indices {
  {true, {
    {{0,0}},
    {{0,0},{1,1},{1,0},{0,1}},
    {{0,0},{1,0},{2,0},{0,1},{1,1},{2,1},{2,0},{2,1},{2,2}}
         }},
  {false, {
    {{0,0}},
    {{0,0},{1,1},{0,1}},
    {{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}}
          }}};
  const std::map<bool, std::vector<unsigned int>> required {
    {true, {1, 4, 9}},
    {false, {1, 3, 6}}};
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

  const bool _ld;
  unsigned int _ndisp;
  unsigned int _ncomps;

  const MaterialProperty<RankTwoTensor> &_stress;
  const MaterialProperty<RankFourTensor> &_material_jacobian;
  const MaterialProperty<RankTwoTensor> &_F;

  std::vector<const Function*> _targets;

  enum class ConstraintType { Stress, Strain };
  std::vector<ConstraintType> _ctypes;

  unsigned int _qp;
  unsigned int _h;
  unsigned int _hh;

  HomogenizationConstants::index_list _indices;

  RankTwoTensor _residual;
  RankFourTensor _jacobian;
};
