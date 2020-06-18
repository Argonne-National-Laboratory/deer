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

#include "HomogenizationConstraintKernel.h"

class HomogenizationConstraintIntegral;

class HomogenizationConstraintIntegral : public ElementUserObject
{
 public:
  static InputParameters validParams();

  HomogenizationConstraintIntegral(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

  virtual Real getResidual(unsigned int h) const;
  virtual RankTwoTensor getJacobian(unsigned int h) const;

 protected:
  virtual Real computeResidual();
  virtual RankTwoTensor computeJacobian();

  const bool _ld;

  unsigned int _ndisp;
  unsigned int _num_hvars;

  const MaterialProperty<RankTwoTensor> &_stress;
  const MaterialProperty<RankFourTensor> &_material_jacobian;
  const MaterialProperty<RankTwoTensor> &_F;

  std::vector<const Function*> _targets;

  enum class ConstraintType { Stress, Strain };
  std::vector<ConstraintType> _ctypes;

  unsigned int _qp;
  unsigned int _h;

  std::vector<Real> _residual;
  std::vector<RankTwoTensor> _jacobian;
  
  HomogenizationConstants::index_list _indices;
  
};
