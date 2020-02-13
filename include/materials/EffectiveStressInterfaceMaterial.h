//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceMaterial.h"
#include "RankTwoTensor.h"
class EffectiveStressInterfaceMaterial;
template <> InputParameters validParams<EffectiveStressInterfaceMaterial>();

/**
 * EffectiveStressInterfaceMaterial uses the namespace EffectiveStressTools to
 * compute scalar values from Rank-2 tensors.
 */
class EffectiveStressInterfaceMaterial : public InterfaceMaterial {
public:
  EffectiveStressInterfaceMaterial(const InputParameters &parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  const MaterialProperty<RankTwoTensor> &_stress_master;
  const MaterialProperty<RankTwoTensor> &_stress_slave;
  MaterialProperty<Real> &_effective_stress;

  /// the type offective stress to be computed
  MooseEnum _effective_stress_type;
  /// the type offective stress to be computed
  MooseEnum _interface_value_type;
  /// vector of paramters
  std::vector<Real> _params_vector;

  const Point _point1;
  const Point _point2;
  Point _input_direction;
  const bool _stateful;
};
