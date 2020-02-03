//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PureElasticCZM.h"

class ViscousSlidingCZM;
template <> InputParameters validParams<ViscousSlidingCZM>();
/**
 * Implementation of the non-stateful exponential traction separation law
 * proposed by Salehani, Mohsen Khajeh and Irani, Nilgoon 2018
 **/
class ViscousSlidingCZM : public PureElasticCZM {
public:
  ViscousSlidingCZM(const InputParameters &parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  virtual void ComputeShearTraction(RealVectorValue &traction) override;
  virtual void
  ComputeShearTractionDerivatives(RankTwoTensor &traction_derivatives) override;

  virtual Real ComputeShearViscosity();
  virtual RankTwoTensor ComputeShearViscosityDerivatives();

  /// the displacement jump associated to the maximum traction
  const Real _shear_viscosity;

  /// the value of the old traction in local coordinates
  ///@{
  const MaterialProperty<RealVectorValue> &_traction_old;
  ///@}

  /// the coupled displacement and neighbor displacement values
  ///@{
  MaterialProperty<RealVectorValue> &_displacement_jump_dot;
  MaterialProperty<RealVectorValue> &_displacement_jump_global_dot;
  ///@}

  /// the coupled displacement and neighbor displacement values
  ///@{
  std::vector<const VariableValue *> _disp_dot;
  std::vector<const VariableValue *> _disp_neighbor_dot;
  ///@}
};
